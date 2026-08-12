#!/usr/bin/env python3
"""[eval · ColabFold AF2-Multimer] 한 복합체의 5 pose -> 각 DockQ + 신뢰도 -> summary.csv.

eval_chai.py를 본떴다. 차이:
  - 입력 = PDB (out_colabfold/<id>_unrelaxed_rank_*_model_*_seed_*.pdb, 보통 5개)
  - 신뢰도 = ColabFold scores json (iptm=global 스칼라, ptm, plddt 리스트)
  - AF2-Multimer는 per-chain-pair ipTM이 없음 -> pair_iptm_abag = "" (NA)
  - pose 식별 = model 번호(seed 열에 기록), rank는 notes에.
chains.json(antigen/antibody/order)로 역할·사슬순서. ColabFold는 입력 order대로 A/B/C 라벨.
run-dir 안: chains.json, <ID>.pdb/.cif(정답), out_colabfold/*.pdb + *_scores_*.json
사용(pipeline 루트에서, DockQ 설치된 env):
  python scripts/eval_colabfold.py --target 10ZO --run-dir runs/10zo \
      --phase Benchmark --model AF2-Multimer --csv results/summary.csv
필요: DockQ(CLI). lib_log(같은 scripts/ 폴더).
"""
import argparse, glob, json, os, re, subprocess, sys
import lib_log


def find_one(patterns):
    for p in patterns:
        h = sorted(glob.glob(p, recursive=True))
        if h:
            return h[0]
    return None


def _interfaces(out):
    inter, cur = {}, None
    for line in out.splitlines():
        m = re.search(r"Native chains:\s*([A-Za-z0-9]),\s*([A-Za-z0-9])", line)
        if m:
            cur = frozenset([m.group(1), m.group(2)])
            continue
        m = re.search(r"^\s*DockQ:\s*([\d.]+)", line)
        if m and cur is not None:
            inter[cur] = float(m.group(1))
            cur = None
    return inter


def parse_dockq(model_pdb, native, mapping=None):
    def run(mp):
        cmd = ["DockQ", model_pdb, native]
        if mp:
            cmd += ["--mapping", mp]
        return subprocess.run(cmd, capture_output=True, text=True).stdout
    out = run(mapping)                          # 예측 사슬 -> native 같은 사슬 고정
    inter = _interfaces(out)
    if not inter and mapping:                   # 다중카피 비-cognate -> auto-map 재시도
        out2 = run(None)
        i2 = _interfaces(out2)
        if i2:
            return i2, out2
    return inter, out


def ab_ag_dockq(inter, antigen, antibody):
    abag = {ab: inter.get(frozenset([antigen, ab])) for ab in antibody}
    have = [v for v in abag.values() if v is not None]
    if not have and len(antibody) == 1 and len(inter) == 1:
        # 나노바디 + 다중카피: 유일 계면이 곧 Ab-Ag
        only = next(iter(inter.values()))
        abag = {antibody[0]: only}
        have = [only]
    return abag, have


def parse_conf(scores_json):
    """ColabFold scores json: iptm(global 스칼라) · ptm · plddt(리스트). pair-ipTM 없음."""
    with open(scores_json) as f:
        c = json.load(f)
    giptm = c.get("iptm")
    ptm = c.get("ptm")
    pl = c.get("plddt")
    mplddt = (sum(pl) / len(pl)) if isinstance(pl, list) and pl else None
    return giptm, ptm, mplddt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--run-dir", default=".")
    ap.add_argument("--results", default="out_colabfold")   # run-dir 하위 출력 폴더
    ap.add_argument("--phase", default="Benchmark")
    ap.add_argument("--model", default="AF2-Multimer")
    ap.add_argument("--csv", default="results/summary.csv")
    a = ap.parse_args()
    rd, tid = a.run_dir, a.target

    cj = os.path.join(rd, "chains.json")
    if not os.path.exists(cj):
        sys.exit(f"[에러] chains.json 없음: {cj} (extract_chains 먼저)")
    d = json.load(open(cj))
    antigen, antibody, order = d["antigen"], d["antibody"], d["order"]

    native = find_one([f"{rd}/{tid}.pdb", f"{rd}/{tid}.cif",
                       f"{rd}/{tid.upper()}.pdb", f"{rd}/{tid.upper()}.cif",
                       f"{rd}/{tid.lower()}.pdb", f"{rd}/{tid.lower()}.cif"])
    if not native:
        sys.exit(f"[에러] native 구조 없음 (run-dir={rd})")

    resdir = os.path.join(rd, a.results)
    pdbs = sorted(glob.glob(f"{resdir}/*_unrelaxed_rank_*.pdb"))
    if not pdbs:
        sys.exit(f"[에러] pose pdb 없음: {resdir}/*_unrelaxed_rank_*.pdb")

    _m = "".join(order)
    n_ok, best = 0, None
    for pdb in pdbs:
        base = os.path.basename(pdb)
        mR = re.search(r"rank_(\d+)", base)
        mM = re.search(r"model_(\d+)", base)
        rank = int(mR.group(1)) if mR else 0
        model_no = int(mM.group(1)) if mM else rank
        scores = pdb.replace("_unrelaxed_", "_scores_").rsplit(".pdb", 1)[0] + ".json"
        giptm = ptm = mplddt = None
        if os.path.exists(scores):
            giptm, ptm, mplddt = parse_conf(scores)

        inter, dockq_out = parse_dockq(pdb, native, mapping=f"{_m}:{_m}")
        abag, have = ab_ag_dockq(inter, antigen, antibody)
        if not have:
            print(f"[경고] {base}: Ab-Ag DockQ 파싱 실패 "
                  f"(interfaces={[tuple(k) for k in inter]}) -> 건너뜀")
            continue
        dockq_min = min(have)
        success = dockq_min >= 0.23
        n_ok += 1
        best = dockq_min if best is None else max(best, dockq_min)

        detail = " / ".join(f"{antigen}-{ab} {abag[ab]:.3f}" if abag[ab] is not None
                            else f"{antigen}-{ab} NA" for ab in antibody)
        notes = (f"rank{rank} model{model_no}; DockQ(Ab-Ag) {detail}; "
                 f"pair ipTM NA(AF2-M); ptm {round(ptm,3) if ptm is not None else 'NA'}; "
                 f"mean pLDDT {round(mplddt,1) if mplddt is not None else 'NA'}")
        row = {"phase": a.phase, "model": a.model, "target": tid, "seed": model_no,
               "dockq_abag": round(dockq_min, 3),
               "global_iptm": round(giptm, 3) if giptm is not None else "",
               "pair_iptm_abag": "",            # AF2-Multimer는 pair-ipTM 없음
               "success": success, "notes": notes}
        lib_log.append_row(a.csv, row)
        print(f"  {base}: DockQ={dockq_min:.3f}  ipTM={giptm}  success={success}")

    print(f"=== {tid}: {n_ok} pose 기록 -> {a.csv} | best DockQ = "
          f"{round(best,3) if best is not None else 'NA'} "
          f"({'성공' if (best is not None and best >= 0.23) else '실패'})")


if __name__ == "__main__":
    main()
