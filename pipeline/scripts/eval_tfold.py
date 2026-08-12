#!/usr/bin/env python3
"""[eval · tFold-Ag] pose(PDB) -> DockQ + confidence(REMARK 250) -> summary.csv.

tFold 사슬 = H/L/A. native(chains.json order = 항체..+항원)와 위치 대응 -> DockQ --mapping "HLA:<order>".
Confidence는 별도 파일 없이 PDB REMARK 250 헤더(ipTM/pTM/lDDT-Ca). pair-ipTM 없음(global만).
pose = seed별 PDB 여러 개(out_tfold/*_seed*.pdb).
사용(pipeline 루트, DockQ 설치된 env):
  python scripts/eval_tfold.py --target 10ZO --run-dir runs/10zo --model tFold-Ag --csv results/summary.csv
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
    out = run(mapping)
    inter = _interfaces(out)
    if not inter and mapping:
        out2 = run(None)
        i2 = _interfaces(out2)
        if i2:
            return i2, out2
    return inter, out


def ab_ag_dockq(inter, antigen, antibody):
    abag = {ab: inter.get(frozenset([antigen, ab])) for ab in antibody}
    have = [v for v in abag.values() if v is not None]
    if not have and len(antibody) == 1 and len(inter) == 1:
        only = next(iter(inter.values()))
        abag = {antibody[0]: only}
        have = [only]
    return abag, have


def parse_conf(pdb_path):
    """tFold confidence = PDB REMARK 250 헤더 (별도 json/npz 없음)."""
    txt = open(pdb_path).read()
    def grab(label):
        m = re.search(rf"REMARK 250 {re.escape(label)}:\s*([\d.]+)", txt)
        return float(m.group(1)) if m else None
    return {"iptm": grab("Predicted ipTM score"),
            "ptm": grab("Predicted pTM score"),
            "plddt": grab("Predicted lDDT-Ca score"),
            "ag_plddt": grab("Antigen Predicted lDDT-Ca score")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--run-dir", default=".")
    ap.add_argument("--results", default="out_tfold")   # run-dir 하위 pose 폴더
    ap.add_argument("--phase", default="Benchmark")
    ap.add_argument("--model", default="tFold-Ag")
    ap.add_argument("--csv", default="results/summary.csv")
    a = ap.parse_args()
    rd, tid = a.run_dir, a.target

    cj = os.path.join(rd, "chains.json")
    if not os.path.exists(cj):
        sys.exit(f"[에러] chains.json 없음: {cj}")
    d = json.load(open(cj))
    antigen, antibody, order = d["antigen"], d["antibody"], d["order"]

    native = find_one([f"{rd}/{tid}.pdb", f"{rd}/{tid}.cif",
                       f"{rd}/{tid.upper()}.pdb", f"{rd}/{tid.upper()}.cif",
                       f"{rd}/{tid.lower()}.pdb", f"{rd}/{tid.lower()}.cif"])
    if not native:
        sys.exit(f"[에러] native 없음 (run-dir={rd})")

    pdbs = sorted(glob.glob(f"{rd}/{a.results}/*.pdb"))
    if not pdbs:
        sys.exit(f"[에러] pose pdb 없음: {rd}/{a.results}/*.pdb")

    # tFold 사슬 = H(+L)+A, native = order (heavy..light,antigen 순) -> 위치 대응
    model_str = "H" + ("L" if len(antibody) >= 2 else "") + "A"
    native_str = "".join(order)
    n_ok, best = 0, None
    for pdb in pdbs:
        base = os.path.basename(pdb)
        ms = re.search(r"seed(\d+)", base)
        seed = int(ms.group(1)) if ms else 0
        conf = parse_conf(pdb)
        inter, dockq_out = parse_dockq(pdb, native, mapping=f"{model_str}:{native_str}")
        abag, have = ab_ag_dockq(inter, antigen, antibody)
        if not have:
            print(f"[경고] {base}: Ab-Ag DockQ 파싱 실패 "
                  f"(interfaces={[tuple(k) for k in inter]}) -> 건너뜀")
            continue
        dq = min(have)
        ok = dq >= 0.23
        n_ok += 1
        best = dq if best is None else max(best, dq)
        detail = " / ".join(f"{antigen}-{ab} {abag[ab]:.3f}" if abag[ab] is not None
                            else f"{antigen}-{ab} NA" for ab in antibody)
        notes = (f"seed{seed}; DockQ(Ab-Ag) {detail}; pair ipTM NA(tFold global만); "
                 f"ptm {conf['ptm']}; lDDT-Ca {conf['plddt']}; antigen lDDT {conf['ag_plddt']}; "
                 f"MSA=항체free/항원only")
        row = {"phase": a.phase, "model": a.model, "target": tid, "seed": seed,
               "dockq_abag": round(dq, 3),
               "global_iptm": round(conf["iptm"], 3) if conf["iptm"] is not None else "",
               "pair_iptm_abag": "",              # tFold는 pair-ipTM 없음
               "success": ok, "notes": notes}
        lib_log.append_row(a.csv, row)
        print(f"  {base}: DockQ={dq:.3f}  ipTM={conf['iptm']}  success={ok}")

    print(f"=== {tid}: {n_ok} pose 기록 -> {a.csv} | best DockQ = "
          f"{round(best,3) if best is not None else 'NA'} "
          f"({'성공' if (best is not None and best >= 0.23) else '실패'})")


if __name__ == "__main__":
    main()
