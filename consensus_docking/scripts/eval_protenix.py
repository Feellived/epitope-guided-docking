#!/usr/bin/env python3
"""[eval · Protenix] 예측 1건 -> DockQ + 신뢰도(summary JSON) -> summary.csv 기록.

chains.json에서 antigen/antibody/순서. Protenix 출력 구조:
  out_protenix/<id>/seed_*/predictions/<id>_sample_<seed>.cif
  out_protenix/<id>/seed_*/predictions/<id>_summary_confidence_sample_<seed>.json
사용:
  python eval_protenix.py --target 10ZO --run-dir runs/10zo --seed 0 \
      --phase Benchmark --model Protenix --csv results/summary.csv
필요: DockQ(CLI).  (boltz env에서 실행)
"""
import argparse, glob, json, os, re, subprocess, sys
import log_lib


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


def parse_dockq(model_cif, native, mapping=None):
    def run(mp):
        cmd = ["DockQ", model_cif, native]
        if mp:
            cmd += ["--mapping", mp]
        return subprocess.run(cmd, capture_output=True, text=True).stdout
    out = run(mapping)
    inter = _interfaces(out)
    if not inter and mapping:                # 다중카피 비-cognate -> auto-map 재시도
        out2 = run(None)
        i2 = _interfaces(out2)
        if i2:
            return i2, out2
    return inter, out


def parse_conf(json_path, order, antigen, antibody):
    d = json.load(open(json_path))
    giptm = d.get("iptm")
    if giptm is None:                        # 상위 iptm 없으면 per-chain iptm 평균으로 대체
        ci = d.get("chain_iptm")
        if ci:
            giptm = sum(ci) / len(ci)
    mat = d.get("chain_pair_iptm") or d.get("chain_pair_iptm_global")
    idx = {c: i for i, c in enumerate(order)}
    ai = idx.get(antigen)
    vals = []
    if mat is not None:
        for ab in antibody:
            bi = idx.get(ab)
            if ai is None or bi is None:
                continue
            for a, b in [(ai, bi), (bi, ai)]:
                try:
                    vals.append(float(mat[a][b]))
                except Exception:
                    pass
    return (float(giptm) if giptm is not None else None), (min(vals) if vals else None), vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--run-dir", default=".")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--phase", default="Benchmark")
    ap.add_argument("--model", default="Protenix")
    ap.add_argument("--csv", default="results/summary.csv")
    a = ap.parse_args()
    rd, tid = a.run_dir, a.target
    low = tid.lower()

    cj = os.path.join(rd, "chains.json")
    if not os.path.exists(cj):
        sys.exit(f"[에러] chains.json 없음: {cj} (extract_chains 먼저)")
    d = json.load(open(cj))
    antigen, antibody, order = d["antigen"], d["antibody"], d["order"]

    native = find_one([f"{rd}/{tid}.pdb", f"{rd}/{tid}.cif",
                       f"{rd}/{tid.upper()}.pdb", f"{rd}/{tid.upper()}.cif",
                       f"{rd}/{low}.pdb", f"{rd}/{low}.cif"])
    cif = find_one([f"{rd}/out_protenix/**/{low}_sample_{a.seed}.cif"])
    conf = find_one([f"{rd}/out_protenix/**/{low}_summary_confidence_sample_{a.seed}.json"])
    miss = [n for n, v in [("native", native), ("cif", cif), ("conf", conf)] if not v]
    if miss:
        sys.exit(f"[에러] 파일 못 찾음: {miss} (run-dir={rd})")

    _m = "".join(order)
    inter, dockq_out = parse_dockq(cif, native, mapping=f"{_m}:{_m}")
    giptm, pair, pair_vals = parse_conf(conf, order, antigen, antibody)

    abag = {ab: inter.get(frozenset([antigen, ab])) for ab in antibody}
    have = [v for v in abag.values() if v is not None]
    if not have and len(antibody) == 1 and len(inter) == 1:
        only = next(iter(inter.values()))
        abag = {antibody[0]: only}
        have = [only]
    if not have:
        sys.exit(f"[에러] Ab-Ag interface DockQ 파싱 실패. antigen={antigen} antibody={antibody}\n"
                 f"DockQ가 본 interface: {[tuple(k) for k in inter]}\n--- DockQ 원문 ---\n{dockq_out}")
    dockq_min = min(have)
    success = dockq_min >= 0.23

    detail = " / ".join(f"{antigen}-{ab} {abag[ab]:.3f}" if abag[ab] is not None else f"{antigen}-{ab} NA"
                        for ab in antibody)
    pr = f"{min(pair_vals):.2f}-{max(pair_vals):.2f}" if pair_vals else "NA"
    notes = f"DockQ(Ab-Ag) {detail}; pair ipTM {pr}; interfaces={len(inter)}; MSA=protenix"

    row = {"phase": a.phase, "model": a.model, "target": tid, "seed": a.seed,
           "dockq_abag": round(dockq_min, 3),
           "global_iptm": round(giptm, 3) if giptm is not None else "",
           "pair_iptm_abag": round(pair, 3) if pair is not None else "",
           "success": success, "notes": notes}
    saved = log_lib.append_row(a.csv, row)
    print("=== 기록됨 ->", a.csv, "===")
    for k in log_lib.FIELDS:
        print(f"  {k}: {saved[k]}")


if __name__ == "__main__":
    main()
