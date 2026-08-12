#!/usr/bin/env python3
"""[eval · Chai-1] 예측 1건 -> DockQ + 신뢰도(npz) -> summary.csv 기록.

chains.json(추출 기록)을 읽어 antigen/antibody/사슬순서를 얻음.
run-dir 안: chains.json, <ID>.pdb/.cif(정답),
            out_chai/pred.model_idx_<seed>.cif, out_chai/scores.model_idx_<seed>.npz
사용:
  python eval_chai.py --target 10ZO --run-dir runs/10zo --seed 0 \
      --phase Benchmark --model Chai-1 --csv results/summary.csv
필요: DockQ(CLI), numpy.  (eval을 chai env에서 돌리면 `pip install DockQ` 먼저)
"""
import argparse, glob, json, os, re, subprocess, sys
import numpy as np
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


def parse_conf(npz_path, order, antigen, antibody):
    d = np.load(npz_path)
    giptm = float(d["iptm"][0])
    mat = d["per_chain_pair_iptm"][0]        # (n,n), 인덱스 = 사슬 등장 순서(order)
    idx = {c: i for i, c in enumerate(order)}
    ai = idx.get(antigen)
    vals = []
    for ab in antibody:
        bi = idx.get(ab)
        if ai is None or bi is None:
            continue
        for a, b in [(ai, bi), (bi, ai)]:
            try:
                vals.append(float(mat[a][b]))
            except Exception:
                pass
    return giptm, (min(vals) if vals else None), vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--run-dir", default=".")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--phase", default="Benchmark")
    ap.add_argument("--model", default="Chai-1")
    ap.add_argument("--out-dir", default="out_chai", help="Chai 출력 폴더 (no-MSA는 out_chai_nomsa)")
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
    cif = find_one([f"{rd}/{a.out_dir}/pred.model_idx_{a.seed}.cif"])
    npz = find_one([f"{rd}/{a.out_dir}/scores.model_idx_{a.seed}.npz"])
    miss = [n for n, v in [("native", native), ("cif", cif), ("npz", npz)] if not v]
    if miss:
        sys.exit(f"[에러] 파일 못 찾음: {miss} (run-dir={rd})")

    _m = "".join(order)
    inter, dockq_out = parse_dockq(cif, native, mapping=f"{_m}:{_m}")
    giptm, pair, pair_vals = parse_conf(npz, order, antigen, antibody)

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
    notes = f"DockQ(Ab-Ag) {detail}; pair ipTM {pr}; interfaces={len(inter)}"

    row = {"phase": a.phase, "model": a.model, "target": tid, "seed": a.seed,
           "dockq_abag": round(dockq_min, 3),
           "global_iptm": round(giptm, 3) if giptm is not None else "",
           "pair_iptm_abag": round(pair, 3) if pair is not None else "",
           "success": success, "notes": notes}
    saved = lib_log.append_row(a.csv, row)
    print("=== 기록됨 ->", a.csv, "===")
    for k in lib_log.FIELDS:
        print(f"  {k}: {saved[k]}")


if __name__ == "__main__":
    main()
