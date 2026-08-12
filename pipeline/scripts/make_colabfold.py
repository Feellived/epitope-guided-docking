#!/usr/bin/env python3
"""[Layer 2 · AF2-Multimer/ColabFold] chains.json -> colabfold_<id>.fasta.
ColabFold multimer 입력 = 단일 레코드, 사슬을 ':'로 이어붙임 (chains.json order 순서).
검증: 9XF9 = 항원(A,240) : H(219) : L(218) = 677 (colabfold.log query length과 일치).
사용: python make_colabfold.py --run-dir runs [--only 10ZO,9K6J]
"""
import argparse, glob, json, os


def iter_runs(run_dir, only):
    for p in sorted(glob.glob(os.path.join(run_dir, "*", "chains.json"))):
        with open(p) as f:
            d = json.load(f)
        if only and d["pdb_id"].upper() not in only:
            continue
        yield os.path.dirname(p), d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", default="runs")
    ap.add_argument("--only")
    a = ap.parse_args()
    only = set(x.upper() for x in a.only.split(",")) if a.only else None
    for rd, d in iter_runs(a.run_dir, only):
        pid = d["pdb_id"].lower()
        joined = ":".join(c["seq"] for c in d["chains"])   # 복합체 = 사슬을 ':'로
        out = os.path.join(rd, f"colabfold_{pid}.fasta")
        with open(out, "w") as f:
            f.write(f">{pid}\n{joined}\n")
        # 항원만 FASTA (tFold 항원 MSA용; make_tfold의 antigen_<id>.fasta와 별개로 colabfold 입력)
        ag = next((c for c in d["chains"] if c["role"] == "antigen"), None)
        if ag:
            with open(os.path.join(rd, f"antigen_cf_{pid}.fasta"), "w") as f:
                f.write(f">A\n{ag['seq']}\n")
        n = sum(len(c["seq"]) for c in d["chains"])
        print(f"[colabfold] {d['pdb_id']} ({d['order']}) len={n} -> {out}")


if __name__ == "__main__":
    main()
