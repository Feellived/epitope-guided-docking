#!/usr/bin/env python3
"""[Layer 2 · Chai-1] chains.json -> chai_<id>.fasta.
Chai는 등장 순서대로 사슬을 A,B,C.. 로 배정 (name은 라벨).
사용: python make_chai.py --run-dir runs [--only 10ZO,9K6J]
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
        out = os.path.join(rd, f"chai_{d['pdb_id'].lower()}.fasta")
        with open(out, "w") as f:
            for c in d["chains"]:
                f.write(f">protein|name={c['id']}\n{c['seq']}\n")
        print(f"[chai] {d['pdb_id']} ({d['order']}) -> {out}")


if __name__ == "__main__":
    main()
