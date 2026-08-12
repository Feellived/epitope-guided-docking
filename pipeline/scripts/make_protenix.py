#!/usr/bin/env python3
"""[Layer 2 · Protenix-v2] chains.json -> protenix_<id>.json.
사용: python make_protenix.py --run-dir runs [--only 10ZO,9K6J]
⚠️ 정확한 입력 스키마는 첫 실행 시 `protenix`문서로 확인 후 조정.
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
        job = {"name": d["pdb_id"].lower(),
               "sequences": [{"proteinChain": {"sequence": c["seq"], "count": 1}} for c in d["chains"]]}
        out = os.path.join(rd, f"protenix_{d['pdb_id'].lower()}.json")
        with open(out, "w") as f:
            json.dump([job], f, indent=2)
        print(f"[protenix] {d['pdb_id']} ({d['order']}) -> {out}")


if __name__ == "__main__":
    main()
