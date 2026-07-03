#!/usr/bin/env python3
"""[Layer 2 · Boltz-2] chains.json -> boltz_<id>.yaml.
사용: python make_boltz.py --run-dir runs [--only 10ZO,9K6J]
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
        out = os.path.join(rd, f"boltz_{d['pdb_id'].lower()}.yaml")
        with open(out, "w") as f:
            f.write("version: 1\nsequences:\n")
            for c in d["chains"]:
                f.write(f'  - protein: {{ id: {c["id"]}, sequence: "{c["seq"]}" }}\n')
        print(f"[boltz] {d['pdb_id']} ({d['order']}) -> {out}")


if __name__ == "__main__":
    main()
