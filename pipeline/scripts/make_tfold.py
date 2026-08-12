#!/usr/bin/env python3
"""[Layer 2 · tFold-Ag] chains.json -> tfold_<id>.fasta(>H,>L,>A) + antigen_<id>.fasta(>A).
tFold-Ag는 JSON 모드가 icf_path KeyError 버그라 FASTA 모드 사용.
매핑: 항체 heavy=H, light=L(있으면), 항원=A. antigen_<id>.fasta는 항원 a3m 생성용(원격 MMseqs2).
사용: python make_tfold.py --run-dir runs [--only 10ZO,9K6J]
"""
import argparse, glob, json, os


def iter_runs(run_dir, only):
    for p in sorted(glob.glob(os.path.join(run_dir, "*", "chains.json"))):
        d = json.load(open(p))
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
        by_id = {c["id"]: c["seq"] for c in d["chains"]}
        ab, ag = d["antibody"], d["antigen"]
        # 항체 heavy=H, light=L(Fab만), 항원=A  (native order와 위치 대응)
        recs = [("H", by_id[ab[0]])]
        if len(ab) >= 2:
            recs.append(("L", by_id[ab[1]]))
        recs.append(("A", by_id[ag]))
        fa = os.path.join(rd, f"tfold_{pid}.fasta")
        with open(fa, "w") as f:
            for cid, seq in recs:
                f.write(f">{cid}\n{seq}\n")
        agfa = os.path.join(rd, f"antigen_{pid}.fasta")   # 항원 단독(a3m 생성용)
        with open(agfa, "w") as f:
            f.write(f">A\n{by_id[ag]}\n")
        print(f"[tfold] {d['pdb_id']} ({[c[0] for c in recs]}) -> {fa} + {agfa}")


if __name__ == "__main__":
    main()
