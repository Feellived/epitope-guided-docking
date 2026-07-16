#!/usr/bin/env python3
"""항원 표면을 k-means로 K개 패치 분할(전체 커버). SASA로 표면잔기 선별 → CA k-means → 패치 잔기리스트.
   → haddock/<t>/patches/patch_N.txt (chain B 원번호). 사용(boltz): python scripts/partition_surface.py [--k 6]"""
import argparse, glob, os
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import build_msafree_summary as B

def kmeans(X, k, seed=0, iters=100):
    rng = np.random.RandomState(seed); C = X[rng.choice(len(X), k, replace=False)].copy()
    for _ in range(iters):
        lab = ((X[:, None, :] - C[None, :, :]) ** 2).sum(-1).argmin(1)
        newC = np.array([X[lab == j].mean(0) if (lab == j).any() else C[j] for j in range(k)])
        if np.allclose(newC, C): break
        C = newC
    return lab

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--haddock", default="haddock"); ap.add_argument("--k", type=int, default=6)
    ap.add_argument("--sasa", type=float, default=15.0); ap.add_argument("--targets", default=" ".join(B.TARGETS))
    a = ap.parse_args()
    for t in a.targets.split():
        t = t.lower()
        cand = sorted(glob.glob(os.path.join(a.haddock, t, "antigen*.pdb")))
        if not cand: print(f"{t}: antigen*.pdb 없음"); continue
        m = PDBParser(QUIET=True).get_structure("x", cand[0])[0]
        ShrakeRupley().compute(m, level="R")
        res = []; ca = []
        for ch in m:
            for r in ch:
                if "CA" in r and getattr(r, "sasa", None) and r.sasa >= a.sasa:
                    res.append(r.id[1]); ca.append(r["CA"].get_coord())
        if len(ca) < a.k: print(f"{t}: 표면잔기 {len(ca)} < k → skip"); continue
        lab = kmeans(np.array(ca), a.k)
        od = os.path.join(a.haddock, t, "patches"); os.makedirs(od, exist_ok=True)
        sizes = []
        for j in range(a.k):
            rr = sorted(res[i] for i in range(len(res)) if lab[i] == j)
            open(os.path.join(od, f"patch_{j}.txt"), "w").write(" ".join(map(str, rr)) + "\n"); sizes.append(len(rr))
        print(f"{t}: 표면 {len(res)}잔기 → {a.k}패치 크기{sizes}  ({os.path.basename(cand[0])})")

if __name__ == "__main__": main()
