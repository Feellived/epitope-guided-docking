#!/usr/bin/env python3
"""HADDOCK 물리 도킹 리포트(per-cluster): 타깃별 클러스터 수·크기·HADDOCK점수·에피토프(recall/RBM%)
   + co-folder(server) 대비 divergence(Q1). build_msafree_summary 로직 재사용.
   ⚠️ repo 루트(consensus_docking)에서 boltz env로 실행. 사용: python scripts/report_haddock.py"""
import argparse, glob, gzip, os, csv
from collections import Counter, defaultdict
import numpy as np
from Bio.PDB import PDBParser
import build_msafree_summary as B

def load_gz(path):
    h = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    return PDBParser(QUIET=True).get_structure("x", h)[0]

def find_score_tsv(rundir):
    cands = glob.glob(os.path.join(rundir, "*clustfcc*", "*.tsv")) + \
            glob.glob(os.path.join(rundir, "**", "*.tsv"), recursive=True)
    for f in cands:
        try:
            hdr = next(csv.reader(open(f), delimiter="\t"))
            if {"model_name", "score", "cluster_id"} <= set(hdr):
                return f
        except Exception:
            pass
    return None

def pose_epitope(pose_paths, ni):
    votes = Counter(); n = 0
    for pose in pose_paths:
        if not os.path.exists(pose): continue
        try:
            pm = load_gz(pose); n += 1
            pid, _, ps, pca, _ = B.best(pm, ni["ag_seq"]); used = {pid}; pab = []
            for s in ni["ab_seqs"]:
                r = B.best(pm, s, exclude=used)
                if r[0]: used.add(r[0]); pab.append(r[3])
            pab = np.vstack(pab) if pab else np.zeros((0, 3))
            mp = B.amap(ps, ni["ag_seq"])
            for i in B.cset(pca, pab):
                if i in mp and mp[i] < len(ni["auth"]): votes[ni["auth"][mp[i]]] += 1
        except Exception:
            continue
    if n == 0: return None
    thr = max(1, (n + 1) // 2); pred = {a for a, v in votes.items() if v >= thr}
    return pred, n

def scores(pred, ni):
    recall = len(pred & ni["true"]) / len(ni["true"]) if ni["true"] else 0.0
    rbm = len([x for x in pred if B.is_rbm(x)]) / len(pred) if pred else 0.0
    return rbm, recall

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--haddock", default="haddock"); ap.add_argument("--emref", default="4_emref")
    ap.add_argument("--csv", default="results/msafree_summary.csv")
    a = ap.parse_args()
    cof = defaultdict(dict)
    if os.path.exists(a.csv):
        for r in csv.DictReader(open(a.csv)):
            if r["experiment"] == "cofolder":
                cof[r["target"].upper()][r["model"]] = (r["rbm_frac"], r["epitope_recall"])
    for t in B.TARGETS:
        print(f"\n===== {t} ({B.EPICLASS.get(t,'')}) =====")
        ni = B.native_info(t)
        rundir = os.path.join(a.haddock, t.lower(), "run")
        if ni is None: print("  native 없음"); continue
        tsv = find_score_tsv(rundir)
        if not tsv: print("  clustfcc tsv 없음"); continue
        rows = list(csv.DictReader(open(tsv), delimiter="\t"))
        emdir = os.path.join(rundir, a.emref)
        clusters = defaultdict(list)
        for r in rows:
            cid = (r.get("cluster_id") or "-").strip()
            if cid and cid != "-":
                clusters[cid].append((float(r["score"]), r["model_name"]))
        truelist = sorted(ni["true"])
        print(f"  클러스터 {len(clusters)}개 | native 에피토프({len(truelist)}): {truelist[:10]}{'...' if len(truelist)>10 else ''}")
        # 전체 최고 단일 모델(rank1, 미클러스터 포함)
        allm = sorted((float(r["score"]), r["model_name"]) for r in rows)
        if allm:
            ep = pose_epitope([os.path.join(emdir, allm[0][1] + ".gz")], ni)
            if ep:
                rbm, rec = scores(ep[0], ni)
                print(f"  best_single {allm[0][1]:14} score={allm[0][0]:7.1f}  RBM {rbm*100:3.0f}%  recall {rec:.2f}")
        # 클러스터별
        clist = []
        for cid, mem in clusters.items():
            mem.sort()
            poses = [os.path.join(emdir, m + ".gz") for _, m in mem]
            ep = pose_epitope(poses, ni)
            if ep:
                rbm, rec = scores(ep[0], ni)
                clist.append((mem[0][0], cid, len(mem), rbm, rec, sorted(ep[0])))
        clist.sort()
        print(f"  {'clust':7}{'size':5}{'score':>9}{'RBM%':>7}{'recall':>8}   epitope(auth, 일부)")
        for sc, cid, size, rbm, rec, pred in clist:
            print(f"  c{cid:<6}{size:<5}{sc:>9.1f}{rbm*100:>6.0f}%{rec:>8.2f}   {pred[:8]}")
        if t.upper() in cof:
            cfstr = "  ".join(f"{m.split('-')[0]}:{v[0]}/{v[1]}" for m, v in cof[t.upper()].items())
            print(f"  [co-folder server RBM%/recall] {cfstr}")
    print("\n주: HADDOCK RBM% 낮고 recall 높은 클러스터가 상위점수면=물리가 진짜자리 찾음(co-folder와 divergence=가치).")
    print("    진짜자리가 어느 클러스터에도 없으면=물리도 못 만듦→패치-스캐닝 필요(Q2).")

if __name__ == "__main__":
    main()
