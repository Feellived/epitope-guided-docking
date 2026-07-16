#!/usr/bin/env python3
"""HADDOCK emref 풀 오라클: 타깃별 전체 emref 포즈를 개별 스캔(클러스터/최고점 무관) →
   진짜자리 포즈가 풀에 있나 + 있으면 recall·RBM% + 점수순위(HADDOCK score가 얼마나 묻었나).
   → 진짜자리 풀에 있고 점수순위 나쁨=재랭커 구제가능 / 아예 없음=패치스캔(강제생성) 필요.
   사용(boltz env, repo루트): python scripts/haddock_pool_oracle.py"""
import glob, gzip, os, csv
import numpy as np
from Bio.PDB import PDBParser
import build_msafree_summary as B

def load_gz(p):
    h = gzip.open(p, "rt") if p.endswith(".gz") else open(p)
    return PDBParser(QUIET=True).get_structure("x", h)[0]

def pose_rec(pose, ni):
    try:
        pm = load_gz(pose)
        pid, _, ps, pca, _ = B.best(pm, ni["ag_seq"]); used = {pid}; pab = []
        for s in ni["ab_seqs"]:
            r = B.best(pm, s, exclude=used)
            if r[0]: used.add(r[0]); pab.append(r[3])
        pab = np.vstack(pab) if pab else np.zeros((0, 3))
        mp = B.amap(ps, ni["ag_seq"]); pred = set()
        for i in B.cset(pca, pab):
            if i in mp and mp[i] < len(ni["auth"]): pred.add(ni["auth"][mp[i]])
        rec = len(pred & ni["true"]) / len(ni["true"]) if ni["true"] else 0
        rbm = len([x for x in pred if B.is_rbm(x)]) / len(pred) if pred else 0
        return rec, rbm
    except Exception:
        return None

def find_score_tsv(rundir):
    for f in glob.glob(os.path.join(rundir, "**", "*.tsv"), recursive=True):
        try:
            hdr = next(csv.reader(open(f), delimiter="\t"))
            if {"model_name", "score", "cluster_id"} <= set(hdr): return f
        except Exception:
            pass
    return None

print(f"{'target':7}{'class':9}{'pool':6}{'oracle_rec':11}{'#rec>=.3':9}  best-recall 포즈: rec/RBM%/HADDOCKscore(점수순위)")
print("-"*100)
for t in B.TARGETS:
    ni = B.native_info(t)
    if ni is None: print(f"{t:7} native 없음"); continue
    rundir = os.path.join("haddock", t.lower(), "run")
    poses = sorted(glob.glob(os.path.join(rundir, "4_emref", "*.pdb.gz")))
    tsv = find_score_tsv(rundir); score = {}; srank = {}
    if tsv:
        rows = sorted((float(r["score"]), r["model_name"]) for r in csv.DictReader(open(tsv), delimiter="\t"))
        for i, (sc, mn) in enumerate(rows): score[mn] = sc; srank[mn] = i + 1
    res = []
    for p in poses:
        rr = pose_rec(p, ni)
        if rr is None: continue
        res.append((rr[0], rr[1], os.path.basename(p)[:-3]))
    if not res: print(f"{t:7} 포즈 없음"); continue
    res.sort(reverse=True)
    top = res[0]; nhi = sum(1 for r in res if r[0] >= 0.3)
    mn = top[2]; sc = score.get(mn); rk = srank.get(mn)
    scstr = f"{sc:.0f}(#{rk}/{len(res)})" if sc is not None else "?"
    print(f"{t:7}{B.EPICLASS.get(t,''):9}{len(res):<6}{top[0]:<11.2f}{nhi:<9}  {mn}: {top[0]:.2f}/{top[1]*100:.0f}%/{scstr}")
print("\n주: oracle_rec 높은데 점수순위(#)가 나쁨 = 진짜자리 포즈가 풀에 있으나 HADDOCK점수가 묻음 → 재랭커 구제가능.")
print("    oracle_rec≈0 = 진짜자리 포즈 자체가 풀에 없음(RBM이 샘플링 독식) → 패치스캔(강제 생성) 필요.")
