#!/usr/bin/env python3
"""HADDOCK 포즈 풀 오라클: 특정 stage 전체 포즈 개별 스캔 → 진짜자리 포즈 있나 + (emref면)점수순위.
   --stage 4_emref(기본,정제풀) / 1_rigidbody(넓은샘플, 탐색갭 확인).
   사용(boltz, repo루트): python scripts/analyze_haddock_pool.py [--stage 1_rigidbody] [--targets "9ZDU 8SIT"]"""
import argparse, glob, gzip, os, csv
import numpy as np
from Bio.PDB import PDBParser
import lib_structure as B

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
        except Exception: pass
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="4_emref"); ap.add_argument("--targets", default=" ".join(B.TARGETS))
    a = ap.parse_args()
    print(f"[stage={a.stage}]  target class  pool  oracle_rec #rec>=.3  best포즈 rec/RBM%/점수(순위)")
    print("-"*92)
    for t in a.targets.split():
        t = t.upper(); ni = B.native_info(t)
        if ni is None: print(f"{t:7} native 없음"); continue
        rundir = os.path.join("haddock", t.lower(), "run")
        poses = sorted(glob.glob(os.path.join(rundir, a.stage, "*.pdb.gz")))
        score = {}; srank = {}
        if a.stage.endswith("emref"):
            tsv = find_score_tsv(rundir)
            if tsv:
                rows = sorted((float(r["score"]), r["model_name"]) for r in csv.DictReader(open(tsv), delimiter="\t"))
                for i, (sc, mn) in enumerate(rows): score[mn] = sc; srank[mn] = i + 1
        res = []
        for p in poses:
            rr = pose_rec(p, ni)
            if rr is not None: res.append((rr[0], rr[1], os.path.basename(p)[:-3]))
        if not res: print(f"{t:7} 포즈 없음(stage {a.stage} 확인)"); continue
        res.sort(reverse=True); top = res[0]; nhi = sum(1 for r in res if r[0] >= 0.3)
        mn = top[2]; sc = score.get(mn); rk = srank.get(mn)
        scstr = f"{sc:.0f}(#{rk}/{len(res)})" if sc is not None else "-"
        print(f"{t:7}{B.EPICLASS.get(t,''):9}{len(res):<6}{top[0]:<11.2f}{nhi:<9}{mn}: {top[0]:.2f}/{top[1]*100:.0f}%/{scstr}")
    print("\n주: 1_rigidbody서 진짜자리 있음(oracle 높음)+emref서 낮았으면=점수컷에 잘림(seletop↑로 부활).")
    print("    1_rigidbody서도 oracle≈0=진짜 탐색 갭 → 패치스캔 필요.")

if __name__ == "__main__": main()
