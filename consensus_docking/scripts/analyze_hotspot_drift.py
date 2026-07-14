#!/usr/bin/env python3
"""[Study B/E2 분석] 각 모델이 항원의 '과대표집 hotspot'으로 쏠렸나 vs 진짜 자리 찾았나.

기존 예측물(runs_*/<id>) + 항원별 over-rep 지도(results/epitope_overrep_<id>.csv) 사용 (재예측 X).
모든 좌표를 chains.json 항원서열 index로 통일(정렬) → hotspot·pred·true 직접 비교.

각 (타깃,포맷,모델): majority-vote 예측 에피토프 →
  true_recall     = pred∩true / true      (진짜 자리 회복)
  pred_hot        = pred∩hot / pred        (예측이 hotspot에 몰림)
  true_hot        = true∩hot / true        (진짜 자리의 hotspot 비율)
  DRIFT           = pred_hot - true_hot     (+면 진짜보다 hotspot으로 더 쏠림 = 편향)
가설: MSA모델 DRIFT>0(hotspot 쏠림)·true_recall↓, tFold DRIFT≈0·true_recall↑.

사용: conda activate boltz
  python scripts/analyze_hotspot_drift.py --run-dir runs_diverse --targets 9K6J,9Y0A,9W43,...
  python scripts/analyze_hotspot_drift.py --run-dir runs_rbd        # 새 세트
⚠️ 먼저 9K6J로 검증(MSA 쏠림·tFold 회피 나와야 = 정상). hotspot = over-rep map 상위 --topk 잔기.
"""
import argparse, csv, glob, json, os
from collections import Counter
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_3to1

POSE_GLOB = {"Boltz-2":"*boltz/**/*_model_*.cif","Chai-1":"out_chai/pred.model_idx_*.cif",
    "Protenix-base":"out_protenix/**/*_sample_*.cif","AF2-Multimer":"*colabfold/*_unrelaxed_rank_*.pdb",
    "tFold-Ag":"out_tfold/*seed*.pdb"}
MODELS=list(POSE_GLOB); EPI=8.0
T3TO1={k.upper():v for k,v in protein_letters_3to1.items()}
_al=PairwiseAligner();_al.mode="global";_al.match_score=1;_al.mismatch_score=-1;_al.open_gap_score=-3;_al.extend_gap_score=-0.5

def load(p): return (MMCIFParser(QUIET=True) if p.endswith(".cif") else PDBParser(QUIET=True)).get_structure("x",p)[0]
def seq_ca(ch):
    s,c=[],[]
    for r in ch:
        if "CA" not in r: continue
        aa=T3TO1.get(r.get_resname().upper())
        if aa: s.append(aa); c.append(r["CA"].get_coord())
    return "".join(s),(np.array(c) if c else np.zeros((0,3)))
def best(model,tgt,exclude=()):
    b=(None,-1e9,"",np.zeros((0,3)))
    for ch in model:
        if ch.id in exclude: continue
        s,ca=seq_ca(ch)
        if len(s)<5: continue
        sc=_al.score(s,tgt)/max(len(s),len(tgt),1)
        if sc>b[1]: b=(ch.id,sc,s,ca)
    return b
def amap(sf,st):
    if not sf or not st: return {}
    aln=_al.align(sf,st)[0]; m={}
    for (a0,a1),(b0,b1) in zip(aln.aligned[0],aln.aligned[1]):
        for k in range(a1-a0): m[a0+k]=b0+k
    return m
def contacts(ag,ab):
    if len(ag)==0 or len(ab)==0: return set()
    d=np.sqrt(((ag[:,None,:]-ab[None,:,:])**2).sum(-1))
    return set(np.where(d.min(1)<=EPI)[0])
def find_native(rd,T):
    tags=("boltz","chai","protenix","colabfold","tfold")
    for f in os.listdir(rd):
        fl=f.lower()
        if fl.endswith((".pdb",".cif")) and T.lower() in fl and not any(m in fl for m in tags): return os.path.join(rd,f)
    return None

def load_hotspot(tag,topk,mincount):
    """over-rep map csv 상위 잔기 → chains.json-seq index(=ref_pos-1) 집합."""
    p=f"results/epitope_overrep_{tag}.csv"
    if not os.path.exists(p): return None,0
    rows=list(csv.DictReader(open(p)))
    # 컬럼: ref_pos, ... , (n_complexes | clustered_n | freq)
    key="clustered_n" if rows and "clustered_n" in rows[0] else ("n_complexes" if rows and "n_complexes" in rows[0] else None)
    if not key: return None,0
    vals=[(int(r["ref_pos"])-1,int(float(r[key]))) for r in rows if r.get(key,"").strip() not in ("","0")]
    total=max((v for _,v in vals),default=0)
    vals.sort(key=lambda x:-x[1])
    hot={i for i,v in vals[:topk] if v>=mincount}
    return hot,total

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--run-dir",default="runs_diverse")
    ap.add_argument("--targets",default="9K6J,9XF9,9Y0A,9VXL,9LY5,9MQR,9KKJ,9W43,9YC6,9J4C")
    ap.add_argument("--formats",default="Fab"); ap.add_argument("--topk",type=int,default=15)
    ap.add_argument("--minhot",type=int,default=3,help="hotspot 잔기 최소 결합복합체수(과대표집 근거)")
    a=ap.parse_args()
    RUN=a.run_dir
    print(f"{'target':6} {'ovrep':>6} {'model':13} {'pred_hot':>8} {'true_hot':>8} {'DRIFT':>6} {'true_recall':>11}")
    print("-"*70)
    for T in [t.strip().upper() for t in a.targets.split(",")]:
        hot,ovrep=load_hotspot(T,a.topk,a.minhot)
        if hot is None:
            print(f"# {T}: over-rep map 없음(results/epitope_overrep_{T}.csv) → 스킵"); continue
        for fmt in [f.strip() for f in a.formats.split(",")]:
            low=T.lower()+("_fv" if fmt=="Fv" else "")
            rd=os.path.join(RUN,low); cj=os.path.join(rd,"chains.json")
            if not os.path.exists(cj): continue
            d=json.load(open(cj)); ag=str(d["antigen"]); ab=d["antibody"]; ab=[ab] if isinstance(ab,str) else [str(x) for x in ab]
            sm={c["id"]:c["seq"] for c in d.get("chains",[])}; ag_seq=sm.get(ag,""); ab_seqs=[sm.get(i,"") for i in ab]
            natf=find_native(rd,T)
            if not natf or not ag_seq: continue
            nat=load(natf); nids=[c.id for c in nat]
            nseq,nca=(seq_ca(nat[ag]) if ag in nids else best(nat,ag_seq)[2:])
            nab=np.vstack([seq_ca(nat[i])[1] for i in ab if i in nids]) if any(i in nids for i in ab) else np.zeros((0,3))
            # true epitope → chains.json-seq index
            tmap=amap(nseq,ag_seq); true={tmap[i] for i in contacts(nca,nab) if i in tmap}
            if not true: continue
            true_hot=len(true&hot)/len(true)
            for model in MODELS:
                votes=Counter(); npose=0
                for pose in sorted(glob.glob(os.path.join(rd,POSE_GLOB[model]),recursive=True)):
                    try:
                        pm=load(pose); npose+=1
                        pid,_,ps,pca=best(pm,ag_seq); used={pid}; pab=[]
                        for s in ab_seqs:
                            r=best(pm,s,exclude=used)
                            if r[0]: used.add(r[0]); pab.append(r[3])
                        pab=np.vstack(pab) if pab else np.zeros((0,3))
                        m=amap(ps,ag_seq)
                        for i in contacts(pca,pab):
                            if i in m: votes[m[i]]+=1
                    except Exception: continue
                if npose==0: continue
                thr=max(1,(npose+1)//2); pred={r for r,v in votes.items() if v>=thr}
                if not pred:
                    print(f"{T:6} {ovrep:>6} {model:13} {'(no pred)':>8}"); continue
                ph=len(pred&hot)/len(pred); tr=len(pred&true)/len(true); drift=ph-true_hot
                flag=" <-- 쏠림" if drift>0.15 and tr<0.4 else ("  회피?" if tr>=0.5 else "")
                print(f"{T:6} {ovrep:>6} {model:13} {ph:>8.2f} {true_hot:>8.2f} {drift:>+6.2f} {tr:>11.2f}{flag}")
        print()

if __name__=="__main__":
    main()
