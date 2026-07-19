#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""재랭커 타당성 파일럿 (Exp3, RBD-10 pool) — make-or-break go/no-go.
질문: 이미 디스크에 있는 pool(co-folder blind + 감소-MSA + HADDOCK raw)에서, native 없이 좌표만으로 뽑은
값싼 계면 피처로 blind top-1 pose를 고르면 ipTM/단일모델보다 낫고 oracle에 근접하는가?
GPU 0. Biopython+numpy+scipy. PLM/PAE 없음(1차 기하·물리 피처만).

실행(cwd = repo/consensus_docking, run_phase0.sh와 동일):
  # 1) pool 피처+DockQ 라벨 캐시 생성(느림, resumable) — 먼저 --fast로 빠르게 훑고, 이후 전체
  python scripts/rerank_pilot.py build --fast --limit-per-source 20
  python scripts/rerank_pilot.py build                    # BSA 포함 전체
  # 2) 분석(즉시)
  python scripts/rerank_pilot.py analyze
주의: DockQ는 '평가 전용 라벨'이며 피처/선택에는 절대 안 씀. n=10 사례기반(통계 아님), 무학습 heuristic.
"""
import argparse, glob, gzip, json, os, re, csv, tempfile, sys, math
import numpy as np
from scipy.spatial import cKDTree
import build_msafree_summary as B
import haddock_dockq as HD           # load_gz, write_merged, dockq, native_merged (runs_rbd 기반)

CACHE = "results/rerank_pool.csv"
OUT   = "results/rerank_pilot_results.csv"

# ---- pool 소스 (research 스펙 §1; glob은 이 파일이 source of truth) ----
POOL = [
 ("blind","boltz",    "runs_rbd/{t}",        "*boltz/**/*_model_*.cif",           "cif"),
 ("blind","chai",     "runs_rbd/{t}",        "out_chai/pred.model_idx_*.cif",     "cif"),
 ("blind","protenix", "runs_rbd/{t}",        "out_protenix/**/*_sample_*.cif",    "cif"),
 ("blind","af2m",     "runs_rbd/{t}",        "*colabfold/*_unrelaxed_rank_*.pdb", "pdb"),
 ("blind","tfold",    "runs_rbd/{t}",        "out_tfold/*seed*.pdb",              "pdb"),
 ("msad143","boltz",  "runs_msad_143/{t}",   "*boltz/**/*_model_*.cif",           "cif"),
 ("msad8","boltz",    "runs_msad_8/{t}",     "*boltz/**/*_model_*.cif",           "cif"),
 ("msad1","boltz",    "runs_msad_1/{t}",     "*boltz/**/*_model_*.cif",           "cif"),
 ("msad143","protenix","runs_msad_143/{t}",  "out_protenix/**/*sample*.cif",      "cif"),
 ("haddock_emref","haddock","haddock/{t}",   "run/4_emref/*.pdb.gz",              "pdbgz"),
 ("haddock_rigid","haddock","haddock/{t}",   "run/1_rigidbody/*.pdb.gz",          "pdbgz"),
]

VDW = {"C":1.70,"N":1.55,"O":1.52,"S":1.80,"P":1.80,"F":1.47,"CL":1.75,"SE":1.90}
HYD = set("AVLIMFWY"); POS = set("KRH"); NEG = set("DE")
from Bio.Data.IUPACData import protein_letters_3to1
T3 = {k.upper():v for k,v in protein_letters_3to1.items()}
CUT = 5.0   # heavy-atom 접촉 컷오프(프로젝트 KNN12 정합)

def seed_of(path):
    for pat in (r"_model_(\d+)", r"model_idx_(\d+)", r"_sample_(\d+)", r"rank_(\d+)", r"seed[_]?(\d+)"):
        m = re.search(pat, os.path.basename(path))
        if m: return int(m.group(1))
    return -1

def enumerate_pool(t):
    """(family, model, seed, path, kind) 리스트. 존재하는 파일만."""
    out = []
    for fam, model, rd_t, gl, kind in POOL:
        rd = rd_t.format(t=t.lower())
        if not os.path.isdir(rd): continue
        for p in sorted(glob.glob(os.path.join(rd, gl), recursive=True)):
            out.append((fam, model, seed_of(p), p, kind))
    return out

def load_pose(path, kind):
    return HD.load_gz(path) if kind == "pdbgz" else B.load(path)

def residues(chain, model):
    """chain의 (CA, CB, resname1) per residue + heavy atom 좌표/원소."""
    ca={}; cb={}; rname={}; atoms=[]; ael=[]; ares=[]
    for r in model[chain]:
        if r.id[0] != ' ': continue                       # hetero 제외
        rid = r.id[1]
        aa = T3.get(r.get_resname().upper())
        if aa is None: continue
        rname[rid]=aa
        if "CA" in r: ca[rid]=r["CA"].coord
        cb[rid]= (r["CB"].coord if "CB" in r else (r["CA"].coord if "CA" in r else None))
        for a in r:
            el = (a.element or a.get_name()[0]).strip().upper()
            if el == "H" or a.get_name().strip().startswith("H"): continue
            atoms.append(a.coord); ael.append(el); ares.append(rid)
    return ca, cb, rname, np.array(atoms) if atoms else np.zeros((0,3)), ael, ares

def side(model, ids):
    """여러 사슬 합쳐 per-residue(사슬,잔기)키 dict + heavy atom 배열."""
    CA={}; CB={}; RN={}; X=[]; EL=[]; RES=[]
    for cid in ids:
        if cid not in [c.id for c in model]: continue
        ca,cb,rn,at,el,ares = residues(cid, model)
        for k in ca: CA[(cid,k)]=ca[k]
        for k in cb:
            if cb[k] is not None: CB[(cid,k)]=cb[k]
        for k in rn: RN[(cid,k)]=rn[k]
        for a,e,rr in zip(at,el,ares): X.append(a); EL.append(e); RES.append((cid,rr))
    return CA, CB, RN, (np.array(X) if X else np.zeros((0,3))), EL, RES

def features(model, ag_id, ab_ids, fast):
    agCA,agCB,agRN,agX,agEL,agRES = side(model,[ag_id])
    abCA,abCB,abRN,abX,abEL,abRES = side(model, ab_ids)
    f = {}
    if len(agX)==0 or len(abX)==0: return None
    tree = cKDTree(abX)
    # 접촉(5Å)
    nb = tree.query_ball_point(agX, CUT)
    n_atom=0; cres=set(); ag_if=set(); ab_if=set(); pair_res=set()
    for i,pl in enumerate(nb):
        for j in pl:
            n_atom+=1; ag_if.add(agRES[i]); ab_if.add(abRES[j]); pair_res.add((agRES[i],abRES[j]))
    f["n_contact_atoms"]=n_atom
    f["n_contact_res"]=len(pair_res)
    f["n_iface_ag"]=len(ag_if); f["n_iface_ab"]=len(ab_if)
    f["contact_density"]= n_atom/max(1,len(pair_res))
    # clash (vdw sum - 0.4)
    nb2 = tree.query_ball_point(agX, 3.8)
    ncl=0; sev=0.0
    for i,pl in enumerate(nb2):
        vi=VDW.get(agEL[i],1.7)
        for j in pl:
            d=float(np.linalg.norm(agX[i]-abX[j])); thr=vi+VDW.get(abEL[j],1.7)-0.4
            if d<thr: ncl+=1; sev+=(thr-d)
    f["n_clash"]=ncl; f["clash_severity"]=sev
    # 조성/상보성 (계면 잔기·접촉쌍 resname)
    if_res = list(ag_if|ab_if)
    RNall={**agRN,**abRN}
    aas=[RNall.get(r,"X") for r in if_res]
    nres=max(1,len(aas))
    f["frac_hydrophobic"]=sum(a in HYD for a in aas)/nres
    f["frac_charged"]=sum(a in (POS|NEG) for a in aas)/nres
    hyd_pair=0; chg_pair=0; opp=0; same=0
    for (ar,br) in pair_res:
        a=agRN.get(ar,"X"); b=abRN.get(br,"X")
        if a in HYD and b in HYD: hyd_pair+=1
        ac = 1 if a in POS else (-1 if a in NEG else 0)
        bc = 1 if b in POS else (-1 if b in NEG else 0)
        if ac and bc:
            chg_pair+=1
            if ac*bc<0: opp+=1
            else: same+=1
    f["hydrophobic_complementarity"]= hyd_pair/max(1,len(pair_res))
    f["charge_complementarity"]= (opp-same)/chg_pair if chg_pair else 0.0
    # 기하: 계면 잔기 CB로 planarity, CA로 Rg/clustering/KNN
    CBall={**agCB,**abCB}; CAall={**agCA,**abCA}
    ifCB=np.array([CBall[r] for r in if_res if r in CBall])
    ifCA=np.array([CAall[r] for r in if_res if r in CAall])
    if len(ifCB)>=4:
        c=ifCB-ifCB.mean(0); _,_,vt=np.linalg.svd(c,full_matrices=False)
        f["planarity_rmsd"]=float(np.sqrt((c@vt[-1])**2).mean())
    else: f["planarity_rmsd"]=0.0
    if len(ifCA)>=2:
        f["iface_Rg"]=float(np.sqrt(((ifCA-ifCA.mean(0))**2).sum(1).mean()))
        # clustering: CB-CB<8Å 그래프
        pts=ifCB if len(ifCB)>=2 else ifCA
        D=np.linalg.norm(pts[:,None,:]-pts[None,:,:],axis=-1); A=(D<8.0)&(D>0)
        deg=A.sum(1); cc=[]
        for i in range(len(pts)):
            nb_i=np.where(A[i])[0]
            if len(nb_i)<2: cc.append(0.0); continue
            sub=A[np.ix_(nb_i,nb_i)]; cc.append(sub.sum()/(len(nb_i)*(len(nb_i)-1)))
        f["iface_clustering"]=float(np.mean(cc)) if cc else 0.0
        # KNN(12) edge len
        k=min(12,len(ifCA)-1)
        if k>=1:
            Dca=np.linalg.norm(ifCA[:,None,:]-ifCA[None,:,:],axis=-1)
            np.fill_diagonal(Dca,1e9); knn=np.sort(Dca,1)[:,:k]
            f["knn_edge_len_mean"]=float(knn.mean())
        else: f["knn_edge_len_mean"]=0.0
    else:
        f["iface_Rg"]=0.0; f["iface_clustering"]=0.0; f["knn_edge_len_mean"]=0.0
    # BSA (Tier-2, --fast면 skip)
    if not fast:
        try:
            from Bio.PDB.SASA import ShrakeRupley
            sr=ShrakeRupley(probe_radius=1.4, n_points=100)
            import copy
            # 복합체(항원+항체 사슬만) SASA
            def sasa_sum(chs):
                from Bio.PDB.Structure import Structure; from Bio.PDB.Model import Model as M; from Bio.PDB.Chain import Chain as C
                s=Structure("x"); m=M(0); s.add(m)
                for cid in chs:
                    if cid in [c.id for c in model]:
                        m.add(copy.deepcopy(model[cid]))
                sr.compute(s, level="A")
                return sum(a.sasa for a in s.get_atoms())
            s_cplx=sasa_sum([ag_id]+list(ab_ids)); s_ag=sasa_sum([ag_id]); s_ab=sasa_sum(list(ab_ids))
            f["bsa_total"]=max(0.0, s_ag+s_ab-s_cplx)
        except Exception:
            f["bsa_total"]=""
    else:
        f["bsa_total"]=""
    return f

FEATCOLS=["n_contact_atoms","n_contact_res","n_iface_ag","n_iface_ab","contact_density","n_clash",
 "clash_severity","frac_hydrophobic","frac_charged","hydrophobic_complementarity","charge_complementarity",
 "planarity_rmsd","iface_Rg","iface_clustering","knn_edge_len_mean","bsa_total"]

def dockq_label(model, ag_id, ab_ids, natm, td):
    mp=os.path.join(td,"pose.pdb"); HD.write_merged(model, ag_id, ab_ids, mp)
    return HD.dockq(mp, natm)

# =============================== BUILD ===============================
def cmd_build(a):
    os.makedirs("results", exist_ok=True)
    done=set()
    if os.path.exists(CACHE):
        for r in csv.DictReader(open(CACHE)):
            done.add((r["target"], r["path"]))
    hdr=["target","epitope_class","family","model","seed","path","dockq"]+FEATCOLS
    write_hdr = not os.path.exists(CACHE)
    fh=open(CACHE,"a",newline=""); w=csv.DictWriter(fh,fieldnames=hdr)
    if write_hdr: w.writeheader()
    targets=a.targets.split() if a.targets else B.TARGETS
    for t in targets:
        t=t.upper(); ni=B.native_info(t)
        if ni is None: print(f"{t}: native 없음 skip"); continue
        poses=enumerate_pool(t)
        # 소스별 상한
        if a.limit_per_source:
            from collections import defaultdict
            cnt=defaultdict(int); kept=[]
            for fam,model,sd,p,kind in poses:
                if cnt[(fam,model)]<a.limit_per_source: kept.append((fam,model,sd,p,kind)); cnt[(fam,model)]+=1
            poses=kept
        print(f"[{t}] {len(poses)} pose (캐시된 것 skip)")
        with tempfile.TemporaryDirectory() as td:
            natm=HD.native_merged(t, td)
            if natm is None: print(f"  {t}: native merge 실패 skip"); continue
            n_new=0
            for fam,model,sd,p,kind in poses:
                if (t,p) in done: continue
                try:
                    m=load_pose(p,kind)
                    aid,_,_,_,_=B.best(m, ni["ag_seq"]); abids=[c.id for c in m if c.id!=aid]
                    if not abids: continue
                    feat=features(m, aid, abids, a.fast)
                    if feat is None: continue
                    q=dockq_label(m, aid, abids, natm, td)
                    row={"target":t,"epitope_class":B.EPICLASS.get(t,""),"family":fam,"model":model,
                         "seed":sd,"path":p,"dockq":(f"{q:.4f}" if q is not None else "")}
                    row.update({k:(feat.get(k,"")) for k in FEATCOLS})
                    w.writerow(row); fh.flush(); n_new+=1; done.add((t,p))
                except Exception as e:
                    continue
            print(f"  {t}: 신규 {n_new}행")
    fh.close(); print(f"→ {CACHE}")

# =============================== ANALYZE ===============================
def zscore(v):
    v=np.array(v,float); m=v[~np.isnan(v)]
    if len(m)<2 or m.std()==0: return np.zeros(len(v))
    return np.nan_to_num((v-m.mean())/m.std())

def load_summary_iptm():
    """results/summary.csv → {(target,model,seed): global_iptm}. 없거나 커버 안되면 {}."""
    p="results/summary.csv"
    if not os.path.exists(p): return {}
    d={}
    for r in csv.DictReader(open(p)):
        try:
            key=(r.get("target","").upper(), r.get("model","").lower(), int(float(r.get("seed",-1))))
            g=r.get("global_iptm","")
            if g not in ("","nan","NaN"): d[key]=float(g)
        except Exception: continue
    return d

MODEL_ALIAS={"boltz":"boltz-2","chai":"chai-1","protenix":"protenix-base","af2m":"af2-multimer","tfold":"tfold-ag"}

def cmd_analyze(a):
    if not os.path.exists(CACHE): sys.exit(f"{CACHE} 없음 — 먼저 build 실행")
    rows=[r for r in csv.DictReader(open(CACHE)) if r["dockq"] not in ("","nan")]
    by={}
    for r in rows: by.setdefault(r["target"],[]).append(r)
    iptm=load_summary_iptm()
    # 무학습 스코어 정의 (물리 prior; 값싼 조합). z는 타깃 내 정규화.
    def combo_score(rs):
        def col(k): return [float(x[k]) if x.get(k,"")!="" else np.nan for x in rs]
        z=lambda k: zscore(col(k))
        s = z("n_contact_res") + z("hydrophobic_complementarity") + z("charge_complementarity") \
            + z("iface_clustering") - 2.0*z("clash_severity") - 0.5*z("knn_edge_len_mean")
        if any(x.get("bsa_total","")!="" for x in rs): s = s + z("bsa_total")
        return s
    SINGLE={"contacts":("n_contact_res",+1),"neg_clash":("clash_severity",-1),
            "bsa":("bsa_total",+1),"clustering":("iface_clustering",+1),"hyd_compl":("hydrophobic_complementarity",+1)}
    tiers=[0.23,0.49,0.80]
    methods=["reranker_combo"]+list(SINGLE)+["iptm_pick","protenix_ceiling","random_mean","ORACLE"]
    agg={mth:{"top1":[], "win_iptm":[]} for mth in methods}
    perrow=[]; iptm_cov=0
    print(f"\n{'target':7}{'class':10}{'#pose':6}{'oracle':7}{'combo':7}{'iptm':7}{'protenix':9}{'rand':6}")
    print("-"*60)
    for t,rs in sorted(by.items()):
        dq=np.array([float(x["dockq"]) for x in rs])
        orc=float(dq.max())
        # combo
        cs=combo_score(rs); combo_dq=float(dq[int(np.argmax(cs))])
        # single-feature
        sing={}
        for name,(col,sign) in SINGLE.items():
            vals=np.array([float(x[col]) if x.get(col,"")!="" else np.nan for x in rs])
            if np.all(np.isnan(vals)): sing[name]=np.nan; continue
            idx=int(np.nanargmax(sign*vals)); sing[name]=float(dq[idx])
        # iptm-pick (co-folder만, summary 커버 시)
        ip=[]
        for x in rs:
            k=(t, MODEL_ALIAS.get(x["model"],x["model"]), int(float(x["seed"])) if x["seed"] not in("","-1") else -1)
            ip.append(iptm.get(k, np.nan))
        ip=np.array(ip); iptm_dq=float(dq[int(np.nanargmax(ip))]) if not np.all(np.isnan(ip)) else np.nan
        # best single model (protenix) = protenix pose 중 최고 DockQ (모델 천장)
        prot=[float(x["dockq"]) for x in rs if x["model"]=="protenix"]
        prot_dq=max(prot) if prot else np.nan
        rand=float(dq.mean())
        vals={"reranker_combo":combo_dq,"iptm_pick":iptm_dq,"protenix_ceiling":prot_dq,
              "random_mean":rand,"ORACLE":orc}; vals.update(sing)
        if not math.isnan(iptm_dq): iptm_cov+=1
        for mth in methods:
            v=vals.get(mth,np.nan); agg[mth]["top1"].append(v)
            if mth!="ORACLE" and not math.isnan(iptm_dq) and not (isinstance(v,float) and math.isnan(v)):
                agg[mth]["win_iptm"].append(1 if v>iptm_dq+1e-6 else 0)
        print(f"{t:7}{B.EPICLASS.get(t,''):10}{len(rs):<6}{orc:<7.2f}{combo_dq:<7.2f}"
              f"{(f'{iptm_dq:.2f}' if not math.isnan(iptm_dq) else '-'):<7}{(f'{prot_dq:.2f}' if not math.isnan(prot_dq) else '-'):<9}{rand:<6.2f}")
        perrow.append({"target":t,"n_pose":len(rs),"oracle":f"{orc:.3f}","combo":f"{combo_dq:.3f}",
                       "iptm":("" if math.isnan(iptm_dq) else f"{iptm_dq:.3f}"),
                       "protenix":("" if math.isnan(prot_dq) else f"{prot_dq:.3f}")})
    # 집계
    if iptm_cov==0:
        print("\n⚠️ ipTM baseline 미커버 — results/summary.csv에 이 RBD-10 타깃(seed 매칭)이 없음.")
        print("   → 판정은 ORACLE·단일모델·random 대비로 하고, ipTM 비교는 별도(summary.csv 세트/키 확인 필요).")
    print(f"\n=== 집계 (n={len(by)} 타깃) ===")
    orc_arr=np.array(agg["ORACLE"]["top1"],float)
    n_orc23=int(np.sum(orc_arr>=0.23))            # oracle이 acceptable 가능한 타깃 수(=선택 상한)
    print(f"{'method':20}{'meanTop1':10}{'@.23':6}{'@.49':6}{'@.80':6}{'Regret':8}{'HitOrc@.23':12}{'WinVsipTM'}")
    for mth in methods:
        arr=np.array(agg[mth]["top1"],float); valid=~np.isnan(arr)
        if valid.sum()==0: continue
        t1=arr[valid]; succ={tt:int(np.sum(t1>=tt)) for tt in tiers}
        if mth=="ORACLE": rg="0.000"
        else:
            m=valid & ~np.isnan(orc_arr)
            rg=f"{float(np.mean(orc_arr[m]-arr[m])):.3f}" if m.sum() else "-"     # NaN 타깃 제외한 평균 regret
        hit=int(np.sum((arr>=0.23)&(orc_arr>=0.23)))                              # Hit-Oracle@.23 분자
        win=agg[mth]["win_iptm"]; winr=f"{sum(win)}/{len(win)}" if win else "-"
        print(f"{mth:20}{t1.mean():<10.3f}{succ[0.23]:<6}{succ[0.49]:<6}{succ[0.80]:<6}{rg:<8}{f'{hit}/{n_orc23}':<12}{winr}")
    with open(OUT,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["target","n_pose","oracle","combo","iptm","protenix"]); w.writeheader(); w.writerows(perrow)
    print(f"\n→ {OUT}")
    print("판정: combo(무학습 물리 재랭커) succ@.49·HitOrc·Regret가 random/단일모델을 이기고 ORACLE에 근접하면 GO(재랭커 값어치).")
    print("      ipTM 커버되면 Win-vs-ipTM이 '망가진 선택기를 이긴다'의 직접 증거. 단일피처(contacts/neg_clash/bsa)가 이미 이기면 그 신호가 핵심.")
    print("      못 이기면 기하피처 부족 → PLM 임베딩·min-edge PAE로 확장 필요(2차).")
    print("정직: DockQ=평가전용 라벨(피처·선택 미사용). n=10 사례(통계 아님). combo 가중치=물리 prior(학습 아님). protenix_ceiling=Protenix pose 중 최댓값(천장).")

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    sub=ap.add_subparsers(dest="cmd", required=True)
    b=sub.add_parser("build"); b.add_argument("--targets",default=""); b.add_argument("--fast",action="store_true")
    b.add_argument("--limit-per-source",type=int,default=0)
    an=sub.add_parser("analyze")
    a=ap.parse_args()
    if a.cmd=="build": cmd_build(a)
    else: cmd_analyze(a)
