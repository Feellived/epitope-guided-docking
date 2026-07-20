#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""재랭커 타당성 파일럿 (Exp3, RBD-10 pool) — make-or-break go/no-go.
질문: 이미 디스크에 있는 pool(co-folder blind + 감소-MSA + HADDOCK raw)에서, native 없이 좌표만으로 뽑은
값싼 계면 피처로 blind top-1 pose를 고르면 ipTM/단일모델보다 낫고 oracle에 근접하는가?
GPU 0. Biopython+numpy(+scipy 선택). 기하·물리 피처 + min-edge PAE(Boltz·AF2, PAE가 디스크에 있을 때). PLM 임베딩은 2차.

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
try:
    from scipy.spatial import cKDTree; _SCIPY=True
except Exception:
    _SCIPY=False                      # scipy 없으면 numpy 폴백(설치 불필요)
import build_msafree_summary as B
import haddock_dockq as HD           # load_gz, write_merged, dockq, native_merged (runs_rbd 기반)

def neigh(agX, abX, r):
    """각 ag atom에 대해 반경 r 내 ab atom index 리스트. scipy 있으면 cKDTree, 없으면 numpy."""
    if len(abX)==0: return [[] for _ in agX]
    if _SCIPY:
        return cKDTree(abX).query_ball_point(agX, r)
    out=[]
    for i in range(len(agX)):
        out.append(list(np.where(np.linalg.norm(abX-agX[i],axis=1)<r)[0]))
    return out

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

def load_pae(model_name, path):
    """모델별 PAE 행렬 로드(N_tok×N_tok). Boltz npz / ColabFold scores json. 실패 시 None."""
    d=os.path.dirname(path)
    try:
        if model_name=="af2m":
            m=re.search(r"rank_(\d+)", os.path.basename(path)); rk=m.group(1) if m else ""
            cands=sorted(glob.glob(os.path.join(d,"*scores*rank*"+rk+"*.json"))) or sorted(glob.glob(os.path.join(d,"*scores*.json")))
            for c in cands:
                j=json.load(open(c))
                for k in ("pae","predicted_aligned_error"):
                    if k in j:
                        a=np.asarray(j[k],float)
                        if a.ndim==2 and a.shape[0]==a.shape[1]: return a
        elif model_name=="boltz":
            m=re.search(r"model_(\d+)", os.path.basename(path)); sd=m.group(1) if m else ""
            pats=[f"pae*model_{sd}*.npz",f"*pae*model_{sd}*.npz",f"pae*model_{sd}*.npy","*pae*.npz"]
            for base in (d, os.path.dirname(d)):
                for pat in pats:
                    for c in sorted(glob.glob(os.path.join(base,pat))):
                        obj=np.load(c, allow_pickle=True)
                        arrs=[np.asarray(obj[k],float) for k in obj.files] if hasattr(obj,"files") else [np.asarray(obj,float)]
                        for a in arrs:
                            if a.ndim==2 and a.shape[0]==a.shape[1]: return a
    except Exception:
        return None
    return None

def load_iptm(model_name, path):
    """per-pose global ipTM (라벨 아님, 모델 자기신뢰도) — baseline·공통집합 정의용. 실패 시 None."""
    d=os.path.dirname(path); bn=os.path.basename(path)
    try:
        if model_name=="boltz":
            m=re.search(r"model_(\d+)",bn); sd=m.group(1) if m else ""
            for c in sorted(glob.glob(os.path.join(d,f"confidence*model_{sd}*.json")))+sorted(glob.glob(os.path.join(d,"confidence*.json"))):
                j=json.load(open(c))
                for k in ("iptm","complex_iptm"):
                    if k in j: return float(j[k])
        elif model_name=="af2m":
            m=re.search(r"rank_(\d+)",bn); rk=m.group(1) if m else ""
            for c in sorted(glob.glob(os.path.join(d,f"*scores*rank*{rk}*.json")))+sorted(glob.glob(os.path.join(d,"*scores*.json"))):
                j=json.load(open(c))
                if "iptm" in j: return float(j["iptm"])
        elif model_name=="chai":
            m=re.search(r"idx_(\d+)",bn); sd=m.group(1) if m else ""
            for c in sorted(glob.glob(os.path.join(d,f"scores*idx_{sd}*.npz")))+sorted(glob.glob(os.path.join(d,"scores*.npz"))):
                o=np.load(c,allow_pickle=True)
                if "iptm" in getattr(o,"files",[]): return float(np.asarray(o["iptm"]).mean())
        elif model_name=="protenix":
            m=re.search(r"sample_(\d+)",bn); sd=m.group(1) if m else ""
            for c in sorted(glob.glob(os.path.join(d,f"*summary_confidence*sample_{sd}*.json")))+sorted(glob.glob(os.path.join(d,"*summary_confidence*.json"))):
                j=json.load(open(c))
                if "iptm" in j: return float(j["iptm"])
        elif model_name=="tfold":
            with open(path) as fh:
                for line in fh:
                    if not (line.startswith("REMARK") or line.startswith("HEADER") or line.startswith("TITLE")): break
                    if "iptm" in line.lower():
                        mm=re.search(r"[-+]?\d*\.\d+",line)
                        if mm: return float(mm.group())
    except Exception:
        return None
    return None

def token_ranges(model, ag_id, ab_ids):
    """structure 사슬 순서대로 CA 잔기에 토큰 인덱스 부여 → (ag_tok, ab_tok, n_total)."""
    ag=[]; ab=[]; idx=0
    for ch in model:
        n=sum(1 for r in ch if r.id[0]==' ' and 'CA' in r and T3.get(r.get_resname().upper()))
        rng=list(range(idx,idx+n)); idx+=n
        if ch.id==ag_id: ag+=rng
        elif ch.id in ab_ids: ab+=rng
    return ag, ab, idx

def min_edge_pae(pae, ag_tok, ab_tok, topn=10):
    """계면 cross block(ag×ab, ab×ag)에서 가장 작은 topn PAE 평균(낮을수록 신뢰). = min-edge PAE 근사."""
    if pae is None or not ag_tok or not ab_tok: return None
    if max(max(ag_tok),max(ab_tok))>=pae.shape[0]: return None
    blk=np.concatenate([pae[np.ix_(ag_tok,ab_tok)].ravel(), pae[np.ix_(ab_tok,ag_tok)].ravel()])
    if blk.size==0: return None
    return float(np.sort(blk)[:min(topn,blk.size)].mean())

def features(model, ag_id, ab_ids, fast, model_name="", path="", with_pae=True):
    agCA,agCB,agRN,agX,agEL,agRES = side(model,[ag_id])
    abCA,abCB,abRN,abX,abEL,abRES = side(model, ab_ids)
    f = {}
    if len(agX)==0 or len(abX)==0: return None
    # 접촉(5Å)
    nb = neigh(agX, abX, CUT)
    n_atom=0; cres=set(); ag_if=set(); ab_if=set(); pair_res=set()
    for i,pl in enumerate(nb):
        for j in pl:
            n_atom+=1; ag_if.add(agRES[i]); ab_if.add(abRES[j]); pair_res.add((agRES[i],abRES[j]))
    f["n_contact_atoms"]=n_atom
    f["n_contact_res"]=len(pair_res)
    f["n_iface_ag"]=len(ag_if); f["n_iface_ab"]=len(ab_if)
    f["contact_density"]= n_atom/max(1,len(pair_res))
    # clash (vdw sum - 0.4)
    nb2 = neigh(agX, abX, 3.8)
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
    # min-edge PAE (Boltz·AF2만; PAE 디스크에 있을 때. 실패=빈값 → 기하 파일럿 무영향)
    f["min_edge_pae"]=""
    if with_pae and model_name in ("boltz","af2m"):
        pae=load_pae(model_name, path)
        if pae is not None:
            agt,abt,ntok=token_ranges(model, ag_id, ab_ids)
            if pae.shape[0]==ntok:                          # 토큰수=잔기수 일치할 때만(순서 가정 가드)
                v=min_edge_pae(pae, agt, abt)
                if v is not None: f["min_edge_pae"]=v
    return f

FEATCOLS=["n_contact_atoms","n_contact_res","n_iface_ag","n_iface_ab","contact_density","n_clash",
 "clash_severity","frac_hydrophobic","frac_charged","hydrophobic_complementarity","charge_complementarity",
 "planarity_rmsd","iface_Rg","iface_clustering","knn_edge_len_mean","bsa_total","min_edge_pae"]

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
    hdr=["target","epitope_class","family","model","seed","path","dockq","iptm"]+FEATCOLS
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
                    feat=features(m, aid, abids, a.fast, model, p, not a.no_pae)
                    if feat is None: continue
                    q=dockq_label(m, aid, abids, natm, td)
                    ip=load_iptm(model, p)
                    row={"target":t,"epitope_class":B.EPICLASS.get(t,""),"family":fam,"model":model,
                         "seed":sd,"path":p,"dockq":(f"{q:.4f}" if q is not None else ""),
                         "iptm":(f"{ip:.4f}" if ip is not None else "")}
                    row.update({k:(feat.get(k,"")) for k in FEATCOLS})
                    w.writerow(row); fh.flush(); n_new+=1; done.add((t,p))
                except Exception as e:
                    continue
            print(f"  {t}: 신규 {n_new}행")
    fh.close(); print(f"→ {CACHE}")

# =============================== ANALYZE ===============================
def _fv(x,k):
    try: return float(x[k]) if x.get(k,"")!="" else float('nan')
    except: return float('nan')
def _ci(diffs):
    """paired diff bootstrap 95% CI of mean(oracle-method) (deterministic seed)."""
    d=np.array([x for x in diffs if isinstance(x,float) and not math.isnan(x)],float)
    if len(d)<3: return "-"
    rng=np.random.RandomState(0); mu=[rng.choice(d,len(d),replace=True).mean() for _ in range(2000)]
    return f"[{np.percentile(mu,2.5):.2f},{np.percentile(mu,97.5):.2f}]"
# 단일피처 (표시명, 컬럼, 부호) — 부호=물리 prior 고정(튜닝 DOF 0). PAE 최우선.
SINGLE=[("min_pae","min_edge_pae",-1),("contacts","n_contact_res",+1),("bsa","bsa_total",+1),
        ("neg_clash","clash_severity",-1),("clustering","iface_clustering",+1),("hyd_compl","hydrophobic_complementarity",+1)]
def _pick(rs,col,sign):
    v=np.array([_fv(x,col) for x in rs])
    if np.all(np.isnan(v)): return float('nan')
    dq=np.array([float(x["dockq"]) for x in rs]); return float(dq[int(np.nanargmax(sign*v))])

def _report(title, rec, methods, ref, orc_key):
    print(f"\n{title}")
    orc=np.array(rec[orc_key],float); ov=orc[~np.isnan(orc)]
    n_o49=int(np.sum(ov>=0.49)); n_o23=int(np.sum(ov>=0.23))
    refarr=np.array(rec[ref],float) if ref in rec else None
    print(f"{'method':18}{'mTop1':8}{'@.49':6}{'@.23':6}{'Regret':8}{'Reg95CI':14}{'Win@.49':9}{'HitOrc@.49'}")
    for m in methods:
        arr=np.array(rec[m],float); valid=~np.isnan(arr)
        if valid.sum()==0: continue
        av=arr[valid]; s49=int(np.sum(av>=0.49)); s23=int(np.sum(av>=0.23))
        mo=valid & ~np.isnan(orc); diffs=list(orc[mo]-arr[mo])
        rg=f"{np.mean(diffs):.3f}" if mo.sum() else "-"; ci=_ci(diffs) if m!=orc_key else "-"
        hit=int(np.sum((arr[mo]>=0.49)&(orc[mo]>=0.49)))
        if refarr is not None and m not in (ref,orc_key):
            w=valid & ~np.isnan(refarr); winr=f"{int(np.sum(arr[w]>refarr[w]+1e-6))}/{int(w.sum())}" if w.sum() else "-"
        else: winr="-"
        print(f"{m:18}{np.nanmean(arr):<8.3f}{s49:<6}{s23:<6}{rg:<8}{ci:<14}{winr:<9}{hit}/{n_o49}")

def cmd_analyze(a):
    if not os.path.exists(CACHE): sys.exit(f"{CACHE} 없음 — 먼저 build 실행")
    rows=[r for r in csv.DictReader(open(CACHE)) if r["dockq"] not in ("","nan")]
    by={}
    for r in rows: by.setdefault(r["target"],[]).append(r)
    from collections import Counter
    print("=== pool 구성(타깃별 total pose · ipTM 있는 pose) ===")
    for t in sorted(by):
        rs=by[t]; nip=sum(1 for x in rs if not math.isnan(_fv(x,"iptm")))
        fam=Counter(x["family"] for x in rs)
        print(f"  {t:6} total {len(rs):4}  ipTM {nip:4} | "+" ".join(f"{k}:{v}" for k,v in sorted(fam.items())))

    # ===== 분석 A: 공통 후보집합(ipTM 존재 pose=co-folder) — PRIMARY =====
    methsA=[nm for nm,_,_ in SINGLE]+["iptm_pick","protenix_ownIptm","ORACLE"]
    recA={m:[] for m in methsA}
    for t in sorted(by):
        common=[x for x in by[t] if not math.isnan(_fv(x,"iptm"))]
        if not common:
            for m in methsA: recA[m].append(float('nan'))
            continue
        dqc=np.array([float(x["dockq"]) for x in common])
        recA["ORACLE"].append(float(dqc.max()))
        recA["iptm_pick"].append(_pick(common,"iptm",+1))
        for nm,col,sg in SINGLE: recA[nm].append(_pick(common,col,sg))
        prot=[x for x in common if x["model"]=="protenix"]
        recA["protenix_ownIptm"].append(_pick(prot,"iptm",+1) if prot else float('nan'))
    _report(f"### 분석 A — 공통집합(ipTM 존재 co-folder pose) head-to-head [PRIMARY · n={len(by)}] ###",
            recA, methsA, ref="iptm_pick", orc_key="ORACLE")
    best=None
    for nm,_,_ in SINGLE:
        arr=np.array(recA[nm],float); ref=np.array(recA["iptm_pick"],float)
        w=~np.isnan(arr)&~np.isnan(ref); wins=int(np.sum(arr[w]>ref[w]+1e-6))
        if best is None or wins>best[1]: best=(nm,wins,int(w.sum()))
    print(f"\n★ PRIMARY 판정 = 최고 단일피처 '{best[0]}' 의 Win-vs-ipTM@.49 = {best[1]}/{best[2]} (부호 고정=튜닝 DOF 0)")
    print("   → 크고 robust(Regret CI가 0 밑, LOTO 안정)하면 GO: 값싼 계면신호가 망가진 ipTM을 이김.")

    # ===== 분석 B: 확장 pool(전체) — SECONDARY(다른 질문: pool 풍부화 이득) =====
    def combo_unit(rs):   # 단위가중(±1), 물리 prior 부호만 = exploratory(판정 아님)
        def z(col):
            v=np.array([_fv(x,col) for x in rs]); m=v[~np.isnan(v)]
            return np.zeros(len(v)) if len(m)<2 or m.std()==0 else np.nan_to_num((v-m.mean())/m.std())
        s=z("n_contact_res")+z("hydrophobic_complementarity")+z("charge_complementarity")+z("iface_clustering")-z("clash_severity")
        if any(x.get("bsa_total","")!="" for x in rs): s=s+z("bsa_total")
        if any(x.get("min_edge_pae","")!="" for x in rs): s=s-z("min_edge_pae")
        return s
    methsB=["combo_unit","min_pae","strat_random","ORACLE"]
    recB={m:[] for m in methsB}; rngB=np.random.RandomState(0)
    for t in sorted(by):
        rs=by[t]; dq=np.array([float(x["dockq"]) for x in rs])
        recB["ORACLE"].append(float(dq.max()))
        recB["combo_unit"].append(float(dq[int(np.argmax(combo_unit(rs)))]))
        recB["min_pae"].append(_pick(rs,"min_edge_pae",-1))
        fams={}
        for i,x in enumerate(rs): fams.setdefault(x["family"],[]).append(i)
        fk=list(fams); picks=[dq[rngB.choice(fams[fk[rngB.randint(len(fk))]])] for _ in range(200)]
        recB["strat_random"].append(float(np.mean(picks)))
    _report(f"### 분석 B — 확장 pool 전체(co-folder+감소MSA+HADDOCK) [SECONDARY] ###",
            recB, methsB, ref="strat_random", orc_key="ORACLE")

    with open(OUT,"w",newline="") as f:
        w=csv.writer(f); cols=sorted(by); w.writerow(["set","method"]+cols)
        for m in methsA: w.writerow(["A",m]+[("" if math.isnan(v) else f"{v:.3f}") for v in recA[m]])
        for m in methsB: w.writerow(["B",m]+[("" if math.isnan(v) else f"{v:.3f}") for v in recB[m]])
    print(f"\n→ {OUT}")
    print("판정 규율(정직):")
    print(" • PRIMARY = 분석A 공통집합서 '최고 단일피처'(고정부호, 튜닝0)가 ipTM을 이기나(Win@.49 + Regret CI). combo 아님.")
    print(" • min_pae(Boltz·AF2 PAE 있을 때)가 최강 예상(초안 근거). 물리피처는 PAE 없는 pool(tFold/HADDOCK) 보완.")
    print(" • n=10 → 큰·robust 효과만 GO. NO-GO는 underpowered(음성 무정보) — 사전등록 더 큰 셋으로 진행이 GO의 뜻.")
    print(" • HADDOCK=ab-initio(CDR active+항원 표면 passive, native 에피토프 restraint 없음)=blind. min-DockQ=선택 상한.")
    print(" • DockQ/oracle=평가 라벨(피처·선택 미사용). ipTM=모델 자기신뢰도. protenix_ownIptm=실전 단일모델 baseline / ORACLE=상한.")

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    sub=ap.add_subparsers(dest="cmd", required=True)
    b=sub.add_parser("build"); b.add_argument("--targets",default=""); b.add_argument("--fast",action="store_true")
    b.add_argument("--limit-per-source",type=int,default=0); b.add_argument("--no-pae",action="store_true")
    an=sub.add_parser("analyze")
    a=ap.parse_args()
    if a.cmd=="build": cmd_build(a)
    else: cmd_analyze(a)
