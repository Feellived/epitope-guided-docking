#!/usr/bin/env python3
"""[Study B/E2] 과대표집 편향 테스트용 RBD 복합체 선별 = 9K6J류(진짜 에피토프가 과대표집 hotspot에서 벗어난 항체).

과대표집 hotspot = RBD 지도 상위 잔기(E484/F486/Q493 등 escape·ACE2 코어).
9K6J처럼 '진짜 자리가 hotspot을 안 덮는' 항체여야 → 모델이 (틀리게) hotspot으로 쏠리나 볼 수 있음.
(⚠️ 'non-RBM'이 아님 — RBM 안이어도 hotspot만 아니면 됨. non-RBM은 희귀, off-hotspot은 흔함.)

각 최근 RBD 복합체의 에피토프 → hotspot 겹침% 계산 → 낮은 순 = 9K6J류 후보. 항원·H/L 식별해 yaml 출력.
사용: python select_rbd_subdominant.py --max 150 --since 2024-01-01 --hot-max 0.6 --out-yaml configs/targets_rbd_sub.yaml
  ⚠️ --since: 높을수록 leakage-clean(모델 컷오프 이후)이나 후보 적음. 2025-07+면 전 모델 clean이나 매우 적음.
"""
import argparse, io, json, os, re, urllib.request
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_3to1

REF=("NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRG"
     "DEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGST"
     "PCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF")
OFF=330; T3TO1={k.upper():v for k,v in protein_letters_3to1.items()}
# 과대표집 hotspot (우리 RBD 지도 상위 잔기; auth 번호)
HOT=set([449,452,453,455,456,475,476,477,484,485,486,487,489,490,492,493,494,498,501,502,505])
_al=PairwiseAligner();_al.mode="global";_al.match_score=1;_al.mismatch_score=-1;_al.open_gap_score=-3;_al.extend_gap_score=-0.5

def ident(a,b): return _al.score(a,b)/max(len(a),len(b),1) if a and b else -1.0
def amap(sf,st):
    if not sf or not st: return {}
    aln=_al.align(sf,st)[0]; m={}
    for (a0,a1),(b0,b1) in zip(aln.aligned[0],aln.aligned[1]):
        for k in range(a1-a0): m[a0+k]=b0+k
    return m
def cseq(ch): return "".join(T3TO1.get(r.get_resname().upper(),"") for r in ch if "CA" in r)

def fasta_roles(pdb):
    try: txt=urllib.request.urlopen(f"https://www.rcsb.org/fasta/entry/{pdb}/display",timeout=30).read().decode()
    except Exception: return {}
    out={}
    for r in txt.split(">")[1:]:
        head=r.split("\n")[0]; chunk=head.split("|")
        if len(chunk)<3: continue
        chains=set(re.findall(r'auth ([A-Za-z0-9]+)',chunk[1])) or set(re.findall(r'\b([A-Za-z0-9])\b',chunk[1].replace("Chains","").replace("Chain","")))
        name=chunk[2].lower()
        role="H" if ("heavy" in name or re.search(r'\bh chain',name)) else ("L" if ("light" in name or "kappa" in name or "lambda" in name or re.search(r'\bl chain',name)) else None)
        if role:
            for c in chains: out[c]=role
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--max",type=int,default=150); ap.add_argument("--since",default="2024-01-01")
    ap.add_argument("--res",type=float,default=3.2); ap.add_argument("--hot-max",type=float,default=0.6,help="hotspot겹침 이하만 9K6J류 후보")
    ap.add_argument("--out-yaml",default="configs/targets_rbd_sub.yaml"); ap.add_argument("--take",type=int,default=12)
    a=ap.parse_args()
    q={"query":{"type":"group","logical_operator":"and","nodes":[
        {"type":"terminal","service":"sequence","parameters":{"evalue_cutoff":1,"identity_cutoff":0.9,"sequence_type":"protein","value":REF}},
        {"type":"terminal","service":"text","parameters":{"attribute":"rcsb_assembly_info.polymer_entity_instance_count","operator":"greater_or_equal","value":2}},
        {"type":"terminal","service":"text","parameters":{"attribute":"rcsb_entry_info.resolution_combined","operator":"less_or_equal","value":a.res}},
        {"type":"terminal","service":"text","parameters":{"attribute":"rcsb_accession_info.initial_release_date","operator":"greater","value":a.since}}]},
       "return_type":"entry","request_options":{"paginate":{"start":0,"rows":a.max}}}
    req=urllib.request.Request("https://search.rcsb.org/rcsbsearch/v2/query",data=json.dumps(q).encode(),headers={"Content-Type":"application/json"})
    pdbs=[r["identifier"] for r in json.load(urllib.request.urlopen(req,timeout=60)).get("result_set",[])]
    print(f"최근({a.since}+) RBD 복합체 {len(pdbs)}개 분류 중 (hotspot 겹침%)...")
    cand=[]
    for i,pdb in enumerate(pdbs,1):
        try:
            m=MMCIFParser(QUIET=True).get_structure(pdb,io.StringIO(urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb}.cif",timeout=40).read().decode("utf-8","replace")))[0]
            seqs={c.id:cseq(c) for c in m if len([r for r in c if "CA" in r])>=8}
            ag=max(seqs,key=lambda c:ident(seqs[c],REF))
            if ident(seqs[ag],REF)<0.85: continue
            bid=[c for c in seqs if c!=ag and len(seqs[c])>=30 and ident(seqs[c],REF)<0.5]
            if not bid: continue
            bat=np.array([at.get_coord() for c in bid for r in m[c] for at in r])
            A=[r for r in m[ag] if "CA" in r]; cm=amap(seqs[ag],REF); hit=set()
            for idx,r in enumerate(A):
                if idx not in cm: continue
                aa=np.array([at.get_coord() for at in r])
                if np.sqrt(((aa[:,None,:]-bat[None,:,:])**2).sum(-1)).min()<=5: hit.add(cm[idx]+1+OFF)
            if len(hit)<4: continue
            hf=len(hit&HOT)/len(hit); roles=fasta_roles(pdb); nano=len(bid)==1
            H=[c for c in bid if roles.get(c)=="H"]; L=[c for c in bid if roles.get(c)=="L"]
            hc = bid[0] if nano else (H[0] if H else (bid[0] if bid else None))  # 나노바디=단일사슬을 H로
            lc = None if nano else (L[0] if L else (bid[1] if len(bid)>1 else None))
            cand.append({"pdb":pdb,"ag":ag,"H":hc,"L":lc,
                         "n":len(hit),"hot":hf,"nano":nano,"off":sorted(x for x in hit if x not in HOT)})
        except Exception: continue
        if i%15==0: print(f"  {i}/{len(pdbs)} (후보 {len(cand)})")
    cand.sort(key=lambda x:x["hot"])
    print(f"\n=== hotspot 겹침 낮은 순 (=9K6J류: 진짜자리가 과대표집 hotspot 밖) ===")
    for c in cand[:30]:
        hl=f"H={c['H']} L={c['L']}" if not c['nano'] else f"nanobody({c['H']})"
        print(f"  {c['pdb']} hotspot={c['hot']:.0%} epi={c['n']:2} {hl}  off={c['off'][:6]}")
    sub=[c for c in cand if c["hot"]<=a.hot_max and (c["nano"] or (c["H"] and c["L"]))][:a.take]
    print(f"\nhotspot<={a.hot_max:.0%} & 사슬OK = {len(sub)}개 선택 → {a.out_yaml}")
    os.makedirs(os.path.dirname(a.out_yaml) or ".",exist_ok=True)
    with open(a.out_yaml,"w") as f:
        f.write("# Study B/E2: RBD 9K6J류(off-hotspot) 세트 (select_rbd_subdominant.py 생성)\nrbd_subdominant:\n")
        for c in sub:
            if c["nano"]:
                f.write(f"  - id: {c['pdb']}\n    format: nanobody\n    antibody: {{H: {c['H']}}}\n    antigen: {{chain: {c['ag']}, name: SARS-CoV-2 RBD}}\n    keep_chains: [{c['ag']}, {c['H']}]\n")
            else:
                f.write(f"  - id: {c['pdb']}\n    format: Fab\n    antibody: {{H: {c['H']}, L: {c['L']}}}\n    antigen: {{chain: {c['ag']}, name: SARS-CoV-2 RBD}}\n    keep_chains: [{c['ag']}, {c['H']}, {c['L']}]\n")
    print("확인: 위 표 hotspot% 낮은 것들 = 모델이 hotspot으로 쏠리는지 볼 테스트셋. 예측 후 where_docked로 검증.")

if __name__=="__main__":
    main()
