#!/usr/bin/env python3
"""HADDOCK 항체 준비 (ANARCI): CDR active 식별 + 단일사슬 재넘버 → clean PDB + actpass.
IMGT CDR = 27-38/56-65/105-117. 어떤 넘버링이든 처리(순차·Kabat·Chothia·나노바디).
실행: anarci env에서  python prep_haddock_antibody.py --dir haddock/8sit
필요: anarci, biopython."""
import argparse, os
from Bio.PDB import PDBParser, PDBIO, Chain, Model, Structure
from Bio.Data.IUPACData import protein_letters_3to1
from anarci import run_anarci
T3TO1={k.upper():v for k,v in protein_letters_3to1.items()}
IMGT_CDR=[(27,38),(56,65),(105,117)]

def cdr_positions(seq):
    """ANARCI IMGT → 입력서열 위치(0-based) 중 CDR 집합 + 사슬타입."""
    res=run_anarci([("q",seq)], scheme="imgt")
    numbering=res[1]; details=res[2] if len(res)>2 else None
    if not numbering or numbering[0] is None: return set(),"?"
    dom,start,end=numbering[0][0]
    ctype=(details[0][0].get("chain_type","?") if details and details[0] else "?")
    cdr=set(); pos=start
    for (num,ins),aa in dom:
        if aa=="-": continue
        if any(lo<=num<=hi for lo,hi in IMGT_CDR): cdr.add(pos)
        pos+=1
    return cdr,ctype

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--dir",required=True); ap.add_argument("--segid",default="A"); a=ap.parse_args()
    m=PDBParser(QUIET=True).get_structure("x",os.path.join(a.dir,"antibody.pdb"))[0]
    nc=Chain.Chain(a.segid); new=1; active=[]; info=[]
    for ch in m:
        res=[r for r in ch if "CA" in r]
        if len(res)<50: continue
        seq="".join(T3TO1.get(r.get_resname().upper(),"X") for r in res)
        try: cdr,ctype=cdr_positions(seq)
        except Exception as e: print(f"  {ch.id} ANARCI 실패: {e}"); cdr,ctype=set(),"ERR"
        info.append(f"{ch.id}({ctype},{len(res)}res,CDR{len(cdr)})")
        for i,r in enumerate(res):
            if i in cdr: active.append(new)
            r.detach_parent(); r.id=(' ',new,' '); nc.add(r); new+=1
    st=Structure.Structure("ab"); md=Model.Model(0); md.add(nc); st.add(md)
    io=PDBIO(); io.set_structure(st); io.save(os.path.join(a.dir,"antibody_clean.pdb"))
    act=sorted(set(active))
    open(os.path.join(a.dir,"antibody.actpass"),"w").write(" ".join(map(str,act))+"\n\n")
    print(f"{os.path.basename(a.dir)}: {' '.join(info)} | 단일사슬 {a.segid} {new-1}res | CDR active {len(act)} → antibody_clean.pdb, antibody.actpass")
if __name__=="__main__": main()
