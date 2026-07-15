#!/usr/bin/env python3
"""HADDOCK 입력: native → 항체(H+L 원본사슬) + 항원(단일사슬 B로 rename). segid 충돌 방지(항체=A, 항원=B).
biopython만(설치 0). 사용: python prep_haddock.py --target 8SIT --out-dir haddock/8sit"""
import argparse, json, os
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Chain, Model, Structure
def load(p):
    return (MMCIFParser(QUIET=True) if p.endswith(".cif") else PDBParser(QUIET=True)).get_structure("x",p)[0]
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--run-dir",default="runs_rbd"); ap.add_argument("--target",required=True)
    ap.add_argument("--out-dir",required=True); a=ap.parse_args()
    rd=os.path.join(a.run_dir,a.target.lower())
    d=json.load(open(os.path.join(rd,"chains.json")))
    ag=str(d["antigen"]); ab=d["antibody"]; ab=[ab] if isinstance(ab,str) else [str(x) for x in ab]
    nat=[f for f in os.listdir(rd) if a.target.lower() in f.lower() and f.endswith((".pdb",".cif"))
         and not any(mm in f.lower() for mm in ("boltz","chai","protenix","colabfold","tfold"))][0]
    m=load(os.path.join(rd,nat)); ids=[c.id for c in m]; os.makedirs(a.out_dir,exist_ok=True)
    # 항원 → 단일사슬 B (원래 잔기번호 유지)
    st=Structure.Structure("ag"); md=Model.Model(0); bc=Chain.Chain("B")
    if ag in ids:
        for r in list(m[ag]):
            if r.id[0]==' ' and "CA" in r: r.detach_parent(); bc.add(r)
    md.add(bc); st.add(md); io=PDBIO(); io.set_structure(st); io.save(os.path.join(a.out_dir,"antigen.pdb"))
    # 항체 → 원본 사슬(H/L; haddock_ab_prep이 A로 병합)
    st2=Structure.Structure("ab"); md2=Model.Model(0)
    for cid in ab:
        if cid not in ids: continue
        ch=Chain.Chain(cid)
        for r in list(m[cid]):
            if r.id[0]==' ' and "CA" in r: r.detach_parent(); ch.add(r)
        md2.add(ch)
    st2.add(md2); io2=PDBIO(); io2.set_structure(st2); io2.save(os.path.join(a.out_dir,"antibody.pdb"))
    print(f"{a.target}: 항원 {ag}->B({len(list(bc))}res), 항체 {ab} -> antigen.pdb·antibody.pdb")
if __name__=="__main__": main()
