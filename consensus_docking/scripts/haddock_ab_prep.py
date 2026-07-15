#!/usr/bin/env python3
"""HADDOCK 항체 준비: Kabat CDR active 식별 + 단일사슬 재넘버 → clean PDB + actpass.
Kabat: H 31-35/50-65/95-102, L 24-34/50-56/89-97 (삽입코드 자동 포함). ANARCI 불필요.
사용: python haddock_ab_prep.py --target 8SIT --dir haddock/8sit  (biopython)"""
import argparse, json, os
from Bio.PDB import PDBParser, PDBIO, Chain, Model, Structure
KABAT={"H":[(31,35),(50,65),(95,102)], "L":[(24,34),(50,56),(89,97)]}
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--run-dir",default="runs_rbd"); ap.add_argument("--target",required=True)
    ap.add_argument("--dir",required=True); ap.add_argument("--segid",default="A"); a=ap.parse_args()
    d=json.load(open(os.path.join(a.run_dir,a.target.lower(),"chains.json")))
    ab=d["antibody"]; ab=[ab] if isinstance(ab,str) else [str(x) for x in ab]
    heavy=ab[0]; light=ab[1] if len(ab)>1 else None
    m=PDBParser(QUIET=True).get_structure("x",os.path.join(a.dir,"antibody.pdb"))[0]
    ids=[c.id for c in m]
    nc=Chain.Chain(a.segid); new=1; active=[]
    for cid in ([heavy]+([light] if light else [])):
        if cid not in ids: continue
        ranges=KABAT["H"] if cid==heavy else KABAT["L"]
        for r in list(m[cid]):
            if "CA" not in r: continue
            if any(lo<=r.id[1]<=hi for lo,hi in ranges): active.append(new)
            r.detach_parent(); r.id=(' ',new,' '); nc.add(r); new+=1
    st=Structure.Structure("ab"); md=Model.Model(0); md.add(nc); st.add(md)
    io=PDBIO(); io.set_structure(st); io.save(os.path.join(a.dir,"antibody_clean.pdb"))
    act=sorted(set(active))
    open(os.path.join(a.dir,"antibody.actpass"),"w").write(" ".join(map(str,act))+"\n\n")
    print(f"{a.target}: CDR active {len(act)}개, 단일사슬 {a.segid} 재넘버(총 {new-1}res) → antibody_clean.pdb, antibody.actpass")
    print("  active:", act)
if __name__=="__main__": main()
