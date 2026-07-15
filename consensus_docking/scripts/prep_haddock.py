#!/usr/bin/env python3
"""HADDOCK 입력: native PDB에서 항체(H+L)·항원 분리 → 클린 PDB 2개. biopython만(설치 0).
HETATM/water/altloc 제거. 사용: python prep_haddock.py --target 8SIT --out-dir haddock/8sit"""
import argparse, json, os
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
def load(p):
    par = MMCIFParser(QUIET=True) if p.endswith(".cif") else PDBParser(QUIET=True)
    return par.get_structure("x", p)[0]
class Sel(Select):
    def __init__(self, chains): self.ch=set(chains)
    def accept_chain(self, c): return c.id in self.ch
    def accept_residue(self, r): return r.id[0]==' '                 # 표준 잔기만
    def accept_atom(self, a): return (not a.is_disordered()) or a.get_altloc() in ('A',' ')  # altloc A만
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--run-dir",default="runs_rbd"); ap.add_argument("--target",required=True)
    ap.add_argument("--out-dir",required=True); a=ap.parse_args()
    rd=os.path.join(a.run_dir,a.target.lower())
    d=json.load(open(os.path.join(rd,"chains.json")))
    ag=str(d["antigen"]); ab=d["antibody"]; ab=[ab] if isinstance(ab,str) else [str(x) for x in ab]
    nat=[f for f in os.listdir(rd) if a.target.lower() in f.lower() and f.endswith((".pdb",".cif"))
         and not any(m in f.lower() for m in ("boltz","chai","protenix","colabfold","tfold"))][0]
    m=load(os.path.join(rd,nat)); os.makedirs(a.out_dir,exist_ok=True)
    io=PDBIO(); io.set_structure(m)
    io.save(os.path.join(a.out_dir,"antigen.pdb"), Sel([ag]))
    io.save(os.path.join(a.out_dir,"antibody.pdb"), Sel(ab))
    # 잔기 수 리포트
    ag_n=len([r for c in m if c.id==ag for r in c if r.id[0]==' '])
    ab_n=len([r for c in m if c.id in ab for r in c if r.id[0]==' '])
    print(f"{a.target}: 항원={ag}({ag_n}res) 항체={ab}({ab_n}res) native={nat} → {a.out_dir}/antigen.pdb·antibody.pdb")
if __name__=="__main__": main()
