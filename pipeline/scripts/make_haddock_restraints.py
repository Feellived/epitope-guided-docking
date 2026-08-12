#!/usr/bin/env python3
"""HADDOCK 제약+config: 항원 표면 passive + ambig.tbl + run.cfg. (항체 segid A, 항원 segid B 고정)
실행(haddock3 env): python make_haddock_restraints.py --target 8SIT --dir haddock/8sit --sampling 1000"""
import argparse, glob, os, subprocess
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--run-dir",default="runs_rbd"); ap.add_argument("--target",required=True)
    ap.add_argument("--dir",required=True); ap.add_argument("--sampling",type=int,default=1000); a=ap.parse_args()
    dd=a.dir
    for f in glob.glob(os.path.join(dd,"*_passive_*.actpass")): os.remove(f)   # 재실행 위해 정리
    subprocess.run(["haddock3-restraints","calc_accessibility","antigen.pdb","-c","0.4","--export_to_actpass"],check=True,cwd=dd)
    cand=glob.glob(os.path.join(dd,"*_passive_*.actpass"))
    if not cand: raise SystemExit("항원 passive actpass 생성 실패")
    ag_ap=os.path.basename(cand[0]); print("  항원 actpass:", ag_ap)
    with open(os.path.join(dd,"ambig.tbl"),"w") as tbl:
        subprocess.run(["haddock3-restraints","active_passive_to_ambig","antibody.actpass",ag_ap,
                        "--segid-one","A","--segid-two","B"],check=True,cwd=dd,stdout=tbl)
    print("  ambig.tbl (항체 A <-> 항원 B)")
    dabs=os.path.abspath(dd)
    open(os.path.join(dd,"run.cfg"),"w").write(
      f'run_dir = "{dabs}/run"\nmolecules = ["{dabs}/antibody_clean.pdb", "{dabs}/antigen.pdb"]\n'
      f'[topoaa]\n[rigidbody]\nambig_fname = "{dabs}/ambig.tbl"\nsampling = {a.sampling}\n'
      f'[seletop]\nselect = 200\n[flexref]\n[emref]\n[clustfcc]\n[seletopclusts]\ntop_models = 4\n')
    print(f"  run.cfg (sampling {a.sampling}) -> haddock3 {dd}/run.cfg")
if __name__=="__main__": main()
