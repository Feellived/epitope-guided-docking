#!/usr/bin/env python3
"""HADDOCK3 ab-ag ab-initio config 생성 (CDR active + 항원 surface passive).
사용: python make_haddock_cfg.py --dir haddock/8sit --sampling 1000 --out haddock/8sit/run.cfg
(caprieval 생략 — pose는 우리 where_docked/merge_dockq로 분석)"""
import argparse, os
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--dir",required=True); ap.add_argument("--sampling",type=int,default=1000)
    ap.add_argument("--out",required=True); a=ap.parse_args()
    dd=os.path.abspath(a.dir)
    open(a.out,"w").write(f'''# HADDOCK3 ab-ag ab-initio (항체 CDR active + 항원 surface passive = epitope-agnostic)
run_dir = "{dd}/run"
molecules = ["{dd}/antibody_clean.pdb", "{dd}/antigen.pdb"]

[topoaa]

[rigidbody]
ambig_fname = "{dd}/ambig.tbl"
sampling = {a.sampling}

[seletop]
select = 200

[flexref]

[emref]

[clustfcc]

[seletopclusts]
top_models = 4
''')
    print(f"config → {a.out} (sampling {a.sampling})")
if __name__=="__main__": main()
