#!/usr/bin/env python3
"""MSA-depth 실험용 Boltz YAML: 항원 chain=subsampled a3m, 항체=single-seq(msa: empty).
항원 MSA 깊이만 변화 → drift 변화 격리. 사용:
  python make_boltz_msa_depth.py --run-dir runs_rbd --target 8SDF --antigen-a3m /abs/A_d8.a3m --out out.yaml"""
import argparse, json, os
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--run-dir",required=True); ap.add_argument("--target",required=True)
    ap.add_argument("--antigen-a3m",required=True); ap.add_argument("--out",required=True)
    a=ap.parse_args()
    d=json.load(open(os.path.join(a.run_dir,a.target.lower(),"chains.json")))
    ag=str(d["antigen"])
    L=["version: 1","sequences:"]
    for c in d["chains"]:
        msa=a.antigen_a3m if c["id"]==ag else "empty"
        L+=["  - protein:",f"      id: {c['id']}",f'      sequence: "{c["seq"]}"',f"      msa: {msa}"]
    open(a.out,"w").write("\n".join(L)+"\n")
    print(f"[boltz-msa] {a.target} 항원={ag} msa={a.antigen_a3m}")
if __name__=="__main__": main()
