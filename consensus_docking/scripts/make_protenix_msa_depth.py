#!/usr/bin/env python3
"""Protenix MSA-depth 입력 JSON: 항원 proteinChain.unpairedMsaPath=subsampled a3m, 항체=query-only a3m(single-seq).
사용: python make_protenix_msa_depth.py --run-dir runs_rbd --target 8SIT --antigen-a3m /abs/A_d8.a3m --dir runs_msad_8/8sit"""
import argparse, json, os
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--run-dir",default="runs_rbd"); ap.add_argument("--target",required=True)
    ap.add_argument("--antigen-a3m",required=True); ap.add_argument("--dir",required=True); a=ap.parse_args()
    d=json.load(open(os.path.join(a.run_dir,a.target.lower(),"chains.json")))
    ag=str(d["antigen"]); ab=d["antibody"]; ab=[ab] if isinstance(ab,str) else [str(x) for x in ab]
    sm={c["id"]:c["seq"] for c in d["chains"]}; os.makedirs(a.dir,exist_ok=True)
    seqs=[{"proteinChain":{"sequence":sm[ag],"count":1,"unpairedMsaPath":os.path.abspath(a.antigen_a3m)}}]
    for i in ab:
        qa=os.path.join(a.dir,f"ab_{i}_query.a3m"); open(qa,"w").write(f">{i}\n{sm[i]}\n")
        seqs.append({"proteinChain":{"sequence":sm[i],"count":1,"unpairedMsaPath":os.path.abspath(qa)}})
    out=os.path.join(a.dir,f"protenix_{a.target.lower()}.json")
    json.dump([{"name":a.target.lower(),"sequences":seqs}], open(out,"w"), indent=2)
    print(f"{a.target}: 항원={ag}(msa={os.path.basename(a.antigen_a3m)}), 항체={ab}(query-only) -> {out}")
if __name__=="__main__": main()
