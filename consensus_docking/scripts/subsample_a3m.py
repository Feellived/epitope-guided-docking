#!/usr/bin/env python3
"""a3m를 depth N으로 subsample (query 유지 + 나머지 무작위 N-1). MSA-depth 인과실험용.
사용: python subsample_a3m.py in.a3m out.a3m --depth 8 --seed 0"""
import argparse, random
def read_a3m(p):
    recs=[]; hdr=None; buf=[]
    for line in open(p, encoding="utf-8", errors="replace"):
        line=line.rstrip("\n")
        if line.startswith(">"):
            if hdr is not None: recs.append((hdr,"".join(buf)))
            hdr=line; buf=[]
        elif hdr is not None: buf.append(line)
    if hdr is not None: recs.append((hdr,"".join(buf)))
    return recs
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("inp"); ap.add_argument("out")
    ap.add_argument("--depth",type=int,required=True,help="query 포함 총 서열 수")
    ap.add_argument("--seed",type=int,default=0)
    a=ap.parse_args()
    recs=read_a3m(a.inp)
    if not recs: raise SystemExit("빈 a3m")
    query=recs[0]; rest=recs[1:]; k=max(0,a.depth-1)
    if k < len(rest): rest=random.Random(a.seed).sample(rest,k)
    keep=[query]+rest
    with open(a.out,"w") as f:
        for h,s in keep: f.write(h+"\n"+s+"\n")
    print(f"{a.inp}: {len(recs)}→{len(keep)} (depth {a.depth}, seed {a.seed}) → {a.out}")
if __name__=="__main__": main()
