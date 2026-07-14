#!/usr/bin/env python3
"""[tFold 픽스] colabfold/mmseqs a3m의 불량 행 제거 = tFold MSA 파서 KeyError('\x00') 방지.

원인: a3m의 각 행은 query와 같은 match-컬럼 수여야 하는데, mmseqs가 가끔 정렬길이가 다른 파편 행을
      1개 섞음(예: 길이 194 정상 143행 + 길이 6 불량 1행). tFold가 짧은 행을 194로 null(\x00) 패딩 →
      HHBLITS_AA_TO_ID['\x00'] KeyError. → query와 match-컬럼 수가 다르거나 비표준문자 행을 버림.

순수 stdlib(env 무관). 사용:
  python3 sanitize_a3m.py runs_rbd                     # runs_rbd/*/msa_*/*.a3m 전부
  python3 sanitize_a3m.py runs_diverse runs_rbd        # 여러 run-dir
  python3 sanitize_a3m.py path/to/A.a3m                # 개별 파일
백업: 첫 정리 시 <file>.orig 보존(있으면 건너뜀).
"""
import glob, os, re, shutil, sys

TFOLD_AA = set("ACDEFGHIKLMNPQRSTVWYXBZJOU-")    # tFold가 매핑하는 match-컬럼 알파벳(대문자+gap)
def match_cols(seq): return re.sub("[a-z]", "", seq)   # 소문자(insertion) 제거 → tFold가 실제 보는 컬럼

def clean(path):
    recs=[]; hdr=None; buf=[]
    for line in open(path, encoding="utf-8", errors="replace"):
        line=line.rstrip("\n")
        if line.startswith(">"):
            if hdr is not None: recs.append((hdr,"".join(buf)))
            hdr=line; buf=[]
        elif hdr is not None: buf.append(line)
    if hdr is not None: recs.append((hdr,"".join(buf)))
    if not recs: return 0,0,0
    qlen=len(match_cols(recs[0][1])); keep=[recs[0]]   # query는 항상 유지
    for h,s in recs[1:]:
        mc=match_cols(s)                                # tFold 관점: 길이==query & 알파벳 정상이어야
        if len(mc)==qlen and set(mc) <= TFOLD_AA: keep.append((h,s))
    if len(keep)!=len(recs):
        if not os.path.exists(path+".orig"): shutil.copy2(path, path+".orig")
        with open(path,"w") as f:
            for h,s in keep: f.write(h+"\n"+s+"\n")
    return len(recs), len(keep), qlen

def main():
    args=sys.argv[1:] or ["runs_diverse"]
    files=[]
    for a in args:
        files += [a] if a.endswith(".a3m") else glob.glob(os.path.join(a,"*","msa_*","*.a3m"))
    bad=0
    for f in sorted(files):
        n,k,q=clean(f)
        tag="OK" if n==k else f"불량 {n-k}행 제거"
        if n!=k: bad+=1
        print(f"  {f}: {n}→{k}행 (query len {q}) {tag}")
    print(f"완료: {len(files)}파일 중 {bad}개 정리")

if __name__=="__main__": main()
