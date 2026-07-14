#!/usr/bin/env python3
"""[Study B/E1-Foldseek] Foldseek multimersearch 결과로 복합체를 결합모드별 클러스터링 →
클러스터(=비중복 결합모드) 기준으로 에피토프 과대표집 지도. (easy-multimercluster가 macOS 크래시라 우회)

입력: complexes/ (prep_complexes.py 산출, chain A=항원)  +  msearch_report (foldseek easy-multimersearch)
방법: 복합체쌍 TM-score(msearch_report 5열) >= --tm 이면 같은 모드 → 연결성분 클러스터.
      각 클러스터 대표의 에피토프(A와 5Å내 binder 잔기) → 기준서열 위치.
출력: raw(복합체별) vs clustered(클러스터별) 잔기 빈도 → 중복이 얼마나 부풀렸나 대비.
"""
import argparse, csv, glob, os
import numpy as np
from Bio.PDB import PDBParser
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_3to1

DEFAULT_REF = ("NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRG"
               "DEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGST"
               "PCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF")
OFFSET = 330
T3TO1 = {k.upper(): v for k, v in protein_letters_3to1.items()}
_al = PairwiseAligner(); _al.mode = "global"; _al.match_score = 1; _al.mismatch_score = -1
_al.open_gap_score = -3; _al.extend_gap_score = -0.5


def region(auth):
    if 437 <= auth <= 508: return "RBM"
    if 369 <= auth <= 386: return "core"
    if (337 <= auth <= 346) or auth in (356, 440, 441, 444): return "class3"
    return "other"


def aln_map(sf, st):
    if not sf or not st: return {}
    aln = _al.align(sf, st)[0]; m = {}
    for (a0, a1), (b0, b1) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(a1 - a0): m[a0 + k] = b0 + k
    return m


def epitope_refpos(pdb_path, ref):
    """chain A(항원) 잔기 중 binder(그외 사슬)와 5Å내 → ref 위치 집합."""
    m = PDBParser(QUIET=True).get_structure("x", pdb_path)[0]
    if "A" not in [c.id for c in m]: return set()
    A = [(r, "".join(T3TO1.get(r.get_resname().upper(), ""))) for r in m["A"] if "CA" in r]
    aseq = "".join(aa for _, aa in A)
    binder = np.array([at.get_coord() for c in m if c.id != "A" for r in c for at in r])
    if len(binder) == 0 or not aseq: return set()
    cmap = aln_map(aseq, ref); hit = set()
    for idx, (r, _) in enumerate(A):
        if idx not in cmap: continue
        ag = np.array([at.get_coord() for at in r])
        if np.sqrt(((ag[:, None, :] - binder[None, :, :]) ** 2).sum(-1)).min() <= 5.0:
            hit.add(cmap[idx])
    return hit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", default="msearch_report")
    ap.add_argument("--complexes", default="complexes")
    ap.add_argument("--tm", type=float, default=0.5, help="복합체 TM>=이면 같은 결합모드")
    ap.add_argument("--ref-fasta"); ap.add_argument("--out", default="results/epitope_overrep_clustered.csv")
    a = ap.parse_args()
    ref = DEFAULT_REF
    if a.ref_fasta: ref = "".join(l.strip() for l in open(a.ref_fasta) if not l.startswith(">"))
    L = len(ref)

    ids = sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(a.complexes, "*.pdb")))
    idx = {x: i for i, x in enumerate(ids)}
    parent = list(range(len(ids)))
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(x, y): parent[find(x)] = find(y)
    for line in open(a.report):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6: continue
        q, t = f[0], f[1]
        if q == t or q not in idx or t not in idx: continue
        try: tm = float(f[4])
        except ValueError: continue
        if tm >= a.tm: union(idx[q], idx[t])
    clusters = {}
    for x in ids: clusters.setdefault(find(idx[x]), []).append(x)
    clist = list(clusters.values())
    print(f"복합체 {len(ids)}개 → 클러스터 {len(clist)}개 (TM>={a.tm})  = 중복 {len(ids)-len(clist)}개 합쳐짐")

    raw = np.zeros(L, int); clu = np.zeros(L, int)
    for x in ids:
        for rp in epitope_refpos(os.path.join(a.complexes, x + ".pdb"), ref): raw[rp] += 1
    for members in clist:
        rep = members[0]
        for rp in epitope_refpos(os.path.join(a.complexes, rep + ".pdb"), ref): clu[rp] += 1

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["auth", "aa", "raw_n", "clustered_n", "region"])
        for p in range(L):
            auth = p + 1 + OFFSET
            w.writerow([auth, ref[p], int(raw[p]), int(clu[p]), region(auth)])
    print(f"[저장] {a.out}")
    print(f"\n=== 상위 15 (auth aa | raw→clustered) ===")
    for p in np.argsort(-clu)[:15]:
        auth = p + 1 + OFFSET
        print(f"  {auth:4} {ref[p]}  raw={int(raw[p]):3} → clustered={int(clu[p]):3}  [{region(auth)}]")


if __name__ == "__main__":
    main()
