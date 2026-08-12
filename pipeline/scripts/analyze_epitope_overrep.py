#!/usr/bin/env python3
"""[Study B · E1] 에피토프 과대표집 지도 = 한 항원의 어느 잔기에 항체가 몇 번 붙었나.

한 항원(기본 SARS-CoV-2 RBD)의 PDB 복합체를 전부 모아 → 각 복합체의 에피토프(항체/binder와
5Å 내 항원 잔기) 추출 → 항원 기준서열 위치로 정렬·집계 → '잔기별 결합 빈도' = 과대표집 패치.

검증 목표: RBD면 E484/F486·RBM(437-508)이 hot으로 떠야 함(= Phase 0에서 실패 모델들이 쏠린 곳,
Barnes class1). 뜨면 방법 검증 완료 → 다른 항원에 --ref-fasta·--rcsb-seq로 재사용.

방법: 항원 사슬 = 기준서열 최고 매칭 / binder = 그 외 단백질 사슬(항원 복사본 제외, len>=30).
      에피토프 = binder 원자와 5Å 내 항원 잔기 → 기준서열 위치로 alignment 매핑.
      ⚠️ 1차 버전은 binder에 ACE2 등 비-항체도 포함(RBM 신호 강화 = 과대표집과 같은 방향, 무해).
         항체-only가 필요하면 --antibody-only(ANARCI 필요) 추후.

설치: biopython (+ 인터넷). 사용:
  python analyze_epitope_overrep.py                      # RBD, RCSB에서 자동 목록(res<=3.0), max 200
  python analyze_epitope_overrep.py --max 400
  python analyze_epitope_overrep.py --pdb-list ids.txt   # 직접 목록
결과: results/epitope_overrep_<tag>.csv  (auth_resnum, ref_pos, aa, n_complexes, freq) + 상위 출력
"""
import argparse, csv, io, json, os, urllib.request
import numpy as np
from Bio.PDB import MMCIFParser
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_3to1

# 기본 기준서열 = SARS-CoV-2 RBD (Spike 331~ ; ref_pos 1 = Spike 331 → auth = ref_pos + OFFSET)
DEFAULT_REF = ("NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRG"
               "DEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGST"
               "PCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF")
DEFAULT_OFFSET = 330       # ref_pos + 330 = Spike auth number (E484 = ref_pos 154)
CONTACT = 5.0
T3TO1 = {k.upper(): v for k, v in protein_letters_3to1.items()}
_al = PairwiseAligner()
_al.mode = "global"; _al.match_score = 1; _al.mismatch_score = -1
_al.open_gap_score = -3; _al.extend_gap_score = -0.5


def region_rbd(auth):
    if 437 <= auth <= 508: return "RBM"
    if 369 <= auth <= 386: return "core(class4)"
    if (337 <= auth <= 346) or auth in (356, 440, 441, 444): return "class3"
    return "other"


def rcsb_seq_search(ref, max_res, rows, identity=0.9):
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "sequence",
             "parameters": {"evalue_cutoff": 1, "identity_cutoff": identity, "sequence_type": "protein", "value": ref}},
            {"type": "terminal", "service": "text",
             "parameters": {"attribute": "rcsb_assembly_info.polymer_entity_instance_count", "operator": "greater_or_equal", "value": 2}},
            {"type": "terminal", "service": "text",
             "parameters": {"attribute": "rcsb_entry_info.resolution_combined", "operator": "less_or_equal", "value": max_res}},
         ]},
         "return_type": "entry",
         "request_options": {"paginate": {"start": 0, "rows": rows}}}
    req = urllib.request.Request("https://search.rcsb.org/rcsbsearch/v2/query",
                                 data=json.dumps(q).encode(), headers={"Content-Type": "application/json"})
    d = json.load(urllib.request.urlopen(req, timeout=60))
    return [r["identifier"] for r in d.get("result_set", [])]


def fetch_cif(pdb):
    url = f"https://files.rcsb.org/download/{pdb}.cif"
    txt = urllib.request.urlopen(url, timeout=60).read().decode("utf-8", "replace")
    return MMCIFParser(QUIET=True).get_structure(pdb, io.StringIO(txt))[0]


def res_list(chain):
    """[(residue_obj, one-letter)] CA 있는 잔기만, 순서대로."""
    out = []
    for r in chain:
        if "CA" in r:
            aa = T3TO1.get(r.get_resname().upper())
            if aa: out.append((r, aa))
    return out


def ident(a, b):
    if not a or not b: return -1.0
    return _al.score(a, b) / max(len(a), len(b), 1)


def aln_map(s_from, s_to):
    if not s_from or not s_to: return {}
    aln = _al.align(s_from, s_to)[0]
    m = {}
    for (a0, a1), (b0, b1) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(a1 - a0): m[a0 + k] = b0 + k
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-fasta", help="기준 항원서열 FASTA(없으면 SARS-CoV-2 RBD)")
    ap.add_argument("--offset", type=int, default=DEFAULT_OFFSET, help="ref_pos+offset=auth 번호(RBD=330)")
    ap.add_argument("--pdb-list", help="PDB id 목록 파일(없으면 RCSB 자동 검색)")
    ap.add_argument("--max", type=int, default=200)
    ap.add_argument("--max-res", type=float, default=3.0)
    ap.add_argument("--identity", type=float, default=0.9, help="RCSB 서열검색 identity(비-RBD 항원=0.7 권장, isoform 포함)")
    ap.add_argument("--min-ag-ident", type=float, default=0.8, help="항원 사슬 매칭 최소 identity(다도메인/fragment 항원은 낮춰야)")
    ap.add_argument("--tag", default="rbd")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    ref = DEFAULT_REF
    if a.ref_fasta:
        ref = "".join(l.strip() for l in open(a.ref_fasta) if not l.startswith(">"))
    L = len(ref)

    if a.pdb_list:
        pdbs = [x.strip().upper() for x in open(a.pdb_list) if x.strip()]
    else:
        print(f"[RCSB] 기준서열({L}aa) 복합체 검색 (res<={a.max_res}, identity>={a.identity}) ...")
        try:
            pdbs = rcsb_seq_search(ref, a.max_res, a.max, a.identity)
        except Exception as e:                                   # 검색 실패해도 안 죽고 빈 지도 저장
            print(f"!! RCSB 검색 실패({type(e).__name__}: {e}) → 빈 지도. --identity 낮추거나 --pdb-list 사용")
            pdbs = []
    pdbs = pdbs[:a.max]
    print(f"대상 PDB {len(pdbs)}개 (max {a.max})")
    if not pdbs:
        print("!! 검색 0개 — 원인: 긴서열/isoform/네트워크. --identity 0.7 또는 --pdb-list 권장")

    freq = np.zeros(L, dtype=int)      # ref position별 몇 개 복합체가 에피토프로
    used = 0; skipped = 0
    for i, pdb in enumerate(pdbs, 1):
        try:
            m = fetch_cif(pdb)
            chains = {ch.id: res_list(ch) for ch in m}
            chains = {k: v for k, v in chains.items() if len(v) >= 20}
            if not chains: skipped += 1; continue
            seqs = {k: "".join(aa for _, aa in v) for k, v in chains.items()}
            ag = max(seqs, key=lambda c: ident(seqs[c], ref))
            if ident(seqs[ag], ref) < a.min_ag_ident: skipped += 1; continue   # 항원 사슬 못 찾음
            # binder = 항원 아님 + RBD 복사본 아님(ident<0.5) + len>=30
            binder_atoms = []
            for c, s in seqs.items():
                if c == ag: continue
                if len(chains[c]) < 30: continue
                if ident(s, ref) >= 0.5: continue                    # 다른 RBD/spike 복사본 제외
                for r, _ in chains[c]:
                    for at in r: binder_atoms.append(at.get_coord())
            if not binder_atoms: skipped += 1; continue
            B = np.array(binder_atoms)
            # 항원 잔기별 min 거리 → 에피토프
            cmap = aln_map(seqs[ag], ref)                            # 항원 사슬 idx → ref idx
            hit_ref = set()
            for idx, (r, _) in enumerate(chains[ag]):
                if idx not in cmap: continue
                ag_at = np.array([at.get_coord() for at in r])
                d = np.sqrt(((ag_at[:, None, :] - B[None, :, :]) ** 2).sum(-1)).min()
                if d <= CONTACT: hit_ref.add(cmap[idx])
            for rp in hit_ref: freq[rp] += 1
            used += 1
        except Exception as e:
            skipped += 1
        if i % 25 == 0: print(f"  {i}/{len(pdbs)}  (사용 {used}, 스킵 {skipped})")

    out = a.out or f"results/epitope_overrep_{a.tag}.csv"
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["ref_pos", "auth_resnum", "aa", "n_complexes", "freq", "region"])
        for p in range(L):
            auth = p + 1 + a.offset
            w.writerow([p + 1, auth, ref[p], int(freq[p]), round(freq[p] / max(used, 1), 3),
                        region_rbd(auth) if not a.ref_fasta else ""])
    print(f"\n[저장] {out}  (복합체 사용 {used} / 스킵 {skipped})")
    order = np.argsort(-freq)
    print(f"\n=== 과대표집 상위 25 잔기 (auth번호, aa, 결합복합체수, 빈도) ===")
    for p in order[:25]:
        auth = p + 1 + a.offset
        reg = f"  [{region_rbd(auth)}]" if not a.ref_fasta else ""
        print(f"  {auth:4} {ref[p]}   n={int(freq[p]):4}  freq={freq[p]/max(used,1):.2f}{reg}")
    if not a.ref_fasta:
        from collections import Counter
        top = [p + 1 + a.offset for p in order[:40]]
        c = Counter(region_rbd(x) for x in top)
        print(f"\n상위40 영역분포: RBM {c['RBM']} · core {c['core(class4)']} · class3 {c['class3']} · other {c['other']}")
        print("검증: RBM/E484(484)·F486(486) 상위면 = 과대표집 지도 방법 OK.")


if __name__ == "__main__":
    main()
