#!/usr/bin/env python3
"""[ZDOCK 1/2] 입력 준비 — native 쪼갠 bound 단량체 + ANARCI CDR 블록목록.
  runs_diverse/<id>/{chains.json, <ID>.pdb} 를 읽어 zdock/<id>/ 에:
    antigen.pdb   — 항원 단일사슬(chain A)
    antibody.pdb  — 항체 H+L 병합 단일사슬(chain B), 1..N 재넘버
    block_noncdr.txt — 비-CDR 잔기 번호(한 줄에 하나) = block.pl 입력
  ⚠️ block.pl은 '단일사슬 + 잔기번호 매칭'만 → H+L을 한 사슬로 병합·유일번호화해야 CDR만 열 수 있음.
     CDR = ANARCI IMGT(27-38/56-65/105-117). haddock_ab_prep.py 로직 재사용.
실행(anarci env):  python scripts/zdock_prep.py [--only 9K6J] [--runs runs_diverse] [--out zdock]
필요: anarci, biopython.
"""
import argparse, glob, json, os
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Chain, Model, Structure
from Bio.Data.IUPACData import protein_letters_3to1
from anarci import run_anarci

T3 = {k.upper(): v for k, v in protein_letters_3to1.items()}
IMGT_CDR = [(27, 38), (56, 65), (105, 117)]

def load(p):
    return (MMCIFParser(QUIET=True) if p.endswith(".cif") else PDBParser(QUIET=True)).get_structure("x", p)[0]

def residues(chain):
    return [r for r in chain if "CA" in r]

def cdr_index(seq):
    """ANARCI IMGT → 입력서열 내 CDR 위치(0-based 인덱스 집합)."""
    res = run_anarci([("q", seq)], scheme="imgt")
    numbering = res[1]
    if not numbering or numbering[0] is None:
        return set()
    dom, start, end = numbering[0][0]
    cdr = set(); pos = start
    for (num, ins), aa in dom:
        if aa == "-":
            continue
        if any(lo <= num <= hi for lo, hi in IMGT_CDR):
            cdr.add(pos)
        pos += 1
    return cdr

def write_single(res_lists, cid, out):
    """여러 residue-list를 한 사슬 cid로 병합, 1..N 순차 재넘버 → PDB 저장. 총 잔기수 반환."""
    st = Structure.Structure("x"); md = Model.Model(0); nc = Chain.Chain(cid); n = 1
    for rl in res_lists:
        for r in rl:
            rr = r.copy(); rr.detach_parent(); rr.id = (" ", n, " "); nc.add(rr); n += 1
    md.add(nc); st.add(md); io = PDBIO(); io.set_structure(st); io.save(out)
    return n - 1

def find_native(rd, tid):
    for f in sorted(os.listdir(rd)):
        low = f.lower()
        if (tid.lower() in low and low.endswith((".pdb", ".cif"))
                and not any(m in low for m in ("boltz", "chai", "protenix", "colabfold", "tfold", "antigen", "antibody", "_fv"))):
            return os.path.join(rd, f)
    return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs_diverse")
    ap.add_argument("--out", default="zdock")
    ap.add_argument("--only", help="쉼표구분 PDB id (예: 9K6J,9Y0A)")
    a = ap.parse_args()
    only = set(x.upper() for x in a.only.split(",")) if a.only else None

    for rd in sorted(glob.glob(os.path.join(a.runs, "*"))):
        tid = os.path.basename(rd)
        if not os.path.isdir(rd) or tid.endswith("_fv"):
            continue
        if only and tid.upper() not in only:
            continue
        cjp = os.path.join(rd, "chains.json")
        if not os.path.exists(cjp):
            continue
        cj = json.load(open(cjp))
        ag = str(cj["antigen"]); ab = cj["antibody"]; ab = [ab] if isinstance(ab, str) else [str(x) for x in ab]
        nat = find_native(rd, tid)
        if not nat:
            print(f"[skip] {tid}: native 없음"); continue
        m = load(nat); ids = [c.id for c in m]
        if ag not in ids or not all(x in ids for x in ab):
            print(f"[warn] {tid}: 사슬 {ag}/{ab} 중 없음 (있는 것 {ids})"); continue

        od = os.path.join(a.out, tid); os.makedirs(od, exist_ok=True)
        n_ag = write_single([residues(m[ag])], "A", os.path.join(od, "antigen.pdb"))
        ab_lists = [residues(m[c]) for c in ab]
        n_ab = write_single(ab_lists, "B", os.path.join(od, "antibody.pdb"))

        # CDR: 각 항체 사슬 → CDR 인덱스 → 병합 오프셋으로 1..N 번호화
        cdr_nums = set(); off = 0; info = []
        for c in ab:
            rl = residues(m[c]); seq = "".join(T3.get(r.get_resname().upper(), "X") for r in rl)
            try:
                cs = cdr_index(seq)
            except Exception as e:
                cs = set(); print(f"  {tid} {c} ANARCI 실패: {e}")
            for i in cs:
                if i < len(rl):
                    cdr_nums.add(off + i + 1)     # 병합 후 1-based 번호
            info.append(f"{c}({len(rl)}res,CDR{len(cs)})"); off += len(rl)

        noncdr = [n for n in range(1, n_ab + 1) if n not in cdr_nums]
        open(os.path.join(od, "block_noncdr.txt"), "w").write("\n".join(map(str, noncdr)) + "\n")
        print(f"[ok] {tid}: 항원 {ag}({n_ag}res)->A | 항체 {'+'.join(ab)} {info}->B {n_ab}res | "
              f"CDR {len(cdr_nums)} / 비-CDR blocked {len(noncdr)} -> {od}/")

if __name__ == "__main__":
    main()
