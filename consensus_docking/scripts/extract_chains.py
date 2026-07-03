#!/usr/bin/env python3
"""[Layer 1] 사슬 서열 추출 — 모델 무관, 복합체당 한 번.
  native PDB 다운(pdb→cif) + keep_chains 서열 추출 → runs/<id>/chains.json (+ <id>.fasta).
  chains.json = 이 실험의 "정답 서열·역할" 단일 기록 (Layer 2가 이걸 읽음).

사용: python extract_chains.py --targets targets.yaml --out-dir runs [--only 10ZO,9K6J]
필요: biopython, pyyaml
"""
import argparse, json, os, urllib.request, urllib.error


def fetch_native(pdb_id, run_dir):
    for ext in ("pdb", "cif"):
        dst = os.path.join(run_dir, f"{pdb_id}.{ext}")
        if os.path.exists(dst):
            return dst
        try:
            urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id}.{ext}", dst)
            return dst
        except urllib.error.HTTPError:
            continue
    return None


def chain_seqs(path):
    from Bio.PDB import PDBParser, MMCIFParser, PPBuilder
    parser = MMCIFParser(QUIET=True) if path.endswith(".cif") else PDBParser(QUIET=True)
    st = parser.get_structure("x", path)
    ppb = PPBuilder()
    out = {}
    for ch in st[0]:
        seq = "".join(str(pp.get_sequence()) for pp in ppb.build_peptides(ch))
        if seq:
            out[ch.id] = seq
    return out


def load_targets(path):
    import yaml
    with open(path) as f:
        data = yaml.safe_load(f)
    items = []
    for lst in (data or {}).values():
        items.extend(lst)
    return items


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="targets.yaml")
    ap.add_argument("--out-dir", default="runs")
    ap.add_argument("--only", help="쉼표구분 PDB id (예: 10ZO,9K6J)")
    a = ap.parse_args()
    only = set(x.upper() for x in a.only.split(",")) if a.only else None

    for t in load_targets(a.targets):
        tid = t["id"]
        if only and tid.upper() not in only:
            continue
        keep = [str(c) for c in t.get("keep_chains", [])]
        antigen = str(t["antigen"]["chain"])
        antibody = [str(v) for v in t["antibody"].values()]
        rd = os.path.join(a.out_dir, tid.lower())
        os.makedirs(rd, exist_ok=True)
        native = fetch_native(tid, rd)
        if not native:
            print(f"[skip] {tid}: native 다운 실패")
            continue
        allseqs = chain_seqs(native)
        miss = [c for c in keep if c not in allseqs]
        if miss:
            print(f"[warn] {tid}: 서열 없는 사슬 {miss} (있는 것: {list(allseqs)})")
        chains = [{"id": c, "role": ("antigen" if c == antigen else "antibody"), "seq": allseqs[c]}
                  for c in keep if c in allseqs]
        if not chains:
            print(f"[skip] {tid}: 쓸 사슬 없음")
            continue
        data = {"pdb_id": tid, "antigen": antigen, "antibody": antibody,
                "order": [c["id"] for c in chains], "chains": chains}
        with open(os.path.join(rd, "chains.json"), "w") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        with open(os.path.join(rd, f"{tid.lower()}.fasta"), "w") as f:
            for c in chains:
                f.write(f">{c['id']}|{c['role']}\n{c['seq']}\n")
        print(f"[ok] {tid}: {[c['id'] for c in chains]} (항원 {antigen}) -> {rd}/chains.json")


if __name__ == "__main__":
    main()
