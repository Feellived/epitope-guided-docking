#!/usr/bin/env python3
"""targets.yaml의 모든 복합체를 한 번에 프렙:
  native 다운로드(pdb→cif fallback) + keep_chains 서열 추출(1:1 한 쌍) +
  모델별 입력 생성: Boltz YAML / Chai FASTA / Protenix JSON + chains.json(역할·순서).

사용:
  python prep_all.py --targets targets.yaml --out-dir runs           # 전부
  python prep_all.py --targets targets.yaml --out-dir runs --only 10ZO,9K6J
필요: biopython, pyyaml.  각 복합체는 runs/<id>/ 폴더로.
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
    for chain in st[0]:
        seq = "".join(str(pp.get_sequence()) for pp in ppb.build_peptides(chain))
        if seq:
            out[chain.id] = seq
    return out


def write_boltz(path, chains):
    with open(path, "w") as f:
        f.write("version: 1\nsequences:\n")
        for cid, seq in chains:
            f.write(f'  - protein: {{ id: {cid}, sequence: "{seq}" }}\n')


def write_chai(path, chains):
    with open(path, "w") as f:
        for cid, seq in chains:
            f.write(f">protein|name={cid}\n{seq}\n")


def write_protenix(path, chains, name):
    job = {"name": name, "sequences": [{"proteinChain": {"sequence": s, "count": 1}} for _, s in chains]}
    with open(path, "w") as f:
        json.dump([job], f, indent=2)


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
        rd = os.path.join(a.out_dir, tid.lower())
        os.makedirs(rd, exist_ok=True)
        native = fetch_native(tid, rd)
        if not native:
            print(f"[skip] {tid}: native 다운로드 실패")
            continue
        allseqs = chain_seqs(native)
        chains = [(c, allseqs[c]) for c in keep if c in allseqs]
        miss = [c for c in keep if c not in allseqs]
        if miss:
            print(f"[warn] {tid}: keep_chains 중 서열 없음 {miss} (있는 사슬: {list(allseqs)})")
        if not chains:
            print(f"[skip] {tid}: 쓸 사슬 없음")
            continue
        low = tid.lower()
        write_boltz(os.path.join(rd, f"boltz_{low}.yaml"), chains)
        write_chai(os.path.join(rd, f"chai_{low}.fasta"), chains)
        write_protenix(os.path.join(rd, f"protenix_{low}.json"), chains, low)
        roles = {"antigen": str(t["antigen"]["chain"]),
                 "antibody": [str(v) for v in t["antibody"].values()],
                 "order": [c for c, _ in chains]}
        with open(os.path.join(rd, "chains.json"), "w") as f:
            json.dump(roles, f, indent=2)
        print(f"[ok] {tid}: chains {[c for c, _ in chains]} -> {rd}")


if __name__ == "__main__":
    main()
