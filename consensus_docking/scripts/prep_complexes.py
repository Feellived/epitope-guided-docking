#!/usr/bin/env python3
"""[Study B/E1-Foldseek] PDB 복합체 → (항원 1사슬 + 결합 항체사슬)만 추출·표준화 → Foldseek 입력.
각 복합체를 chain A=항원(RBD), B/C..=binder 로 재라벨해 out-dir에 저장. Foldseek-Multimer 클러스터용.
"""
import argparse, io, json, os, urllib.request
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model as BModel
from Bio.PDB.Chain import Chain as BChain
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_3to1

DEFAULT_REF = ("NITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRG"
               "DEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGST"
               "PCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNF")
T3TO1 = {k.upper(): v for k, v in protein_letters_3to1.items()}
_al = PairwiseAligner(); _al.mode = "global"; _al.match_score = 1; _al.mismatch_score = -1
_al.open_gap_score = -3; _al.extend_gap_score = -0.5


def fetch(pdb):
    txt = urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb}.cif", timeout=60).read().decode("utf-8", "replace")
    return MMCIFParser(QUIET=True).get_structure(pdb, io.StringIO(txt))[0]


def seq_of(chain):
    return "".join(T3TO1.get(r.get_resname().upper(), "") for r in chain if "CA" in r)


def ident(a, b):
    if not a or not b: return -1.0
    return _al.score(a, b) / max(len(a), len(b), 1)


def contacts(agc, others, thr=5.0):
    ag_at = np.array([at.get_coord() for r in agc if "CA" in r for at in r])
    ot_at = np.array([at.get_coord() for ch in others for r in ch if "CA" in r for at in r])
    if len(ag_at) == 0 or len(ot_at) == 0: return 0
    d = np.sqrt(((ag_at[:, None, :] - ot_at[None, :, :]) ** 2).sum(-1))
    return int((d.min(1) <= thr).sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb-list"); ap.add_argument("--max", type=int, default=400)
    ap.add_argument("--ref-fasta"); ap.add_argument("--out-dir", default="complexes")
    a = ap.parse_args()
    ref = DEFAULT_REF
    if a.ref_fasta: ref = "".join(l.strip() for l in open(a.ref_fasta) if not l.startswith(">"))

    if a.pdb_list:
        pdbs = [x.strip().upper() for x in open(a.pdb_list) if x.strip()]
    else:
        q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
                {"type": "terminal", "service": "sequence", "parameters": {"evalue_cutoff": 1, "identity_cutoff": 0.9, "sequence_type": "protein", "value": ref}},
                {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_assembly_info.polymer_entity_instance_count", "operator": "greater_or_equal", "value": 2}},
                {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_entry_info.resolution_combined", "operator": "less_or_equal", "value": 3.0}}]},
             "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": a.max}}}
        req = urllib.request.Request("https://search.rcsb.org/rcsbsearch/v2/query", data=json.dumps(q).encode(), headers={"Content-Type": "application/json"})
        pdbs = [r["identifier"] for r in json.load(urllib.request.urlopen(req, timeout=60)).get("result_set", [])]
    pdbs = pdbs[:a.max]
    os.makedirs(a.out_dir, exist_ok=True)
    print(f"대상 {len(pdbs)}개 → {a.out_dir}/")
    made = skip = 0
    for i, pdb in enumerate(pdbs, 1):
        try:
            m = fetch(pdb)
            chains = {ch.id: ch for ch in m if len([r for r in ch if "CA" in r]) >= 20}
            seqs = {k: seq_of(v) for k, v in chains.items()}
            ag = max(seqs, key=lambda c: ident(seqs[c], ref))
            if ident(seqs[ag], ref) < 0.8: skip += 1; continue
            binders = [c for c in chains if c != ag and len([r for r in chains[c] if "CA" in r]) >= 30 and ident(seqs[c], ref) < 0.5]
            binders = [c for c in binders if contacts(chains[ag], [chains[c]]) >= 3]   # 항원에 실제 접촉
            if not binders: skip += 1; continue
            # 표준 재라벨: A=항원, B,C..=binder
            s = Structure("x"); mm = BModel(0); s.add(mm)
            newlab = ["A", "B", "C", "D", "E", "F"]
            A = BChain("A")
            for r in chains[ag]:
                if "CA" in r: A.add(r.copy())
            mm.add(A)
            for bi, c in enumerate(binders[:5]):
                B = BChain(newlab[bi + 1])
                for r in chains[c]:
                    if "CA" in r: B.add(r.copy())
                mm.add(B)
            io_ = PDBIO(); io_.set_structure(s); io_.save(os.path.join(a.out_dir, f"{pdb}.pdb"))
            made += 1
        except Exception:
            skip += 1
        if i % 25 == 0: print(f"  {i}/{len(pdbs)} (made {made}, skip {skip})")
    print(f"완료: {made}개 저장, {skip} 스킵 → {a.out_dir}/")


if __name__ == "__main__":
    main()
