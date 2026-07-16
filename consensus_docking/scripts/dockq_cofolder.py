#!/usr/bin/env python3
"""co-folder(+MSA-depth) 전체 포즈 DockQ: (target,experiment,model,condition)별 best DockQ. 항체 H+L 병합.
   recall만 있던 build_msafree_summary의 빈 dockq 열을 채움. HADDOCK은 haddock_dockq로 별도.
   사용(boltz env, DockQ 설치필요): python scripts/dockq_cofolder.py [--out results/dockq_cofolder.csv]"""
import argparse, csv, glob, json, os, re, subprocess, tempfile
from Bio.PDB import PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model as BModel
from Bio.PDB.Chain import Chain as BChain
import build_msafree_summary as B

def write_merged(model, ag_id, ab_ids, out_pdb):
    s = Structure("m"); m = BModel(0); s.add(m); A = BChain("A"); Bc = BChain("B")
    i = 1
    for r in model[ag_id]:
        if "CA" not in r: continue
        rr = r.copy(); rr.id = (" ", i, " "); A.add(rr); i += 1
    j = 1
    for aid in ab_ids:
        if aid not in [c.id for c in model]: continue
        for r in model[aid]:
            if "CA" not in r: continue
            rr = r.copy(); rr.id = (" ", j, " "); Bc.add(rr); j += 1
    m.add(A); m.add(Bc); io = PDBIO(); io.set_structure(s); io.save(out_pdb)

def dockq(mp, nat):
    try:
        r = subprocess.run(["DockQ", mp, nat], capture_output=True, text=True, timeout=180)
    except Exception:
        return None
    out = (r.stdout or "") + (r.stderr or "")
    vals = re.findall(r"\bDockQ\b[:\s]+([0-9]*\.[0-9]+)", out)
    return max(float(v) for v in vals) if vals else None

def native_merged(t, td):
    rd = os.path.join("runs_rbd", t.lower()); d = json.load(open(os.path.join(rd, "chains.json")))
    ag = str(d["antigen"]); ab = d["antibody"]; ab = [ab] if isinstance(ab, str) else [str(x) for x in ab]
    nat = next((os.path.join(rd, f) for f in os.listdir(rd) if t.lower() in f.lower() and f.endswith((".pdb", ".cif"))
                and not any(x in f.lower() for x in ("boltz", "chai", "protenix", "colabfold", "tfold"))), None)
    if not nat: return None
    nm = B.load(nat); ids = [c.id for c in nm]; sm = {c["id"]: c["seq"] for c in d["chains"]}
    ag_id = ag if ag in ids else B.best(nm, sm[ag])[0]
    ab_ids = [i for i in ab if i in ids]
    if not ab_ids: return None
    out = os.path.join(td, "nat.pdb"); write_merged(nm, ag_id, ab_ids, out); return out

def id_chains(pm, ni):
    pid, _, _, _, _ = B.best(pm, ni["ag_seq"]); used = {pid}; abids = []
    for s in ni["ab_seqs"]:
        r = B.best(pm, s, exclude=used)
        if r[0]: used.add(r[0]); abids.append(r[0])
    return pid, abids

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default="results/dockq_cofolder.csv"); a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    ni_all = {t: B.native_info(t) for t in B.TARGETS}
    rows = []
    for exp, model, cond, rundir, gl in B.SOURCES:
        if exp == "haddock": continue
        for t in B.TARGETS:
            ni = ni_all.get(t)
            if not ni: continue
            poses = sorted(glob.glob(os.path.join(rundir, t.lower(), gl), recursive=True))
            if not poses: continue
            with tempfile.TemporaryDirectory() as td:
                natm = native_merged(t, td)
                if natm is None: continue
                best = None
                for pi, pose in enumerate(poses):
                    try:
                        pm = B.load(pose); pid, abids = id_chains(pm, ni)
                        if not abids: continue
                        mp = os.path.join(td, f"m{pi}.pdb"); write_merged(pm, pid, abids, mp)
                        q = dockq(mp, natm)
                        if q is not None and (best is None or q > best): best = q
                    except Exception:
                        continue
                rows.append({"target": t, "epitope_class": B.EPICLASS.get(t, ""), "experiment": exp, "model": model,
                             "condition": cond, "best_dockq": f"{best:.3f}" if best is not None else "", "n_poses": len(poses)})
                print(f"{t:6} {model:14} {cond:7} DockQ {best if best is not None else '-'}  ({len(poses)}p)")
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["target", "epitope_class", "experiment", "model", "condition", "best_dockq", "n_poses"])
        w.writeheader(); w.writerows(rows)
    print(f"\n{len(rows)}행 → {a.out}")

if __name__ == "__main__": main()
