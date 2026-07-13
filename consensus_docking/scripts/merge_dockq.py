#!/usr/bin/env python3
"""[min-DockQ 컨벤션 검증] 항체 H+L을 단일 receptor로 병합한 표준 Ab-Ag DockQ 재평가.
기존 eval은 DockQ(항원-H)·DockQ(항원-L)의 min을 씀 → 두 사슬 다 맞아야 성공(과엄격).
여기선 H+L을 한 사슬(B)로 묶어 항원(A) vs 항체(B) DockQ 1개를 산출 = CAPRI 표준 Ab-Ag DockQ.
목적: 9XF9·9YC6가 min 컨벤션에 가려진 진짜 성공인지 확정.

사용: conda activate boltz   # DockQ + biopython 있는 env
      python scripts/merge_dockq.py                       # 기본 9XF9,9YC6 × Fab/Fv × 5모델
      python scripts/merge_dockq.py --targets 9XF9,9KKJ   # 다른 대상
각 (타깃,포맷,모델)의 모든 pose에 대해 병합 DockQ를 구해 best를 보고. min·max(H,L)와 나란히 비교.
"""
import argparse, glob, json, os, re, subprocess, tempfile
from Bio.PDB import MMCIFParser, PDBParser, PDBIO

RUN = "runs_diverse"
POSE_GLOB = {
    "Boltz-2":       "results_boltz/**/*_model_*.cif",
    "Chai-1":        "out_chai/pred.model_idx_*.cif",
    "Protenix-base": "out_protenix/**/*_sample_*.cif",
    "AF2-Multimer":  "results_colabfold/*_unrelaxed_rank_*.pdb",
    "tFold-Ag":      "out_tfold/*seed*.pdb",
}
MODELS = list(POSE_GLOB)


def load(path):
    p = MMCIFParser(QUIET=True) if path.endswith(".cif") else PDBParser(QUIET=True)
    return p.get_structure("x", path)[0]


def ca(chain):
    return sum(1 for r in chain if "CA" in r)


def order_ab(model, ab_ids, n_heavy):
    """항체 사슬들을 heavy(먼저)·light 순으로 정렬 (native heavy CA수에 가까운 것=heavy)."""
    present = [c for c in ab_ids if c in [ch.id for ch in model]]
    if len(present) <= 1:
        return present
    present.sort(key=lambda cid: abs(ca(model[cid]) - n_heavy))
    return present  # 첫 원소가 heavy에 가장 가까움


def write_merged(model, ag_id, ab_ids, out_pdb):
    """새 PDB: chain A=antigen, chain B=항체(ab_ids 순서대로 이어붙임), resnum 1..N 재넘버링."""
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model as BModel
    from Bio.PDB.Chain import Chain as BChain
    s = Structure("m"); m = BModel(0); s.add(m)
    A = BChain("A"); B = BChain("B")
    i = 1
    for r in model[ag_id]:
        if "CA" not in r: continue
        rr = r.copy(); rr.id = (" ", i, " "); A.add(rr); i += 1
    j = 1
    for aid in ab_ids:
        if aid not in [c.id for c in model]: continue
        for r in model[aid]:
            if "CA" not in r: continue
            rr = r.copy(); rr.id = (" ", j, " "); B.add(rr); j += 1
    m.add(A); m.add(B)
    io = PDBIO(); io.set_structure(s); io.save(out_pdb)


_DBG = {"n": 0}


def dockq(model_pdb, native_pdb):
    try:
        r = subprocess.run(["DockQ", model_pdb, native_pdb],
                           capture_output=True, text=True, timeout=180)
    except Exception as e:
        if _DBG["n"] < 1:
            print("[DEBUG] DockQ 실행 예외:", e); _DBG["n"] += 1
        return None
    out = (r.stdout or "") + (r.stderr or "")
    if _DBG["n"] < 1:                       # 첫 호출 원문 출력 (파서 확인용)
        print("[DEBUG] ===== DockQ raw output =====\n" + out[:800] + "\n[DEBUG] ====="); _DBG["n"] += 1
    vals = re.findall(r"\bDockQ\b[:\s]+([0-9]*\.[0-9]+)", out)
    return max(float(v) for v in vals) if vals else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="9XF9,9YC6")
    a = ap.parse_args()
    targets = [t.strip().upper() for t in a.targets.split(",")]

    print(f"{'target':6} {'fmt':4} {'model':14} {'merged':>8} {'poses':>6}")
    for T in targets:
        for fmt in ("Fab", "Fv"):
            low = T.lower() + ("_fv" if fmt == "Fv" else "")
            rd = os.path.join(RUN, low)
            cj = os.path.join(rd, "chains.json")
            if not os.path.exists(cj):
                continue
            d = json.load(open(cj))
            ag = str(d["antigen"])
            ab = d["antibody"]; ab = [ab] if isinstance(ab, str) else [str(x) for x in ab]
            nat = next((os.path.join(rd, f) for f in os.listdir(rd)
                        if f.lower().endswith((".pdb", ".cif")) and T.lower() in f.lower()
                        and "boltz" not in f and "chai" not in f and "protenix" not in f
                        and "colabfold" not in f and "tfold" not in f), None)
            if not nat:
                nat = next((os.path.join(rd, f) for f in os.listdir(rd) if f.lower().endswith((".pdb", ".cif"))), None)
            if not nat:
                continue
            natm = load(nat)
            n_heavy = len(d["chains"][0]["seq"]) if d.get("chains") else 200
            with tempfile.TemporaryDirectory() as td:
                nat_m = os.path.join(td, "nat.pdb")
                write_merged(natm, ag, ab, nat_m)   # native: 선언된 사슬(heavy first)
                for model in MODELS:
                    poses = sorted(glob.glob(os.path.join(rd, POSE_GLOB[model]), recursive=True))
                    best = None
                    for pi, pose in enumerate(poses):
                        try:
                            pm = load(pose)
                            n_ag = sum(1 for c in d["chains"] if c["id"] == ag for _ in c["seq"])
                            pag = min(pm, key=lambda ch: abs(ca(ch) - n_ag)).id
                            pab = order_ab(pm, [c.id for c in pm if c.id != pag], n_heavy)
                            mp = os.path.join(td, f"m{pi}.pdb")
                            write_merged(pm, pag, pab, mp)
                            q = dockq(mp, nat_m)
                            if q is not None and (best is None or q > best):
                                best = q
                        except Exception:
                            continue
                    bs = f"{best:.3f}" if best is not None else "  -  "
                    print(f"{T:6} {fmt:4} {model:14} {bs:>8} {len(poses):>6}")
        print()


if __name__ == "__main__":
    main()
