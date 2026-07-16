#!/usr/bin/env python3
"""패치-강제 HADDOCK DockQ: 타깃별 전 패치(patch_*/run/<stage>)의 포즈에서 best DockQ.
   ab-initio(haddock_dockq) 대비 '패치강제가 seletop 편향을 막았나' 확인.
   사용(boltz env, DockQ 설치): python scripts/haddock_patch_dockq.py [--stage 4_emref] [--topk 25]"""
import argparse, glob, os, csv, tempfile
import build_msafree_summary as B
from haddock_dockq import load_gz, write_merged, dockq, id_and_recall, native_merged

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="4_emref"); ap.add_argument("--topk", type=int, default=25)
    ap.add_argument("--targets", default=" ".join(B.TARGETS)); a = ap.parse_args()
    print(f"[패치-강제, stage={a.stage}, DockQ on top{a.topk} by recall]")
    print(f"{'target':7}{'class':9}{'best_DockQ':11}{'recall':8}{'patch':7}{'#>=.23'}")
    print("-" * 50)
    rows = []
    for t in a.targets.split():
        t = t.upper(); ni = B.native_info(t)
        if ni is None: print(f"{t:7} native 없음"); continue
        poses = sorted(glob.glob(os.path.join("haddock", t.lower(), "patch_*", "run", a.stage, "*.pdb.gz")))
        if not poses:
            print(f"{t:7}{B.EPICLASS.get(t,''):9} (패치 포즈 없음 — 아직 안 돌았거나 stage 확인)"); continue
        scored = []
        for p in poses:
            try:
                rec, pid, abids = id_and_recall(load_gz(p), ni)
                pn = p.split("patch_")[1].split("/")[0]
                scored.append((rec, p, pid, abids, pn))
            except Exception: continue
        scored.sort(reverse=True)
        with tempfile.TemporaryDirectory() as td:
            natm = native_merged(t, td)
            if natm is None: print(f"{t:7} native merge 실패"); continue
            bq = brec = bpn = None; nhit = nd = 0
            for rec, p, pid, abids, pn in scored[:a.topk]:
                try:
                    mp = os.path.join(td, "m.pdb"); write_merged(load_gz(p), pid, abids, mp)
                    q = dockq(mp, natm); nd += 1
                    if q is not None:
                        if q >= 0.23: nhit += 1
                        if bq is None or q > bq: bq, brec, bpn = q, rec, pn
                except Exception: continue
            bs = f"{bq:.3f}" if bq is not None else "-"
            rs = f"{brec:.2f}" if brec is not None else "-"
            print(f"{t:7}{B.EPICLASS.get(t,''):9}{bs:>10} {rs:>7}  p{str(bpn):>4}  {nhit}/{nd}")
            rows.append({"target": t, "epitope_class": B.EPICLASS.get(t, ""), "stage": "patch_" + a.stage,
                         "best_dockq": (f"{bq:.3f}" if bq is not None else ""), "best_recall": rs, "n_hit": nhit, "n_scored": nd})
    os.makedirs("results", exist_ok=True)
    with open("results/haddock_dockq_patch.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["target", "epitope_class", "stage", "best_dockq", "best_recall", "n_hit", "n_scored"])
        w.writeheader(); w.writerows(rows)
    print("\n→ results/haddock_dockq_patch.csv")
    print("주: 패치 best_DockQ vs ab-initio(rigidbody) 비교 → 패치강제가 seletop 편향을 막고 near-native를 살렸나.")

if __name__ == "__main__": main()
