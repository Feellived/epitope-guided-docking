#!/usr/bin/env python3
"""MSA-depth 통합 리포트(5모델): build_msafree_summary CSV → 모델별 depth 스윕(RBM%+recall).
HADDOCK은 report_haddock.py로 분리. 사용: python scripts/report_msa_depth.py [csv]"""
import sys, csv, glob
from collections import defaultdict

cands = [sys.argv[1]] if len(sys.argv) > 1 else sorted(glob.glob("*.csv") + glob.glob("results/*.csv"))
path = None
for c in cands:
    try:
        hdr = next(csv.reader(open(c)))
        if "model" in hdr and "condition" in hdr:
            path = c
    except Exception:
        pass
if not path:
    sys.exit("summary CSV 못 찾음 — build_msafree_summary.py 먼저 실행")
print(f"[리포트] {path}\n")
rows = list(csv.DictReader(open(path)))

DORD = {"full": 999, "d143": 143, "d8": 8, "d1": 1}
def dnum(c): return DORD.get(str(c), -1)
def pct(x):
    try:
        v = float(x); return f"{v*100:.0f}" if v <= 1.0 else f"{v:.0f}"
    except Exception:
        return (str(x) or "-")
def r2(x):
    try:
        return f"{float(x):.2f}"
    except Exception:
        return (str(x) or "-")

data = defaultdict(lambda: defaultdict(dict)); eclass = {}
for r in rows:
    if r.get("experiment") not in ("cofolder", "msa_depth"):
        continue
    m, t, c = r.get("model", ""), r.get("target", "").upper(), r.get("condition", "")
    data[m][t][c] = (r.get("rbm_frac", ""), r.get("epitope_recall", ""), r.get("n_poses", ""))
    eclass[t] = r.get("epitope_class", "")

MODEL_ORD = ["Boltz-2", "Protenix-base", "tFold", "Chai-1", "AF2-Multimer"]
TARGET_ORD = ["8SIT", "9ML9", "8SIQ", "9ML8", "8P5M", "9ZDU", "8SDF", "8XSI", "8SIS", "9SBB"]
models = [m for m in MODEL_ORD if m in data] + [m for m in data if m not in MODEL_ORD]

for m in models:
    conds = sorted({c for t in data[m] for c in data[m][t]}, key=dnum, reverse=True)
    print(f"===== {m}   (셀 = RBM% / recall) =====")
    print(f"{'target':7}{'class':11}" + "".join(f"{str(c):>13}" for c in conds))
    tgts = [t for t in TARGET_ORD if t in data[m]] + [t for t in sorted(data[m]) if t not in TARGET_ORD]
    for t in tgts:
        cells = []
        for c in conds:
            if c in data[m][t]:
                rbm, rec, n = data[m][t][c]; cells.append(f"{pct(rbm):>3}/{r2(rec):>4}")
            else:
                cells.append(f"{'-':>8}")
        print(f"{t:7}{eclass.get(t,''):11}" + "".join(f"{x:>13}" for x in cells))
    print()
print("주: RBM%↓ + recall↑ = 항원 MSA 감소로 off-hotspot rescue. full(양사슬 MSA)/d143(항체 single)/d8/d1.")
