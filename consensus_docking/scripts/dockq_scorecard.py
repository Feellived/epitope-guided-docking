#!/usr/bin/env python3
"""DockQ 통합 스코어카드 (방법론별 분리 저장 뷰). 파이프라인=HADDOCK(전반부 후보)→co-folder(정교 refine), 병렬 union 아님.
   ① co-folder(dockq_cofolder.csv, 모델×조건) ② HADDOCK(haddock_dockq_*.csv) 별도 표시.
사용: python scripts/dockq_scorecard.py"""
import csv, glob, os
from collections import defaultdict

TORDER = ["8SIT", "9ML9", "9ZDU", "8SIQ", "8SIS", "9ML8", "8SDF", "8XSI", "8P5M", "9SBB"]
def dq(x):
    try: return float(x)
    except Exception: return None
def mark(v):
    if v is None: return "-"
    return f"{v:.2f}" + ("*" if v >= 0.49 else ("." if v >= 0.23 else " "))

cof = defaultdict(dict); cls = {}
p = "results/dockq_cofolder.csv"
if os.path.exists(p):
    for r in csv.DictReader(open(p)):
        cof[r["target"]][(r["model"], r["condition"])] = dq(r["best_dockq"]); cls[r["target"]] = r["epitope_class"]
conds = []
for m in ["Protenix-base", "Boltz-2", "Chai-1", "AF2-Multimer", "tFold", "tFold-Ag"]:
    for c in ["server", "d143", "d8", "d1"]:
        if any((m, c) in cof[t] for t in cof): conds.append((m, c))

print("=== ① co-folder DockQ (모델×조건 분리 | *=medium+ .=accept) ===")
print("tgt    class     " + " ".join(f"{m.split('-')[0][:4]}/{c[:4]:<4}" for m, c in conds))
for t in [x for x in TORDER if x in cof]:
    print(f"{t:6} {cls.get(t,''):9} " + " ".join(f"{mark(cof[t].get((m,c))):>9}" for m, c in conds))

had = defaultdict(dict)
for f in glob.glob("results/haddock_dockq_*.csv"):
    stage = os.path.basename(f).replace("haddock_dockq_", "").replace(".csv", "")
    for r in csv.DictReader(open(f)):
        had[r["target"]][stage] = dq(r["best_dockq"]); cls.setdefault(r["target"], r.get("epitope_class", ""))
if had:
    print("\n=== ② HADDOCK DockQ (전반부 후보 생성기, 별도) ===")
    stages = sorted({s for t in had for s in had[t]})
    print("tgt    class     " + " ".join(f"{s:>11}" for s in stages))
    for t in [x for x in TORDER if x in had]:
        print(f"{t:6} {cls.get(t,''):9} " + " ".join(f"{mark(had[t].get(s)):>11}" for s in stages))
else:
    print("\n(HADDOCK csv 없음 — haddock_dockq.py 실행 후 나타남)")
print("\n주: 파이프라인=HADDOCK(후보)→co-folder(refine). co-folder server=blind(편향), depth=MSA감소 bias회피(blind).")
print("    진짜 파이프라인(HADDOCK 후보→co-folder guided refine)은 아직 미실험 = 다음 단계.")
