#!/usr/bin/env python3
"""HADDOCK DockQ 스코어카드: results/haddock_dockq_*.csv(emref·rigidbody) 결합 표.
사용: python scripts/analyze_scorecard_haddock.py"""
import csv, glob, os
from collections import defaultdict
data = defaultdict(dict); cls = {}
files = glob.glob("results/haddock_dockq_*.csv")
if not files:
    raise SystemExit("results/haddock_dockq_*.csv 없음 — eval_haddock.py 먼저 실행")
for f in files:
    stage = os.path.basename(f).replace("haddock_dockq_", "").replace(".csv", "")
    for r in csv.DictReader(open(f)):
        data[r["target"]][stage] = (r.get("best_dockq", ""), r.get("best_recall", ""), r.get("n_hit", ""), r.get("n_scored", ""))
        cls[r["target"]] = r.get("epitope_class", "")
stages = sorted({s for t in data for s in data[t]}, key=lambda s: (s != "emref", s))
ORDER = ["8XSI", "9SBB", "8SIS", "9ML9", "8P5M", "9ML8", "8SDF", "8SIQ", "8SIT", "9ZDU"]
print(f"{'target':7}{'class':9}" + "".join(f"{'DockQ/'+s:>22}" for s in stages) + "   판정 (셀=DockQ/recall(hit/n))")
print("-" * (16 + 22 * len(stages) + 20))
for t in [x for x in ORDER if x in data] + [x for x in data if x not in ORDER]:
    cells = []
    for s in stages:
        if s in data[t]:
            dq, rec, nh, ns = data[t][s]; cells.append(f"{dq or '-'}/{rec or '-'}({nh}/{ns})")
        else:
            cells.append("-")
    best = max([float(data[t][s][0]) for s in data[t] if data[t][s][0]] + [0])
    v = "✅ near-native" if best >= 0.49 else ("🔶 acceptable" if best >= 0.23 else "❌ recall착시/hard")
    print(f"{t:7}{cls.get(t,''):9}" + "".join(f"{c:>22}" for c in cells) + f"   {v}")
print("\n주: emref=정제 생존자, rigidbody=raw 생성물. rigidbody≥0.49 & emref낮음=점수컷이 죽임(재랭크 회수).")
