#!/usr/bin/env python3
"""Consensus 재분석: summary.csv -> 모델별 성공 · Oracle · ipTM-pick · regret.

데이터가 있는 서버에서 돌린다(짧은 콘솔 출력만 공유). stdlib만 사용(pandas 불필요).
사용(pipeline 루트에서):
  python scripts/analyze_consensus.py --csv results/summary.csv --phase Benchmark
컬럼 가정: model, target, seed, dockq_abag, global_iptm, pair_iptm_abag, phase.
"""
import argparse, csv, collections


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="results/summary.csv")
    ap.add_argument("--phase", default="Benchmark", help="빈 문자열이면 전체")
    ap.add_argument("--thr", type=float, default=0.23)
    a = ap.parse_args()
    thr = a.thr

    rows = []
    with open(a.csv) as f:
        for r in csv.DictReader(f):
            if a.phase and r.get("phase") != a.phase:
                continue
            dq = fnum(r.get("dockq_abag"))
            if dq is None:
                continue
            rows.append({"model": r["model"], "target": r["target"], "seed": r.get("seed"),
                         "dockq": dq, "giptm": fnum(r.get("global_iptm")),
                         "pair": fnum(r.get("pair_iptm_abag"))})
    if not rows:
        print(f"[에러] {a.csv} (phase={a.phase!r}) 에 유효 행 없음")
        return

    models = sorted({r["model"] for r in rows})
    targets = sorted({r["target"] for r in rows})
    by_mt = collections.defaultdict(list)
    by_t = collections.defaultdict(list)
    for r in rows:
        by_mt[(r["model"], r["target"])].append(r)
        by_t[r["target"]].append(r)
    print(f"로드: {len(rows)}행 | 모델 {models} | 타깃 {targets}")

    # 1) 모델 x 타깃 best DockQ
    print(f"\n=== best DockQ  (성공 = >= {thr}, '*'=성공) ===")
    print("target".ljust(8) + "".join(m[:12].ljust(14) for m in models))
    msucc = {m: 0 for m in models}
    mn = {m: 0 for m in models}
    for t in targets:
        line = t.ljust(8)
        for m in models:
            rs = by_mt.get((m, t))
            if rs:
                best = max(x["dockq"] for x in rs)
                ok = best >= thr
                msucc[m] += int(ok); mn[m] += 1
                line += (f"{best:.2f}{'*' if ok else ' '}").ljust(14)
            else:
                line += "-".ljust(14)
        print(line)
    print("성공".ljust(8) + "".join(f"{msucc[m]}/{mn[m]}".ljust(14) for m in models))

    # 2) Oracle (전 모델·pose 최고)
    print("\n=== Oracle (전 모델·pose 중 최고 = 합의 천장) ===")
    orc_succ = 0; orc_dq = {}
    for t in targets:
        best = max(x["dockq"] for x in by_t[t])
        orc_dq[t] = best
        ok = best >= thr; orc_succ += int(ok)
        # 어느 모델이 그 best를 냈나
        who = max(by_t[t], key=lambda x: x["dockq"])["model"]
        print(f"  {t}: {best:.3f} ({who}) {'성공' if ok else '실패'}")
    print(f"Oracle 성공 = {orc_succ}/{len(targets)}")

    # 3) Cross-model global ipTM-pick + 평균 regret
    print("\n=== Cross-model global ipTM-pick (배포 가능한 선택 baseline) ===")
    pick_succ = 0; regrets = []
    for t in targets:
        cand = [x for x in by_t[t] if x["giptm"] is not None]
        if not cand:
            print(f"  {t}: global ipTM 없음"); continue
        p = max(cand, key=lambda x: x["giptm"])
        ok = p["dockq"] >= thr; pick_succ += int(ok)
        reg = orc_dq[t] - p["dockq"]; regrets.append(reg)
        print(f"  {t}: pick={p['model']} ipTM={p['giptm']:.2f} DockQ={p['dockq']:.3f} "
              f"{'성공' if ok else '실패'} (regret {reg:.3f})")
    mr = sum(regrets) / len(regrets) if regrets else 0.0
    print(f"ipTM-pick 성공 = {pick_succ}/{len(targets)} | 평균 regret = {mr:.3f}")

    # 4) 모델 내부 ipTM-pick (각 모델의 top-ipTM pose가 성공인가 = 랭킹 품질)
    print("\n=== 모델 내부 ipTM-pick (각 모델 자기 pose 중 top-ipTM 선택) ===")
    for m in models:
        s = n = 0
        for t in targets:
            rs = [x for x in by_mt.get((m, t), []) if x["giptm"] is not None]
            if not rs:
                continue
            n += 1
            s += int(max(rs, key=lambda x: x["giptm"])["dockq"] >= thr)
        print(f"  {m}: {s}/{n}")

    print("\n[요약] 단일최고 vs Oracle vs ipTM-pick 비교로 'consensus 천장'과 "
          "'배포 신호(ipTM)의 한계'를 읽는다.")


if __name__ == "__main__":
    main()
