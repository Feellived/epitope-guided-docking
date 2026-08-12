#!/usr/bin/env python3
"""[선택기 검증] 임의의 자세 묶음을 AbEpiScore-1.0 으로 채점하고, ipTM·천장과 맞대어
   "정답을 모르는 상태에서 좋은 자세를 고를 수 있는가"를 정량화한다.

왜 필요한가 — 데모 표의 최고값은 채점할 때 정답을 보고 뽑은 **천장**이다. 배포에서 쓰려면
정답 없이 그 자세를 지목해야 하는데, 지금까지 후보였던 ipTM 은 타깃 내 Spearman 중앙값이
+0.04 로 동전 던지기와 구별되지 않는다(pose_features.csv 5,760 자세). AbEpiScore-1.0 은
**모델링된 계면의 정확도**를 채점하므로(항체 특이적, Science Advances 2025) 같은 자리를
대신할 수 있는지 재야 한다. **이미 정답을 아는 자세를 쓰므로 새 예측이 필요 없다.**

두 모드
  ① score  구조 파일들을 채점            → <out>.csv  (path,target,run,pose,score)
  ② eval   ①의 결과 + 정답 CSV 를 붙여   → 타깃 내 상관 · 천장 대 선택 비교

핵심 — 세 값을 나란히 낸다: **천장**(정답 보고 최고) · **AbEpiScore 선택** · **ipTM 선택**.
최고값만 보고하는 관행을 버리고 "고른 것"의 성적을 보고한다.

사용 (pipeline/ 에서. AbEpiTope venv 는 ~/tools/AbEpiTope-1.0/.venv)
  # ① 8ulr_HL 의 조성 실행 자세 전부 채점
  python scripts/eval_abepiscore.py score \\
      --glob '/mnt/data/admuser/msadepth/compreps/seedrep_cand/protenix/8ulr_HL/**/*.cif' \\
      --out results/abepiscore_8ulr_HL.csv

  # ② 정답과 맞대기 (compreps CSV 또는 pose_features.csv)
  python scripts/eval_abepiscore.py eval \\
      --scores results/abepiscore_8ulr_HL.csv \\
      --truth ~/projects/bk21-msa-depth-bias/pipeline/results/compreps_8ulr_HL.csv

  # 데모 자세에 적용 (자세 5개 중 하나를 고르는 층)
  python scripts/eval_abepiscore.py score --glob 'runs_sites_guided/8k3k_d/**/*.cif' \\
      --out results/abepiscore_demo_8k3k_D.csv
"""
import argparse, csv, glob, os, re, shutil, statistics as st, subprocess, sys, tempfile
from collections import defaultdict

DRIVER = '''
import sys
from pathlib import Path
from abepitope.main import StructureData, EvalAbAgs
inp, enc, tmp, out = (Path(x) for x in sys.argv[1:5])
data = StructureData()
data.encode_proteins(inp, enc, tmp)
EvalAbAgs(data).predict(out)
print("ABEPITOPE_OK")
'''

# 경로에서 (타깃, 실행, 자세)를 뽑는다. compreps · 데모 두 배치를 모두 안다.
PATS = [
    # /compreps/seedrep_cand/<model>/<target>/<depth>/<run>/results/**/<pose>.cif
    re.compile(r"/(?P<model>boltz|protenix|chai)/(?P<target>[^/]+)/(?P<depth>d[^/]+)/"
               r"(?P<run>seed[^/]+)/.*?(?P<pose>[^/]+\.(?:cif|pdb))$"),
    # runs_sites_guided/<target>/<arm>/**/<pose>.cif
    re.compile(r"/(?P<target>[^/]+)/(?P<run>blind|noconstraint|fullmsa|sizematch|ours|all)/"
               r".*?(?P<pose>[^/]+\.(?:cif|pdb))$"),
]


def parse_key(p):
    for rx in PATS:
        m = rx.search(p)
        if m:
            d = m.groupdict()
            return d.get("target", "?"), d.get("run", "?"), d["pose"]
    return "?", "?", os.path.basename(p)


def spearman(xs, ys):
    if len(xs) < 3:
        return None
    def rank(v):
        s = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v)
        i = 0
        while i < len(s):                      # 동순위 평균
            j = i
            while j + 1 < len(s) and v[s[j + 1]] == v[s[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[s[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    mx, my = st.mean(rx), st.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** .5
    return num / den if den else None


# ─────────────────────────── ① score ───────────────────────────
def cmd_score(a):
    files = []
    for g in a.glob:
        files += glob.glob(g, recursive=True)
    files = sorted(set(f for f in files if f.endswith((".cif", ".pdb"))))
    if not files:
        sys.exit(f"!! 구조 파일이 없다: {a.glob}")
    if os.path.exists(a.out) and not a.force:
        sys.exit(f"!! {a.out} 이 이미 있다. 덮어쓰려면 --force (원자료를 잃지 않도록 기본 거부)")
    # 실행 방식 결정 — venv 가 있으면 그걸 쓰고, 없으면 현재 파이썬 + PYTHONPATH 로 돈다.
    act = os.path.join(a.repo, ".venv", "bin", "activate")
    env = os.environ.copy()
    env["PYTHONPATH"] = a.repo + os.pathsep + env.get("PYTHONPATH", "")
    use_venv = os.path.exists(act) and not a.no_venv
    if not use_venv:
        chk = subprocess.run([a.python, "-c", "import abepitope; print('OK')"],
                             capture_output=True, text=True, env=env)
        if "OK" not in chk.stdout:
            print(chk.stderr[-800:])
            sys.exit(f"!! abepitope 를 import 하지 못했다.\n"
                     f"   현재 파이썬: {a.python}\n"
                     f"   찾은 경로  : {a.repo}\n"
                     f"   해결 ① 패키지가 있는 conda env 에서 실행 (예: conda activate boltz)\n"
                     f"   해결 ② venv 설치: bash scripts/run_setup_tools.sh --only abepitope --apply")
        print(f"■ venv 없이 현재 파이썬으로 실행한다 ({a.python})")
    print(f"■ 자세 {len(files)}개를 AbEpiScore-1.0 으로 채점 (묶음 {a.batch}개씩)")

    rows = []
    for b0 in range(0, len(files), a.batch):
        chunk = files[b0:b0 + a.batch]
        work = tempfile.mkdtemp(prefix="abepi_")
        inp = os.path.join(work, "in"); os.makedirs(inp)
        alias = {}
        for i, p in enumerate(chunk):
            nm = f"p{b0 + i:05d}{os.path.splitext(p)[1]}"
            alias[f"p{b0 + i:05d}"] = p
            shutil.copy(p, os.path.join(inp, nm))
        drv = os.path.join(work, "driver.py"); open(drv, "w").write(DRIVER)
        enc, tmp, out = (os.path.join(work, x) for x in ("enc", "tmp", "out"))
        if use_venv:
            r = subprocess.run(["bash", "-c",
                                f'source "{act}" && python "{drv}" "{inp}" "{enc}" "{tmp}" "{out}"'],
                               capture_output=True, text=True)
        else:
            r = subprocess.run([a.python, drv, inp, enc, tmp, out],
                               capture_output=True, text=True, env=env)
        if "ABEPITOPE_OK" not in r.stdout:
            print(r.stdout[-1200:]); print(r.stderr[-1200:])
            sys.exit(f"!! AbEpiTope 실행 실패 (묶음 {b0}~). 반쪽 결과를 남기지 않고 멈춘다.\n"
                     f"   흔한 원인 = GPU 메모리. --batch 를 줄여 볼 것.")
        res = next(iter(glob.glob(os.path.join(out, "**", "output.csv"), recursive=True)), None)
        if not res:
            sys.exit(f"!! output.csv 를 못 찾았다: {out}")
        got = 0
        for r0 in csv.DictReader(open(res)):
            keys = {k.lower(): k for k in r0}
            kf = next((keys[k] for k in keys if "file" in k), None)
            ks = next((keys[k] for k in keys if "abepiscore" in k), None)
            if not (kf and ks):
                sys.exit(f"!! output.csv 열 이름을 못 알아봤다: {list(r0)}")
            m = re.search(r"(p\d{5})", str(r0[kf]))
            if not m or m.group(1) not in alias:
                continue
            p = alias[m.group(1)]
            t, run, pose = parse_key(p)
            rows.append(dict(path=p, target=t, run=run, pose=pose,
                             score=round(float(r0[ks]), 5)))
            got += 1
        print(f"  묶음 {b0//a.batch + 1}: {got}/{len(chunk)} 채점")
        shutil.rmtree(work, ignore_errors=True)

    if len(rows) < len(files) * 0.5:
        sys.exit(f"!! 채점된 자세가 {len(rows)}/{len(files)} 로 절반 미만이다. 파일을 쓰지 않는다.")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["path", "target", "run", "pose", "score"])
        w.writeheader(); w.writerows(rows)
    print(f"\n→ {a.out}  ({len(rows)}행. 점수가 높을수록 계면이 정확하다는 뜻)")
    print(f"   다음: python scripts/eval_abepiscore.py eval --scores {a.out} --truth <정답 CSV>")


# ─────────────────────────── ② eval ───────────────────────────
def cmd_eval(a):
    sc = {}
    for r in csv.DictReader(open(a.scores)):
        sc[(r["target"], r["run"], r["pose"])] = float(r["score"])
    if not sc:
        sys.exit(f"!! {a.scores} 가 비어 있다")

    def num(r, k):
        try:
            return float(r[k])
        except (KeyError, TypeError, ValueError):
            return None

    joined = []
    miss = 0
    for tf in a.truth:
        for r in csv.DictReader(open(tf)):
            run = r.get("seed") or r.get("arm") or r.get("run") or "?"
            key = (r["target"], run, r.get("pose", "?"))
            if key not in sc:
                miss += 1; continue
            joined.append(dict(target=r["target"], run=run, score=sc[key],
                               dockq=num(r, "dockq"), recall=num(r, "recall"),
                               iptm=num(r, "iptm")))
    if not joined:
        sys.exit("!! 정답 CSV 와 붙는 행이 없다. run 키(seed/arm)와 pose 파일명이 맞는지 볼 것.")
    print(f"■ 정답 CSV {len(a.truth)}개 · 붙은 자세 {len(joined)}개 (못 붙은 행 {miss}개)\n")

    by_t = defaultdict(list)
    for r in joined:
        by_t[r["target"]].append(r)

    print(f"{'타깃':<12}{'자세':>5}{'ρ(점수,DockQ)':>14}{'ρ(점수,자리)':>13}"
          f"{'ρ(ipTM,DockQ)':>14}   실행 내 자세 선택: 천장 대 AbEpi 대 ipTM")
    print("-" * 108)
    agg = defaultdict(list)
    for t, ps in sorted(by_t.items()):
        dq = [p["dockq"] for p in ps if p["dockq"] is not None]
        s_dq = spearman([p["score"] for p in ps if p["dockq"] is not None], dq) if dq else None
        rc = [p["recall"] for p in ps if p["recall"] is not None]
        s_rc = spearman([p["score"] for p in ps if p["recall"] is not None], rc) if rc else None
        ip = [p for p in ps if p["iptm"] is not None and p["dockq"] is not None]
        s_ip = spearman([p["iptm"] for p in ip], [p["dockq"] for p in ip]) if len(ip) >= 3 else None
        # 실행(=자세 5개) 안에서 하나를 고르는 층
        orc = ab = it = n = 0.0
        for run in sorted({x["run"] for x in ps}):
            gg = [x for x in ps if x["run"] == run and x["dockq"] is not None]
            if not gg:
                continue
            orc += max(x["dockq"] for x in gg)
            ab += max(gg, key=lambda x: x["score"])["dockq"]
            gi = [x for x in gg if x["iptm"] is not None]
            it += max(gi, key=lambda x: x["iptm"])["dockq"] if gi else max(x["dockq"] for x in gg)
            n += 1
        f = lambda v: "  n/a" if v is None else f"{v:+.3f}"
        pick = f"{orc/n:.3f} / {ab/n:.3f} / {it/n:.3f}" if n else "n/a"
        print(f"{t:<12}{len(ps):>5}{f(s_dq):>14}{f(s_rc):>13}{f(s_ip):>14}   {pick}")
        for k, v in (("dq", s_dq), ("rc", s_rc), ("ip", s_ip)):
            if v is not None:
                agg[k].append(v)
        if n:
            agg["orc"].append(orc/n); agg["ab"].append(ab/n); agg["it"].append(it/n)

    # ── 조성(=실행) 하나를 고르는 층 ─────────────────────────────
    # 한 조성 안의 자세 5개는 거의 같으므로(조성 내 산포 ≪ 조성 간 산포) 위의 "실행 내
    # 자세 선택"은 고를 것이 없는 층이다. 이 연구가 실제로 필요한 것은 조성 여러 개 중
    # 어느 실행의 출력을 믿을지 고르는 일이다. 그 층을 여기서 잰다.
    print("\n■ 조성(=실행) 하나를 고를 때 — 이 연구가 실제로 필요한 층")
    print(f"{'타깃':<12}{'지표':>8}{'조성수':>7}{'천장':>8}{'AbEpi':>8}{'무작위':>8}"
          f"{'AbEpi 순위%':>12}   성공기준 도달")
    print("-" * 96)
    across = defaultdict(list)
    for t, ps in sorted(by_t.items()):
        for key, lab, thr in (("dockq", "DockQ", 0.49), ("recall", "자리겹침", 0.40)):
            runs = defaultdict(list)
            for r in ps:
                if r[key] is not None:
                    runs[r["run"]].append(r)
            if len(runs) < 3:
                continue
            # 실행 점수 = 그 실행 자세들의 최대 (정답을 보지 않는 값)
            agg = {k: (max(x["score"] for x in v), max(x[key] for x in v))
                   for k, v in runs.items()}
            vals = sorted((v[1] for v in agg.values()), reverse=True)
            picked = agg[max(agg, key=lambda k: agg[k][0])][1]
            ceil, rand = vals[0], st.mean(vals)
            # 고른 실행이 전체에서 몇 등인가 (100% = 최고, 0% = 최악)
            worse = sum(1 for v in vals if v < picked)
            pct = 100.0 * worse / max(len(vals) - 1, 1)
            hit = sum(1 for v in vals if v >= thr)
            ok = "○ 도달" if picked >= thr else "× 미달"
            print(f"{t:<12}{lab:>8}{len(agg):>7}{ceil:>8.3f}{picked:>8.3f}{rand:>8.3f}"
                  f"{pct:>11.0f}%   {ok} (기준 충족 조성 {hit}/{len(vals)})")
            across[key].append((ceil, picked, rand, pct, ceil >= thr, picked >= thr))

    # 전체 집계 — 이 검증의 결론이 나오는 자리
    for key, lab, thr in (("dockq", "DockQ", 0.49), ("recall", "자리겹침", 0.40)):
        v = across[key]
        if not v:
            continue
        n = len(v)
        reach = sum(1 for x in v if x[4])          # 천장이 기준을 넘는 타깃 = 건질 것이 있는 곳
        got = sum(1 for x in v if x[5])            # 그중 실제로 집어낸 곳
        print(f"\n■ 전체 {lab} (타깃 {n}개, 기준 {thr})")
        print(f"  평균   천장 {st.mean([x[0] for x in v]):.3f}  ·  "
              f"AbEpi 선택 {st.mean([x[1] for x in v]):.3f}  ·  "
              f"무작위 {st.mean([x[2] for x in v]):.3f}")
        print(f"  고른 조성의 순위 백분위 중앙값 {st.median([x[3] for x in v]):.0f}% "
              f"(무작위 선택의 기대값 = 50%)")
        print(f"  ★ 정답이 존재하는 타깃 {reach}개 중 실제로 집어낸 것 {got}개")

    print("\n■ 요약 (타깃 내 상관의 중앙값)")
    for k, lab in (("dq", "AbEpiScore ↔ DockQ"), ("rc", "AbEpiScore ↔ 자리 겹침"),
                   ("ip", "ipTM ↔ DockQ")):
        if agg[k]:
            v = agg[k]
            print(f"  {lab:<22} 중앙값 {st.median(v):+.3f}  (양수 {sum(1 for x in v if x>0)}/{len(v)})")
    if agg["orc"]:
        print(f"\n■ 실행 안에서 자세 하나를 고를 때의 평균 DockQ")
        print(f"  천장(정답 보고 최고) {st.mean(agg['orc']):.3f}  ·  "
              f"AbEpiScore 선택 {st.mean(agg['ab']):.3f}  ·  ipTM 선택 {st.mean(agg['it']):.3f}")
        print("  → AbEpiScore 값이 천장에 가까울수록 '정답 없이 고를 수 있다'는 뜻이다.")


def main():
    ap = argparse.ArgumentParser(description="AbEpiScore-1.0 로 자세를 채점하고 선택 성능을 재다")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("score")
    s.add_argument("--glob", nargs="+", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--repo", default=os.path.expanduser("~/tools/AbEpiTope-1.0"))
    s.add_argument("--batch", type=int, default=100)
    s.add_argument("--force", action="store_true")
    s.add_argument("--no-venv", action="store_true",
                   help="venv 를 무시하고 현재 파이썬 + PYTHONPATH=repo 로 실행")
    s.add_argument("--python", default=sys.executable)
    s.set_defaults(fn=cmd_score)
    e = sub.add_parser("eval")
    e.add_argument("--scores", required=True)
    e.add_argument("--truth", nargs="+", required=True,
                   help="정답 CSV 여럿을 한 번에 받아 표 하나로 낸다")
    e.set_defaults(fn=cmd_eval)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()
