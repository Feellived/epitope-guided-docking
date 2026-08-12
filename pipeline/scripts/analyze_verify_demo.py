#!/usr/bin/env python3
"""[적대적 점검] 데모 표를 믿기 전에, 비교가 성립하는지 자체를 검사한다.

왜 필요한가 — 2026-07-29 하루에 이 파이프라인에서 '멈추지 않고 조용히 잘못되는' 실패를
다섯 건 잡았다. 숫자가 나왔다는 것이 비교가 성립한다는 뜻이 아니다. 아래는 실제로 겪은
실패를 그대로 검사 항목으로 옮긴 것이다.

  1) 팔마다 자세 수가 다르다        → 예산이 안 맞아 best-of-N 이 한쪽에 유리
     (실제: ours 1개 vs sizematch 5개로 표가 나왔었다)
  2) 두 팔의 제약 잔기가 같다       → 다른 팔인 척하지만 같은 실험
     (고른 후보가 원래 MSA 자리와 같으면 fullmsa 와 ours 가 같아진다)
  3) sizematch 크기가 ours 와 다르다 → 크기 교란을 배제하지 못한다
  4) noconstraint 입력에 제약이 있다 → 기준선이 기준선이 아니다
  5) 모델·샘플 수가 팔마다 다르다    → 모델 차이가 섞인다
     (실제: PROT_MODEL 이 blind 와 달라 비교가 깨질 뻔했다)
  6) run.log 에 조용한 경고가 있다   → MSA 를 버리고 단일서열로 예측한 사고와 같은 계열
  7) blind 깊이가 sites JSON 과 다르다 → 다른 실험의 자료를 기준선으로 쓴 것

문제를 찾으면 종료 코드 1. 발표에 쓸 표는 이 검사를 통과한 뒤에 쓴다.

사용:
  python scripts/analyze_verify_demo.py --targets "8k3k_D 8sit_HL 8ume_HL"
  python scripts/analyze_verify_demo.py --targets 8k3k_D --models boltz --arms "noconstraint ours"
"""
import argparse, glob, json, os, re, sys

ARMS = ["noconstraint", "fullmsa", "sizematch", "ours"]


def boltz_contacts(p):
    """boltz YAML 에서 pocket 접촉 잔기 집합을 읽는다."""
    if not os.path.exists(p):
        return None
    s = open(p).read()
    if "pocket:" not in s:
        return set()
    return {int(m) for m in re.findall(r"^\s*-\s*\[\s*\w+\s*,\s*(\d+)\s*\]", s, re.M)}


def protenix_contacts(p):
    if not os.path.exists(p):
        return None
    d = json.load(open(p))
    c = (d[0] or {}).get("constraint") or {}
    pk = c.get("pocket") or {}
    return {int(x["position"]) for x in pk.get("contact_residues", [])}


def poses(d):
    return glob.glob(os.path.join(d, "**", "*.cif"), recursive=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", required=True)
    ap.add_argument("--outdir", default="runs_sites_guided")
    ap.add_argument("--models", default="boltz protenix")
    ap.add_argument("--arms", default=" ".join(ARMS))
    ap.add_argument("--msad", default=os.path.expanduser("~/projects/msa-composition-bias/pipeline"))
    ap.add_argument("--msad-data", default="/mnt/data/msadepth")
    a = ap.parse_args()
    arms = a.arms.split()
    models = a.models.replace(",", " ").split()
    bad = warn = 0

    def fail(msg):
        nonlocal bad
        bad += 1
        print("  ❌ " + msg)

    def note(msg):
        nonlocal warn
        warn += 1
        print("  ⚠️  " + msg)

    for t in a.targets.replace(",", " ").split():
        print("\n■ %s" % t)
        base = os.path.join(a.outdir, t.lower())
        if not os.path.isdir(base):
            fail("%s 없음 — 데모를 돌리지 않았다" % base); continue

        # ── 제약 잔기 비교 (2,3,4번) ──────────────────────────────────────────
        con = {}
        for arm in arms:
            d = os.path.join(base, arm)
            for m in models:
                p = (os.path.join(d, "boltz_%s_patch0.yaml" % t.lower()) if m == "boltz"
                     else os.path.join(d, "protenix_%s_patch0.json" % t.lower()))
                c = boltz_contacts(p) if m == "boltz" else protenix_contacts(p)
                if c is not None:
                    con[(arm, m)] = c
        for m in models:
            got = {arm: con[(arm, m)] for arm in arms if (arm, m) in con}
            if not got:
                continue
            if "noconstraint" in got and got["noconstraint"]:
                fail("[%s] noconstraint 입력에 제약 잔기가 %d개 있다 — 기준선이 아니다"
                     % (m, len(got["noconstraint"])))
            if "ours" in got and "fullmsa" in got:
                if got["ours"] == got["fullmsa"]:
                    fail("[%s] ours 와 fullmsa 의 제약 잔기가 완전히 같다 — 두 팔이 같은 실험이다" % m)
                else:
                    ov = len(got["ours"] & got["fullmsa"]) / max(1, len(got["ours"] | got["fullmsa"]))
                    if ov > 0.8:
                        note("[%s] ours 와 fullmsa 가 %.0f%% 겹친다 — 대조가 약하다" % (m, ov * 100))
            if "ours" in got and "sizematch" in got:
                if len(got["ours"]) != len(got["sizematch"]):
                    fail("[%s] sizematch %d잔기 ≠ ours %d잔기 — 크기 교란을 배제하지 못한다"
                         % (m, len(got["sizematch"]), len(got["ours"])))
                if got["sizematch"] == got["ours"]:
                    fail("[%s] sizematch 가 ours 와 같은 잔기다 — 대조가 아니다" % m)
            sizes = ", ".join("%s %d" % (k, len(v)) for k, v in sorted(got.items()))
            print("  제약 잔기 수 [%s]: %s" % (m, sizes))

        # ── 자세 수·로그 (1,5,6번) ────────────────────────────────────────────
        for m in models:
            counts = {}
            for arm in arms:
                o = os.path.join(base, arm, "out_%s" % m)
                if not os.path.isdir(o):
                    continue
                counts[arm] = len(poses(o))
                lg = os.path.join(o, "run.log")
                if os.path.exists(lg):
                    txt = open(lg, errors="ignore").read().lower()
                    for pat, why in (("no msa", "MSA 없이 예측했을 수 있다"),
                                     ("could not", "실행 중 오류"),
                                     ("traceback", "예외 발생"),
                                     ("out of memory", "메모리 부족"),
                                     ("warning: ", "경고")):
                        if pat in txt:
                            note("[%s/%s] run.log 에 '%s' — %s" % (m, arm, pat.strip(), why))
                            break
                elif counts[arm]:
                    note("[%s/%s] run.log 이 없다 — 어떻게 돌았는지 확인 불가" % (m, arm))
            if not counts:
                continue
            uniq = set(counts.values())
            if len(uniq) > 1:
                fail("[%s] 팔마다 자세 수가 다르다 %s — 예산이 안 맞아 비교가 성립하지 않는다"
                     % (m, counts))
            else:
                print("  자세 수 [%s]: 모든 팔 %d개 (예산 일치)" % (m, uniq.pop()))

        # ── blind 기준선의 깊이 (7번) ────────────────────────────────────────
        sp = os.path.join(a.msad, "results", "sites_%s.json" % t)
        if os.path.exists(sp):
            sj = json.load(open(sp))
            bd = os.path.join(a.msad_data, "compreps", "seedrep_cand",
                              sj["model"], t, sj["depth"])
            n = len(glob.glob(os.path.join(bd, "seedfull_r*")))
            if n == 0:
                note("blind 기준선 폴더가 비었다: %s" % bd)
            else:
                print("  blind 기준선: %s %s · 실행 %d개" % (sj["model"], sj["depth"], n))
            others = glob.glob(os.path.join(os.path.dirname(bd), "d*"))
            if len(others) > 1:
                fail("blind 쪽에 깊이 폴더가 %d개다 (%s) — 다른 실험 자료가 섞일 수 있다"
                     % (len(others), ", ".join(sorted(os.path.basename(x) for x in others))))
        else:
            note("%s 없음 — blind 기준선을 확인할 수 없다" % sp)

    print("\n" + "=" * 60)
    if bad:
        print("❌ 문제 %d건 · 경고 %d건 — 이 표를 발표에 쓰지 말 것" % (bad, warn))
        sys.exit(1)
    print("✅ 문제 없음 · 경고 %d건" % warn)
    print("   (경고는 판단 재료다. 통과가 '결과가 옳다'는 뜻은 아니고, '비교가 성립한다'는 뜻이다.)")


if __name__ == "__main__":
    main()
