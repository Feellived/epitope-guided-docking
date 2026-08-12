#!/usr/bin/env python3
"""[선택] MSA 조성이 만든 결합자리 후보 중 하나를 **정답 구조를 보지 않고** 고른다.

생성(조성 다양화)과 선택은 다른 문제다. bk21-msa-depth-bias 본 검정에서 확인된 것:
  · 생성 ○ — 조성을 다양화하면 정답 자리가 후보 안에 들어온다(천장 16/18)
  · 선택 ✗ — '가장 많은 조성이 모인 후보'를 고르면 13/18 로, 아무것도 안 한 것(13/18)과 같다
따라서 이 파일의 목적은 **더 나은 선택 규칙을 붙일 자리**를 만들어 두는 것이다.

═══ 정답을 보지 않는다는 규약 ═══
입력 JSON(sites_<타깃>.json)에는 보고용으로 true_covered·precision 이 들어 있다. 이 둘은
native 구조를 보고 계산한 값이라 **고르는 데 쓰면 데모가 무효**가 된다. 그래서 여기서는
불러오는 즉시 버리고, 버렸다는 사실을 화면에 찍는다.

═══ 규칙(--method) ═══
  ncomp      조성이 가장 많이 모인 후보. 추가 설치 불필요. 위에 적었듯 기준선을 못 이긴다(대조용).
  discotope  DiscoTope-3.0 잔기별 점수의 후보 평균. 항원 구조만 보고 'ω여기가 에피토프다운가'를
             재므로 MSA를 쓰지 않는다 = 이 프로젝트 설계 철학과 일치. 별도 설치 필요.
             입력 = --scores CSV (열: chain,pos,score 또는 residue,score). 1-based 서열 위치 기준.
  abepitope  AbEpiTope-1.0 등 '계면 타당성' 점수. 후보별 점수 CSV(cand,score)를 받는다.
  manual     --pick 으로 사람이 지정(발표에서 구두 고지하는 경우).

사용:
  python scripts/select_site.py --sites results/sites_8k3k_D.json --method ncomp
  python scripts/select_site.py --sites ... --method discotope --scores discotope_8k3k.csv
  python scripts/select_site.py --sites ... --method manual --pick 1
출력: 고른 후보 번호를 stdout 마지막 줄에 'PICK=<n>' 으로 찍는다(스크립트에서 받아쓰기 좋게).
"""
import argparse, csv, json, os, sys

BLIND_DROP = ("true_covered", "precision")   # native 를 본 값 — 선택에 쓰지 않는다


def load_sites(path):
    d = json.load(open(path))
    dropped = 0
    for c in d["candidates"]:
        for k in BLIND_DROP:
            if k in c:
                c.pop(k); dropped += 1
    if dropped:
        print(f"[눈먼 선택] 정답 구조를 본 값 {dropped}개(true_covered·precision)를 버리고 시작한다.")
    return d


def score_ncomp(d, _a):
    """조성이 많이 모일수록 높다. 동점이면 잔기 수가 적은 쪽(더 좁게 특정한 쪽)."""
    return {c["cand"]: (c["n_comp"], -len(c["residues"])) for c in d["candidates"]}


def _read_residue_scores(path):
    """DiscoTope 류 잔기별 점수 → {(사슬순번, 0-based 위치): 점수}.
    사슬 열이 없으면 항원이 하나라고 보고 순번 0 으로 둔다. 위치는 1-based 로 받아 0-based 로 낮춘다."""
    out = {}
    for r in csv.DictReader(open(path)):
        k = {x.lower(): v for x, v in r.items()}
        pos = k.get("pos") or k.get("position") or k.get("residue") or k.get("resid")
        sc = k.get("score") or k.get("discotope") or k.get("value")
        if pos is None or sc is None:
            continue
        ch = k.get("chain_ord") or k.get("chain_index")
        try:
            ci = int(ch) if ch not in (None, "") and str(ch).lstrip("-").isdigit() else 0
            out[(ci, int(float(pos)) - 1)] = float(sc)
        except Exception:
            continue
    if not out:
        sys.exit(f"!! {path} 에서 잔기 점수를 못 읽었다. 열 이름 = chain_ord,pos,score 를 기대한다.")
    return out


def score_discotope(d, a):
    if not a.scores:
        sys.exit("!! --method discotope 에는 --scores <잔기별 점수 CSV> 가 필요하다.")
    rs = _read_residue_scores(a.scores)
    out = {}
    for c in d["candidates"]:
        v = [rs[(k[0], k[1])] for k in map(tuple, c["residues"]) if (k[0], k[1]) in rs]
        if not v:
            print(f"  ! 후보 {c['cand']}: 점수와 겹치는 잔기가 없다 — 최하위로 둔다")
            out[c["cand"]] = (float("-inf"), 0)
        else:
            # 평균 = '이 자리 전체가 얼마나 에피토프다운가'. 잔기 수가 적은 후보가 우연히 높게
            # 나오는 것을 막으려면 --min-res 로 하한을 준다.
            out[c["cand"]] = (sum(v) / len(v), len(v))
    return out


def score_abepitope(d, a):
    if not a.scores:
        sys.exit("!! --method abepitope 에는 --scores <후보별 점수 CSV(cand,score)> 가 필요하다.")
    m = {}
    for r in csv.DictReader(open(a.scores)):
        k = {x.lower(): v for x, v in r.items()}
        if k.get("cand") is not None and k.get("score") is not None:
            m[int(k["cand"])] = float(k["score"])
    return {c["cand"]: (m.get(c["cand"], float("-inf")), 0) for c in d["candidates"]}


METHODS = {"ncomp": score_ncomp, "discotope": score_discotope, "abepitope": score_abepitope}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", required=True, help="msa-depth 의 results/sites_<타깃>.json")
    ap.add_argument("--method", default="ncomp", choices=list(METHODS) + ["manual"])
    ap.add_argument("--scores", default="", help="discotope/abepitope 점수 파일")
    ap.add_argument("--pick", type=int, default=0, help="--method manual 일 때 고를 후보 번호")
    ap.add_argument("--min-res", type=int, default=0, help="이 잔기 수보다 작은 후보는 제외")
    ap.add_argument("--out", default="", help="선택 결과를 JSON 으로도 남긴다")
    a = ap.parse_args()

    d = load_sites(a.sites)
    cands = [c for c in d["candidates"] if len(c["residues"]) >= a.min_res]
    if not cands:
        sys.exit("!! 조건을 만족하는 후보가 없다")
    print(f"■ {d['target']} · {d['model']} · {d['depth']}   후보 {len(d['candidates'])}개"
          f"{f' (잔기 {a.min_res}개 미만 제외 후 {len(cands)}개)' if a.min_res else ''}")

    if a.method == "manual":
        if not a.pick:
            sys.exit("!! --method manual 에는 --pick <후보번호> 가 필요하다")
        pick = a.pick
        print(f"  [사람이 지정] 후보 {pick}")
        print("  ⚠️ 이 선택은 눈먼 선택이 아니다. 발표·문서에서 '사람이 골랐다'를 반드시 함께 말할 것.")
    else:
        sc = METHODS[a.method](d, a)
        print(f"\n  {'후보':5}{'조성수':7}{'잔기':6}{'원래MSA':>8}{'점수':>12}")
        print("  " + "-" * 42)
        for c in sorted(cands, key=lambda c: -sc[c["cand"]][0] if sc[c["cand"]][0] == sc[c["cand"]][0] else 0):
            s = sc[c["cand"]][0]
            print(f"  {c['cand']:<5}{c['n_comp']:<7}{len(c['residues']):<6}"
                  f"{'예' if c['from_full_msa'] else '아니오':>8}"
                  f"{(f'{s:.3f}' if s == s and s != float('-inf') else '  -  '):>12}")
        pick = max(cands, key=lambda c: sc[c["cand"]])["cand"]
        print(f"\n  [{a.method}] 고른 후보 = {pick}")

    chosen = next(c for c in cands if c["cand"] == pick)
    full = next((c for c in d["candidates"] if c["from_full_msa"]), None)
    print(f"  잔기 {len(chosen['residues'])}개 · 모인 조성 {chosen['n_comp']}가지 ({','.join(chosen['comps'])})")
    if full:
        same = full["cand"] == pick
        print(f"  원래 MSA 가 간 자리 = 후보 {full['cand']}"
              + ("  ← 고른 것과 같다(대조군과 구별이 안 되므로 데모 효과 없음)" if same
                 else "  ← 고른 것과 다르다(대조군으로 쓸 수 있다)"))
    else:
        print("  원래 MSA 가 속한 후보가 없다 — 대조군은 별도로 정해야 한다")

    if a.out:
        json.dump(dict(target=d["target"], method=a.method, pick=pick,
                       full_msa_cand=(full["cand"] if full else None)),
                  open(a.out, "w"), indent=1)
        print(f"  → {a.out}")
    print(f"PICK={pick}")


if __name__ == "__main__":
    main()
