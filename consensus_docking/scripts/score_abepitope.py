#!/usr/bin/env python3
"""[선택] AbEpiTope-1.0 으로 후보별 예측 복합체의 **계면 타당성**을 채점 → rank_sites.py 입력.

왜 이 도구인가 — DiscoTope 는 항원 표면만 보고 '일반적으로 항원성이 높은 자리'를 선호한다.
그런데 이 프로젝트의 명제가 바로 "MSA 조성이 모델을 **인기 있는 자리**로 끌고 간다"이므로,
인기 자리를 선호하는 도구로 후보를 고르면 빠져나오려던 편향을 선택 단계에서 다시 넣는 꼴이 된다.
AbEpiTope 은 **복합체**를 보고 "이 항체가 여기 붙는 게 말이 되나"를 채점하므로 항체 특이적이고,
설계 목적 자체가 예측된 항체-항원 구조를 ipTM 보다 잘 랭킹하는 것이라 우리 명제와 같은 자리에 있다.

═══ 흐름 (2패스) ═══
  ① sites_to_pocket.py --arm all      후보 전부를 제약 입력으로
  ② 후보마다 값싸게 1샘플 예측         prescreen_abepitope.sh
  ③ 이 파일: 후보별 복합체 → 점수      → results/abepitope_<타깃>.csv (cand,score)
  ④ rank_sites.py --method abepitope --scores ...   → 이긴 후보
  ⑤ 그 후보로 본 유도 예측(run_demo_guided.sh METHOD=manual --pick)

⚠️ 이 도구의 파이썬 API·CLI 는 판올림마다 달라진다. 아래 순서로 시도하고, 전부 실패하면
   무엇을 만들어야 하는지(cand,score 두 열 CSV)를 찍고 멈춘다 — 조용히 틀린 값을 쓰지 않는다.
     (1) --cmd 로 사용자가 준 커맨드 틀
     (2) 설치된 파이썬 패키지의 흔한 진입점 몇 가지
   설치 후 README 를 확인해 (1) 로 고정해 두는 것을 권한다.

사용:
  python scripts/score_abepitope.py --target 8k3k_D --dir runs_sites_guided/8k3k_d/all
  python scripts/score_abepitope.py --target 8k3k_D --dir ... \\
      --cmd "conda run -n abepitope abepitope score --pdb {structure} --out {out}"
"""
import argparse, csv, glob, json, os, re, subprocess, sys, tempfile


def find_structures(d, models=("boltz", "protenix")):
    """후보 번호 → 그 후보의 예측 구조 하나. 파일명/폴더명의 patch<번호> 로 후보를 알아낸다."""
    out = {}
    for p in sorted(glob.glob(os.path.join(d, "**", "*.cif"), recursive=True)
                    + glob.glob(os.path.join(d, "**", "*.pdb"), recursive=True)):
        m = re.search(r"patch(\d+)", p)
        if not m:
            continue
        k = int(m.group(1))
        out.setdefault(k, p)          # 후보당 하나면 충분(사전 스크리닝이라 1샘플)
    return out


def score_by_cmd(cmd, struct):
    """--cmd 틀 실행 → stdout 또는 --out 파일에서 숫자 하나를 뽑는다."""
    with tempfile.TemporaryDirectory() as td:
        outf = os.path.join(td, "score.txt")
        c = cmd.format(structure=struct, out=outf)
        r = subprocess.run(c, shell=True, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"    ! 실패: {c}\n      {r.stderr.strip()[-300:]}")
            return None
        txt = open(outf).read() if os.path.exists(outf) else r.stdout
        nums = re.findall(r"-?\d+\.\d+", txt)
        if not nums:
            print(f"    ! 출력에서 점수를 못 찾았다: {txt.strip()[:200]}")
            return None
        return float(nums[-1])


def score_by_api(struct):
    """설치된 패키지의 흔한 진입점들을 시도. 하나도 안 맞으면 None."""
    try:
        import abepitope                                          # noqa: F401
    except Exception:
        return None
    for path in ("abepitope.AbEpiTope", "abepitope.model.AbEpiTope", "abepitope.api.score"):
        mod, _, attr = path.rpartition(".")
        try:
            m = __import__(mod, fromlist=[attr])
            fn = getattr(m, attr)
            r = fn(struct) if not isinstance(fn, type) else fn().score(struct)
            if isinstance(r, (int, float)):
                return float(r)
            if isinstance(r, dict):
                for k in ("score", "AbEpiScore", "abepiscore", "interface_score"):
                    if k in r:
                        return float(r[k])
        except Exception:
            continue
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--dir", required=True, help="후보별 예측물이 있는 폴더(--arm all 실행 결과)")
    ap.add_argument("--cmd", default=os.environ.get("ABEPITOPE_CMD", ""),
                    help="커맨드 틀. {structure} 와 {out} 을 채워 실행한다")
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    st = find_structures(a.dir)
    if not st:
        sys.exit(f"!! {a.dir} 에서 patch<번호> 가 붙은 예측 구조를 못 찾았다.\n"
                 f"   먼저: bash scripts/prescreen_abepitope.sh {a.target}")
    print(f"■ {a.target}   후보 {len(st)}개의 예측 구조를 채점한다")

    rows, failed = [], 0
    for k in sorted(st):
        s = None
        if a.cmd:
            s = score_by_cmd(a.cmd, st[k])
        if s is None:
            s = score_by_api(st[k])
        if s is None:
            failed += 1
            print(f"  후보 {k}: 채점 실패 ({os.path.basename(st[k])})")
        else:
            print(f"  후보 {k}: {s:.4f}   ({os.path.basename(st[k])})")
            rows.append(dict(cand=k, score=round(s, 5)))

    if not rows:
        print("\n!! 한 개도 채점하지 못했다. AbEpiTope 의 실제 사용법을 확인하고 --cmd 로 지정할 것.")
        print("   예: --cmd \"conda run -n abepitope python -m abepitope --pdb {structure} --out {out}\"")
        print(f"   또는 직접 다음 형식의 CSV 를 만들어 rank_sites.py 에 넘겨도 된다:")
        print(f"     cand,score")
        print(f"     1,0.87")
        print(f"     2,0.31")
        raise SystemExit(2)
    if failed:
        print(f"\n  ⚠️ {failed}개 후보를 채점하지 못했다 — 그 후보는 선택에서 빠진다(최하위 아님, 아예 제외).")

    out = a.out or f"results/abepitope_{a.target}.csv"
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["cand", "score"])
        w.writeheader(); w.writerows(rows)
    print(f"\n→ {out}")
    print(f"   다음: python scripts/rank_sites.py --sites <sites_{a.target}.json> "
          f"--method abepitope --scores {out}")


if __name__ == "__main__":
    main()
