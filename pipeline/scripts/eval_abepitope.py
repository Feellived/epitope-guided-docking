#!/usr/bin/env python3
"""[선택] AbEpiTope-1.0 (AbEpiScore-1.0) 로 후보별 예측 복합체의 계면 정확도를 채점.

왜 이 도구인가 — DiscoTope 는 항원 표면만 보고 '일반적으로 항원성이 높은 자리'를 선호한다.
이 프로젝트의 명제가 "MSA 조성이 모델을 **인기 있는 자리**로 끌고 간다"이므로, 인기 자리를
선호하는 도구로 후보를 고르면 빠져나오려던 편향을 선택 단계에서 다시 넣는 순환 논리가 된다.
AbEpiScore-1.0 은 **모델링된 항체-항원 계면이 얼마나 정확한가**를 채점하므로 항체 특이적이고,
논문(Science Advances 2025)에서 AlphaFold 신뢰도보다 계면 정확도를 잘 맞힌다고 보고됐다 —
"망가진 ipTM 대신 계면-국소 재랭커"라는 우리 명제와 같은 자리에 있는 도구다.

═══ 흐름 (2패스) ═══
  ① make_pocket_from_sites.py --arm all   후보 전부를 제약 입력으로
  ② 후보마다 값싸게 1샘플 예측       run_prescreen_abepitope.sh
  ③ 이 파일                        후보별 복합체 → results/abepitope_<타깃>.csv (cand,score)
  ④ select_site.py --method abepitope --scores ...   → 이긴 후보
  ⑤ 그 후보로 본 유도 예측

AbEpiTope 은 uv 로 만든 **자체 venv**(레포 안 .venv)에서 돈다. 그 venv 의 activate 에
CUDA 라이브러리 경로가 들어 있으므로, 파이썬을 직접 부르지 않고 activate 를 거쳐 부른다.

사용:
  python scripts/eval_abepitope.py --target 8k3k_D --dir runs_sites_guided/8k3k_d/all
  python scripts/eval_abepitope.py --target 8k3k_D --dir ... --repo ~/tools/AbEpiTope-1.0
"""
import argparse, csv, glob, os, re, shutil, subprocess, sys, tempfile

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


def find_structures(d):
    """후보 번호 → 그 후보의 예측 구조 하나. 경로의 patch<번호> 로 후보를 알아낸다."""
    out = {}
    for p in sorted(glob.glob(os.path.join(d, "**", "*.cif"), recursive=True)
                    + glob.glob(os.path.join(d, "**", "*.pdb"), recursive=True)):
        m = re.search(r"patch(\d+)", p)
        if m:
            out.setdefault(int(m.group(1)), p)     # 후보당 하나면 충분(사전 스크리닝 1샘플)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--dir", required=True, help="후보별 예측물 폴더(--arm all 실행 결과)")
    ap.add_argument("--repo", default=os.path.expanduser("~/tools/AbEpiTope-1.0"))
    ap.add_argument("--keep", default="", help="중간 산출물(인코딩·output.csv)을 남길 폴더")
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    act = os.path.join(a.repo, ".venv", "bin", "activate")
    if not os.path.exists(act):
        sys.exit(f"!! AbEpiTope venv 를 못 찾았다: {act}\n"
                 f"   설치: bash scripts/run_setup_tools.sh --only abepitope --apply")

    st = find_structures(a.dir)
    if not st:
        sys.exit(f"!! {a.dir} 에서 patch<번호> 가 붙은 예측 구조를 못 찾았다.\n"
                 f"   먼저: bash scripts/run_prescreen_abepitope.sh {a.target}")
    print(f"■ {a.target}   후보 {len(st)}개의 예측 구조를 AbEpiScore-1.0 으로 채점")

    work = a.keep or tempfile.mkdtemp(prefix="abepi_")
    inp = os.path.join(work, "in"); os.makedirs(inp, exist_ok=True)
    # 파일 이름에 후보 번호를 박아 둔다 — 출력 CSV 의 filename 으로 되짚기 위해.
    for k, p in sorted(st.items()):
        shutil.copy(p, os.path.join(inp, f"cand{k}{os.path.splitext(p)[1]}"))
        print(f"  후보 {k} ← {os.path.basename(p)}")
    drv = os.path.join(work, "driver.py"); open(drv, "w").write(DRIVER)
    enc, tmp, out = (os.path.join(work, x) for x in ("enc", "tmp", "out"))

    cmd = f'source "{act}" && python "{drv}" "{inp}" "{enc}" "{tmp}" "{out}"'
    print(f"  $ bash -c '{cmd[:90]}...'")
    r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    if "ABEPITOPE_OK" not in r.stdout:
        print(r.stdout[-1500:]); print(r.stderr[-1500:])
        sys.exit("!! AbEpiTope 실행 실패 — 위 메시지를 볼 것 (GPU 메모리·CUDA 경로가 흔한 원인)")

    res = None
    for c in glob.glob(os.path.join(out, "**", "output.csv"), recursive=True):
        res = c; break
    if not res:
        sys.exit(f"!! output.csv 를 못 찾았다: {out}")

    rows, miss = [], 0
    for r0 in csv.DictReader(open(res)):
        keys = {k.lower(): k for k in r0}
        kf = next((keys[k] for k in keys if "file" in k), None)
        ks = next((keys[k] for k in keys if "abepiscore" in k), None)
        if not (kf and ks):
            sys.exit(f"!! output.csv 열 이름을 못 알아봤다: {list(r0)}")
        m = re.search(r"cand(\d+)", str(r0[kf]))
        if not m:
            miss += 1; continue
        rows.append(dict(cand=int(m.group(1)), score=round(float(r0[ks]), 5)))
    if not rows:
        sys.exit(f"!! {res} 에서 후보 점수를 못 읽었다")
    rows.sort(key=lambda x: -x["score"])
    print(f"\n  {'후보':6}{'AbEpiScore-1.0':>16}")
    for x in rows:
        print(f"  {x['cand']:<6}{x['score']:>16.4f}")
    if miss:
        print(f"  ⚠️ 후보 번호를 못 읽은 행 {miss}개")

    o = a.out or f"results/abepitope_{a.target}.csv"
    os.makedirs(os.path.dirname(o) or ".", exist_ok=True)
    with open(o, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["cand", "score"]); w.writeheader(); w.writerows(rows)
    print(f"\n→ {o}   (점수가 높을수록 계면이 정확하다는 뜻)")
    if a.keep:
        print(f"   중간 산출물: {work}  (interface.csv 에 에피토프·파라토프 잔기가 있다)")


if __name__ == "__main__":
    main()
