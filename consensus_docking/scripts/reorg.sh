#!/usr/bin/env bash
# consensus_docking 폴더 재구성 (2026-07-14, 사용자 승인: 전체정리 + 레거시 archive/ + .pse git rm)
# 방식: git mv/git rm (history 보존). ★ 커밋·push는 하지 않는다 —
#       실행 후 `git status`로 검토하고 직접 커밋/푸시.
# 실행: repo 클론에서
#   cd <repo>/consensus_docking && bash scripts/reorg.sh
# 런타임 출력 폴더(runs_diverse·results_boltz·out_*·results_colabfold)는 절대 안 건드림.
set -euo pipefail
cd "$(dirname "$0")/.."                       # -> consensus_docking
[ -d ../.git ] || { echo "!! git repo 루트(../.git) 없음. repo 클론에서 실행하세요."; exit 1; }
G(){ echo "  + $*"; "$@"; }

echo "== 0) 대상 폴더 생성 =="
mkdir -p docs configs archive figures

echo "== 1) .pse 바이너리 제거 (git rm; history엔 남고 재생성 가능) =="
G git rm -q --ignore-unmatch fig_diverse/*.pse fig_5model/*.pse

echo "== 2) 그림 폴더 통합 -> figures/ (PNG만) =="
[ -d fig ]         && G git mv fig figures/analysis
[ -d fig_diverse ] && G git mv fig_diverse figures/cases_diverse
[ -d fig_5model ]  && G git mv fig_5model figures/cases_5model

echo "== 3) 문서 -> docs/ (이름 통일) =="
[ -f benchmark_targets_2026.md ] && G git mv benchmark_targets_2026.md docs/benchmark_2026.md
[ -f RUN_ALL.md ]                && G git mv RUN_ALL.md docs/run_guide.md
[ -f AUTOMATION.md ]             && G git mv AUTOMATION.md docs/automation.md

echo "== 4) 설정 통합 -> configs/ (옛 9VP0 config는 archive) =="
[ -f configs/targets.yaml ]          && G git mv configs/targets.yaml archive/configs_targets_9vp0.yaml
[ -f scripts/targets_diverse.yaml ]  && G git mv scripts/targets_diverse.yaml configs/targets_diverse.yaml
[ -f scripts/targets.yaml ]          && G git mv scripts/targets.yaml configs/targets_hackathon.yaml

echo "== 5) 레거시 스크립트 -> archive/ =="
for f in 01_prep_9vp0.py 02_run_boltz.sh 03_eval.sh render_cases.sh rerun_boltz_failed.sh; do
  [ -f "scripts/$f" ] && G git mv "scripts/$f" "archive/$f"
done

echo "== 6) 참조 수정 (이동으로 깨지는 4곳; perl -i = mac/linux 공통) =="
perl -pi -e 's{OUT=fig_diverse\b}{OUT=figures/cases_diverse}g' scripts/render_cases_diverse.sh
perl -pi -e 's{OUT=fig_5model\b}{OUT=figures/cases_5model}g'   scripts/render_5model.sh
perl -pi -e 's{default="targets\.yaml"}{default="configs/targets_diverse.yaml"}g' scripts/extract_chains.py
for d in docs/run_guide.md docs/automation.md; do
  [ -f "$d" ] && perl -pi -e 's{scripts/targets\.yaml}{configs/targets_hackathon.yaml}g' "$d"
done
echo "  (수정: render_cases_diverse OUT · render_5model OUT · extract_chains 기본값 · docs targets 경로)"

echo "== 7) README 재작성 (현재 diverse 흐름 + Phase 0) =="
cat > README.md <<'READMEEOF'
# Consensus Docking — 항체-항원 pose 재랭킹 / 오케스트레이션

여러 공개 co-folder(Boltz-2 · Chai-1 · Protenix · AF2-Multimer · tFold-Ag)를 한 파이프라인으로
돌려, 단일 ipTM 대신 **모델 간 합의(consensus)로 native-like pose를 고른다.**
성공 기준 = **DockQ ≥ 0.49 (CAPRI Medium; 2026-07-14 상향, 3-tier 0.23/0.49/0.80 병기)**
+ **RMSD 이원화**(항체 내부 fold RMSD / Ligand-RMSD 배치) + epitope recall.

> 메인 프로젝트(사수 논문 라인)와 해커톤이 개념 수렴 — 모델 풀·방법 공유, 산출물은 각각.

## 폴더 구조
- `docs/`     — 방법론·가이드: `benchmark_2026.md` · `run_guide.md` · `automation.md`
- `configs/`  — 타깃 세트: `targets_diverse.yaml`(현재 메인 10-타깃 Fab) · `targets_hackathon.yaml`(9VP0 + 5-타깃)
- `scripts/`  — 코드(아래 파이프라인)
- `results/`  — 커밋되는 `.txt` 산출(diagnose_* · where_docked_* · merged_*). `summary.csv`는 gitignore
- `figures/`  — `analysis/`(집계 그림) · `cases_diverse/` · `cases_5model/`  (PNG만; `.pse`는 gitignore)
- `archive/`  — 레거시(9VP0 브링업 01/02/03 · render_cases · rerun_boltz_failed · 옛 configs)

## 파이프라인
1. `extract_chains.py` — native 다운 + keep_chains 서열 → `runs_diverse/<id>/chains.json`
2. `make_fv.py` — 항체 사슬 ANARCI(IMGT)로 Fv(VH+VL) 절단 → `runs_diverse/<id>_fv/`
3. `make_{boltz,chai,protenix,colabfold,tfold}.py` — 모델별 입력 생성
4. 예측: `run_diverse.sh`(chai·protenix·colabfold) · `run_tfold.sh` · `run_5model.sh`(=colabfold+tfold 묶음)
   → 각 `runs_diverse/<id>/{results_boltz,out_chai,out_protenix,results_colabfold,out_tfold}`
5. eval: `eval_{boltz,chai,colabfold,protenix,tfold}.py` → `results/summary.csv`(source of truth)
   → `sync_notion.py`(개인 Mac, 토큰 env)
6. 분석: `consensus_analysis.py`·`plot_iptm_dockq.py`(집계) /
   `diagnose_antibody_side.py`·`where_docked.py`·`merge_dockq.py`(정밀 진단)
7. 렌더: `pymol_pose.py` + `render_{cases_diverse,5model}.sh` → `figures/`

## Phase 0 (go/no-go 게이트)
`run_phase0.sh` = diagnose + where_docked + merge_dockq 를 10-타깃 × Fab/Fv 전부에 실행.
새 예측 없이 기존 pose만 재분석 → 0.49 기준 oracle·placement 여부·ipTM 품질로 프로젝트 go/no-go 판정.

## 규칙
- 런타임 출력 폴더(`runs_diverse`·`results_boltz`·`out_*`·`results_colabfold`)는 스크립트에 하드코딩 → **이름 변경 금지.**
- 결과 구조·MSA·weight·CSV·`.pse`는 git에 안 올림(.gitignore).
- Notion 토큰은 서버·깃·문서에 저장 금지(개인 Mac env에서만).
READMEEOF
echo "  README.md 재작성 완료"

echo "== 8) .gitignore 보강 (repo 루트) =="
GI=../.gitignore
grep -qxF '*.pse' "$GI" || { printf '\n# PyMOL 세션(재생성 가능, 대용량)\n*.pse\n' >> "$GI"; echo "  + *.pse"; }
grep -qxF 'consensus_docking/runs_diverse/' "$GI" || { echo 'consensus_docking/runs_diverse/' >> "$GI"; echo "  + consensus_docking/runs_diverse/"; }

echo "== 9) 스테이징 =="
G git add -A ../.gitignore .

echo
echo "== 10) 자기점검: 이동한 옛 경로가 (archive 제외) 아직 참조되나 =="
LEFT=$(grep -rn -E 'fig_diverse|fig_5model|scripts/targets\.yaml|configs/targets\.yaml|01_prep_9vp0|RUN_ALL\.md|AUTOMATION\.md|benchmark_targets_2026|(^|[^a-z])src/' \
        --include='*.sh' --include='*.py' --include='*.md' . 2>/dev/null \
        | grep -v '^\./archive/' || true)
if [ -n "$LEFT" ]; then echo "  ⚠️ 남은 참조(확인 요):"; echo "$LEFT" | sed 's/^/    /'; else echo "  ✓ 깨진 참조 없음"; fi

echo
echo "======================================================================"
echo " 완료(커밋 안 함). 아래로 검토 후 직접 커밋:"
echo "   git status"
echo "   git commit -m 'consensus_docking 재구성: docs/configs/figures/archive 분리, .pse 제거, README 갱신'"
echo "   git push"
echo "======================================================================"
