#!/usr/bin/env bash
# consensus_docking 정리 v2 (2026-07-14).
# v1은 수정된 .pse에서 git rm 실패(set -e) → 중단됐음. 이번엔 그걸 수정 + 대폭 단순화:
#   · 문서 4개 -> README 1개 (RUN_ALL·AUTOMATION·benchmark_targets 삭제, 핵심은 README에 통합)
#   · 런타임 아티팩트(예측·로그·렌더·클론레포)는 조직화 대신 .gitignore로 숨김 (git status 정리)
#   · 출력 폴더명 out_* 통일 (reader는 suffix 호환) · 설정 configs/ · 레거시 archive/
#   · figures 폴더 통합은 이번엔 안 함(fig/ 가 untracked 렌더 투성이 → gitignore가 맞음)
# git mv/rm 기반. ★ 커밋 안 함 — 실행 후 git status 검토하고 직접 커밋.
# ⚠️ 트래킹 안 된 네 코드(deeprank_batch.py·immunebuilder_rmsd.py·make_colabfold.exp1.py)는 안 건드림.
# 실행: repo 클론에서  cd <repo>/consensus_docking && bash scripts/reorg.sh
set -uo pipefail                       # -e 안 씀 — 한 단계 실패해도 계속(단계별로 처리)
cd "$(dirname "$0")/.." || exit 1
[ -d ../.git ] || { echo "!! git repo 루트(../.git) 없음. repo 클론에서 실행하세요."; exit 1; }
G(){ echo "  + $*"; "$@" || echo "    (경고: 위 명령 실패 — 계속 진행)"; }

echo "== 1) .gitignore 보강: 런타임 아티팩트 전부 숨김 (git status 정리) =="
GI=../.gitignore
add_ignore(){ grep -qxF "$1" "$GI" 2>/dev/null || { echo "$1" >> "$GI"; echo "  + $1"; }; }
add_ignore '*.pse'
add_ignore '*.log'
add_ignore 'consensus_docking/runs/'
add_ignore 'consensus_docking/runs_diverse/'
add_ignore 'consensus_docking/tfold/'
add_ignore 'consensus_docking/DeepRank-GNN-esm/'
add_ignore 'consensus_docking/fig/'
add_ignore 'consensus_docking/results/*.png'

echo "== 2) 트래킹된 .pse 강제 제거 (수정본 포함; pymol_pose.py로 재생성 가능) =="
G git rm -f -q --ignore-unmatch fig_diverse/*.pse fig_5model/*.pse

echo "== 3) 문서 통합: README 하나만 남기고 3개 삭제 (핵심은 아래 README로) =="
G git rm -q --ignore-unmatch RUN_ALL.md AUTOMATION.md benchmark_targets_2026.md

echo "== 3b) 승인된 불필요 코드 삭제 (사용자 승인 2026-07-14; untracked → rm) =="
for f in deeprank_batch.py immunebuilder_rmsd.py scripts/make_colabfold.exp1.py; do
  if [ -e "$f" ]; then G rm -f "$f"; echo "    삭제: $f"; fi
done

echo "== 4) 설정 통합 -> configs/ (옛 9VP0 config는 archive) =="
mkdir -p archive
[ -f configs/targets.yaml ]         && G git mv configs/targets.yaml archive/configs_targets_9vp0.yaml
[ -f scripts/targets_diverse.yaml ] && G git mv scripts/targets_diverse.yaml configs/targets_diverse.yaml
[ -f scripts/targets.yaml ]         && G git mv scripts/targets.yaml configs/targets_hackathon.yaml

echo "== 5) 레거시 스크립트 -> archive/ =="
for f in 01_prep_9vp0.py 02_run_boltz.sh 03_eval.sh render_cases.sh rerun_boltz_failed.sh; do
  [ -f "scripts/$f" ] && G git mv "scripts/$f" "archive/$f"
done

echo "== 6) 출력 폴더명 통일: results_boltz->out_boltz · results_colabfold->out_colabfold =="
# reader(Phase 0 핵심 + eval_boltz) = suffix glob(*boltz/*colabfold)로 서버 rename 전/후 모두 호환.
for f in scripts/diagnose_antibody_side.py scripts/where_docked.py scripts/merge_dockq.py scripts/eval_boltz.py; do
  [ -f "$f" ] && perl -pi -e 's{results_boltz/}{*boltz/}g; s{results_colabfold/}{*colabfold/}g' "$f"
done
# writer(run_diverse)·eval_colabfold = 새 이름(out_*)으로 확정.
[ -f scripts/run_diverse.sh ]    && perl -pi -e 's{results_colabfold}{out_colabfold}g' scripts/run_diverse.sh
[ -f scripts/eval_colabfold.py ] && perl -pi -e 's{results_colabfold}{out_colabfold}g' scripts/eval_colabfold.py
echo "  (reader=suffix 호환 / writer=out_*. 서버 실물 폴더는 migrate_server_outdirs.sh)"

echo "== 7) extract_chains 기본 타깃 경로 -> configs/ =="
[ -f scripts/extract_chains.py ] && perl -pi -e 's{default="targets\.yaml"}{default="configs/targets_diverse.yaml"}g' scripts/extract_chains.py

echo "== 8) README 단일 문서로 재작성 =="
cat > README.md <<'READMEEOF'
# Consensus Docking — 항체-항원 pose 재랭킹 / 오케스트레이션

여러 공개 co-folder(Boltz-2 · Chai-1 · Protenix · AF2-Multimer · tFold-Ag)를 한 파이프라인으로
돌려, 단일 ipTM 대신 **모델 간 합의(consensus)로 native-like pose를 고른다.**
> 이 폴더 문서는 이 README **하나로 통합**했다(옛 RUN_ALL/AUTOMATION/benchmark 내용은 여기로).

## 성공 기준 (2026-07-14 갱신)
- **DockQ ≥ 0.49 (CAPRI Medium)** 주기준 — 3-tier 병기(0.23 Acceptable / 0.49 Medium / 0.80 High).
- **RMSD 이원화**: 항체 내부 CA-RMSD(fold 품질) / Ligand-RMSD(배치 품질) + epitope recall.
- single seed는 노이지 → multi-seed 분포 + **ipTM-vs-oracle**(Win Rate / Hit-Oracle / Mean Regret).

## 벤치마크 세트 (configs/)
- **`targets_diverse.yaml` — 현재 메인**: 2024+ 공개 Fab 10종(다양 항원, leakage-free).
- **`targets_hackathon.yaml`**: 9VP0(골든/dev) + 2026 공개 5종(10ZO·9Y0A·9OBN·9K6J·9XF9; 나노바디 2 포함).
- 교체 후보(필요시): 9PIU(Fab+mTOR FRB) · 9UOI(나노바디+MERS ternary) · 9SAT/9MNT · 큰 항원 9VDY/9TH6.

## leakage / 사슬 역할 규칙 (사수 보고용 핵심)
- leakage 판정 = RCSB `initial_release_date`. **SAbDab 날짜는 deposit 미러라 부정확 → 쓰지 말 것.**
- 사슬 역할(H/L/항원)·형식(Fab/나노바디) = SAbDab 표기. 다중카피(10ZO·9OBN 등)는 1:1 한 쌍만 사용.

## 폴더 구조
- `README.md` — 이 문서(단일 소스)
- `configs/`  — 타깃 세트 yaml (diverse=메인 · hackathon)
- `scripts/`  — 코드(아래 파이프라인)
- `results/`  — 커밋되는 `.txt` 산출(diagnose_* · where_docked_* · merged_*). png·csv·log는 gitignore
- `archive/`  — 레거시(9VP0 브링업 01/02/03 · render_cases · rerun_boltz_failed · 옛 config)
- 아티팩트(`runs/`·`runs_diverse/`·`tfold/`·`DeepRank-GNN-esm/`·`fig*/`·`*.pse`·`*.log`) = 전부 gitignore

## 파이프라인
1. `extract_chains.py` — native 다운 + keep_chains 서열 → `runs_diverse/<id>/chains.json`
2. `make_fv.py` — 항체 ANARCI(IMGT)로 Fv(VH+VL) 절단 → `runs_diverse/<id>_fv/`
3. `make_{boltz,chai,protenix,colabfold,tfold}.py` — 모델별 입력 생성
4. 예측: `run_diverse.sh`(chai·protenix·colabfold) · `run_tfold.sh` · `run_5model.sh`
   → 각 `runs_diverse/<id>/{out_boltz,out_chai,out_protenix,out_colabfold,out_tfold}`
5. eval: `eval_*.py` → `results/summary.csv`(source of truth) → `sync_notion.py`(개인 Mac, 토큰 env)
6. 분석: `consensus_analysis.py`·`plot_iptm_dockq.py`(집계) /
   `diagnose_antibody_side.py`·`where_docked.py`·`merge_dockq.py`(정밀 진단)
7. 렌더: `pymol_pose.py` + `render_*.sh` → `fig*/`

## Phase 0 (go/no-go 게이트)
`run_phase0.sh` = diagnose + where_docked + merge_dockq 를 10-타깃 × Fab/Fv 전부에.
새 예측 없이 기존 pose 재분석 → 0.49 기준 oracle·placement·ipTM 품질로 프로젝트 go/no-go 판정.

## 규칙
- 출력 폴더명(`runs_diverse`·`out_*`)은 스크립트에 하드코딩(reader는 `*boltz`/`*colabfold` 호환).
- 예측 구조·MSA·weight·CSV·`.pse`·클론 툴(tfold·DeepRank)은 git에 안 올림.
- Notion 토큰은 서버·깃·문서 저장 금지(개인 Mac env만).
READMEEOF
echo "  README.md 재작성 완료"

echo "== 9) v1이 만든 빈 폴더 정리 =="
rmdir docs figures 2>/dev/null || true

echo "== 10) 스테이징 =="
G git add -A ../.gitignore .

echo
echo "== 11) 자기점검: 옛 경로가 (archive 제외) 아직 참조되나 =="
LEFT=$(grep -rn -E 'scripts/targets\.yaml|configs/targets\.yaml|RUN_ALL\.md|AUTOMATION\.md|benchmark_targets_2026\.md|results_boltz|results_colabfold' \
        --include='*.sh' --include='*.py' --include='*.md' . 2>/dev/null | grep -v '^\./archive/' || true)
[ -n "$LEFT" ] && { echo "  ⚠️ 남은 참조(확인 요):"; echo "$LEFT" | sed 's/^/    /'; } || echo "  ✓ 깨진 참조 없음"

echo
echo "== 12) 참고: 아직 남은 대형 아티팩트 (gitignore됨; 서버 디스크 정리는 별도) =="
for d in DeepRank-GNN-esm tfold runs runs_diverse fig fig_diverse fig_5model; do
  [ -e "$d" ] && echo "  · $d/  (git엔 안 올라감. fig* 삭제/정리는 ls 확인 후 별도 안내)"
done

echo
echo "======================================================================"
echo " 완료(커밋 안 함). 검토 후 직접 커밋:"
echo "   git status"
echo "   git commit -m 'consensus_docking 정리: 문서 1개 통합, 아티팩트 gitignore, out_* 통일, 레거시 archive'"
echo "   git push"
echo
echo " ★ 서버 실물 예측 폴더 rename(선택): bash scripts/migrate_server_outdirs.sh [--apply]"
echo "   (reader가 양쪽 호환이라 안 해도 Phase 0 정상 동작)"
echo "======================================================================"
