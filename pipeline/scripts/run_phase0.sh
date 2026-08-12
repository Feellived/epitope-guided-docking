#!/usr/bin/env bash
# Phase 0 재분석 — 새 예측 없이 '기존 예측물'만 재분석 (diagnose/where_docked/merge_dockq).
#   (1) diagnose_antibody_side -> results/diagnose_<T>.txt   epi_recall + 항체 CA-RMSD (fold vs 배치)
#   (2) where_docked           -> results/where_docked_<T>.txt  majority-vote 예측 에피토프(+RBD 영역)
#   (3) merge_dockq            -> results/merged_all.txt      CAPRI merged DockQ (H+L 병합)
# 전제: <RUN>/<id>{,_fv}/ 에 모델 pose + native + chains.json 존재. GPU 불필요, 읽기+DockQ만.
#
# 사용:
#   bash scripts/run_phase0.sh                                              # runs_diverse 10 (Fab,Fv)
#   ONLY=8P5M,8SDF,8SDH,8SIQ,8SIS,8SIT,8XSI,9ML8,9ML9,9SBB,9ZDU RUN=runs_rbd FORMATS=Fab bash scripts/run_phase0.sh
#   SKIP_MERGE=1 ... (느린 merge_dockq 빼고 진단 2개만)
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
ENV="${ENV:-boltz}"
conda activate "$ENV" 2>/dev/null
RUN="${RUN:-runs_diverse}"; export RUN_DIR="$RUN"     # 스크립트들이 RUN_DIR로 run-dir 인식
FMT="${FORMATS:-Fab,Fv}"

T="${ONLY:-9K6J,9XF9,9Y0A,9VXL,9LY5,9MQR,9KKJ,9W43,9YC6,9J4C}"
mkdir -p results
LOG="results/phase0.log"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

say "Phase 0 — run=$RUN | targets: $T | env=$ENV | formats=$FMT"

say "[1/3] diagnose_antibody_side (fold vs 배치)"
python scripts/analyze_antibody_side.py --targets "$T" --formats "$FMT" 2>&1 | tee -a "$LOG"

say "[2/3] where_docked (어디에 놓나 / 예측 에피토프)"
python scripts/analyze_where_docked.py --targets "$T" --formats "$FMT" 2>&1 | tee -a "$LOG"

if [ -z "${SKIP_MERGE:-}" ]; then
  say "[3/3] merge_dockq (CAPRI merged DockQ; 느림)"
  python scripts/analyze_merged_dockq.py --targets "$T" 2>&1 | tee results/merged_all.txt
  cat results/merged_all.txt >> "$LOG"
else
  say "[3/3] merge_dockq SKIP (SKIP_MERGE set)"
fi

say "Phase 0 완료 → results/diagnose_*.txt · where_docked_*.txt · merged_all.txt"
echo
echo "==== 회수(아래를 실행해 출력을 붙여넣어 주세요) ===="
echo "cat results/merged_all.txt"
echo "for t in $(echo $T | tr ',' ' '); do echo \"== where \$t ==\"; cat results/where_docked_\$t.txt 2>/dev/null; done"
echo "for t in $(echo $T | tr ',' ' '); do echo \"== diagnose \$t ==\"; cat results/diagnose_\$t.txt 2>/dev/null; done"
