#!/usr/bin/env bash
# Phase 0 — go/no-go 게이트: 새 co-folder 예측 없이 '기존 5모델 예측물'만 재분석한다.
# 3개 진단을 10 타깃 × Fab/Fv 전부에 실행:
#   (1) diagnose_antibody_side.py -> results/diagnose_<T>.txt      epi_recall + 항체 CA-RMSD
#                                                                  = 실패가 'fold 문제'냐 '배치 문제'냐 분해
#   (2) where_docked.py           -> results/where_docked_<T>.txt  majority-vote 예측 에피토프
#                                                                  = 실패 모델들이 공통 자리에 쏠리나(prior) 흩어지나
#   (3) merge_dockq.py            -> results/merged_all.txt         CAPRI merged DockQ (H+L 병합)
#                                                                  = 0.49(Medium) 기준 oracle 재판정용
#
# 전제: runs_diverse/<id>{,_fv}/ 에 5모델 pose(.cif/.pdb) + native + chains.json 이미 존재(예측 끝난 상태).
#       ★ 새 예측/설치/삭제 없음. 파일 읽기 + DockQ 서브프로세스만. GPU 불필요, 서버 안전.
#
# 사용(서버; merge_dockq는 DockQ 수백 회 호출 → 수십 분, tmux 권장):
#   tmux new -s phase0
#   conda activate boltz          # Bio.PDB + DockQ 있는 env
#   cd ~/projects/bk21-antibody-ml/consensus_docking
#   bash scripts/run_phase0.sh                    # 10 타깃 전부
#   ONLY=9MQR,9W43 bash scripts/run_phase0.sh     # 일부만
#   SKIP_MERGE=1 bash scripts/run_phase0.sh       # 느린 merge_dockq 빼고 진단 2개만 먼저
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
ENV="${ENV:-boltz}"
conda activate "$ENV" 2>/dev/null

T="${ONLY:-9K6J,9XF9,9Y0A,9VXL,9LY5,9MQR,9KKJ,9W43,9YC6,9J4C}"
mkdir -p results
LOG="results/phase0.log"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

say "Phase 0 시작 — targets: $T (env=$ENV)"

say "[1/3] diagnose_antibody_side (fold vs 배치)  -> results/diagnose_<T>.txt"
python scripts/diagnose_antibody_side.py --targets "$T" --formats Fab,Fv 2>&1 | tee -a "$LOG"

say "[2/3] where_docked (어디에 놓나 / majority-vote 에피토프)  -> results/where_docked_<T>.txt"
python scripts/where_docked.py --targets "$T" --formats Fab,Fv 2>&1 | tee -a "$LOG"

if [ -z "${SKIP_MERGE:-}" ]; then
  say "[3/3] merge_dockq (CAPRI merged DockQ; 느림)  -> results/merged_all.txt"
  python scripts/merge_dockq.py --targets "$T" 2>&1 | tee results/merged_all.txt
  cat results/merged_all.txt >> "$LOG"
else
  say "[3/3] merge_dockq SKIP (SKIP_MERGE set)"
fi

say "Phase 0 완료. 산출: results/diagnose_*.txt · where_docked_*.txt · merged_all.txt"
echo
echo "==== 회수(아래 3줄을 실행해 출력 전체를 붙여넣어 주세요) ===="
echo "cat results/merged_all.txt"
echo "for t in 9K6J 9XF9 9Y0A 9VXL 9LY5 9MQR 9KKJ 9W43 9YC6 9J4C; do echo \"== diagnose \$t ==\"; cat results/diagnose_\$t.txt 2>/dev/null; done"
echo "for t in 9K6J 9XF9 9Y0A 9VXL 9LY5 9MQR 9KKJ 9W43 9YC6 9J4C; do echo \"== where \$t ==\"; cat results/where_docked_\$t.txt 2>/dev/null; done"
