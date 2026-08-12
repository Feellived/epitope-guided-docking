#!/usr/bin/env bash
# 지금 돌고 있는 run_diverse(chai/protenix/colabfold)가 끝나면 → Boltz → tFold 순차 실행.
# GPU 안 겹치게 대기(run_diverse.sh 프로세스 사라질 때까지) 후 시작. 이미 된 건 스킵(이어받기).
# 새 tmux 창에서:  RUN=runs_rbd bash scripts/queue_boltz_tfold.sh
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
RUN="${RUN:-runs_rbd}"
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }

# --- run_diverse 끝날 때까지 대기 ([r] = pgrep 자기자신 매칭 방지) ---
say "run_diverse(chai/protenix/colabfold) 끝나길 대기 ... (60초 간격 확인)"
while pgrep -f '[r]un_diverse.sh' >/dev/null; do sleep 60; done
say "run_diverse 완료 감지 → Boltz 시작"

# --- Boltz ---
conda activate boltz 2>/dev/null
for d in "$RUN"/*/; do
  t=$(basename "$d")
  [ -f "$d/boltz_$t.yaml" ] || { say "skip boltz $t (yaml 없음)"; continue; }
  if find "$d/results_boltz" -name '*_model_*.cif' 2>/dev/null | grep -q .; then say "skip boltz $t (이미 있음)"; continue; fi
  say "boltz $t ..."
  ( cd "$d" && boltz predict "boltz_$t.yaml" --use_msa_server --out_dir results_boltz --no_kernels --diffusion_samples 5 > boltz.log 2>&1 ) \
    && say "  OK boltz $t" || say "  !! boltz $t 실패 → $d/boltz.log"
done
say "Boltz 완료 → tFold 시작"

# --- tFold (Phase1 a3m + sanitize + Phase2 predict; 자체 순차) ---
RUN="$RUN" bash scripts/run_tfold.sh
say "======== 전부 완료: chai/protenix/colabfold + Boltz + tFold ========"
