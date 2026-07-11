#!/usr/bin/env bash
# 다양-항원 Fab+Fv 세트에 co-folder(Chai-1·Protenix) 순차 예측 (unattended·resumable·logged).
# 예측만 함 (eval은 아침에 별도 — Chai/Protenix DockQ 사슬매핑 확인 필요).
#
# 사용:
#   tmux new -s night
#   bash scripts/run_diverse.sh                 # chai+protenix 전체(Fab+Fv 20세트)
#   ONLY=9k6j bash scripts/run_diverse.sh       # 스모크: 9k6j Fab+Fv만
#   bash scripts/run_diverse.sh chai            # chai만
#   CHAI_ENV=mychai PROT_ENV=myprot bash scripts/run_diverse.sh   # env 이름 다르면 override
#
# GPU 1개 가정 → 순차. 출력(cif) 있으면 스킵하므로 중단 후 다시 실행하면 이어받음.
cd "$(dirname "$0")/.." || exit 1
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null

RUN=runs_diverse
MODELS="${*:-chai protenix}"
CHAI_ENV="${CHAI_ENV:-chai}"
PROT_ENV="${PROT_ENV:-protenix}"
CF_ENV="${CF_ENV:-colabfold}"          # AF2-Multimer(ColabFold) env — 실제 이름 확인 필요
PROT_MODEL="${PROT_MODEL:-protenix_base_20250630_v1.0.0}"
LOG="run_diverse.log"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

if [ -n "${ONLY:-}" ]; then
  dirs="$RUN/${ONLY} $RUN/${ONLY}_fv"
else
  dirs=$(ls -d $RUN/*/ 2>/dev/null | sed 's:/$::')
fi

# --- 입력 생성 (없으면 덮어씀, 무해) ---
say "입력 생성 (make_chai / make_protenix)"
conda activate boltz 2>/dev/null
python scripts/make_chai.py      --run-dir "$RUN" >> "$LOG" 2>&1
python scripts/make_protenix.py  --run-dir "$RUN" >> "$LOG" 2>&1
python scripts/make_colabfold.py --run-dir "$RUN" >> "$LOG" 2>&1

for model in $MODELS; do
  say "================ MODEL: $model ================"
  case $model in
    chai)
      conda activate "$CHAI_ENV" 2>/dev/null
      export CHAI_DOWNLOADS_DIR=/mnt/data/admuser/chai_downloads
      command -v chai-lab >/dev/null || { say "!! chai-lab 없음 (env '$CHAI_ENV' 확인) → 스킵"; continue; }
      out=out_chai ;;
    protenix)
      conda activate "$PROT_ENV" 2>/dev/null
      export PROTENIX_ROOT_DIR=/mnt/data/admuser/protenix_weights
      export LAYERNORM_TYPE=torch
      command -v protenix >/dev/null || { say "!! protenix 없음 (env '$PROT_ENV' 확인) → 스킵"; continue; }
      out=out_protenix ;;
    colabfold)   # AF2-Multimer — ⚠️ 초안(플래그 미검증), 내일 스모크 후 사용
      conda activate "$CF_ENV" 2>/dev/null
      command -v colabfold_batch >/dev/null || { say "!! colabfold_batch 없음 (env '$CF_ENV' 확인) → 스킵"; continue; }
      out=results_colabfold ;;
    *) say "unknown model $model"; continue ;;
  esac

  for d in $dirs; do
    [ -d "$d" ] || continue
    base=$(basename "$d"); t=$(echo "${base%_fv}" | tr 'A-Z' 'a-z')
    if find "$d/$out" \( -name '*.cif' -o -name '*.pdb' \) 2>/dev/null | grep -q .; then
      say "skip $base ($model 이미 있음)"; continue
    fi
    say "run  $base ($model) ..."
    (
      cd "$d" || exit 1
      case "$model" in
        chai)
          chai-lab fold "chai_$t.fasta" out_chai --use-msa-server > chai.log 2>&1 ;;
        protenix)
          rm -rf out_protenix
          protenix pred -i "protenix_$t.json" -o out_protenix \
            -n "$PROT_MODEL" -s 101 -e 5 \
            --trimul_kernel torch --triatt_kernel torch > protenix.log 2>&1 ;;
        colabfold)   # AF2-Multimer (colabfold 1.6.1) — 플래그 --help로 확정
          env -u LD_LIBRARY_PATH colabfold_batch "colabfold_$t.fasta" results_colabfold \
            --num-models 5 --model-type alphafold2_multimer_v3 --num-recycle 3 > colabfold.log 2>&1 ;;
      esac
    )
    if find "$d/$out" \( -name '*.cif' -o -name '*.pdb' \) 2>/dev/null | grep -q .; then
      say "  OK  $base ($model)"
    else
      say "  FAIL $base ($model) — $d/$model.log 확인"
    fi
  done
done
say "======== ALL DONE ========  (eval은 아침에: 각 dir의 out_chai/out_protenix)"
