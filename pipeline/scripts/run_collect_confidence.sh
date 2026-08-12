#!/usr/bin/env bash
# [나중에 돌릴 배치] 모델이 낼 수 있는데 우리 기존 예측물엔 아직 없는 confidence 재수집.
#   · AF2-M  pair-ipTM  : colabfold_batch --calc-extra-ptm  -> results_colabfold_xtra/
#   · Protenix PAE       : protenix pred (PAE 플래그 탐지)   -> out_protenix_pae/
#   · Chai    PAE        : chai-lab fold (PAE 옵션 탐지)      -> out_chai_pae/
#   (Boltz = pair+PAE 이미 있음 / tFold = 더 낼 것 없음 → 스킵)
#
# ★ 안전: Phase 0 출력(out_chai/out_protenix/results_colabfold)은 절대 안 건드림 — 전부 별도 companion 디렉토리.
# ★ 견고: 플래그가 버전마다 달라, 실행 전 `--help`로 탐지 → 있으면 그 플래그로 실행, 없으면 명확히 SKIP+TODO(헛돌지 않음).
# ★ resumable: companion 디렉토리에 출력 있으면 스킵.
#
# 사용(집 갈 때, GPU 큼 → tmux):
#   tmux new -s xtra
#   ONLY=9k6j bash scripts/run_collect_confidence.sh      # ← 먼저 스모크(어떤 플래그가 있나 확인)
#   bash scripts/run_collect_confidence.sh                # 전체(diverse 10 × Fab/Fv)
#   MODELS="protenix" bash scripts/run_collect_confidence.sh   # 특정 모델만
cd "$(dirname "$0")/.." || exit 1
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
RUN=runs_diverse
MODELS="${MODELS:-colabfold protenix chai}"
CHAI_ENV="${CHAI_ENV:-chai}"; PROT_ENV="${PROT_ENV:-protenix}"; CF_ENV="${CF_ENV:-colabfold}"
PROT_MODEL="${PROT_MODEL:-protenix_base_20250630_v1.0.0}"
LOG="results/collect_full_confidence.log"; mkdir -p results
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

if [ -n "${ONLY:-}" ]; then dirs="$RUN/${ONLY} $RUN/${ONLY}_fv"
else dirs=$(ls -d $RUN/*/ 2>/dev/null | sed 's:/$::'); fi

# help에서 정규식으로 플래그 1개 탐지(첫 매치 반환, 없으면 빈 문자열)
detect(){ "$@" --help 2>&1 | grep -oiE -- "$FLAG_RE" | head -1; }

for model in $MODELS; do
  say "================ MODEL: $model ================"
  case $model in
    colabfold)
      conda activate "$CF_ENV" 2>/dev/null
      command -v colabfold_batch >/dev/null || { say "!! colabfold_batch 없음 → 스킵"; continue; }
      FLAG_RE='--calc-extra-ptm|--calc_extra_ptm'
      flag=$(colabfold_batch --help 2>&1 | grep -oiE -- "$FLAG_RE" | head -1)
      if [ -z "$flag" ]; then say "!! colabfold에 pair-ipTM 플래그(--calc-extra-ptm) 없음 → SKIP(버전업 필요)"; continue; fi
      say "  플래그 감지: $flag  (pair-ipTM: pairwise_iptm 출력)"
      for d in $dirs; do [ -d "$d" ] || continue
        base=$(basename "$d"); t=$(echo "${base%_fv}" | tr 'A-Z' 'a-z'); out="$d/results_colabfold_xtra"
        find "$out" -name '*_scores_*.json' 2>/dev/null | grep -q . && { say "skip $base (있음)"; continue; }
        say "run $base ..."
        ( cd "$d" && env -u LD_LIBRARY_PATH colabfold_batch "colabfold_$t.fasta" results_colabfold_xtra \
            --num-models 5 --model-type alphafold2_multimer_v3 --num-recycle 3 $flag > colabfold_xtra.log 2>&1 )
      done ;;
    protenix)
      conda activate "$PROT_ENV" 2>/dev/null
      command -v protenix >/dev/null || { say "!! protenix 없음 → 스킵"; continue; }
      export PROTENIX_ROOT_DIR="${PROTENIX_ROOT_DIR:-/mnt/data/admuser/protenix_weights}" LAYERNORM_TYPE=torch
      flag=$(protenix pred --help 2>&1 | grep -oiE -- '--write_pae|--save_pae|--output_pae|--write-pae' | head -1)
      if [ -z "$flag" ]; then say "!! protenix에 PAE 저장 플래그 없음 → SKIP (protenix pred --help 확인 필요; base 미지원일 수 있음)"; continue; fi
      say "  플래그 감지: $flag  (PAE 매트릭스 출력)"
      for d in $dirs; do [ -d "$d" ] || continue
        base=$(basename "$d"); t=$(echo "${base%_fv}" | tr 'A-Z' 'a-z'); out="$d/out_protenix_pae"
        find "$out" -name '*pae*' 2>/dev/null | grep -q . && { say "skip $base (있음)"; continue; }
        say "run $base ..."
        ( cd "$d" && protenix pred -i "protenix_$t.json" -o out_protenix_pae -n "$PROT_MODEL" -s 101 -e 5 \
            --trimul_kernel torch --triatt_kernel torch $flag > protenix_pae.log 2>&1 )
      done ;;
    chai)
      conda activate "$CHAI_ENV" 2>/dev/null
      command -v chai-lab >/dev/null || { say "!! chai-lab 없음 → 스킵"; continue; }
      export CHAI_DOWNLOADS_DIR="${CHAI_DOWNLOADS_DIR:-/mnt/data/admuser/chai_downloads}"
      flag=$(chai-lab fold --help 2>&1 | grep -oiE -- '--save-pae|--write-pae|--output-pae|--pae' | head -1)
      if [ -z "$flag" ]; then
        say "!! chai-lab fold에 PAE 저장 옵션 없음 → SKIP. Chai PAE는 CLI로 안 나옴(issue#348):"
        say "   → run_inference API의 StructureCandidates.pae를 파이썬으로 덤프해야 함(버전 확인 후 chai_pae.py 별도 제공)."
        continue
      fi
      say "  플래그 감지: $flag"
      for d in $dirs; do [ -d "$d" ] || continue
        base=$(basename "$d"); t=$(echo "${base%_fv}" | tr 'A-Z' 'a-z'); out="$d/out_chai_pae"
        find "$out" -name '*.npz' 2>/dev/null | grep -q . && { say "skip $base (있음)"; continue; }
        say "run $base ..."
        ( cd "$d" && chai-lab fold "chai_$t.fasta" out_chai_pae --use-msa-server $flag > chai_pae.log 2>&1 )
      done ;;
    *) say "unknown model $model" ;;
  esac
done
say "======== DONE ========  companion 출력: results_colabfold_xtra / out_protenix_pae / out_chai_pae"
say "각 모델 SKIP 메시지 있으면 그 플래그는 이 버전에 없음 = 별도 처리 필요(알려주세요)."
