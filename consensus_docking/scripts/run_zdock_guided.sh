#!/usr/bin/env bash
# ZDOCK 물리 pose → co-folder guided 재실행 (diverse). HADDOCK-guided(run_*_guided.sh)의 ZDOCK·diverse판.
# 흐름: zdock_to_pocket.py(에피토프 제약 생성) → co-folder 재실행(blind 커맨드 + constraint만 차이) → dockq_guided_diverse.py 채점.
#   blind baseline = runs_diverse/<t>/out_<model> (기존 5모델 예측). 순수 '힌트 효과' 격리.
#
# 사용(consensus_docking/에서):
#   SMOKE=1 bash scripts/run_zdock_guided.sh boltz 9KKJ                      # 스모크(1패치 1샘플, 포맷·플러밍 검증)
#   bash scripts/run_zdock_guided.sh boltz "9J4C 9KKJ 9Y0A 9W43"            # 전체(topn1 압축)
#   GEN_ARGS="--topn 10 --union" bash scripts/run_zdock_guided.sh boltz "9J4C 9KKJ 9Y0A 9W43"   # 넓게
#   bash scripts/run_zdock_guided.sh protenix "..."   /   bash scripts/run_zdock_guided.sh tfold "..."
# env override: BOLTZ_ENV PROT_ENV TFOLD_ENV / PROT_MODEL / SAMP(=diffusion/samples) / SEEDS(tfold) / TFOLD_DIR
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1; ROOT="$(pwd)"
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
MODEL="${1:-}"; TARGETS="${2:-}"
[ -z "$MODEL" ] || [ -z "$TARGETS" ] && { echo "사용: bash scripts/run_zdock_guided.sh <boltz|protenix|tfold> \"9J4C 9KKJ ...\""; exit 1; }
SMOKE="${SMOKE:-0}"; OUTDIR="${OUTDIR:-runs_zdock_guided}"; RUN="${RUN:-runs_diverse}"
SAMP="${SAMP:-5}"; SEEDS="${SEEDS:-42 43 44 45 46}"; GEN_ARGS="${GEN_ARGS:-}"
BOLTZ_ENV="${BOLTZ_ENV:-boltz}"; PROT_ENV="${PROT_ENV:-protenix}"; TFOLD_ENV="${TFOLD_ENV:-tfold}"
PROT_MODEL="${PROT_MODEL:-protenix_base_default_v1.0.0}"   # ⚠️ blind diverse는 leaky 20250630 사용 → 깨끗한 delta는 blind도 v1.0.0 재실행 필요(아래 경고)
SEED="${SEED:-101}"; TFOLD_DIR="${TFOLD_DIR:-$ROOT/tfold}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
if [ "$SMOKE" = "1" ]; then SAMP=1; SEEDS="42"; [ -z "$GEN_ARGS" ] && GEN_ARGS="--topn 1 --max-patches 1"; fi

# 1) ZDOCK pose → guided 입력 생성
say "[gen] zdock_to_pocket.py --emit $MODEL --targets \"$TARGETS\" $GEN_ARGS"
conda activate "$BOLTZ_ENV" 2>/dev/null   # biopython/torch 있는 env(대개 boltz). tfold emit는 torch 필요 → tfold면 아래서 재생성 가능
python scripts/zdock_to_pocket.py --targets "$TARGETS" --emit "$MODEL" --outdir "$OUTDIR" --run-dir "$RUN" $GEN_ARGS \
  || { say "!! 입력 생성 실패 (tfold면 torch 있는 env 필요 — TFOLD emit는 tfold env에서 먼저 생성)"; exit 1; }

[ "$MODEL" = "protenix" ] && say "⚠️ PROT_MODEL=$PROT_MODEL. blind diverse는 protenix_base_20250630_v1.0.0(leaky)이라 모델버전 confound — 깨끗한 비교는 blind도 v1.0.0로 재생성 필요."

# 2) 모델별 실행
case "$MODEL" in
  boltz)
    conda activate "$BOLTZ_ENV" 2>/dev/null; command -v boltz >/dev/null || { say "!! boltz 없음"; exit 1; } ;;
  protenix)
    conda activate "$PROT_ENV" 2>/dev/null
    export PROTENIX_ROOT_DIR=/mnt/data/admuser/protenix_weights LAYERNORM_TYPE=torch
    command -v protenix >/dev/null || { say "!! protenix 없음"; exit 1; } ;;
  tfold)
    conda activate "$TFOLD_ENV" 2>/dev/null
    python -c "import tfold,torch" 2>/dev/null || say "!! tfold/torch import 경고(env $TFOLD_ENV) — 계속" ;;
  *) say "unknown MODEL '$MODEL'"; exit 1 ;;
esac

for T in $TARGETS; do
  t=$(echo "$T" | tr 'A-Z' 'a-z'); d="$OUTDIR/$t"
  [ -d "$d" ] || { say "skip $T (guided 입력 없음)"; continue; }
  case "$MODEL" in
    boltz)
      for y in "$d"/boltz_${t}_patch*.yaml; do
        [ -e "$y" ] || continue; pk=$(basename "$y" .yaml | sed 's/.*_//'); od="results_$pk"
        if find "$d/$od" -name '*_model_*.cif' 2>/dev/null | grep -q .; then say "skip $T $pk (있음)"; continue; fi
        say "boltz-guided $T $pk (samples=$SAMP) ..."
        ( cd "$d" && boltz predict "$(basename "$y")" --out_dir "$od" --use_msa_server --diffusion_samples "$SAMP" > "boltz_$pk.log" 2>&1 ) \
          && say "  OK $T $pk" || { say "  !! $T $pk 실패 → $d/boltz_$pk.log"; tail -5 "$d/boltz_$pk.log"; }
        [ "$SMOKE" = "1" ] && { say "SMOKE 종료. 확인: find $d/$od -name '*_model_*.cif'"; exit 0; }
      done ;;
    protenix)
      for j in "$d"/protenix_${t}_patch*.json; do
        [ -e "$j" ] || continue; pk=$(basename "$j" .json | sed 's/.*_//'); od="results_$pk"
        if find "$d/$od" -name '*sample*.cif' 2>/dev/null | grep -q .; then say "skip $T $pk (있음)"; continue; fi
        say "protenix-guided $T $pk (samples=$SAMP, $PROT_MODEL) ..."
        ( cd "$d" && protenix pred -i "$(basename "$j")" -o "$od" -n "$PROT_MODEL" -s "$SEED" -e "$SAMP" \
            --trimul_kernel torch --triatt_kernel torch > "protenix_$pk.log" 2>&1 ) \
          && say "  OK $T $pk" || { say "  !! $T $pk 실패 → $d/protenix_$pk.log"; tail -5 "$d/protenix_$pk.log"; }
        [ "$SMOKE" = "1" ] && { say "SMOKE 종료. 확인: find $d/$od -name '*sample*.cif'"; exit 0; }
      done ;;
    tfold)
      fasta="$ROOT/$d/tfold_${t}.fasta"
      a3m="$ROOT/$RUN/$t/msa_$t/A.a3m"; [ -f "$a3m" ] || a3m=$(find "$ROOT/$RUN/$t/msa_$t" -maxdepth 1 -name '*.a3m' 2>/dev/null | head -1)
      [ -f "$fasta" ] && [ -n "$a3m" ] || { say "skip $T (fasta/a3m 없음: $a3m)"; continue; }
      for icf in "$ROOT/$d"/tfold_${t}_patch*.pt; do
        [ -e "$icf" ] || continue; pk=$(basename "$icf" .pt | sed 's/.*_//'); od="$ROOT/$d/out_tfold_$pk"; mkdir -p "$od"
        for s in $SEEDS; do
          outpdb="$od/${t}_seed${s}.pdb"; [ -s "$outpdb" ] && { say "skip $T $pk seed$s (있음)"; continue; }
          say "tfold-guided $T $pk seed$s ..."
          ( cd "$TFOLD_DIR" && python projects/tfold_ag/predict.py --fasta "$fasta" --msa "$a3m" --icf "$icf" \
              --output "$outpdb" --model_version ppi --seed "$s" > "$od/${t}_seed${s}.log" 2>&1 ) \
            && say "  OK $T $pk seed$s" || { say "  !! $T $pk seed$s 실패 → $od/${t}_seed${s}.log"; tail -3 "$od/${t}_seed${s}.log"; }
          [ "$SMOKE" = "1" ] && { say "SMOKE 종료. 확인: $outpdb"; exit 0; }
        done
      done ;;
  esac
done
say "완료 → $OUTDIR/<t>/. 채점: python scripts/dockq_guided_diverse.py --model $MODEL --guided-dir $OUTDIR --targets \"$TARGETS\""
