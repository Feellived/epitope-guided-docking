#!/usr/bin/env bash
# HADDOCK 물리 pose → co-folder guided 재실행 (diverse). run_zdock_guided.sh의 HADDOCK판.
# 흐름: make_pocket_from_haddock.py(에피토프 제약 생성) → co-folder 재실행(blind 커맨드 + constraint만 차이) → eval_guided_diverse.py 채점.
#   blind baseline = runs_diverse/<t>/out_<model> (기존 5모델 예측) → 순수 '힌트 효과' 격리.
#
# ⭐ 힌트 생성 방침(왜 이런지는 make_pocket_from_haddock.py 상단 docstring에 근거 상세):
#   기본 = 1_rigidbody(점수 컷 이전) 에서 '서로 다른 자리 전부'. 클러스터 대표 아님.
#   이유 요약: Exp4에서 (a) best pose 9/10이 singleton이라 클러스터 모드가 배제 (b) seletop 점수컷이
#   near-native 학살(9ML9 0.99→0.41). 힌트엔 '신뢰도'가 아니라 '커버리지'가 필요하므로 두 컷 다 우회.
#
# 사용(pipeline/에서):
#   SMOKE=1 bash scripts/run_haddock_guided.sh boltz 9K6J                # 스모크(1패치 1샘플)
#   bash scripts/run_haddock_guided.sh boltz                             # diverse 10 전체
#   GEN_ARGS="--stage 4_emref --max-patches 20" bash scripts/run_haddock_guided.sh boltz   # 정제본에서 뽑기
#   bash scripts/run_haddock_guided.sh protenix   /   bash scripts/run_haddock_guided.sh tfold
# env: BOLTZ_ENV PROT_ENV TFOLD_ENV / PROT_MODEL / SAMP / SEEDS / TFOLD_DIR / OUTDIR / GEN_ARGS
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1; ROOT="$(pwd)"
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
MODEL="${1:-}"; TARGETS="${2:-9K6J 9XF9 9Y0A 9VXL 9LY5 9MQR 9KKJ 9W43 9YC6 9J4C}"   # 생략=diverse 10 전체
[ -z "$MODEL" ] && { echo "사용: bash scripts/run_haddock_guided.sh <boltz|protenix|tfold> [\"9K6J ...\"]  (타깃 생략=diverse 10)"; exit 1; }
SMOKE="${SMOKE:-0}"; OUTDIR="${OUTDIR:-runs_haddock_guided}"; RUN="${RUN:-runs_diverse}"
SAMP="${SAMP:-5}"; SEEDS="${SEEDS:-42 43 44 45 46}"; GEN_ARGS="${GEN_ARGS:-}"
BOLTZ_ENV="${BOLTZ_ENV:-boltz}"; PROT_ENV="${PROT_ENV:-protenix}"; TFOLD_ENV="${TFOLD_ENV:-tfold}"
PROT_MODEL="${PROT_MODEL:-protenix_base_20250630_v1.0.0}"   # diverse=clean(타깃 2025-09~2026 > 컷오프) + blind와 동일 모델
SEED="${SEED:-101}"; TFOLD_DIR="${TFOLD_DIR:-$ROOT/tfold}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
if [ "$SMOKE" = "1" ]; then SAMP=1; SEEDS="42"; [ -z "$GEN_ARGS" ] && GEN_ARGS="--max-patches 1 --pool 50"; fi
[ -z "$GEN_ARGS" ] && GEN_ARGS="--max-patches 20"   # 기본=rigidbody 전체 pose에서 서로 다른 자리(최대 20)

# 1) HADDOCK pose → guided 입력 생성
say "[gen] make_pocket_from_haddock.py --emit $MODEL --targets \"$TARGETS\" $GEN_ARGS"
conda activate "$BOLTZ_ENV" 2>/dev/null   # biopython 있는 env. tfold emit는 torch 필요 → tfold env 권장
python scripts/make_pocket_from_haddock.py --targets "$TARGETS" --emit "$MODEL" --outdir "$OUTDIR" --run-dir "$RUN" $GEN_ARGS \
  || { say "!! 입력 생성 실패 (HADDOCK 도킹 끝났나? haddock_diverse/<t>/run/1_rigidbody 확인)"; exit 1; }

# 2) 모델별 실행 (blind 커맨드와 동일, constraint만 차이)
case "$MODEL" in
  boltz)    conda activate "$BOLTZ_ENV" 2>/dev/null; command -v boltz >/dev/null || { say "!! boltz 없음"; exit 1; } ;;
  protenix) conda activate "$PROT_ENV" 2>/dev/null
            export PROTENIX_ROOT_DIR=/mnt/data/admuser/protenix_weights LAYERNORM_TYPE=torch
            command -v protenix >/dev/null || { say "!! protenix 없음"; exit 1; }
            say "PROT_MODEL=$PROT_MODEL (blind과 동일 모델)" ;;
  tfold)    conda activate "$TFOLD_ENV" 2>/dev/null
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
        say "boltz-guided(HADDOCK) $T $pk (samples=$SAMP) ..."
        ( cd "$d" && boltz predict "$(basename "$y")" --out_dir "$od" --use_msa_server --no_kernels --diffusion_samples "$SAMP" > "boltz_$pk.log" 2>&1 ) \
          && say "  OK $T $pk" || { say "  !! $T $pk 실패 → $d/boltz_$pk.log"; tail -5 "$d/boltz_$pk.log"; }
        [ "$SMOKE" = "1" ] && { say "SMOKE 종료. 확인: find $d/$od -name '*_model_*.cif'"; exit 0; }
      done ;;
    protenix)
      for j in "$d"/protenix_${t}_patch*.json; do
        [ -e "$j" ] || continue; pk=$(basename "$j" .json | sed 's/.*_//'); od="results_$pk"
        if find "$d/$od" -name '*sample*.cif' 2>/dev/null | grep -q .; then say "skip $T $pk (있음)"; continue; fi
        say "protenix-guided(HADDOCK) $T $pk (samples=$SAMP) ..."
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
          say "tfold-guided(HADDOCK) $T $pk seed$s ..."
          ( cd "$TFOLD_DIR" && python projects/tfold_ag/predict.py --fasta "$fasta" --msa "$a3m" --icf "$icf" \
              --output "$outpdb" --model_version ppi --seed "$s" > "$od/${t}_seed${s}.log" 2>&1 ) \
            && say "  OK $T $pk seed$s" || { say "  !! $T $pk seed$s 실패 → $od/${t}_seed${s}.log"; tail -3 "$od/${t}_seed${s}.log"; }
          [ "$SMOKE" = "1" ] && { say "SMOKE 종료. 확인: $outpdb"; exit 0; }
        done
      done ;;
  esac
done
say "완료 → $OUTDIR/<t>/. 채점: python scripts/eval_guided_diverse.py --model $MODEL --guided-dir $OUTDIR --targets \"$TARGETS\""
