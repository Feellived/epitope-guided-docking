#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [데모] MSA 조성 → 결합자리 후보 → (선택) → 유도 재도킹 → DockQ 상승 확인.
#
# 흐름 (한 타깃):
#   ① 생성  bk21-msa-depth-bias 본 검정 결과 = results/sites_<타깃>.json  (이미 있어야 함)
#   ② 선택  rank_sites.py  → 후보 하나 (정답 구조를 보지 않음. --method manual 은 예외)
#   ③ 유도  sites_to_pocket.py → co-folder pocket 제약 입력 (두 팔)
#   ④ 실행  boltz / protenix 를 blind 와 **같은 커맨드**로, 제약만 붙여서
#   ⑤ 채점  dockq_demo.py 로 세 팔 비교
#
# 세 팔:
#   (a) blind   제약 없음 — msa-depth 의 원래 MSA 예측을 그대로 쓴다(여기서 안 돌림)
#   (b) fullmsa 원래 MSA 가 간 자리로 유도  ← 핵심 대조
#   (c) ours    우리가 고른 후보로 유도
#   판정은 (c) > (b). (c) > (a) 만으로는 '제약을 주면 좋아진다'와 구별이 안 된다.
#
# 사용(consensus_docking/ 에서):
#   SMOKE=1 bash scripts/run_demo_guided.sh 8k3k_D                 # 1샘플로 배관 확인
#   bash scripts/run_demo_guided.sh 8k3k_D
#   METHOD=manual PICK=1 bash scripts/run_demo_guided.sh 8k3k_D    # 사람이 고른 경우
#   MODELS=boltz bash scripts/run_demo_guided.sh "8k3k_D 8sit_HL"
# env: MSAD(msa-depth pipeline 경로) MODELS METHOD PICK SCORES SAMP SEEDS BOLTZ_ENV PROT_ENV OUTDIR
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1; ROOT="$(pwd)"
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
MSAD="${MSAD:-$HOME/projects/bk21-msa-depth-bias/pipeline}"
TARGETS="${1:-}"
[ -z "$TARGETS" ] && { echo "사용: bash scripts/run_demo_guided.sh \"<타깃> [타깃...]\"  (예: 8k3k_D)"; exit 1; }
MODELS="${MODELS:-boltz protenix}"
METHOD="${METHOD:-ncomp}"; PICK="${PICK:-0}"; SCORES="${SCORES:-}"
SMOKE="${SMOKE:-0}"; SAMP="${SAMP:-5}"; SEEDS="${SEEDS:-42 43 44 45 46}"
BOLTZ_ENV="${BOLTZ_ENV:-boltz}"; PROT_ENV="${PROT_ENV:-protenix}"
# ⚠️ blind 팔(msa-depth 본 검정)이 쓴 것과 **같은 모델**이어야 세 팔 비교가 성립한다.
#    comp_x_reps.sh 의 기본값과 맞춘다(RBD/본검정 세트 = default).
PROT_MODEL="${PROT_MODEL:-protenix_base_default_v1.0.0}"
MSAD_DATA="${MSAD_DATA:-/mnt/data/admuser/msadepth}"
OUTDIR="${OUTDIR:-runs_sites_guided}"
SEED="${SEED:-42}"
[ "$SMOKE" = "1" ] && { SAMP=1; SEEDS="42"; }
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }

for T in $TARGETS; do
  echo; say "████ $T ████"
  SITES="$MSAD/results/sites_${T}.json"
  CHAINS="$MSAD/targets/${T}/chains.json"
  [ -f "$SITES" ]  || { say "!! $SITES 없음 — 먼저 msa-depth 에서: bash analyze_target.sh $T"; continue; }
  [ -f "$CHAINS" ] || { say "!! $CHAINS 없음"; continue; }

  # ── ②-0 DiscoTope 점수가 필요한데 없으면 먼저 만든다 (GPU 불필요) ──────────
  if [ "$METHOD" = "discotope" ] && [ -z "$SCORES" ]; then
    SCORES="results/discotope_${T}.csv"
    if [ -s "$SCORES" ]; then
      say "[점수] $SCORES 이미 있음 — 재사용"
    else
      say "[점수] DiscoTope-3.0 실행 (항원 단독 구조, MSA 안 씀)"
      python scripts/run_discotope.py --target "$T" ${DISCO_REPO:+--repo "$DISCO_REPO"} \
             ${DISCO_PY:+--python "$DISCO_PY"} --out "$SCORES" \
        || { say "!! DiscoTope 실패 — METHOD=ncomp 로 내려간다"; METHOD=ncomp; SCORES=""; }
    fi
  fi

  # ── ② 선택 ─────────────────────────────────────────────────────────────────
  say "[선택] rank_sites.py --method $METHOD"
  RS_ARGS=(--sites "$SITES" --method "$METHOD" --out "$OUTDIR/${T,,}/pick.json")
  [ "$METHOD" = "manual" ]   && RS_ARGS+=(--pick "$PICK")
  [ -n "$SCORES" ]           && RS_ARGS+=(--scores "$SCORES")
  mkdir -p "$OUTDIR/${T,,}"
  OUT=$(python scripts/rank_sites.py "${RS_ARGS[@]}") || { say "!! 선택 실패"; continue; }
  echo "$OUT"
  CAND=$(echo "$OUT" | sed -n 's/^PICK=\([0-9]*\)$/\1/p' | tail -1)
  [ -n "$CAND" ] || { say "!! 후보 번호를 못 받았다"; continue; }

  # ── ③ 유도 입력 (두 팔) ────────────────────────────────────────────────────
  for ARM in fullmsa ours; do
    A_ARGS=(--sites "$SITES" --chains "$CHAINS" --arm "$ARM" --emit "${MODELS// /,}" --outdir "$OUTDIR")
    [ "$ARM" = "ours" ] && A_ARGS+=(--cand "$CAND")
    say "[유도입력] $ARM"
    python scripts/sites_to_pocket.py "${A_ARGS[@]}" || say "  ! $ARM 입력 생성 실패(건너뜀)"
  done

  # ── ④ 실행 ─────────────────────────────────────────────────────────────────
  for M in $MODELS; do
    case "$M" in
      boltz)    conda activate "$BOLTZ_ENV" 2>/dev/null
                command -v boltz >/dev/null || { say "!! boltz 없음 — 건너뜀"; continue; } ;;
      protenix) conda activate "$PROT_ENV" 2>/dev/null
                export PROTENIX_ROOT_DIR=/mnt/data/admuser/protenix_weights LAYERNORM_TYPE=torch
                command -v protenix >/dev/null || { say "!! protenix 없음 — 건너뜀"; continue; } ;;
    esac
    for ARM in fullmsa ours; do
      D="$OUTDIR/${T,,}/$ARM"
      case "$M" in
        boltz)    IN="$D/boltz_${T,,}_patch0.yaml" ;;
        protenix) IN="$D/protenix_${T,,}_patch0.json" ;;
      esac
      [ -f "$IN" ] || { say "  $M/$ARM 입력 없음 — 건너뜀"; continue; }
      O="$D/out_$M"
      if [ -d "$O" ] && find "$O" -name '*.cif' 2>/dev/null | grep -q .; then
        say "  $M/$ARM 이미 있음 — 건너뜀"; continue
      fi
      mkdir -p "$O"
      # ⚠️ 호출 방식은 comp_x_reps.sh(이 서버에서 검증된 것)와 똑같이 맞춘다.
      #    다르게 부르면 조용히 아무것도 안 만들고 끝난다(2026-07-29 실제로 그랬다).
      #    · protenix 는 'pred -i -o -n -s -e' 이고 커널 플래그가 필요하다
      #    · boltz 는 --cache 와 --no_kernels 가 필요하다
      #    · 둘 다 출력 폴더로 cd 한 뒤 -o/--out_dir 를 'results' 로 준다
      IN_ABS="$(cd "$(dirname "$IN")" && pwd)/$(basename "$IN")"
      say "  실행 $M/$ARM (샘플 $SAMP · seed $SEED)"
      case "$M" in
        boltz)
          export BOLTZ_CACHE="${BOLTZ_CACHE:-$MSAD_DATA/boltz_cache}"; mkdir -p "$BOLTZ_CACHE"
          ( cd "$O" && boltz predict "$IN_ABS" --out_dir results --cache "$BOLTZ_CACHE" \
                 --no_kernels --use_msa_server --diffusion_samples "$SAMP" >run.log 2>&1 ) \
            || say "    ! boltz 실패 — $O/run.log 확인" ;;
        protenix)
          ( cd "$O" && protenix pred -i "$IN_ABS" -o results -n "$PROT_MODEL" -s "$SEED" -e "$SAMP" \
                 --trimul_kernel torch --triatt_kernel torch --enable_fusion False >run.log 2>&1 ) \
            || say "    ! protenix 실패 — $O/run.log 확인" ;;
      esac
      find "$O/results" -name '*.cif' 2>/dev/null | grep -q . \
        || say "    ! 산출물이 없다 — $O/run.log 의 마지막 줄을 볼 것"
    done
  done
  say "$T 끝 → $OUTDIR/${T,,}/"
done

echo
say "채점:  conda activate DockQ && python scripts/dockq_demo.py --targets \"$TARGETS\" --outdir $OUTDIR"
