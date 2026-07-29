#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [사전 스크리닝] 후보 자리마다 값싸게 한 번 모델링 → AbEpiTope 로 계면 채점 → 이긴 후보 고르기.
#
# 왜 2패스인가 — AbEpiTope 은 항원 표면이 아니라 **복합체**를 채점한다. 그래서 "이 자리에 붙였을
# 때 계면이 말이 되나"를 물으려면 후보마다 구조가 하나씩 있어야 한다. 후보는 보통 1~7개고
# 여기서는 1샘플만 뽑으므로 몇 분이면 끝난다.
#
# 흐름:
#   ① sites_to_pocket.py --arm all   후보 전부 → 제약 입력 (patch<후보번호>)
#   ② 후보마다 co-folder 1샘플 예측
#   ③ score_abepitope.py             후보별 점수 CSV
#   ④ rank_sites.py --method abepitope   → PICK
#
# 이 스크립트는 ①~④ 만 한다. 고른 후보로 **본** 유도 예측을 돌리는 것은 run_demo_guided.sh 몫:
#   METHOD=manual PICK=<나온 번호> bash scripts/run_demo_guided.sh <타깃>
#
# 사용(consensus_docking/ 에서):
#   bash scripts/prescreen_abepitope.sh 8k3k_D
#   MODEL=protenix bash scripts/prescreen_abepitope.sh 8k3k_D
#   ABEPITOPE_CMD="conda run -n abepitope ..." bash scripts/prescreen_abepitope.sh 8k3k_D
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
MSAD="${MSAD:-$HOME/projects/bk21-msa-depth-bias/pipeline}"
T="${1:-}"
[ -z "$T" ] && { echo "사용: bash scripts/prescreen_abepitope.sh <타깃>"; exit 1; }
MODEL="${MODEL:-boltz}"                 # 사전 스크리닝은 한 모델이면 충분
BOLTZ_ENV="${BOLTZ_ENV:-boltz}"; PROT_ENV="${PROT_ENV:-protenix}"
MSAD_DATA="${MSAD_DATA:-/mnt/data/admuser/msadepth}"
PROT_MODEL="${PROT_MODEL:-protenix_base_default_v1.0.0}"   # blind 팔과 같은 모델
OUTDIR="${OUTDIR:-runs_sites_guided}"
SEED="${SEED:-42}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }

SITES="$MSAD/results/sites_${T}.json"
CHAINS="$MSAD/targets/${T}/chains.json"
[ -f "$SITES" ]  || { say "!! $SITES 없음 — msa-depth 에서: bash analyze_target.sh $T"; exit 1; }
[ -f "$CHAINS" ] || { say "!! $CHAINS 없음"; exit 1; }

# ── ① 후보 전부 방출 ──────────────────────────────────────────────────────────
say "[①] 후보 전부 → 제약 입력"
python scripts/sites_to_pocket.py --sites "$SITES" --chains "$CHAINS" \
       --arm all --emit "$MODEL" --outdir "$OUTDIR" || exit 1
ALLD="$OUTDIR/${T,,}/all"

# ── ② 후보마다 1샘플 예측 ─────────────────────────────────────────────────────
case "$MODEL" in
  boltz)    conda activate "$BOLTZ_ENV" 2>/dev/null
            command -v boltz >/dev/null || { say "!! boltz 없음"; exit 1; } ;;
  protenix) conda activate "$PROT_ENV" 2>/dev/null
            export PROTENIX_ROOT_DIR=/mnt/data/admuser/protenix_weights LAYERNORM_TYPE=torch
            command -v protenix >/dev/null || { say "!! protenix 없음"; exit 1; } ;;
esac

N=0
for IN in "$ALLD"/${MODEL}_${T,,}_patch*.{yaml,json}; do
  [ -f "$IN" ] || continue
  K=$(basename "$IN" | sed -n 's/.*patch\([0-9]*\)\..*/\1/p')
  O="$ALLD/out_${MODEL}_patch${K}"
  if [ -d "$O" ] && find "$O" -name '*.cif' 2>/dev/null | grep -q .; then
    say "  후보 $K 이미 있음 — 건너뜀"; N=$((N+1)); continue
  fi
  mkdir -p "$O"
  say "  후보 $K 예측 (1샘플)"
  # 호출 방식은 comp_x_reps.sh 와 동일하게(검증된 형태). 다르게 부르면 조용히 빈 채로 끝난다.
  IN_ABS="$(cd "$(dirname "$IN")" && pwd)/$(basename "$IN")"
  case "$MODEL" in
    boltz)    export BOLTZ_CACHE="${BOLTZ_CACHE:-$MSAD_DATA/boltz_cache}"; mkdir -p "$BOLTZ_CACHE"
              ( cd "$O" && boltz predict "$IN_ABS" --out_dir results --cache "$BOLTZ_CACHE" \
                     --no_kernels --use_msa_server --diffusion_samples 1 >run.log 2>&1 ) \
                || say "    ! 후보 $K 실패 ($O/run.log)" ;;
    protenix) ( cd "$O" && protenix pred -i "$IN_ABS" -o results -n "$PROT_MODEL" -s "$SEED" -e 1 \
                     --trimul_kernel torch --triatt_kernel torch --enable_fusion False >run.log 2>&1 ) \
                || say "    ! 후보 $K 실패 ($O/run.log)" ;;
  esac
  N=$((N+1))
done
[ "$N" -gt 0 ] || { say "!! 방출된 후보 입력이 없다"; exit 1; }
# 채점기가 patch<번호>로 후보를 알아보도록, 출력 폴더 이름에 이미 patch<번호>가 들어 있다.

# ── ③ AbEpiTope 채점 ──────────────────────────────────────────────────────────
say "[③] AbEpiTope 채점"
python scripts/score_abepitope.py --target "$T" --dir "$ALLD" \
       --out "results/abepitope_${T}.csv" || {
  say "!! 채점 실패 — AbEpiTope 사용법을 확인하고 ABEPITOPE_CMD 로 지정할 것."
  say "   그동안은 METHOD=ncomp 로 데모를 돌릴 수 있다."; exit 2; }

# ── ④ 선택 ────────────────────────────────────────────────────────────────────
say "[④] 후보 선택"
python scripts/rank_sites.py --sites "$SITES" --method abepitope \
       --scores "results/abepitope_${T}.csv" --out "$OUTDIR/${T,,}/pick.json"

echo
say "본 유도 예측:  METHOD=manual PICK=<위 PICK 값> bash scripts/run_demo_guided.sh $T"
