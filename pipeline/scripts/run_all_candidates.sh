#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [선택기 천장] 후보를 하나만 고르지 말고 **전부** 재도킹한다.
#
# 왜 필요한가:
#   지금 파이프라인은 후보 하나를 골라 그것만 제약으로 준다. 그래서 "다른 후보를 줬으면
#   어땠을까"를 모른다. 후보의 정답 덮음(자리 축)은 이미 알지만, 그 자리를 줬을 때
#   DockQ(자세 축)가 얼마 나오는지는 **돌려보지 않으면 알 수 없다**. 자리가 좋아져도
#   자세가 안 따라오는 경우를 반복해서 봤기 때문이다(9w43_AB 자리 0.58→0.78, DockQ 그대로).
#
# 한 번 돌려 두 가지를 얻는다:
#   ① 천장   — 선택기가 완벽했다면 얼마나 좋아졌을까. 선택기가 잃고 있는 양의 크기
#   ② 우회   — 후보를 다 도킹한 뒤 AbEpiScore 로 고르면 도킹 전에 고르는 것보다 나은가
#              (선택 병목을 도킹 뒤로 미루는 설계가 되는지)
#
# 대상은 **선택 실패 8종**뿐이다. 이득 7종은 이미 최선 후보를 골랐음을 확인했다.
#
# 후보마다 OUTDIR 을 따로 쓴다. 같은 폴더에 쓰면 arm 이름이 'ours' 로 같아 덮어쓴다.
#
# 사용 (pipeline/ 에서):
#   bash scripts/run_all_candidates.sh              # 무엇을 돌릴지만 출력
#   bash scripts/run_all_candidates.sh --apply
#   ONLY="8k3k_D" bash scripts/run_all_candidates.sh --apply
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1; ROOT="$(pwd)"
MSAD="${MSAD:-$HOME/projects/msa-composition-bias/pipeline}"
HONEST="${HONEST:-$MSAD/results/honest}"
SHADOW="${SHADOW:-$MSAD/.honest_view}"
OUTBASE="${OUTBASE:-runs_allcand}"
MODELS="${MODELS:-boltz}"
SEEDS="${SEEDS:-42}"
SAMP="${SAMP:-5}"
# 선택 실패 8종 — diagnose_harm.py 의 '선택실패' 판정 그대로
TARGETS="${TARGETS:-8sit_HL 9azt_HL 8t4a_PR 9azv_HL 8siq_HL 8t4d_OQ 9mqr_DE 8k3k_D}"
ONLY="${ONLY:-}"
APPLY=0; [ "${1:-}" = "--apply" ] && APPLY=1
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }

[ -d "$HONEST" ] || { say "!! $HONEST 없음 — 정직한 판을 먼저 돌릴 것"; exit 1; }

# 그림자 폴더 — 옛 sites_*.json 을 덮지 않고 정직한 것을 같은 이름으로 보이게 한다
mkdir -p "$SHADOW/results" || exit 1
ln -sfn "$MSAD/targets" "$SHADOW/targets"
n=0
for f in "$HONEST"/sites_*.json; do
  [ -e "$f" ] || continue
  ln -sf "$f" "$SHADOW/results/$(basename "$f")" && n=$((n+1))
done
say "그림자 폴더 $SHADOW  (정직한 sites ${n}개)"
[ "$APPLY" = 0 ] && say "[dry-run — 아무것도 실행하지 않음. 돌리려면 --apply]"

total=0
printf '%-12s %-8s %s\n' 타깃 후보 상태
printf '%s\n' "--------------------------------------------------"
for T in $TARGETS; do
  [ -n "$ONLY" ] && ! echo " $ONLY " | grep -q " $T " && continue
  SITES="$SHADOW/results/sites_${T}.json"
  [ -f "$SITES" ] || { printf '%-12s %-8s %s\n' "$T" "-" "sites 없음 — 건너뜀"; continue; }
  NC=$(python -c "import json,sys; print(len(json.load(open(sys.argv[1]))['candidates']))" "$SITES") || continue
  for C in $(seq 1 "$NC"); do
    D="$ROOT/$OUTBASE/c$C/${T,,}/ours"
    done_n=$(find "$D" -name '*model*.cif' 2>/dev/null | wc -l)
    if [ "$done_n" -ge "$SAMP" ]; then
      printf '%-12s %-8s %s\n' "$T" "후보$C" "완료($done_n) — 건너뜀"; continue
    fi
    printf '%-12s %-8s %s\n' "$T" "후보$C" "돌릴 것"
    total=$((total+1))
    [ "$APPLY" = 0 ] && continue
    MSAD="$SHADOW" MODELS="$MODELS" ARMS="ours" SEEDS="$SEEDS" SAMP="$SAMP" \
      OUTDIR="$OUTBASE/c$C" METHOD=manual PICK="$C" \
      bash scripts/run_demo_guided.sh "$T"
  done
done

echo
say "돌릴 실행 $total 개 (자세 $SAMP 개씩 · 실측 2분/회 기준 약 $((total*2))분)"
[ "$APPLY" = 0 ] && say "돌리려면: bash scripts/run_all_candidates.sh --apply"
say "채점: 후보 폴더마다 eval_demo.py --outdir $OUTBASE/c<번호>"
