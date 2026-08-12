#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [정직한 유도 재도킹] 정답을 보지 않고 고른 결합자리로 boltz 도킹을 다시 돌린다.
#
# 왜 껍데기가 필요한가 (2026-08-01):
#   run_demo_guided.sh 는 $MSAD/results/sites_<타깃>.json 을 읽는다. 그런데 그 파일들은
#   **실행마다 DockQ(정답) 최고 자세를 대표로 골라** 만든 옛 판이다. 정직한 판은
#   results/honest/ 에 따로 있고, 후보 선택도 ncomp 가 아니라 AbEpiScore 최대로 한다.
#   → 옛 파일을 덮지 않고 **심볼릭 링크로 그림자 폴더**를 만들어 그쪽을 읽게 한다.
#
# 정답을 쓰는 곳: 없다. 자세 선택 = AbEpiScore · 후보 선택 = AbEpiScore 최대 · 채점만 정답.
#   ⚠️ METHOD=manual 을 쓰지만 사람이 고른 것이 아니라 **eval_selectors 가 고른 번호**를 넘긴다.
#
# 예산(타깃 하나): 네 팔 × 자세 5개 = 20회 ≈ 45분(protenix 실측 2.2분/회 기준).
#   30종이면 하룻밤. ARMS=ours 로 줄이면 4분의 1.
#
# 사용(pipeline/ 에서):
#   bash scripts/run_honest_guided.sh                      # 무엇을 돌릴지만 출력
#   bash scripts/run_honest_guided.sh --apply
#   ONLY="8k3k_D 8tp5_HL" bash scripts/run_honest_guided.sh --apply
#   ARMS=ours bash scripts/run_honest_guided.sh --apply     # 값싼 훑기
#   SMOKE=1 ONLY=8k3k_D bash scripts/run_honest_guided.sh --apply   # 배관 확인
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1; ROOT="$(pwd)"
MSAD="${MSAD:-$HOME/projects/msa-composition-bias/pipeline}"
HONEST="${HONEST:-$MSAD/results/honest}"
PICKCSV="${PICKCSV:-$HONEST/pick_abepi_max.csv}"
SHADOW="${SHADOW:-$MSAD/.honest_view}"
OUTDIR="${OUTDIR:-runs_sites_guided_honest}"
MODELS="${MODELS:-boltz}"
ARMS="${ARMS:-noconstraint fullmsa sizematch ours}"
SEEDS="${SEEDS:-42}"           # 팔당 자세 5개면 충분. 늘리면 비용이 선형으로 는다
SAMP="${SAMP:-5}"
ONLY="${ONLY:-}"
SMOKE="${SMOKE:-0}"
APPLY=0; [ "${1:-}" = "--apply" ] && APPLY=1
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }

[ -d "$HONEST" ] || { say "!! $HONEST 없음 — msa-depth 에서 정직한 판을 먼저 돌릴 것"; exit 1; }
[ -f "$PICKCSV" ] || { say "!! $PICKCSV 없음 — eval_selectors.py --pick abepi_max --out 로 먼저 만들 것"; exit 1; }

# ── 그림자 폴더: 옛 sites_*.json 을 덮지 않고 정직한 것을 같은 이름으로 보이게 한다 ──
mkdir -p "$SHADOW/results" || exit 1
ln -sfn "$MSAD/targets" "$SHADOW/targets"
n_link=0
for f in "$HONEST"/sites_*.json; do
  [ -e "$f" ] || continue
  ln -sf "$f" "$SHADOW/results/$(basename "$f")" && n_link=$((n_link+1))
done
say "그림자 폴더 $SHADOW  (정직한 sites $n_link개 연결 · 옛 파일은 건드리지 않음)"

# ── 타깃마다 '무엇을 넘길지'는 eval_selectors 가 고른 후보 번호 ──
say "팔=[$ARMS] 모델=[$MODELS] seed=[$SEEDS] 자세=$SAMP"
[ "$APPLY" = 0 ] && say "[dry-run — 아무것도 실행하지 않음. 돌리려면 --apply]"
printf '%-12s %-6s %-6s %-8s %s\n' 타깃 후보 조성수 잔기 상태
printf '%s\n' "------------------------------------------------------------"

tail -n +2 "$PICKCSV" | while IFS=, read -r T SEL CAND NCOMP NRES REC PRE F1 FULL; do
  [ -n "${T:-}" ] || continue
  [ -n "$ONLY" ] && ! echo " $ONLY " | grep -q " $T " && continue
  if [ ! -f "$SHADOW/results/sites_${T}.json" ]; then
    printf '%-12s %-6s %-6s %-8s %s\n' "$T" "$CAND" "$NCOMP" "$NRES" "sites 없음 — 건너뜀"; continue
  fi
  done_n=$(find "$ROOT/$OUTDIR/$(echo "$T" | tr 'A-Z' 'a-z')" -name '*model*.cif' 2>/dev/null | wc -l)
  if [ "$done_n" -ge 20 ]; then
    printf '%-12s %-6s %-6s %-8s %s\n' "$T" "$CAND" "$NCOMP" "$NRES" "완료($done_n) — 건너뜀"; continue
  fi
  printf '%-12s %-6s %-6s %-8s %s\n' "$T" "$CAND" "$NCOMP" "$NRES" "돌릴 것(현재 $done_n)"
  [ "$APPLY" = 0 ] && continue
  MSAD="$SHADOW" MODELS="$MODELS" ARMS="$ARMS" SEEDS="$SEEDS" SAMP="$SAMP" \
    SMOKE="$SMOKE" OUTDIR="$OUTDIR" METHOD=manual PICK="$CAND" \
    bash scripts/run_demo_guided.sh "$T"
done

echo
say "채점: python scripts/eval_demo.py --outdir $OUTDIR   (네 팔 비교)"
say "⚠️ 판정은 ours > sizematch > fullmsa. ours > noconstraint 를 못 넘으면 기여 없음."
