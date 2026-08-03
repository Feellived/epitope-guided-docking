#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [안전 채점] 후보 전체 재도킹(runs_allcand/c*/)을 채점한다.
#
# 안전이라고 붙인 이유:
#   ① 30종 원본 표(results/demo_dockq.csv)를 **건드리지 않는다.** 후보 폴더마다
#      --csv 로 다른 파일에 쓴다. 그래도 혹시 몰라 시작할 때 원본을 한 부 복사해 둔다.
#   ② 기본이 dry-run 이다. 무엇을 채점할지 세어 보여주기만 한다.
#   ③ 채점 전에 **실행 수가 설계값과 맞는지 먼저 센다.** 자세가 모자란 폴더는 표시한다.
#      (이 파이프라인의 실패는 멈추지 않고 조용히 빈다.)
#   ④ 'ours' 팔만 채점한다. 후보 폴더에는 그 팔밖에 없고, blind(자세 40개)를 후보마다
#      다시 채점하면 시간만 몇 배로 든다.
#
# 결과: results/allcand/c<번호>.csv  (+ 같은 이름의 .log)
#       합치기·천장 계산은 bk21-msa-depth-bias/pipeline/allcand_ceiling.py 가 한다.
#
# 사용 (consensus_docking/ 에서, DockQ 되는 env = boltz):
#   bash scripts/score_allcand.sh                 # 무엇을 채점할지만 출력
#   bash scripts/score_allcand.sh --apply
#   ONLY="8k3k_D" bash scripts/score_allcand.sh --apply
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1; ROOT="$(pwd)"

MSAD="${MSAD:-$HOME/projects/bk21-msa-depth-bias/pipeline}"
HONEST="${HONEST:-$MSAD/results/honest}"
VIEW="${VIEW:-$MSAD/.score_view}"
OUTBASE="${OUTBASE:-runs_allcand}"
SAVE="${SAVE:-results/allcand}"
MODELS="${MODELS:-boltz}"
SAMP="${SAMP:-5}"                       # 한 실행의 설계 자세 수
MAIN="results/demo_dockq.csv"
ONLY="${ONLY:-}"
APPLY=0; [ "${1:-}" = "--apply" ] && APPLY=1
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
PY="${PY:-python}"; command -v "$PY" >/dev/null 2>&1 || PY=python3
lower(){ echo "$1" | tr '[:upper:]' '[:lower:]'; }

[ -d "$OUTBASE" ] || { say "!! $OUTBASE 없음 — 재도킹부터"; exit 1; }
[ -d "$MSAD" ]    || { say "!! $MSAD 없음"; exit 1; }

TMP="$(mktemp -d)" || exit 1
trap 'rm -rf "$TMP"' EXIT

# ── ① 원본 한 부 떠 두기 ──────────────────────────────────────────────────────
mkdir -p "$SAVE" results
if [ -f "$MAIN" ]; then
  BAK="results/demo_dockq.bak_$(date '+%m%d_%H%M').csv"
  cp "$MAIN" "$BAK" && say "원본 보관 → $BAK  ($(wc -l < "$BAK" | tr -d ' ') 줄)"
else
  say "! $MAIN 이 없다 (백업 생략)"
fi

# ── 채점용 보기 폴더 — 정직한 sites 를 쓰되 pose_features 등 코드도 보이게 ─────
# .honest_view 에는 코드가 없어서 import pose_features 가 깨진다. 그래서 따로 만든다.
USEMSAD="$MSAD"
if [ -d "$HONEST" ] && ls "$HONEST"/sites_*.json >/dev/null 2>&1; then
  mkdir -p "$VIEW/results"
  for f in "$MSAD"/*; do
    b="$(basename "$f")"
    [ "$b" = "results" ] && continue
    ln -sfn "$f" "$VIEW/$b"
  done
  n=0
  for f in "$HONEST"/sites_*.json; do
    ln -sf "$f" "$VIEW/results/$(basename "$f")" && n=$((n+1))
  done
  USEMSAD="$VIEW"
  say "채점 보기 폴더 $VIEW  (정직한 sites ${n}개 + 코드 심볼릭)"
else
  say "! $HONEST 가 없다 — $MSAD 를 그대로 쓴다"
fi
if ! "$PY" -c "import sys;sys.path.insert(0,'$USEMSAD');import pose_features" 2>"$TMP/imp"; then
  say "!! $USEMSAD 에서 pose_features 를 import 못 한다 ($PY). conda activate boltz 했는지 확인"
  sed 's/^/     /' "$TMP/imp" | tail -3
  exit 1
fi

# ── 소문자 폴더명 → 원래 타깃 이름 (연관배열 없이 파일로) ──────────────────────
: > "$TMP/canon"
for d in "$MSAD"/targets/*/; do
  [ -d "$d" ] || continue
  t="$(basename "$d")"
  echo "$(lower "$t") $t" >> "$TMP/canon"
done
[ -s "$TMP/canon" ] || { say "!! $MSAD/targets 가 비었다"; exit 1; }

# ── ③ 먼저 센다 ──────────────────────────────────────────────────────────────
say "실행 수 점검 (설계 자세 $SAMP 개/실행)"
printf '%-7s %-12s %6s  %s\n' 후보 타깃 자세수 상태
printf '%s\n' "-------------------------------------------------------"
: > "$TMP/plan"
short=0; empty=0; ok=0
for D in "$OUTBASE"/c*/; do
  [ -d "$D" ] || continue
  K="$(basename "$D")"
  for TD in "$D"*/; do
    [ -d "$TD" ] || continue
    low="$(basename "$TD")"
    T="$(awk -v k="$low" '$1==k{print $2; exit}' "$TMP/canon")"
    if [ -z "$T" ]; then
      printf '%-7s %-12s %6s  %s\n' "$K" "$low" "-" "targets/ 에 없는 이름 — 건너뜀"; continue
    fi
    if [ -n "$ONLY" ] && ! echo " $ONLY " | grep -q " $T "; then continue; fi
    n="$(find "$TD/ours" -name '*model*.cif' 2>/dev/null | wc -l | tr -d ' ')"
    if [ "$n" -eq 0 ]; then
      printf '%-7s %-12s %6s  %s\n' "$K" "$T" "$n" "비었다 — 채점 제외"; empty=$((empty+1)); continue
    elif [ "$n" -lt "$SAMP" ]; then
      printf '%-7s %-12s %6s  %s\n' "$K" "$T" "$n" "!! 설계($SAMP)보다 적다 — 채점은 하되 표시"
      short=$((short+1))
    else
      printf '%-7s %-12s %6s  %s\n' "$K" "$T" "$n" "정상"; ok=$((ok+1))
    fi
    echo "$K $T" >> "$TMP/plan"
  done
done
NK="$(awk '{print $1}' "$TMP/plan" | sort -u | wc -l | tr -d ' ')"
echo
say "정상 $ok · 모자람 $short · 비었음 $empty · 채점할 후보 폴더 ${NK}개"
[ "$short" -gt 0 ] && say "⚠ 자세가 모자란 실행이 있다. 그 줄의 DockQ최고는 표본이 적어 낮게 나올 수 있다."
[ ! -s "$TMP/plan" ] && { say "채점할 것이 없다"; exit 0; }

# c1, c2, ... c10 순서로 (sort -V 는 mac 에 없어 숫자만 떼어 정렬한다)
for K in $(awk '{print $1}' "$TMP/plan" | sort -u | sed 's/^c//' | sort -n | sed 's/^/c/'); do
  TS="$(awk -v k="$K" '$1==k{printf " %s", $2}' "$TMP/plan")"
  if [ "$APPLY" = 0 ]; then say "[dry-run] $K ←$TS"; continue; fi
  say "채점 $K :$TS"
  "$PY" -u scripts/dockq_demo.py \
    --targets "$TS" --outdir "$OUTBASE/$K" --models "$MODELS" \
    --arms ours --msad "$USEMSAD" --csv "$SAVE/$K.csv" 2>&1 | tee "$SAVE/$K.log"
  if [ -s "$SAVE/$K.csv" ]; then
    say "  → $SAVE/$K.csv  ($(( $(wc -l < "$SAVE/$K.csv") - 1 )) 줄)"
  else
    say "  !! $SAVE/$K.csv 가 비었다 — 로그 확인: $SAVE/$K.log"
  fi
done

if [ "$APPLY" = 0 ]; then
  say "[dry-run — 아무것도 채점하지 않음. 돌리려면 --apply]"
  exit 0
fi

echo
say "원본 확인: $MAIN 은 그대로여야 한다"
ls -l "$MAIN" 2>/dev/null
say "다음:"
echo "  python -u $MSAD/allcand_ceiling.py --allcand $ROOT/$SAVE --dockq $ROOT/$MAIN"
