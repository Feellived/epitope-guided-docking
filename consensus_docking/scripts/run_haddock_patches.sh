#!/usr/bin/env bash
# 패치별 focused HADDOCK: partition_surface(boltz) → patch별 ambig/cfg → haddock3. self-healing(부분실패 정리 후 재실행).
# 사용: HAD_ENV=haddock3 bash scripts/run_haddock_patches.sh "8SIT 9ZDU" [K=6] [SAMPLES=100]
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
HAD_ENV="${HAD_ENV:-haddock3}"; K="${2:-6}"; SAMPLES="${3:-100}"
TARGETS="${1:-8SIT 9ZDU}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
conda activate boltz 2>/dev/null
python scripts/partition_surface.py --k "$K" --targets "$TARGETS"
conda activate "$HAD_ENV" 2>/dev/null
command -v haddock3 >/dev/null || { say "!! haddock3 없음(env $HAD_ENV)"; exit 1; }
for T in $TARGETS; do
  t=$(echo "$T"|tr 'A-Z' 'a-z'); d="haddock/$t"
  [ -d "$d" ] || { say "skip $T (no $d)"; continue; }
  python scripts/make_haddock_patch_restraints.py --dir "$d" --samples "$SAMPLES" || { say "!! $T 제약생성 실패"; continue; }
  for cfg in "$d"/run_patch_*.cfg; do
    [ -e "$cfg" ] || continue
    pi=$(basename "$cfg" | sed 's/run_patch_//;s/.cfg//')
    rundir="$d/patch_$pi/run"
    if [ -d "$rundir/6_seletopclusts" ]; then say "skip $T patch$pi (완료)"; continue; fi
    [ -d "$rundir" ] && { say "  $T patch$pi 부분run 정리 후 재실행"; rm -rf "$rundir"; }
    say "haddock3 $T patch$pi (sampling $SAMPLES) ..."
    haddock3 "$cfg" > "$d/patch_$pi.log" 2>&1 && say "  OK $T patch$pi" || say "  !! $T patch$pi 실패 → $d/patch_$pi.log"
  done
done
say "완료 → haddock/<t>/patch_<N>/run/6_seletopclusts. 분석: python scripts/haddock_patch_dockq.py"
