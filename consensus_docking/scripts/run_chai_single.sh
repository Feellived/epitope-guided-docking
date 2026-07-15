#!/usr/bin/env bash
# Chai single-seq(MSA 없음, 양 사슬) = d1 극단 조건 → runs_msad_1/<t>/out_chai
# 대비 baseline = 기존 runs_rbd/<t>/out_chai (full server MSA, cofolder). 사용: CHAI_ENV=chai bash scripts/run_chai_single.sh "8SIT ..."
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
ROOT="$(pwd)"
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
SRC="${SRC:-runs_rbd}"; CHAI_ENV="${CHAI_ENV:-chai}"
export CHAI_DOWNLOADS_DIR="${CHAI_DOWNLOADS_DIR:-/mnt/data/admuser/chai_downloads}"
TARGETS="${1:-8P5M 8SDF 8SIQ 8SIS 8SIT 8XSI 9ML8 9ML9 9SBB 9ZDU}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
conda activate "$CHAI_ENV" 2>/dev/null
command -v chai-lab >/dev/null || { say "!! chai-lab 없음(env $CHAI_ENV)"; exit 1; }
for T in $TARGETS; do
  t=$(echo "$T"|tr 'A-Z' 'a-z'); src="$SRC/$t"; dst="runs_msad_1/$t"
  fasta="$ROOT/$src/chai_$t.fasta"
  [ -f "$fasta" ] || { say "skip $T (no chai fasta — make_chai 먼저)"; continue; }
  mkdir -p "$dst"; [ -e "$dst/chains.json" ] || ln -sf "$ROOT/$src/chains.json" "$dst/chains.json"
  if find "$dst/out_chai" -name '*.cif' 2>/dev/null | grep -q .; then say "skip $T (있음)"; continue; fi
  say "chai single-seq $T ..."
  ( cd "$dst" && chai-lab fold "$fasta" out_chai > chai_single.log 2>&1 ) && say "  OK $T" || say "  !! $T 실패 → $dst/chai_single.log"
done
say "완료 → runs_msad_1/<t>/out_chai (single-seq). 대비=runs_rbd/<t>/out_chai(full)"
