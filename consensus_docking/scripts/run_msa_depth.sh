#!/usr/bin/env bash
# MSA depth 인과실험: 항원 MSA 깊이만 변화 → Boltz drift 변화. (항체=single-seq 상수)
# 조건별 run-dir runs_msad_<D>/<t>/ 에 chains.json·native symlink + subsampled a3m + Boltz pose.
# 사용: bash scripts/run_msa_depth.sh "8SDF 8XSI 9ML8 8SIT" "143 32 8 1"
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
SRC="${SRC:-runs_rbd}"; SEED="${SEED:-0}"
TARGETS="${1:-8SDF 8XSI 9ML8 8SIT}"; DEPTHS="${2:-143 32 8 1}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
conda activate boltz 2>/dev/null
for D in $DEPTHS; do
  RUN="runs_msad_${D}"
  for T in $TARGETS; do
    t=$(echo "$T" | tr 'A-Z' 'a-z'); src="$SRC/$t"; dst="$RUN/$t"
    [ -f "$src/chains.json" ] || { say "skip $T (no chains.json)"; continue; }
    [ -f "$src/msa_$t/A.a3m" ] || { say "skip $T (no A.a3m)"; continue; }
    mkdir -p "$dst"
    ln -sf "$(cd "$src"; pwd)/chains.json" "$dst/chains.json"
    for f in "$src"/*.pdb "$src"/*.cif; do
      [ -e "$f" ] || continue; b=$(basename "$f")
      case "$b" in *boltz*|*chai*|*protenix*|*colabfold*|*tfold*) continue;; esac
      ln -sf "$(cd "$src"; pwd)/$b" "$dst/$b"
    done
    python scripts/subsample_a3m.py "$src/msa_$t/A.a3m" "$dst/A_d${D}.a3m" --depth "$D" --seed "$SEED" >/dev/null
    python scripts/make_boltz_msa_depth.py --run-dir "$SRC" --target "$T" --antigen-a3m "$(cd "$dst"; pwd)/A_d${D}.a3m" --out "$dst/boltz_${t}.yaml" >/dev/null
    if find "$dst/results_boltz" -name '*_model_*.cif' 2>/dev/null | grep -q .; then say "skip boltz $T d$D (있음)"; continue; fi
    say "boltz $T depth=$D ..."
    ( cd "$dst" && boltz predict "boltz_${t}.yaml" --out_dir results_boltz --no_kernels --diffusion_samples 5 > boltz.log 2>&1 ) \
      && say "  OK $T d$D" || say "  !! $T d$D 실패 → $dst/boltz.log"
  done
done
say "완료. 분석: 각 depth마다  RUN_DIR=runs_msad_<D> python scripts/where_docked.py --targets ... --formats Fab"
