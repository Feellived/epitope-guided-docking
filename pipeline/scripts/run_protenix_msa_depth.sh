#!/usr/bin/env bash
# Protenix MSA-depth: 항원 MSA 깊이만 변화(항체 single-seq) → runs_msad_<D>/<t>/out_protenix (GPU; HADDOCK CPU와 병렬)
# 사용: tmux new -s protmsad; PROT_ENV=protenix bash scripts/run_protenix_msa_depth.sh "8SIT 8SDF 8P5M ..." "143 8 1"
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
SRC="${SRC:-runs_rbd}"; SEED="${SEED:-0}"; PROT_ENV="${PROT_ENV:-protenix}"
PROT_MODEL="${PROT_MODEL:-protenix_base_default_v1.0.0}"
TARGETS="${1:-8P5M 8SDF 8SIQ 8SIS 8SIT 8XSI 9ML8 9ML9 9SBB 9ZDU}"; DEPTHS="${2:-143 8 1}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
conda activate "$PROT_ENV" 2>/dev/null
export PROTENIX_ROOT_DIR=/mnt/data/admuser/protenix_weights LAYERNORM_TYPE=torch
command -v protenix >/dev/null || { say "!! protenix 없음(env $PROT_ENV)"; exit 1; }
for D in $DEPTHS; do
  RUN="runs_msad_${D}"
  for T in $TARGETS; do
    t=$(echo "$T"|tr 'A-Z' 'a-z'); src="$SRC/$t"; dst="$RUN/$t"
    [ -f "$src/msa_$t/A.a3m" ] || { say "skip $T (no A.a3m)"; continue; }
    mkdir -p "$dst"
    [ -f "$dst/A_d${D}.a3m" ] || python scripts/prep_a3m_subsample.py "$src/msa_$t/A.a3m" "$dst/A_d${D}.a3m" --depth "$D" --seed "$SEED" >/dev/null
    [ -e "$dst/chains.json" ] || ln -sf "$(cd "$src";pwd)/chains.json" "$dst/chains.json"
    for f in "$src"/*.pdb; do [ -e "$f" ] || continue; b=$(basename "$f"); case "$b" in *boltz*|*chai*|*protenix*|*colabfold*|*tfold*) continue;; esac; ln -sf "$(cd "$src";pwd)/$b" "$dst/$b" 2>/dev/null; done
    python scripts/make_protenix_msa_depth.py --run-dir "$SRC" --target "$T" --antigen-a3m "$(cd "$dst";pwd)/A_d${D}.a3m" --dir "$dst" >/dev/null
    if find "$dst/out_protenix" -name '*sample*.cif' 2>/dev/null | grep -q .; then say "skip protenix $T d$D (있음)"; continue; fi
    say "protenix $T depth=$D ..."
    ( cd "$dst" && protenix pred -i "protenix_${t}.json" -o out_protenix -n "$PROT_MODEL" -s 101 -e 5 --trimul_kernel torch --triatt_kernel torch > protenix_msad.log 2>&1 ) \
      && say "  OK $T d$D" || say "  !! $T d$D 실패 → $dst/protenix_msad.log"
  done
done
say "완료 → runs_msad_<D>/<t>/out_protenix. 분석: lib_structure(Protenix msa_depth 소스 포함됨)"
