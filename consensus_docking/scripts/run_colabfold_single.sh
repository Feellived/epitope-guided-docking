#!/usr/bin/env bash
# AF2-Multimer(ColabFold) single-seq(--msa-mode single_sequence, 양 사슬 MSA 없음) = d1 극단 → runs_msad_1/<t>/out_colabfold
# 대비 baseline = 기존 runs_rbd/<t>/out_colabfold(full MSA, cofolder). 사용: CF_ENV=colabfold bash scripts/run_colabfold_single.sh "8SIT ..."
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
ROOT="$(pwd)"
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
SRC="${SRC:-runs_rbd}"; CF_ENV="${CF_ENV:-colabfold}"
TARGETS="${1:-8P5M 8SDF 8SIQ 8SIS 8SIT 8XSI 9ML8 9ML9 9SBB 9ZDU}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
conda activate "$CF_ENV" 2>/dev/null
command -v colabfold_batch >/dev/null || { say "!! colabfold_batch 없음(env $CF_ENV)"; exit 1; }
for T in $TARGETS; do
  t=$(echo "$T"|tr 'A-Z' 'a-z'); src="$SRC/$t"; dst="runs_msad_1/$t"
  fasta="$ROOT/$src/colabfold_$t.fasta"
  [ -f "$fasta" ] || { say "skip $T (no colabfold fasta — make_colabfold 먼저)"; continue; }
  mkdir -p "$dst"; [ -e "$dst/chains.json" ] || ln -sf "$ROOT/$src/chains.json" "$dst/chains.json"
  if find "$dst/out_colabfold" \( -name '*.pdb' -o -name '*.cif' \) 2>/dev/null | grep -q .; then say "skip $T (있음)"; continue; fi
  say "colabfold single-seq $T ..."
  ( cd "$dst" && env -u LD_LIBRARY_PATH colabfold_batch "$fasta" out_colabfold \
      --msa-mode single_sequence --num-models 5 --model-type alphafold2_multimer_v3 --num-recycle 3 > colabfold_single.log 2>&1 ) \
    && say "  OK $T" || say "  !! $T 실패 → $dst/colabfold_single.log"
done
say "완료 → runs_msad_1/<t>/out_colabfold (single-seq). 대비=runs_rbd/<t>/out_colabfold(full)"
