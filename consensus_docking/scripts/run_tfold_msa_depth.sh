#!/usr/bin/env bash
# tFold MSA-depth: 항원 a3m 깊이만 변화(항체는 원래 PLM=MSA-free) → runs_msad_<D>/<t>/out_tfold
# d143 baseline = 기존 runs_rbd/<t>/out_tfold (cofolder). 여기선 d8(스위트) · d1(붕괴테스트)만.
# 사용: tmux new -s tfmsad; SRC=runs_rbd bash scripts/run_tfold_msa_depth.sh "8SIT 9ML9 8P5M ..." "8 1"
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
ROOT="$(pwd)"
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null
SRC="${SRC:-runs_rbd}"; SEED0="${SEED0:-0}"
TFOLD_ENV="${TFOLD_ENV:-tfold}"; TFOLD_DIR="$ROOT/tfold"
SEEDS="${SEEDS:-42 43 44 45 46}"
TARGETS="${1:-8P5M 8SDF 8SIQ 8SIS 8SIT 8XSI 9ML8 9ML9 9SBB 9ZDU}"; DEPTHS="${2:-8 1}"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
# tfold fasta 없으면 생성 (stdlib)
python scripts/make_tfold.py --run-dir "$SRC" >/dev/null 2>&1 || true
conda activate "$TFOLD_ENV" 2>/dev/null
python -c "import tfold" 2>/dev/null || say "!! tfold import 실패(env $TFOLD_ENV) — 계속 시도"
for D in $DEPTHS; do
  RUN="runs_msad_${D}"
  for T in $TARGETS; do
    t=$(echo "$T"|tr 'A-Z' 'a-z'); src="$SRC/$t"; dst="$RUN/$t"
    fasta="$ROOT/$src/tfold_$t.fasta"
    [ -f "$fasta" ] || { say "skip $T (no tfold fasta)"; continue; }
    [ -f "$src/msa_$t/A.a3m" ] || { say "skip $T (no A.a3m)"; continue; }
    mkdir -p "$dst/out_tfold"
    [ -f "$dst/A_d${D}.a3m" ] || python scripts/subsample_a3m.py "$src/msa_$t/A.a3m" "$dst/A_d${D}.a3m" --depth "$D" --seed "$SEED0" >/dev/null
    [ -e "$dst/chains.json" ] || ln -sf "$ROOT/$src/chains.json" "$dst/chains.json"
    a3m="$ROOT/$dst/A_d${D}.a3m"
    for s in $SEEDS; do
      outpdb="$ROOT/$dst/out_tfold/${t}_seed${s}.pdb"
      [ -s "$outpdb" ] && { say "skip $T d$D seed$s (있음)"; continue; }
      say "tfold $T d$D seed$s ..."
      ( cd "$TFOLD_DIR" && python projects/tfold_ag/predict.py --fasta "$fasta" --msa "$a3m" --output "$outpdb" --seed "$s" > "$ROOT/$dst/out_tfold/${t}_seed${s}.log" 2>&1 )
    done
    n=$(find "$ROOT/$dst/out_tfold" -name "${t}_seed*.pdb" -size +0c 2>/dev/null | wc -l)
    say "  $T d$D: $n pose"
  done
done
say "완료 → runs_msad_<D>/<t>/out_tfold"
