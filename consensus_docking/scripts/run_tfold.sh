#!/usr/bin/env bash
# tFold-Ag 무인 실행 (Fab+Fv 세트): (1) colabfold로 항원 a3m 생성 → (2) tFold predict × seed 42-46.
# tFold는 항체 MSA-free(PLM), 항원 a3m만 필요. JSON 모드 버그 → FASTA(H/L/A) 모드.
#
# 사용:
#   tmux new -s tfold
#   ONLY=9k6j bash scripts/run_tfold.sh     # 스모크(9k6j Fab+Fv) — 내일 아침 먼저 이거
#   bash scripts/run_tfold.sh               # 전체(runs_diverse/* Fab+Fv)
#   TFOLD_ENV=<이름> CF_ENV=<이름> bash ...  # env 이름 다르면 override
#
# 출력(cif/pdb) 있으면 스킵 = 이어받기. eval은 아침에 eval_tfold.py로.
# ⚠️ 초안: (a) tFold conda env 이름, (b) colabfold --msa-only a3m 경로를 스모크로 확정할 것.
cd "$(dirname "$0")/.." || exit 1
ROOT="$(pwd)"                          # .../consensus_docking
source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null

RUN="${RUN:-runs_diverse}"
CF_ENV="${CF_ENV:-colabfold}"
TFOLD_ENV="${TFOLD_ENV:-tfold}"
TFOLD_DIR="$ROOT/tfold"                 # tFold 레포 루트 (predict.py가 tfold 패키지 import)
SEEDS="${SEEDS:-42 43 44 45 46}"
LOG="run_tfold.log"
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

if [ -n "${ONLY:-}" ]; then
  dirs="$RUN/${ONLY} $RUN/${ONLY}_fv"
else
  dirs=$(ls -d $RUN/*/ 2>/dev/null | sed 's:/$::')
fi

# --- 입력 생성 (tfold_<id>.fasta + antigen_<id>.fasta) ---
say "입력 생성 (make_tfold)"
conda activate boltz 2>/dev/null
python scripts/make_tfold.py --run-dir "$RUN" >> "$LOG" 2>&1

# ========== Phase 1: 항원 a3m (colabfold --msa-only) ==========
say "===== Phase 1: 항원 a3m (colabfold) ====="
conda activate "$CF_ENV" 2>/dev/null
command -v colabfold_batch >/dev/null || { say "!! colabfold_batch 없음 (env '$CF_ENV') → 중단"; exit 1; }
for d in $dirs; do
  [ -d "$d" ] || continue
  base=$(basename "$d"); t=$(echo "${base%_fv}" | tr 'A-Z' 'a-z')
  if find "$d/msa_$t" -name '*.a3m' 2>/dev/null | grep -q .; then
    say "skip a3m $base (있음)"; continue
  fi
  say "a3m $base ..."
  ( cd "$d" && env -u LD_LIBRARY_PATH colabfold_batch "antigen_$t.fasta" "msa_$t" --msa-only > agmsa.log 2>&1 )
done

# ========== Phase 2: tFold predict × seed ==========
say "===== Phase 2: tFold predict (seed $SEEDS) ====="
conda activate "$TFOLD_ENV" 2>/dev/null
python -c "import tfold" 2>/dev/null || say "!! 경고: '$TFOLD_ENV' env에서 tfold import 실패 — env 이름 확인 필요(계속 시도)"
for d in $dirs; do
  [ -d "$d" ] || continue
  base=$(basename "$d"); t=$(echo "${base%_fv}" | tr 'A-Z' 'a-z')
  a3m=$(find "$ROOT/$d/msa_$t" -name '*.a3m' 2>/dev/null | head -1)
  if [ -z "$a3m" ]; then say "!! $base a3m 없음 → 스킵"; continue; fi
  fasta="$ROOT/$d/tfold_$t.fasta"
  mkdir -p "$ROOT/$d/out_tfold"
  for s in $SEEDS; do
    outpdb="$ROOT/$d/out_tfold/${t}_seed${s}.pdb"
    [ -s "$outpdb" ] && { say "skip $base seed$s (있음)"; continue; }
    say "tfold $base seed$s ..."
    ( cd "$TFOLD_DIR" && python projects/tfold_ag/predict.py \
        --fasta "$fasta" --msa "$a3m" --output "$outpdb" --seed "$s" \
        > "$ROOT/$d/out_tfold/${t}_seed${s}.log" 2>&1 )
  done
  n=$(find "$ROOT/$d/out_tfold" -name "${t}_seed*.pdb" -size +0c 2>/dev/null | wc -l)
  say "  $base: ${n} pose"
done
say "======== tFold DONE ========  (eval은 아침: eval_tfold.py)"
