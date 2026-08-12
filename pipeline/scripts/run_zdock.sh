#!/usr/bin/env bash
# ZDOCK 물리 도킹 (diverse Fab, native-split bound) — 2단계 오케스트레이터.
#   1) prep: anarci env (chain 분리 + ANARCI CDR 블록목록)
#   2) run : DockQ env (mark_sur → block.pl → zdock → create.pl → DockQ)
# 사용:  cd pipeline && bash scripts/run_zdock.sh [PDBID]      (인자 없으면 전체 10타깃)
# 환경(필요시 오버라이드):
#   ZDOCK_HOME=~/tools/zdock3.0.2_linux_x64  LIBG2C=~/tools/lib
#   ANARCI_ENV=anarci  DOCKQ_ENV=boltz  ZDOCK_FLAGS="-F -N 100"
set -e
CD=$(cd "$(dirname "$0")/.." && pwd)          # pipeline
export ZDOCK_HOME=${ZDOCK_HOME:-$HOME/tools/zdock3.0.2_linux_x64}
export LIBG2C=${LIBG2C:-$HOME/tools/lib}
export ZDOCK_FLAGS=${ZDOCK_FLAGS:-"-F -N 100"}
ANARCI_ENV=${ANARCI_ENV:-anarci}
DOCKQ_ENV=${DOCKQ_ENV:-boltz}
ONLY="$1"; ARG=${ONLY:+--only $ONLY}
cd "$CD"
echo "[1/2] prep ($ANARCI_ENV env) — chain 분리 + CDR"
conda run -n "$ANARCI_ENV" python scripts/prep_zdock.py $ARG
echo "[2/2] run ($DOCKQ_ENV env) — mark_sur→block→zdock→create→DockQ  (flags=$ZDOCK_FLAGS)"
conda run -n "$DOCKQ_ENV" python scripts/run_zdock_poses.py $ARG
echo "완료 → results/zdock_dockq.csv"
