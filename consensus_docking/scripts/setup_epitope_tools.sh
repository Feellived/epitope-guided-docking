#!/usr/bin/env bash
# ══════════════════════════════════════════════════════════════════════════════
# [설치] 결합자리 후보를 '정답을 보지 않고' 고르기 위한 구조 기반 채점 도구.
#
#   DiscoTope-3.0  항원 구조만 보고 잔기별 '에피토프다움' → 후보 자리 랭킹 (선택 단계 ②)
#   AbEpiTope-1.0  예측된 항체-항원 복합체의 계면 타당성 → 도킹 결과 재랭킹 (선택 단계 ④ 이후)
#
# 둘 다 역접힘(ESM-IF1) 기반이라 **MSA 를 쓰지 않는다** — MSA 조성이 만든 후보를 고르는 데
# 써도 순환 논리가 되지 않는다는 것이 이 조합을 고른 이유다.
#
# ⚠️ 기본은 dry-run(무엇을 할지 보여주기만). 실제 설치는 --apply.
#    별도 conda env 를 만들고 그 안에만 설치한다. 기존 env(boltz/protenix/tfold/DockQ)는 안 건드린다.
#    GPU 없이 CPU 로도 돈다(느릴 뿐).
#
# 사용:
#   bash scripts/setup_epitope_tools.sh                 # 무엇을 할지만 출력
#   bash scripts/setup_epitope_tools.sh --apply
#   TOOLS=~/tools ENVNAME=discotope bash scripts/setup_epitope_tools.sh --apply
#   bash scripts/setup_epitope_tools.sh --apply --only discotope
# ══════════════════════════════════════════════════════════════════════════════
set -uo pipefail
TOOLS="${TOOLS:-$HOME/tools}"
ENVNAME="${ENVNAME:-discotope}"
PYVER="${PYVER:-3.9}"
APPLY=0; ONLY=""
while [ $# -gt 0 ]; do
  case "$1" in
    --apply) APPLY=1 ;;
    --only)  shift; ONLY="${1:-}" ;;
    *) echo "!! 모르는 인자: $1"; exit 1 ;;
  esac; shift
done
say(){ echo "[$(date '+%H:%M:%S')] $*"; }
run(){
  if [ $APPLY -eq 1 ]; then say "\$ $*"; eval "$@" || { say "!! 실패: $*"; return 1; }
  else say "(dry-run) $*"; fi
}
want(){ [ -z "$ONLY" ] || [ "$ONLY" = "$1" ]; }

echo "설치 대상 = $TOOLS · conda env = $ENVNAME (python $PYVER)"
[ $APPLY -eq 1 ] || echo "※ 지금은 dry-run 이다. 실제로 설치하려면 --apply 를 붙일 것."
echo

source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
run "mkdir -p '$TOOLS'"

# ── 공통 env ──────────────────────────────────────────────────────────────────
if conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$ENVNAME"; then
  say "conda env '$ENVNAME' 이미 있음 — 새로 만들지 않는다"
else
  run "conda create -y -n '$ENVNAME' python=$PYVER"
fi
run "conda run -n '$ENVNAME' pip install --upgrade pip"

# ── DiscoTope-3.0 ─────────────────────────────────────────────────────────────
if want discotope; then
  echo; say "── DiscoTope-3.0 ──"
  if [ -d "$TOOLS/discotope3_web/.git" ]; then
    say "이미 받아져 있음: $TOOLS/discotope3_web (git pull 로 갱신만)"
    run "git -C '$TOOLS/discotope3_web' pull --ff-only"
  else
    run "git clone https://github.com/Magnushhoie/discotope3_web '$TOOLS/discotope3_web'"
  fi
  # torch 는 서버 CUDA 에 맞춰 이미 있는 것을 쓰는 게 안전하다. 없으면 CPU 판을 깐다.
  run "conda run -n '$ENVNAME' python -c 'import torch' 2>/dev/null || conda run -n '$ENVNAME' pip install torch --index-url https://download.pytorch.org/whl/cpu"
  run "conda run -n '$ENVNAME' pip install -r '$TOOLS/discotope3_web/requirements.txt'"
  say "확인:  conda run -n $ENVNAME python $TOOLS/discotope3_web/src/predict_webserver.py --help"
  say "쓰기:  python scripts/run_discotope.py --target 8k3k_D --repo $TOOLS/discotope3_web \\"
  say "             --python \$(conda run -n $ENVNAME which python)"
fi

# ── AbEpiTope-1.0 ─────────────────────────────────────────────────────────────
if want abepitope; then
  echo; say "── AbEpiTope-1.0 ──"
  say "⚠️ 이 도구는 '항원 표면'이 아니라 '예측된 복합체'를 채점한다 → 쓰는 자리가 다르다."
  say "   후보 자리를 고르는 데(②)가 아니라, 유도 도킹 결과를 재랭킹하는 데(④ 이후) 쓴다."
  if [ -d "$TOOLS/AbEpiTope-1.0/.git" ]; then
    run "git -C '$TOOLS/AbEpiTope-1.0' pull --ff-only"
  else
    run "git clone https://github.com/mnielLab/AbEpiTope-1.0 '$TOOLS/AbEpiTope-1.0'"
  fi
  run "conda run -n '$ENVNAME' pip install -e '$TOOLS/AbEpiTope-1.0'"
  say "⚠️ 이 도구의 정확한 커맨드는 판올림마다 달라진다 — 설치 후 README 를 확인하고,"
  say "   후보별 점수를 cand,score 두 열 CSV 로 만들어 rank_sites.py --method abepitope --scores 에 넘길 것."
fi

echo
if [ $APPLY -eq 1 ]; then
  say "끝. 설치된 것은 env '$ENVNAME' 과 $TOOLS 아래뿐이다(기존 env 무영향)."
else
  say "dry-run 끝. 실제 설치는 --apply."
fi
