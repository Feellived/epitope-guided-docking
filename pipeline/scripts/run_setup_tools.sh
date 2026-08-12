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
#   bash scripts/run_setup_tools.sh                 # 무엇을 할지만 출력
#   bash scripts/run_setup_tools.sh --apply
#   TOOLS=~/tools ENVNAME=discotope bash scripts/run_setup_tools.sh --apply
#   bash scripts/run_setup_tools.sh --apply --only discotope
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
  say "쓰기:  python scripts/select_score_discotope.py --target 8k3k_D --repo $TOOLS/discotope3_web \\"
  say "             --python \$(conda run -n $ENVNAME which python)"
fi

# ── AbEpiTope-1.0 ─────────────────────────────────────────────────────────────
# 공식 설치법(레포 README, 2026-07-29 확인) = uv 로 레포 안에 자체 .venv 를 만든다.
# conda env 를 쓰지 않으므로 위에서 만든 $ENVNAME 과 무관하다.
# 논문 = Science Advances 2025. 두 점수 중 우리가 쓰는 것은 AbEpiScore-1.0
# (모델링된 항체-항원 계면이 얼마나 정확한가). AbEpiTarget-1.0 은 '어느 항체가 이 항원에
# 붙을까'라 우리 과제와 다르다.
if want abepitope; then
  echo; say "── AbEpiTope-1.0 ──"
  say "⚠️ 이 도구는 '항원 표면'이 아니라 '예측된 복합체'를 채점한다 → 후보마다 구조가 하나씩 필요하다"
  say "   (그래서 run_prescreen_abepitope.sh 가 후보별로 1샘플씩 먼저 예측한다)"
  if [ -d "$TOOLS/AbEpiTope-1.0/.git" ]; then
    say "이미 받아져 있음 — 갱신만"
    run "git -C '$TOOLS/AbEpiTope-1.0' pull --ff-only"
  else
    run "git clone https://github.com/mnielLab/AbEpiTope-1.0.git '$TOOLS/AbEpiTope-1.0'"
  fi
  # uv (파이썬 패키지 관리자). 이미 있으면 다시 깔지 않는다.
  if command -v uv >/dev/null 2>&1; then
    say "uv 이미 있음: $(command -v uv)"
  else
    run "curl -LsSf https://astral.sh/uv/install.sh | sh"
    say "※ uv 는 ~/.local/bin 에 깔린다. PATH 에 없으면 export PATH=\"\$HOME/.local/bin:\$PATH\""
  fi
  run "cd '$TOOLS/AbEpiTope-1.0' && uv sync"
  # CUDA 런타임 경로 — README 가 activate 에 덧붙이라고 한다. 두 번 붙지 않게 확인 후.
  if [ $APPLY -eq 1 ] && [ -f "$TOOLS/AbEpiTope-1.0/.venv/bin/activate" ] \
     && ! grep -q "AbEpiTope CUDA runtime libraries" "$TOOLS/AbEpiTope-1.0/.venv/bin/activate"; then
    say "CUDA 라이브러리 경로를 .venv/bin/activate 에 추가"
    cat >> "$TOOLS/AbEpiTope-1.0/.venv/bin/activate" <<'EOF'

# AbEpiTope CUDA runtime libraries
export LD_LIBRARY_PATH="$VIRTUAL_ENV/lib/python3.9/site-packages/nvidia/cusparse/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$VIRTUAL_ENV/lib/python3.9/site-packages/nvidia/cublas/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$VIRTUAL_ENV/lib/python3.9/site-packages/nvidia/cuda_runtime/lib:$LD_LIBRARY_PATH"
EOF
  else
    say "(dry-run 이거나 이미 들어 있음) .venv/bin/activate 에 CUDA 경로 추가"
  fi
  say "확인:  bash -c 'source $TOOLS/AbEpiTope-1.0/.venv/bin/activate && python -c \"from abepitope.main import StructureData, EvalAbAgs; print(1)\"'"
  say "쓰기:  python scripts/eval_abepitope.py --target 8k3k_D --repo $TOOLS/AbEpiTope-1.0 \\"
  say "             --dir runs_sites_guided/8k3k_d/all"
fi

echo
if [ $APPLY -eq 1 ]; then
  say "끝. 설치된 것은 env '$ENVNAME' 과 $TOOLS 아래뿐이다(기존 env 무영향)."
else
  say "dry-run 끝. 실제 설치는 --apply."
fi
