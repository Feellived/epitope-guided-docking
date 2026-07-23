#!/usr/bin/env bash
# HADDOCK3 ab-initio (blind, CDR active + 항원 surface passive) on diverse 10 Fab (= Exp2 consensus 세트).
# 목적: 물리 arm(HADDOCK)을 co-folder 5종 + ZDOCK과 '같은 세트'에 올려 consensus에 나란히 편입.
#       (기존 HADDOCK은 RBD off-hotspot 세트 = Exp4/Study B 소관. 이건 diverse용, haddock_diverse/로 분리.)
# ⚠️ CPU 전용 — GPU sweep과 병행 안전. 단 무겁다: 타깃당 rigidbody sampling + 200 pose flexref/emref
#    → 코어수·sampling 의존, 타깃당 수십분~수시간. 먼저 1타깃 smoke 권장.
# 사용(consensus_docking/에서):
#   bash scripts/run_haddock_diverse.sh 9K6J                 # smoke 1타깃
#   bash scripts/run_haddock_diverse.sh                      # 전체 10타깃 (기본 sampling 1000)
#   SAMPLING=1000 ANARCI_ENV=anarci HADDOCK_ENV=haddock3 bash scripts/run_haddock_diverse.sh
set -uo pipefail
CD=$(cd "$(dirname "$0")/.." && pwd); cd "$CD"
SAMPLING=${SAMPLING:-1000}
ANARCI_ENV=${ANARCI_ENV:-anarci}
HADDOCK_ENV=${HADDOCK_ENV:-haddock3}
RUNDIR=${RUNDIR:-runs_diverse}
OUTBASE=${OUTBASE:-haddock_diverse}
TARGETS=("$@"); [ ${#TARGETS[@]} -eq 0 ] && TARGETS=(9K6J 9Y0A 9LY5 9KKJ 9J4C 9XF9 9VXL 9MQR 9W43 9YC6)

# haddock3 env 확인
if ! conda run -n "$HADDOCK_ENV" which haddock3 >/dev/null 2>&1; then
  echo "!! HADDOCK_ENV='$HADDOCK_ENV'에 haddock3 없음 — 올바른 env로 HADDOCK_ENV 지정 (예: HADDOCK_ENV=haddock3.5)"; exit 1
fi
echo "HADDOCK diverse: RUNDIR=$RUNDIR OUT=$OUTBASE sampling=$SAMPLING (anarci=$ANARCI_ENV haddock=$HADDOCK_ENV)"
ok=0; fail=0
for T in "${TARGETS[@]}"; do
  t=$(echo "$T" | tr 'A-Z' 'a-z'); D="$OUTBASE/$t"
  echo "==== $T → $D ===="
  mkdir -p "$D"
  # 1) native → 항체(H+L)/항원 분리
  if ! python scripts/prep_haddock.py --run-dir "$RUNDIR" --target "$T" --out-dir "$D"; then echo "  [x] prep_haddock $T"; ((fail++)); continue; fi
  # 2) ANARCI CDR active + 단일사슬 항체 (anarci env)
  if ! conda run -n "$ANARCI_ENV" python scripts/haddock_ab_prep.py --dir "$D"; then echo "  [x] haddock_ab_prep $T"; ((fail++)); continue; fi
  # 3) 항원 passive + ambig.tbl + run.cfg (haddock3 env)
  if ! conda run -n "$HADDOCK_ENV" python scripts/make_haddock_restraints.py --run-dir "$RUNDIR" --target "$T" --dir "$D" --sampling "$SAMPLING"; then echo "  [x] restraints $T"; ((fail++)); continue; fi
  # 4) HADDOCK3 도킹 (rigidbody→seletop→flexref→emref→clustfcc→seletopclusts)
  if ! conda run -n "$HADDOCK_ENV" haddock3 "$D/run.cfg"; then echo "  [x] haddock3 $T"; ((fail++)); continue; fi
  echo "  [ok] $T → $D/run/"; ((ok++))
done
echo "완료: ok=$ok fail=$fail / ${#TARGETS[@]}"
echo "채점: HADDOCK pose = $OUTBASE/<t>/run/*/  → DockQ vs native(runs_diverse). haddock_dockq.py를 --run-dir runs_diverse --haddock-dir $OUTBASE 로 확장 필요(다음)."
