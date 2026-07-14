#!/usr/bin/env bash
# [서버 전용] 실물 예측 폴더명 통일: results_boltz->out_boltz, results_colabfold->out_colabfold.
# runs_diverse/<id>{,_fv}/ 안에서만. mv(rename)뿐 — 같은 부모로의 이름변경, 데이터 삭제/이동 없음(reversible).
# reader(diagnose/where_docked/merge_dockq/eval_boltz)는 *boltz/*colabfold suffix로 양쪽 호환 →
# 이 rename은 '정리'일 뿐, 안 해도 Phase 0는 정상 동작.
#
# ★ 기본 dry-run(미리보기). 실제 적용은 --apply. 공용계정(admuser) — 내 폴더 밖 절대 안 건드림.
# 사용(서버):
#   cd ~/projects/bk21-antibody-ml/consensus_docking
#   bash scripts/migrate_server_outdirs.sh            # 미리보기
#   bash scripts/migrate_server_outdirs.sh --apply    # 실제 rename
#   RUN=runs bash scripts/migrate_server_outdirs.sh   # 옛 해커톤 세트에도 적용하려면
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
RUN="${RUN:-runs_diverse}"
APPLY=0; [ "${1:-}" = "--apply" ] && APPLY=1

# old -> new (순서 고정 위해 배열)
OLD=(results_boltz results_colabfold)
NEW=(out_boltz     out_colabfold)

[ -d "$RUN" ] || { echo "!! $RUN 없음 (consensus_docking에서 실행했는지 확인)"; exit 1; }
echo "대상 루트: $RUN   모드: $([ $APPLY = 1 ] && echo APPLY || echo DRY-RUN)"
n=0; skip=0
for d in "$RUN"/*/; do
  [ -d "$d" ] || continue
  for i in "${!OLD[@]}"; do
    o="$d${OLD[$i]}"; w="$d${NEW[$i]}"
    [ -d "$o" ] || continue
    if [ -e "$w" ]; then
      echo "[skip] ${w} 이미 존재 → 수동 확인"; skip=$((skip+1)); continue
    fi
    if [ $APPLY = 1 ]; then
      mv "$o" "$w"; echo "[mv]   $o -> ${NEW[$i]}"
    else
      echo "[dry]  $o -> ${NEW[$i]}"
    fi
    n=$((n+1))
  done
done
echo "----"
[ $APPLY = 1 ] && echo "적용 완료: $n건 rename ($skip건 skip)." \
               || echo "미리보기 $n건 ($skip건 skip). 실제 적용: --apply"
echo "확인: find $RUN -maxdepth 2 -type d \\( -name 'out_boltz' -o -name 'out_colabfold' \\) | sort | head"
