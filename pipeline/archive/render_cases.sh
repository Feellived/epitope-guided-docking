#!/usr/bin/env bash
# 의미있는 도킹 성패 사례를 PyMOL로 렌더 (모델별 1~2개; best-DockQ pose 자동선택).
# 실행: pipeline/ 에서   conda activate pymol && bash scripts/render_cases.sh
# 필요: pymol-open-source, results/summary.csv, runs/<id>/ 예측 출력.
set -u
CSV=results/summary.csv
OUT=fig
mkdir -p "$OUT"

render () {  # model target run-dir-low
  local model="$1" target="$2" low="$3"
  local safe; safe=$(echo "${model}_${target}" | tr 'A-Z /' 'a-z__')
  echo "=== $model $target ==="
  pymol -cq scripts/plot_pose.py -- \
    --run-dir "runs/$low" --target "$target" --model "$model" \
    --from-csv "$CSV" --out "$OUT/${safe}.png"
}

# A. 성공 typology — "모델은 도킹을 할 수 있다"
render Protenix     10ZO 10zo    # 0.89 최고 성공 (나노바디, 교과서적)
render Boltz-2      10ZO 10zo    # 0.84 성공

# B. consensus의 핵심 — 복합체마다 챔피언이 다르다 (오라클 5/5)
render tFold-Ag     9K6J 9k6j    # 0.49 ⭐ 유일 크랙 (항체 MSA-free가 이김)
render AF2-Multimer 9XF9 9xf9    # 0.23 유일 크랙
render Chai-1       9Y0A 9y0a    # 0.56 유일 Fab 크랙
render Chai-1       9OBN 9obn    # 0.73 나노바디 챔피언

# C. ipTM는 성패를 못 가른다 — confident-wrong (B의 tFold 9K6J와 짝)
render Boltz-2      9K6J 9k6j    # 0.03 실패인데 ipTM 0.81 (자신만만 오답)

# D. MSA 절제실험 — 같은 복합체 on/off
render Chai-1       10ZO 10zo    # 0.78 MSA on
render Chai-1-noMSA 10ZO 10zo    # 0.08 MSA off (붕괴)

# E. 실패 typology
render Protenix     9Y0A 9y0a    # 0.05 완전 빗나감
render tFold-Ag     9OBN 9obn    # 0.07 (남들은 성공한 나노바디를 flop)

echo "완료 → $OUT/*.png (+ .pse). Mac으로 내려받아 확인/발표에 사용."
