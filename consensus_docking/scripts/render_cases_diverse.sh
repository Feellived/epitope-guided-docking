#!/usr/bin/env bash
# 다양-항원 3모델(Boltz·Chai·Protenix) 주요 사례 PyMOL 렌더 → fig_diverse/*.png
# 스토리: (A) 단일-모델 rescue  (B) Fv가 성공을 파괴  (C) ipTM 오선택.
# 실행: consensus_docking/ 에서   conda activate pymol && bash scripts/render_cases_diverse.sh
# 필요: pymol-open-source, results/summary.csv, runs_diverse/<id>/ 예측 출력.
set -u
CSV=results/summary.csv
OUT=fig_diverse
mkdir -p "$OUT"

render () {  # model target phase run-dir out-name
  local model="$1" target="$2" phase="$3" rd="$4" name="$5"
  echo "=== $name : $model $target ($phase) ==="
  pymol -cq scripts/pymol_pose.py -- \
    --run-dir "$rd" --target "$target" --model "$model" \
    --from-csv "$CSV" --phase "$phase" --out "$OUT/${name}.png"
}

# A. 단일-모델 rescue — "다 실패했는데 이 모델만 성공"
render Boltz-2       9K6J Diverse-Fab runs_diverse/9k6j     A1_9k6j_boltz_fab_0652     # Boltz만 성공
render Chai-1        9Y0A Diverse-Fab runs_diverse/9y0a     A2_9y0a_chai_fab_0662      # Chai만(PLM)
render Protenix-base 9LY5 Diverse-Fab runs_diverse/9ly5     A3_9ly5_protenix_fab_0871  # Protenix만, 최고점

# B. Fv가 성공을 파괴 / 챔피언 이동 (A3와 짝)
render Protenix-base 9LY5 Diverse-Fv  runs_diverse/9ly5_fv  B1_9ly5_protenix_fv_0058   # 0.871→0.058 붕괴
render Protenix-base 9K6J Diverse-Fv  runs_diverse/9k6j_fv  B2_9k6j_protenix_fv_0776   # Fv가 챔피언을 Boltz→Protenix로

# C. ipTM 오선택 — 다 실패인데 자신만만
render Boltz-2       9W43 Diverse-Fab runs_diverse/9w43     C1_9w43_boltz_fab_iptm085  # DockQ 0.067 / ipTM 0.854

echo "완료 → $OUT/*.png  (Mac에서 확인/보고서용). 커밋: git add consensus_docking/$OUT && git commit && git push"
