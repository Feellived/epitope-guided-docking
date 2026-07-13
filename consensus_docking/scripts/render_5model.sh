#!/usr/bin/env bash
# 5모델 확장 + 병합-DockQ 주요 사례 PyMOL 렌더 → fig_5model/*.png
# 스토리: (N1) min이 가린 진짜 성공  (N2) 견고한 rescue  (N3) near-miss.
# 실행: consensus_docking/ 에서  conda activate pymol && bash scripts/render_5model.sh
# 필요: pymol-open-source, results/summary.csv, runs_diverse/<id>/ 예측 출력.
set -u
CSV=results/summary.csv
OUT=fig_5model
mkdir -p "$OUT"

r(){  # model target phase run-dir out-name label
  echo "=== $5 : $1 $2 ($3) ==="
  pymol -cq scripts/pymol_pose.py -- \
    --run-dir "$4" --target "$2" --model "$1" --phase "$3" \
    --from-csv "$CSV" --label "$6" --out "$OUT/$5.png"
}

# N1. min이 가린 진짜 성공 — 9YC6 Protenix (병합 DockQ 0.366, min은 0.12)
r Protenix-base 9YC6 Diverse-Fab runs_diverse/9yc6    N1_9yc6_protenix_fab "9YC6 Protenix Fab  merged DockQ 0.366 (SUCCESS, min 0.12 hid it)"

# N2. 견고한 rescue — 9K6J tFold Fv (5 seed 모두 성공)
r tFold-Ag      9K6J Diverse-Fv  runs_diverse/9k6j_fv N2_9k6j_tfold_fv     "9K6J tFold-Ag Fv  DockQ 0.53 (robust, 5/5 seed)"

# N3. near-miss — 9XF9 AF2-Multimer Fv (병합 0.21)
r AF2-Multimer  9XF9 Diverse-Fv  runs_diverse/9xf9_fv N3_9xf9_af2m_fv      "9XF9 AF2-Multimer Fv  merged 0.21 (near-miss)"

echo "완료 → $OUT/*.png"
echo "커밋: git add consensus_docking/$OUT && git commit -m '5모델 병합-DockQ PyMOL 렌더' && git push"
