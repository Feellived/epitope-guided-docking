#!/usr/bin/env bash
# 예측 vs 정답 DockQ + 신뢰도 출력.
# 사용: bash 03_eval.sh [예측.cif]   (인자 없으면 Boltz 결과 자동 탐색)
set -e

CIF=${1:-$(find results_boltz -name "*_model_0.cif" | head -1)}
JSON=$(find results_boltz -name "confidence_*model_0.json" | head -1)

echo "### DockQ (예측 vs 9VP0.pdb) ###"
echo "# 볼 것: A-H, A-L interface (=항체-항원). H-L 은 항체 내부라 무시."
DockQ "$CIF" 9VP0.pdb

echo
echo "### 신뢰도 (global ipTM / pair_chains_iptm) ###"
echo "# 항원 A(보통 인덱스 0) vs H,L 의 pair ipTM 이 진짜 계면 신뢰도."
python -m json.tool "$JSON"
