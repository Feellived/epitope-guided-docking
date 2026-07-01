#!/usr/bin/env bash
# Boltz-2 로 9VP0 예측.  사용: bash 02_run_boltz.sh
# (tmux 세션 + conda activate dock 상태에서 실행)
set -e

boltz predict boltz_9vp0.yaml \
    --use_msa_server \
    --out_dir results_boltz \
    --no_kernels 2>&1 | tee boltz.log

echo "=== 예측 파일 위치 ==="
find results_boltz -name "*_model_0.cif"
find results_boltz -name "confidence_*model_0.json"
