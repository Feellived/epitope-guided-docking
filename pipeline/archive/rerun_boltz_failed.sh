#!/usr/bin/env bash
# 실패한 Boltz 타깃 재실행 (MSA 서버 truncation으로 미완). 기본 = 9ly5 9w43 (Fab).
# 원인은 서열이 아니라 MSA 서버 응답 잘림(같은 서열을 Chai·Protenix는 성공). 아침엔 서버 덜 붐벼 성공률↑.
# 사용: conda activate boltz; bash scripts/rerun_boltz_failed.sh [9ly5 9w43 ...]
# ⚠️ GPU 사용 — AF2-M/tFold 스모크가 GPU 점유 중이면 끝난 뒤 실행(동시=OOM 위험).
cd "$(dirname "$0")/.." || exit 1
TARGETS="${*:-9ly5 9w43}"
CACHE=/mnt/data/admuser/boltz
for t in $TARGETS; do
  d=runs_diverse/$t
  [ -d "$d" ] || { echo "[skip] $d 없음"; continue; }
  echo "=== rerun $t ==="
  rm -rf "$d/results_boltz" "$d/out_boltz"                          # 이전 실패 잔재
  find "$d" \( -name '*.a3m' -o -path '*msa*' -o -name '*.npz' \) -delete 2>/dev/null  # 깨진 MSA 캐시
  ( cd "$d" && boltz predict "boltz_$t.yaml" --use_msa_server --cache "$CACHE" \
      --out_dir results_boltz --no_kernels --diffusion_samples 5 > boltz.log 2>&1 )
  n=$(find "$d/results_boltz" -name '*_model_*.cif' 2>/dev/null | wc -l)
  if [ "$n" -ge 5 ]; then
    echo "  OK  $t: $n cif"
  else
    echo "  FAIL $t: $n cif — tail boltz.log:"; tail -3 "$d/boltz.log"
    echo "  (truncated면 MSA 서버 혼잡 → 1~2분 뒤 다시 실행)"
  fi
done
echo "완료. eval은: python scripts/eval_boltz.py --target <T> --run-dir runs_diverse/<t> --targets scripts/targets_benchmark_v2.yaml --seed 0..4 --phase Diverse-Fab --model Boltz-2 --csv results/summary.csv"
