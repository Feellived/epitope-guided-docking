#!/usr/bin/env bash
# 각 타깃 항원서열(chains.json)로 over-rep 지도 → results/epitope_overrep_<T>.csv (analyze_hotspot_drift hotspot).
# CPU/RCSB만 = GPU 예측과 병렬 OK. 진단형: 타깃별 서열길이·검색개수·usable·실패이유 출력.
# 사용: RUN=runs_diverse bash scripts/run_overrep_maps.sh
#      IDENTITY=0.7 로 완화(비-RBD isoform 포함). 특정 타깃만: ... run_overrep_maps.sh 9W43,9XF9
cd "$(dirname "$0")/.." || exit 1
RUN="${RUN:-runs_diverse}"
IDENTITY="${IDENTITY:-0.7}"; MINAG="${MINAG:-0.6}"; MAX="${MAX:-300}"; MAXRES="${MAXRES:-3.5}"
TARGETS="${1:-9K6J,9XF9,9Y0A,9VXL,9LY5,9MQR,9KKJ,9W43,9YC6,9J4C}"
mkdir -p results
for T in $(echo "$TARGETS" | tr ',' ' '); do
  low=$(echo "$T" | tr 'A-Z' 'a-z'); d="$RUN/$low"
  echo "==================== $T ===================="
  if [ ! -f "$d/chains.json" ]; then echo "  !! chains.json 없음: $d"; continue; fi
  # 항원 서열 추출 (실패 이유 노출)
  L=$(python3 - "$d/chains.json" "/tmp/ag_$T.fasta" <<'PY'
import json,sys
try:
    d=json.load(open(sys.argv[1])); ag=str(d["antigen"]); s={c["id"]:c["seq"] for c in d["chains"]}[ag]
    open(sys.argv[2],"w").write(">ag\n"+s+"\n"); print(len(s))
except Exception as e: print("ERR:"+type(e).__name__+":"+str(e))
PY
)
  case "$L" in ERR:*) echo "  !! 항원서열 추출 실패 ($L) — chains.json의 antigen 필드 확인"; continue;; esac
  echo "  항원 길이 ${L}aa → over-rep 검색 (identity>=$IDENTITY, max_res<=$MAXRES, max=$MAX)"
  python scripts/epitope_overrep_map.py --ref-fasta "/tmp/ag_$T.fasta" --offset 0 --identity "$IDENTITY" \
    --min-ag-ident "$MINAG" --max "$MAX" --max-res "$MAXRES" --tag "$T" --out "results/epitope_overrep_$T.csv" 2>&1 \
    | grep -E '검색 실패|검색 0개|대상 PDB|저장|^  [0-9]+/' | tail -4
  # 결과 요약: usable 잔기 수(n_complexes>=1) + 최대 n
  if [ -f "results/epitope_overrep_$T.csv" ]; then
    python3 - "results/epitope_overrep_$T.csv" <<'PY'
import csv,sys
rows=list(csv.DictReader(open(sys.argv[1])))
ns=[int(r["n_complexes"]) for r in rows if r["n_complexes"].strip()]
pos=[n for n in ns if n>0]; mx=max(ns) if ns else 0
print(f"  → hotspot 잔기 {len(pos)}개, 최대결합 n={mx}  {'(빈 지도!)' if mx==0 else ''}")
PY
  fi
done
echo; echo "완료 → results/epitope_overrep_<T>.csv"
