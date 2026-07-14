#!/usr/bin/env bash
# 각 타깃 항원서열(chains.json)로 over-rep 지도 → results/epitope_overrep_<T>.csv
# (analyze_hotspot_drift.py가 hotspot 정의에 씀). CPU/RCSB만 = GPU 예측과 병렬 OK.
# 사용: RUN=runs_diverse bash scripts/run_overrep_maps.sh
cd "$(dirname "$0")/.." || exit 1
RUN="${RUN:-runs_diverse}"
TARGETS="${1:-9K6J,9XF9,9Y0A,9VXL,9LY5,9MQR,9KKJ,9W43,9YC6,9J4C}"
mkdir -p results
for T in $(echo "$TARGETS" | tr ',' ' '); do
  low=$(echo "$T" | tr 'A-Z' 'a-z'); d="$RUN/$low"
  [ -f "$d/chains.json" ] || { echo "skip $T (chains.json 없음)"; continue; }
  python3 -c "import json;d=json.load(open('$d/chains.json'));ag=str(d['antigen']);s={c['id']:c['seq'] for c in d['chains']}[ag];open('/tmp/ag_$T.fasta','w').write('>ag\n'+s+'\n')" || continue
  echo "=== $T 항원 over-rep 지도 ==="
  python scripts/epitope_overrep_map.py --ref-fasta "/tmp/ag_$T.fasta" --offset 0 --tag "$T" --max 200 --out "results/epitope_overrep_$T.csv" 2>&1 | tail -2
done
echo "완료 → results/epitope_overrep_<T>.csv"
