#!/usr/bin/env bash
# [Study B/E1] 에피토프 과대표집 지도 (Foldseek 결합모드 클러스터링 = 비중복).
#   prep(복합체 추출) → foldseek easy-multimersearch → cluster_and_map(클러스터+지도).
#   raw(복합체별) vs clustered(비중복) 둘 다 출력 → 중복이 얼마나 부풀렸나 대비.
#
# Foldseek = 정적바이너리(시스템 설치 아님). 없으면 프로젝트 폴더에 받기(서버 부담 최소):
#   wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz && tar xzf foldseek-linux-avx2.tar.gz
#   → 바이너리 = ./foldseek/bin/foldseek
#
# 사용(RBD 기본, 400개):
#   conda activate boltz     # biopython + 인터넷
#   bash scripts/run_epitope_overrep.sh                    # RBD, max 400
#   bash scripts/run_epitope_overrep.sh 200                # 개수 조절(스모크)
#   bash scripts/run_epitope_overrep.sh 400 other_ag.fasta # 다른 항원(기준서열 FASTA)
# ⚠️ multimersearch가 easy-multimercluster 대체(후자는 macOS 크래시). 리눅스선 둘 다 되나 이게 OS-robust.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
FOLDSEEK="${FOLDSEEK:-./foldseek/bin/foldseek}"
MAX="${1:-400}"; REF="${2:-}"
REFARG=""; [ -n "$REF" ] && REFARG="--ref-fasta $REF"
[ -x "$FOLDSEEK" ] || { echo "!! foldseek 실행파일 없음: $FOLDSEEK"; echo "   wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz && tar xzf foldseek-linux-avx2.tar.gz"; exit 1; }
mkdir -p results

echo "[1/3] 복합체 추출 (max $MAX) → complexes/"
python scripts/prep_complexes.py --max "$MAX" --out-dir complexes $REFARG || exit 1
n=$(ls complexes/*.pdb 2>/dev/null | wc -l); echo "   추출 $n개"
[ "$n" -ge 3 ] || { echo "!! 복합체 부족($n) → 중단"; exit 1; }

echo "[2/3] Foldseek easy-multimersearch (복합체 all-vs-all)"
rm -rf tmp_fs msearch msearch_report 2>/dev/null
"$FOLDSEEK" easy-multimersearch complexes/ complexes/ msearch tmp_fs > results/foldseek.log 2>&1 || { echo "!! multimersearch 실패 → results/foldseek.log"; exit 1; }
[ -f msearch_report ] || { echo "!! msearch_report 없음"; exit 1; }

echo "[3/3] 결합모드 클러스터 + 과대표집 지도"
python scripts/cluster_and_map.py --report msearch_report --complexes complexes $REFARG \
  --out results/epitope_overrep_clustered.csv || exit 1

echo
echo "완료 → results/epitope_overrep_clustered.csv (auth, aa, raw_n, clustered_n, region)"
echo "raw_n=복합체별(중복포함) / clustered_n=비중복 결합모드. 상위 잔기는 위에 출력됨."
