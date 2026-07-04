# 5개 복합체 × 모든 core 모델 실행 가이드

흐름: **① 프렙(1회) → ② 모델별 설치(허가) → ③ 5개 실행 → ④ eval → CSV → Notion.**
core 3 = Boltz-2 · Chai-1 · Protenix-v2. **각 모델은 별도 conda env.**

벤치 5개 = 10zo 9y0a 9obn 9k6j 9xf9  (9vp0 = dev/preliminary).

## ① 추출 + 입력생성 (계층 분리)

**Layer 1 — 사슬 추출 (한 번, 모델 무관):**
```bash
conda activate dock
cd ~/projects/bk21-antibody-ml/consensus_docking
pip install biopython pyyaml        # 없으면
python scripts/extract_chains.py --targets scripts/targets.yaml --out-dir runs \
    --only 10ZO,9Y0A,9OBN,9K6J,9XF9
```
-> runs/<id>/{<ID>.pdb, **chains.json**(서열·역할 단일 기록), <id>.fasta}

**Layer 2 — 모델 입력 생성 (그 모델 쓸 때만, 독립):**
```bash
python scripts/make_boltz.py    --run-dir runs   # Boltz 쓸 때
python scripts/make_chai.py     --run-dir runs   # Chai 쓸 때
python scripts/make_protenix.py --run-dir runs   # Protenix 쓸 때
```
-> 각 runs/<id>/에 boltz_<id>.yaml / chai_<id>.fasta / protenix_<id>.json 만 생성

## ② + ③ 모델별 실행 (각 별도 env; 설치 전 사용자 허가)

### Boltz-2  (env: dock, 이미 있음)
```bash
conda activate dock
for t in 10zo 9y0a 9obn 9k6j 9xf9; do
  ( cd runs/$t && boltz predict boltz_$t.yaml --use_msa_server \
      --out_dir results_boltz --no_kernels --diffusion_samples 5 > boltz.log 2>&1 )
done
```

### Chai-1  (새 env)
```bash
conda create -n chai python=3.11 -y && conda activate chai
pip install chai_lab==0.6.1
export CHAI_DOWNLOADS_DIR=~/chai_weights      # 본인 폴더
for t in 10zo 9y0a 9obn 9k6j 9xf9; do
  ( cd runs/$t && chai-lab fold chai_$t.fasta out_chai --use-msa-server > chai.log 2>&1 )
done
# ⚠️ 정확한 CLI/플래그는 `chai-lab fold --help`로 확인 후 조정
```

### Protenix-v2  (새 env; 서브커맨드 = `pred`, 커널 컴파일 회피)
```bash
conda create -n protenix python=3.11 -y && conda activate protenix
pip install protenix
export PROTENIX_ROOT_DIR=/mnt/data/admuser/protenix_weights   # weight 큰 디스크
export LAYERNORM_TYPE=torch                                   # CUDA 커널(LayerNorm) 컴파일 회피 = gcc 불필요
for t in 10zo 9y0a 9obn 9k6j 9xf9; do
  ( cd runs/$t && rm -rf out_protenix && protenix pred -i protenix_$t.json -o out_protenix \
      -n protenix-v2 -s 101 -e 5 --trimul_kernel torch --triatt_kernel torch \
      --msa_server_mode colabfold > protenix.log 2>&1 )
  echo "=== done $t ==="
done
# 서브커맨드 = `protenix pred`(predict 아님). -e 5=pose 5개, 모든 커널 torch로(컴파일 X). eval은 boltz env.
```

## ④ eval -> CSV -> Notion
```bash
# Boltz (자동 파서 완성):
for t in 10zo 9y0a 9obn 9k6j 9xf9; do
  for s in 0 1 2 3 4; do
    python scripts/eval_boltz.py --target ${t^^} --run-dir runs/$t \
       --targets scripts/targets.yaml --seed $s --phase Benchmark --model Boltz-2 \
       --csv results/summary.csv
  done
done
# 개인 Mac에서 Notion 반영:
#   export NOTION_API_TOKEN=ntn_xxx ; python scripts/sync_notion.py --csv results/summary.csv
```

## 메모
- **Chai/Protenix eval 파서(eval_chai.py·eval_protenix.py)는 첫 실행 출력 확인 후 작성.**
  (Chai=scores .npz의 per_chain_pair_iptm / Protenix=summary JSON의 chain_pair_iptm; 출력 사슬이
  A,B,C로 리네임될 수 있어 DockQ `--mapping` 필요 → 실제 출력 보고 맞춤.)
- 공정비교: 세 모델 다 `--use_msa_server` 계열(같은 원격 MSA)·같은 seed 예산(5)·같은 DockQ(≥0.23).
- 무거운 작업이니 **tmux 안에서**, 로그는 각 `*.log`로.
