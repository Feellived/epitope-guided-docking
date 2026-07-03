# 기록 자동화 (Ab-Ag docking 실험)

2층 구조: **서버는 CSV 자동 기록(토큰 X)**, **Notion 동기화는 토큰 있는 곳(개인 Mac)에서**.

```
[서버, conda dock]  eval_boltz.py  -->  results/summary.csv   (자동 append, source of truth)
[개인 Mac]          sync_notion.py -->  Notion 실험 로그 DB    (CSV의 새 행만 push)
```

## 파일
- `log_lib.py`      : summary.csv 한 줄 append (import용)
- `eval_boltz.py`   : Boltz 예측 1건 -> DockQ+신뢰도 파싱 -> CSV 기록 (서버)
- `sync_notion.py`  : CSV -> Notion DB 동기화 (개인 Mac)
- `results/summary.csv` : 숫자 원본 (git 커밋). `results/.notion_synced.txt` 는 gitignore.

## 서버에서 (예측 후, 자동 기록)
```bash
conda activate dock
cd ~/projects/bk21-antibody-ml/consensus_docking/runs/9vp0     # 예측 결과가 있는 폴더
python ../../scripts/eval_boltz.py --target 9VP0 \
    --targets ../../scripts/targets.yaml \
    --seed 0 --phase Preliminary --model Boltz-2 \
    --csv ../../results/summary.csv
# targets.yaml 못 읽으면: 위 명령에 --antigen A --antibody H,L 추가
```
-> `results/summary.csv` 에 한 줄 자동 추가. 여러 seed면 `--seed 0..N` 반복.

## 개인 Mac에서 (Notion 동기화)
```bash
export NOTION_API_TOKEN=ntn_xxx      # 채팅/서버/깃 아님. 이 셸에서만.
cd <repo>/consensus_docking
python scripts/sync_notion.py --csv results/summary.csv
```
-> summary.csv 의 아직 안 올린 행만 Notion 실험 로그 DB로 push.

## 주의
- **토큰은 절대 서버·깃·문서에 저장 금지.** sync_notion 은 개인 Mac 환경변수에서만 읽음.
- eval_boltz 는 지금 **Boltz 출력 형식** 기준. Chai-1/Protenix 등은 confidence/DockQ 형식이 달라
  각 모델용 파서(eval_chai.py 등)를 모델 추가할 때 붙인다. (CSV/Notion 파이프는 공용.)
- 첫 실행 시 9VP0 결과가 preliminary 값(DockQ~0.34/0.46, ipTM 0.84)과 맞는지로 파서 검증.
