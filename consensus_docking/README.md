# Consensus Docking (제4회 인공지능 신약개발 해커톤)

항체-항원 도킹/구조예측 툴들을 한 파이프라인으로 묶고, 학습된 consensus 층으로
신뢰할 수 있는 복합체 구조를 골라내는 프레임워크.

> ⚠️ **메인 프로젝트(epitope binning / CD16a AF3 rescue)와 별개.**
> 메인 데이터는 이 레포 밖(`~/projects/bk21-epitope-binning/`)으로 분리함.
> 이 폴더는 consensus docking 해커톤 전용.

## 골든 케이스
leakage-free **9VP0** (Fab + IL-31, X-ray 2.1A, 2025 공개) 하나를 모든 모델의
고정 입력으로 두고 (신뢰도, DockQ)를 뽑아 비교표를 채운다. (configs/targets.yaml)

## 실행 (서버: admuser, ssh 직접실행형)
```bash
tmux new -s hj_dock
conda activate dock
mkdir -p ~/dock/9vp0 && cd ~/dock/9vp0
# 이 폴더 scripts/ 를 여기로 (git clone 또는 복사)

python 01_prep_9vp0.py     # 9VP0.pdb + fasta + boltz yaml 생성
bash   02_run_boltz.sh     # Boltz-2 예측
bash   03_eval.sh          # DockQ + 신뢰도
```

## 사슬 배치 (9VP0)
A = IL-31 항원(144) / H = Fab 중쇄(217) / L = Fab 경쇄(222)
- 항체-항원 interface = **A-H, A-L** (DockQ에서 이것만 성공지표). H-L 은 항체 내부라 무시.
- pair_chains_iptm 에서 **항원 A(보통 인덱스 0) vs H,L** 값을 본다.

## 비교표 (모델마다 한 줄)
| 모델 | global ipTM | Ab-Ag pair ipTM | DockQ(Ab-Ag) | 성공(>=0.23) |
|---|---|---|---|---|
| Boltz-2 | | | | |
| Chai-1 | | | | |
| Protenix | | | | |
| tFold-Ag | | | | |
| AF3 | | | | |

## 폴더
- `scripts/` : 서버 진입점 (01_prep / 02_run / 03_eval)
- `configs/` : 타깃·모델 목록 (targets.yaml)
- `src/`     : (나중) 재사용 모듈 — 모델별 runner, eval, consensus
- 결과물(예측 구조·MSA 등)은 서버 `~/dock/` 에만, git 에는 안 올림.
