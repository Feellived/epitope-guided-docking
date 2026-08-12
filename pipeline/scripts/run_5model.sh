#!/usr/bin/env bash
# AF2-Multimer(ColabFold) + tFold-Ag 통합 무인 실행 (Fab+Fv, GPU 1개 순차).
# 이미 있는 Boltz·Chai·Protenix에 더해 5모델 co-folder를 완성.
#
# 사용:
#   tmux new -s five
#   ONLY=9k6j bash scripts/run_5model.sh     # 스모크(AF2-M 9k6j는 이미 됨→skip, tFold 9k6j만 검증)
#   bash scripts/run_5model.sh               # 전체 (AF2-M 나머지 + tFold 전부)
#   TFOLD_ENV=<이름> CF_ENV=<이름> bash ...   # env 이름 다르면 override (자식 스크립트로 전파됨)
#
# 각 단계 resumable(출력 있으면 스킵). 로그 = run_diverse.log · run_tfold.log · 각 dir.
[ -n "${ONLY:-}" ] && export ONLY
[ -n "${CF_ENV:-}" ] && export CF_ENV
[ -n "${TFOLD_ENV:-}" ] && export TFOLD_ENV
cd "$(dirname "$0")" || exit 1

echo "########## [1/2] AF2-Multimer (ColabFold) ##########"
bash run_diverse.sh colabfold

echo "########## [2/2] tFold-Ag ##########"
bash run_tfold.sh

echo "########## 5모델 co-folder 완료: Boltz·Chai·Protenix + AF2-M·tFold ##########"
echo "다음: eval_colabfold.py · eval_tfold.py 로 채점 → summary.csv"
