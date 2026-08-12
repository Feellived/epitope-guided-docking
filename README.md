# epitope-guided-docking

서울대학교 의과대학 항체·면역학 실험실(이창한 교수 연구실) **BK21 여름 학부연구 2026** 코드 저장소.

---

## 무엇을 하는 저장소인가

**항체–항원 복합체에서 결합 자리(epitope)를 제안하고, 그 자리를 제약으로 주어 다시 도킹하는
파이프라인.** 공개 co-folder(Boltz-2 · Chai-1 · Protenix · AF2-Multimer · tFold-Ag)는 단백질
구조 자체는 잘 접지만 항체–항원 결합 자리와 자세에서는 한계를 보이며, 모델이 내놓는 전역
신뢰도(ipTM)는 계면의 옳고 그름을 가르지 못한다. 이 저장소는 그 한계를 정량화하고, 결합 자리를
후보로 생성·선택·주입하는 네 단계를 구현한다.

```
① 생성   항원 MSA 조성을 다시 뽑아 서로 다른 결합 자리 후보를 만든다        (Protenix)
② 조립   조성 간 투표로 후보 결합 자리 몇 개를 세운다
③ 선택   정답을 보지 않고 후보 하나를 고른다                                (AbEpiScore-1.0)
④ 재도킹 고른 자리를 pocket 제약으로 주고 다시 접는다                       (Boltz-2)
   게이트  무제약 예측의 신뢰도 < 0.75 이면서 후보가 3개 이하일 때만 개입한다
```

네 단계 어디에서도 정답 구조를 참조하지 않는다.

> **상세 인수인계 문서는 Notion 「인수인계서 Ⅰ · Consensus Docking」에 있다.**
> 연구 경과·복합체 명세·모델별 상세 실행법·코드 지도·알려진 함정이 모두 그쪽에 정리돼 있다.
> 이 README는 저장소를 열었을 때 길을 잃지 않기 위한 최소 안내다.

---

## 확정된 결과 (복합체 30종)

| 층위 | 성적 |
|---|---|
| **생성** — 정답 결합 자리의 50% 이상을 포함하는 후보가 만들어진 복합체 | 23종 (77%) |
| **선택** — 실제로 고른 후보가 그 조건을 만족한 복합체 | 17종 (57%) |
| **적용** — 개입 대상 16종 중 DockQ가 오른 복합체 | 5종 (31%) |

개입한 16종 중 예측 자세가 선택한 자리로 이동한 것이 9종이고, DockQ까지 개선된 5종이 모두
여기에 속한다. 나머지 4종은 **자리는 옮겼으나 자세를 맞히지 못한** 경우로, 재도킹 모델 자체의
한계에 해당한다.

**대표 사례** — 8q7s_H는 DockQ 0.014 → 0.327, 결합 자리 겹침 0.150 → 0.867로 올랐고 나머지
세 대조 조건은 0.015 이하에 머물렀다. 8sis_HL(0.034 → 0.315)과 8tp5_HL(0.072 → 0.217)은
CAPRI acceptable(0.23)을 새로 넘었고, 9zdu_HL은 0.577 → 0.863으로 high 등급에 도달했다.

**개입 대상 선별이 결정적이다.** 무제약 예측의 신뢰도가 0.75 미만이면서 후보가 3개 이하인
16종으로 한정하면 평균 +0.052이고 최대 하락은 −0.051에 그치지만, **30종 전체에 개입하면
평균이 0.049 하락한다.**

**AlphaFold 3 대비** — 자리 축에서는 우위다. 8q7s_H에서 AlphaFold 3는 5개 결과 모두 정답
결합 자리를 전혀 회복하지 못했으나(겹침 0.000) 이 파이프라인은 0.824였다. 다만 자세 축에서는
같은 자리를 찾고도 AlphaFold 3가 더 정확한 경우가 있어 4종 중 2종에서 뒤처진다.

### 알려진 한계

- **선택 단계가 실질적으로 작동하지 않는다.** 비교한 선택기 다섯 개 어느 것도 무작위 선택과
  유의하게 구별되지 않았다(순열검정 p ≥ 0.074, Bonferroni 문턱 0.0056). AbEpiScore의 평균
  F1은 0.382로, 최적 후보를 골랐을 때의 0.517에 미치지 못한다
- **성능 상한 자체가 낮다.** 생성된 후보는 정답 결합 자리의 약 2.5배 크기여서 정밀도가
  0.35에 그친다. 후보를 좁히는 정제 실험에서 상한이 0.482 → 0.586으로 올랐다
- 적용 규칙의 두 임계값(신뢰도 0.75 · 후보 3개)은 **새로운 복합체 조합에서 검증되지 않았다**

---

## 디렉토리

| 경로 | 내용 |
|---|---|
| `pipeline/` | **현재 연구의 전부.** 아래 표 참조 |
| `pipeline/scripts/` | 코드 100여 개. 명명 규칙은 다음 절 |
| `pipeline/configs/` | 타깃 세트 yaml (`_v1`=초기 6종, `_v2`=다양-항원 Fab 10종) |
| `pipeline/dataset/` | MSA 조성·깊이 실험용 데이터셋 조립(SAbDab → manifest → 라벨링 → 사다리) |
| `pipeline/results/` | 커밋되는 `.txt` 산출물. `.png`·`.csv`·`.log`는 추적하지 않는다 |
| `pipeline/figures/` | 발표·보고서에 쓴 선별 그림. `figures/adhoc/`은 추적하지 않는다 |
| `pipeline/archive/` | 대체됐으나 과거 결과 재현을 위해 남긴 코드 |
| `legacy_af3/` | **현 파이프라인과 무관한 초기 과제 산출물.** AF3 서버 출력 처리 도구와 프로토콜 문서 |

`legacy_af3/`는 2026년 6월 온보딩 과제(AlphaFold 3 서버 출력에서 계면 피처를 뽑는 작업)의
결과물이다. `pipeline/`과 import 관계가 전혀 없으며, 그 안의 문서는 6월 시점 기준이다.

---

## 스크립트 명명 규칙

파일명 앞부분만 보고 파이프라인의 어느 단계인지 알 수 있도록 여덟 개 접두사로 통일했다.

| 접두사 | 의미 | 예 |
|---|---|---|
| `prep_` | 구조·서열·MSA 전처리 | `prep_chains.py` `prep_fv.py` `prep_a3m_sanitize.py` |
| `make_` | 모델 입력 파일 생성 | `make_boltz.py` `make_pocket_from_haddock.py` |
| `run_` | 실행 오케스트레이션·설치·동기화 | `run_diverse.sh` `run_tfold.sh` `run_sync_notion.py` |
| `eval_` | 예측 결과 채점 | `eval_boltz.py` `eval_guided_boltz.py` |
| `select_` | 정답을 보지 않고 후보 선택 | `select_site.py` `select_report_pick.py` |
| `analyze_` | 집계·진단·통계 | `analyze_consensus.py` `analyze_where_docked.py` |
| `plot_` | 그림 생성 | `plot_pose.py` `plot_epitope.py` |
| `lib_` | 공유 라이브러리(CLI 없음) | `lib_structure.py` `lib_log.py` |

`lib_structure.py`가 사실상의 중앙 라이브러리다. 서열 추출·정렬 매핑·native 정보·타깃 상수를
담고 있으며 10여 개 모듈이 `import lib_structure as B`로 쓴다.

---

## 실행 순서

각 단계의 대표 명령만 적는다. 전체 인자는 `--help`로 확인한다.

```bash
# 1) 전처리 — native 내려받아 사슬 서열 추출
python pipeline/scripts/prep_chains.py --targets pipeline/configs/targets_benchmark_v2.yaml --out-dir runs_diverse
python pipeline/scripts/prep_fv.py --run-dir runs_diverse          # Fv(VH+VL) 파생 세트

# 2) 모델 입력 생성
python pipeline/scripts/make_boltz.py --run-dir runs_diverse
python pipeline/scripts/make_chai.py --run-dir runs_diverse

# 3) 예측 (서버, GPU)
bash pipeline/scripts/run_diverse.sh chai protenix
bash pipeline/scripts/run_tfold.sh

# 4) 채점 — summary.csv 가 숫자 기록의 정본
python pipeline/scripts/eval_boltz.py --target 9K6J --run-dir runs_diverse/9k6j --seed 0 \
    --phase Diverse-Fab --csv pipeline/results/summary.csv

# 5) 분석
python pipeline/scripts/analyze_consensus.py --csv pipeline/results/summary.csv --thr 0.23
python pipeline/scripts/plot_iptm_dockq.py --csv pipeline/results/summary.csv
```

---

## conda 환경

모델마다 의존성이 충돌해 환경을 나눠 쓴다.

| 환경 | 용도 | 비고 |
|---|---|---|
| `boltz` | Boltz-2 예측 · **DockQ 채점** | 채점은 반드시 이 환경에서 |
| `protenix` | Protenix 예측 | `LAYERNORM_TYPE=torch` 필요 |
| `chai` | Chai-1 예측 | |
| `colabfold` | AF2-Multimer 예측 · a3m 생성 | `env -u LD_LIBRARY_PATH` 필요 |
| `tfold` | tFold-Ag 예측 | |
| `anarci` | ANARCI(IMGT) 번호매김 | HMMER 필요 |
| `haddock3` | HADDOCK3 물리 도킹 | CNS 번들 |
| `pymol` | 그림 렌더 | `pymol-open-source` |

---

## ⚠️ 이 파이프라인은 멈추지 않고 조용히 빈다

실패해도 종료 코드가 정상이고 결과 파일도 생기는 사고를 네 유형 겪었다.
**채점 결과를 믿기 전에 실행 수가 설계값과 맞는지 먼저 확인한다.**

| 유형 | 증상 |
|---|---|
| a3m 질의행 오염 | Boltz가 MSA를 통째로 버리고 경고만 남긴 뒤 정상 종료 |
| 깊이·조건 폴더 혼입 | 한 폴더에 두 조건이 섞여 자세 수가 배로 늘고 최고값이 부풀려짐 |
| DockQ 환경 밖 실행 | `dockq` 열만 통째로 비는데 recall은 정상이라 결과가 있어 보인다 |
| 비교군 설정 불일치 | 팔마다 자세 예산·제약이 달라 비교가 성립하지 않음 |

자가 점검 도구는 `pipeline/scripts/analyze_verify_demo.py`(7개 항목, 실패 시 exit 1)다.
자세한 내용과 사전 게이트는 `pipeline/README.md`와 Notion 인수인계서 Ⅰ에 있다.

---

## 외부 의존

파이프라인의 **① 생성**과 **② 조립** 단계는 형제 저장소
**[`msa-composition-bias`](https://github.com/Feellived/msa-composition-bias)** 에 있다.
이 저장소의 데모·honest·allcand 흐름은 그 저장소를 환경변수 `MSAD`(기본값
`~/projects/msa-composition-bias/pipeline`)로 참조한다.

- 공유 모듈 `lib_pose_features`·`lib_epitope_recall`을 `$MSAD/scripts`에서 불러온다
- 결합 자리 후보 `$MSAD/results/sites_<타깃>.json`과 사슬 정보 `$MSAD/targets/<타깃>/chains.json`을 읽는다

그 저장소가 없으면 해당 흐름은 동작하지 않는다. **한쪽의 파일명을 바꾸면 반대쪽도 함께
고쳐야 한다.**

---

## 비고

- 예측 구조·MSA·weight·대용량 CSV·PyMOL 세션은 저장소에 올리지 않는다(`.gitignore`).
- 서버 접속 정보(IP·포트·계정·비밀번호)와 API 토큰은 문서·코드·커밋에 남기지 않는다.
- 코드는 붙여넣기가 아니라 커밋 후 서버에서 `git pull` 로 전달한다.
