# Consensus Docking — 항체-항원 pose 재랭킹 / 오케스트레이션

여러 공개 co-folder(Boltz-2 · Chai-1 · Protenix · AF2-Multimer · tFold-Ag)를 한 파이프라인으로
돌려, 단일 ipTM 대신 **모델 간 합의(consensus)로 native-like pose를 고른다.**
> 이 폴더 문서는 이 README **하나로 통합**했다(옛 RUN_ALL/AUTOMATION/benchmark 내용은 여기로).

## 성공 기준 (2026-07-14 갱신)
- **DockQ ≥ 0.49 (CAPRI Medium)** 주기준 — 3-tier 병기(0.23 Acceptable / 0.49 Medium / 0.80 High).
- **RMSD 이원화**: 항체 내부 CA-RMSD(fold 품질) / Ligand-RMSD(배치 품질) + epitope recall.
- single seed는 노이지 → multi-seed 분포 + **ipTM-vs-oracle**(Win Rate / Hit-Oracle / Mean Regret).

## 벤치마크 세트 (configs/)
- **`targets_benchmark_v2.yaml` — 현재 메인**: 2024+ 공개 Fab 10종(다양 항원, leakage-free).
- **`targets_benchmark_v1.yaml`**: 9VP0(골든/dev) + 2026 공개 5종(10ZO·9Y0A·9OBN·9K6J·9XF9; 나노바디 2 포함).
- 교체 후보(필요시): 9PIU(Fab+mTOR FRB) · 9UOI(나노바디+MERS ternary) · 9SAT/9MNT · 큰 항원 9VDY/9TH6.

## leakage / 사슬 역할 규칙 (사수 보고용 핵심)
- leakage 판정 = RCSB `initial_release_date`. **SAbDab 날짜는 deposit 미러라 부정확 → 쓰지 말 것.**
- 사슬 역할(H/L/항원)·형식(Fab/나노바디) = SAbDab 표기. 다중카피(10ZO·9OBN 등)는 1:1 한 쌍만 사용.

## 폴더 구조
- `README.md` — 이 문서(단일 소스)
- `configs/`  — 타깃 세트 yaml (`targets_benchmark_v2.yaml`=메인 · `targets_benchmark_v1.yaml`=초기)
- `scripts/`  — 코드(아래 파이프라인)
- `results/`  — 커밋되는 `.txt` 산출(diagnose_* · where_docked_* · merged_*). png·csv·log는 gitignore
- `archive/`  — 레거시. 9VP0 브링업(01/02/03) · 옛 config · `reorg.sh`(2026-07-14 실행 완료) ·
  모델별 pocket 변환기 3종. 변환기는 `make_pocket_from_haddock.py`로 통합됐으나
  과거 결과 재현 경로를 위해 남겨 둔다
- `figures/`  — 발표·보고서에 쓴 선별 그림(`cases_*/`). `figures/adhoc/`은 추적하지 않는다
- 아티팩트(`runs*/`·`haddock*/`·`zdock/`·`tfold/`·`DeepRank-GNN-esm/`·`*.pse`·`*.log`) = 전부 gitignore

## 파이프라인
1. `prep_chains.py` — native 다운 + keep_chains 서열 → `runs_diverse/<id>/chains.json`
2. `prep_fv.py` — 항체 ANARCI(IMGT)로 Fv(VH+VL) 절단 → `runs_diverse/<id>_fv/`
3. `make_{boltz,chai,protenix,colabfold,tfold}.py` — 모델별 입력 생성
4. 예측: `run_diverse.sh`(chai·protenix·colabfold) · `run_tfold.sh` · `run_5model.sh`
   → 각 `runs_diverse/<id>/{out_boltz,out_chai,out_protenix,out_colabfold,out_tfold}`
5. eval: `eval_*.py` → `results/summary.csv`(source of truth) → `run_sync_notion.py`(개인 Mac, 토큰 env)
6. 분석: `analyze_consensus.py`·`plot_iptm_dockq.py`(집계) /
   `analyze_antibody_side.py`·`analyze_where_docked.py`·`analyze_merged_dockq.py`(정밀 진단)
7. 렌더: `plot_pose.py`·`plot_epitope.py` + `run_plot_cases_*.sh` → `figures/`

## Phase 0 (go/no-go 게이트)
`run_phase0.sh` = diagnose + where_docked + merge_dockq 를 10-타깃 × Fab/Fv 전부에.
새 예측 없이 기존 pose 재분석 → 0.49 기준 oracle·placement·ipTM 품질로 프로젝트 go/no-go 판정.

## 규칙
- 출력 폴더명(`runs_diverse`·`out_*`)은 스크립트에 하드코딩(reader는 `*boltz`/`*colabfold` 호환).
- 예측 구조·MSA·weight·CSV·`.pse`·클론 툴(tfold·DeepRank)은 git에 안 올림.
- Notion 토큰은 서버·깃·문서 저장 금지(개인 Mac env만).

## ⚠️ 이 파이프라인은 멈추지 않고 조용히 빈다

채점 결과를 믿기 전에 **실행 수가 설계값과 맞는지 먼저 확인한다.** 아래는 실제로 겪은 사고다.

| 유형 | 증상 | 사전 게이트 |
|---|---|---|
| a3m 질의행 오염 | Boltz가 MSA를 통째로 버리고 경고만 남긴 뒤 정상 종료 | 실행 전 a3m 행 길이 검사, 실행 후 `run.log` 경고 grep |
| 깊이 폴더 혼입 | 한 폴더에 두 조건이 섞여 자세 수가 배로 늘고 최고값이 부풀려짐 | 채점 전 pose 개수 확인 |
| DockQ 환경 밖 실행 | `dockq` 열만 통째로 비는데 recall은 정상이라 겉보기엔 결과가 있다 | `which DockQ` 확인 후 채점 |
| 팔(arm) 설정 불일치 | 비교군끼리 자세 예산·제약이 달라 비교가 성립하지 않음 | `analyze_verify_demo.py` (7개 항목 자가 점검, 실패 시 exit 1) |

## 외부 의존

데모·honest·allcand 흐름(`run_demo_guided.sh`·`run_honest_guided.sh`·`run_all_candidates.sh`·
`make_pocket_from_sites.py`·`select_score_discotope.py`)은 형제 레포
`bk21-msa-depth-bias/pipeline`을 하드코딩 경로로 참조한다. 그 레포가 없으면 동작하지 않는다.
