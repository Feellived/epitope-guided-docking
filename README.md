# bk21-antibody-ml

서울대학교 의과대학 항체·면역학 실험실(이창한 교수님 연구실) **BK21 여름 학부연구 2026** 코드 저장소.
항체 엔지니어링 × 머신러닝.

## 메인 프로젝트

**AlphaFold3 항체-항원 도킹 스크리닝 보강 ML 프레임워크**

AF3의 전역 `ipTM` 신뢰도는 항체-항원 *계면의 국소 기하·물리화학 상보성*을 잘 못 잡아
false-negative/false-positive가 크다. 본 프로젝트는 **PLM 임베딩 + 계면 그래프 기하 피처로
저-ipTM 영역의 native-like 포즈를 rescue**하는 사후(post-prediction) 평가 도구를 만든다.

> 자세한 파이프라인은 초안 문서(`Idea/Epitope_binning_6page.docx`, 이 레포 밖) 참조.

---

## 문서 (docs/) — 사수 지시 업무 프로토콜

직접 손으로 한 번씩 실행해보는 게 목적이라, **번호 순서대로** 따라가면 된다.

| 순서 | 문서 | 내용 |
|---|---|---|
| 00 | [VSCode 파이썬 워크플로](docs/00_VSCode_파이썬_워크플로.md) | 코랩→.py 전환, venv, 통합 터미널, 디버거, Remote-SSH |
| 01 | [AF3 서버 프로토콜](docs/01_AF3_서버_프로토콜.md) | alphafoldserver.com: 입력 순서·PTM·MSA·템플릿·seed 고정 |
| 02 | [AF3 출력 처리](docs/02_AF3_출력_처리.md) | cif/json/pdb, full_data 키 의미, 폴더 정리, 피처→CSV |
| 03 | [서버 ssh·scp·tmux](docs/03_서버_ssh_scp_tmux.md) | 서버 접속·파일전송·세션 유지·Remote-SSH |
| 04 | [RFD3 프로토콜](docs/04_RFD3_프로토콜.md) | RFDiffusion3 binder 설계 흐름 |
| 05 | [Alanine scanning](docs/05_alanine_scanning.md) | CDR 잔기 Ala 치환 hotspot 분석 전처리 |

> ⚠️ 문서 속 외부 도구의 UI/명령(AF3 서버, RFD3, VSCode 등)은 버전에 따라 달라질 수 있다.
> 각 문서 끝의 주의사항을 보고, 실제 화면/`--help`로 한 번씩 대조할 것.

## 스크립트 (scripts/)

모두 `python scripts/<파일> --help`로 옵션 확인 가능. (로컬 검증: `--help`·합성데이터 통과,
**실제 AF3 데이터로는 미검증** — 첫 케이스에서 결과를 눈으로 확인할 것.)

| 스크립트 | 역할 |
|---|---|
| [`af3_extract_features.py`](scripts/af3_extract_features.py) | full_data/summary 재귀 탐색 → 계면 피처(contact·PAE 방향별·pLDDT) CSV. **다중체인 항체 지원** |
| [`af3_convert.py`](scripts/af3_convert.py) | cif → pdb 변환(gemmi, biopython fallback), json 요약 추출 |
| [`clean_af3_outputs.py`](scripts/clean_af3_outputs.py) | msas/template 보조 폴더 재귀 삭제 (기본 DRY-RUN, `--apply`로 실제 삭제) |
| [`cif_distance.py`](scripts/cif_distance.py) | 구조에서 잔기 간 거리·계면 접촉 잔기 추출 |
| [`alanine_scan_prep.py`](scripts/alanine_scan_prep.py) | CDR/계면 잔기 Ala 치환 변이체 + AF3 입력 명세 생성 |

## 시작하기

```bash
# 1) 이 폴더에서 가상환경 만들고 켜기
python3 -m venv .venv
source .venv/bin/activate

# 2) 의존성 설치
pip install -r requirements.txt

# 3) 예: 스크립트 옵션 확인
python scripts/af3_extract_features.py --help
```

자세한 환경 셋업(VSCode 포함)은 [docs/00](docs/00_VSCode_파이썬_워크플로.md) 참조.

## 비고

- 미공개 연구실 데이터·문서·구조파일(`data/`, `*.docx`, `*.pdf`, `*.csv`, `*.cif`, `*.pdb` 등)은
  `.gitignore`로 커밋에서 제외. 의도적으로 올리려면 `git add -f`.
- 서버 접속 정보(IP/포트/ID/PW)는 문서에 넣지 않고 placeholder로만 둔다.
