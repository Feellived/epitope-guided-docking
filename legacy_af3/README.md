# legacy_af3 — 초기 과제 산출물 (2026-06)

**현재 파이프라인과 무관하다.** `pipeline/`과 import 관계가 전혀 없으며, 유지보수도 하지 않는다.

2026년 6월 온보딩 과제로 만든 것들이다. 당시 계획은 AlphaFold 3의 전역 신뢰도(ipTM)가
계면의 국소 기하를 반영하지 못하는 문제를 사후 재랭커로 보정하는 것이었다. 이 계획은
이후 폐기됐고, 연구는 결합 자리 제안·선택·유도 재도킹 방향으로 옮겨 갔다(`pipeline/` 참조).

남겨 두는 이유는 두 가지다. AlphaFold 3 서버 출력을 다룰 일이 다시 생길 수 있고,
`docs/`의 프로토콜 문서가 서버 사용법·Unix 작업 흐름 자료로 여전히 쓸 만하기 때문이다.

## scripts/

| 스크립트 | 역할 |
|---|---|
| `af3_extract_features.py` | AF3 출력 폴더를 재귀 탐색해 계면 피처(접촉 확률·방향별 PAE·pLDDT)를 CSV로 |
| `af3_convert.py` | mmCIF → PDB 변환(gemmi 우선, biopython 대체), `--summary`로 신뢰도 스칼라 추출 |
| `clean_af3_outputs.py` | `msas`·`template` 보조 폴더 재귀 삭제. **기본 dry-run**, `--apply`로 실제 삭제 |
| `cif_distance.py` | 두 사슬 그룹 간 잔기 거리 → 계면 접촉 잔기 |
| `alanine_scan_prep.py` | CDR 잔기 Ala 치환 변이체 + AF3 작업 명세 생성 |
| `fetch_uniprot_domain.py` | UniProt 서열 내려받기 / 도메인 구간 잘라 FASTA |

## docs/

번호 순서대로 읽도록 쓴 실습 문서다(00 VSCode 작업 흐름 · 01 AF3 서버 · 02 AF3 출력 처리 ·
03 ssh/scp/tmux · 04 RFdiffusion3 · 05 alanine scanning).

⚠️ **본문은 2026년 6월 시점 기준이며 이후 갱신하지 않았다.** 특히 다음에 주의한다.

- 02·05 문서는 `af3_convert.py`·`af3_extract_features.py`·`cif_distance.py`를 "예정 스크립트"로
  적고 가칭 인자를 제시한다. 세 스크립트는 모두 작성됐고 실제 인자는 다르다. `--help`로 확인할 것.
- 04(RFdiffusion3)는 전부 예시이며 이 저장소에서 실제로 실행한 적이 없다.
- `PROGRESS.md`는 2026-06-29에 멈춰 있다. 그 이후의 모든 작업은 `pipeline/`에 있다.
