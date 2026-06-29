# bk21-antibody-ml

서울대학교 의과대학 항체·면역학 실험실(이창한 교수님 연구실) **BK21 여름 학부연구 2026** 코드 저장소.
항체 엔지니어링 × 머신러닝.

## 메인 프로젝트

**AlphaFold3 항체-항원 도킹 스크리닝 보강 ML 프레임워크**
(*An Integrated Computational Framework to Augment AlphaFold3 for High-Throughput Virtual
Screening of Antibody-Antigen Interactions*)

AF3의 전역 `ipTM` 신뢰도는 항체-항원 *계면의 국소 기하·물리화학 상보성*을 잘 못 잡아
false-negative/false-positive가 크다. 특히 PDB 템플릿이 적은 신규 표적에서 저-ipTM 영역에
좋은 포즈가 묻힌다. 본 프로젝트는 **PLM 임베딩 + 계면 그래프 기하 피처로 저-ipTM 영역의
native-like 포즈를 rescue**하는 사후(post-prediction) 평가 도구를 만든다.

### 파이프라인

1. **데이터** — PDB 복합체당 15 seed × 5 model = 75 도킹 모델. (성공 기준 DockQ ≥ 0.23)
2. **계면 피처 추출** — AF3 출력 JSON(`full_data`, `summary_confidences`)에서
   contact_probs / PAE / pLDDT 계면 피처를 추출. (`scripts/`)
3. **그래프화** — 항원(A)–항체(B) 계면을 heavy-atom 5.0 Å 컷오프 + KNN(K=12)로 접촉 그래프화.
4. **노드 임베딩** — ESM-2 / ProstT5 / SaProt.
5. **모델** — PCA + penalized LogReg, GroupKFold(PDB 단위) 누수 차단.

> ⚠️ 구체 모델·수치는 초안 기준이며 갱신 중.

## 구조

```
scripts/   AF3 출력 파싱·피처 추출 스크립트
data/      (gitignore — 로컬 전용)
notebooks/ 탐색·분석 노트북
```

## 비고

- 미공개 연구실 데이터·아이디어는 커밋하지 않는다 (`data/`, `*.docx`, `*.pdf` 등 gitignore).
