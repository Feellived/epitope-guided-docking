# 진행 체크리스트 (BK21 항체 ML)

최종 갱신: 2026-06-29

## ✅ 한 것
- [x] GitHub 레포 + VSCode + venv + 통합터미널 워크플로 전환
- [x] PyMOL 설치 (brew), DockQ 설치 (pip)
- [x] docs/ 프로토콜 6종 (00 VSCode / 01 AF3서버 / 02 AF3출력 / 03 ssh·scp·tmux / 04 RFD3 / 05 alanine)
- [x] scripts/ 도구 7종 (af3_extract_features, clean_af3_outputs, af3_convert, cif_distance, alanine_scan_prep, fetch_uniprot_domain)
- [x] **CD16a + anti-CD16a Fab(7SEG) AF3 1사이클 완주:**
  서열수집(UniProt ECD 17–208 + PDB Fab, 신호펩타이드 제거) → 제출(223→224 복붙 교정)
  → clean → extract → CSV → 계면 pLDDT 재추출 → PyMOL 시각비교(정답 7SEG) → DockQ 채점
- [x] **전체 파이프라인 end-to-end 검증 완료**

## 📊 CD16a 케이스 결과 (결과 로그)
- AF3 **항원 fold은 정확**: 정답(7SEG) 대비 항원 RMSD **0.73 Å**.
- 그러나 **도킹은 5모델 모두 실패**: 항원-항체 계면 **DockQ ≈ 0.02–0.03** (≪ 0.23 성공기준),
  fnat 0.000, LRMSD ~33 Å. (heavy-light 내부 계면은 0.74 — trivial, 무시.)
- ipTM ~0.53(낮음)·PAE 15–23Å(높음) = AF3 자가 신호가 "불안"이라 했고, DockQ가 "정말 실패"로 확인.
  → **자가신호↔실제실패 일치 = 프로젝트 가설(저신뢰 신호로 나쁜 포즈 가려내기)의 증거.**
- 결론: 이 **단일 seed로는 native 못 찾음** → 확정된 "어려운/실패" 예시 1개 확보.
- 교훈: 멀티체인 항체 DockQ는 **항원 낀 계면만** 봐야 함(Total은 Fab 내부 때문에 부풀려짐).

## 📊 라이소자임(3HFM, HyHEL-10 Fab) 케이스 결과 (성공 대조)
- **AF3가 도킹 정확히 맞춤**: 항원-항체 DockQ **~0.83–0.89 (high quality)**, 5모델 전부 일관 성공.
- 항원 RMSD 0.56Å + Fab도 정답과 겹침(PyMOL 확인).
- → **CD16a(실패 ~0.02) vs 라이소자임(성공 ~0.85) = 좋은/나쁜 포즈 대조쌍 확보.** 프로젝트 ML의 출발점.
- 이유: 라이소자임 복합체는 PDB에 풍부 + MSA 깊음 = AF3가 쉽게 맞히는 표적.

## ⬜ 다음
- [x] **C-2: PTM(글리칸) v2 완료** — CD16a N-당화 5자리(NAG, ECD 40·47·76·164·171)에 추가.
  결과: 글리칸 넣어도 **도킹 안 바뀜** — 항원-항체 DockQ ~0.02–0.04 (baseline과 동일, 모두 ≪0.23).
  Fab 여전히 반대 면. → 이 epitope은 글리칸 비의존(이 표적/seed 기준). 첫 비교실험(baseline vs +glycan) 완료.
- [ ] seed 늘리기 (서버 15 job vs 로컬 AF3 배치 — 사수 확인)
- [x] 대조용 성공 복합체 완료 — 3HFM 라이소자임 (항원 DockQ ~0.83–0.89, high). CD16a와 대조쌍 확보.
- [ ] 미실행 도구: cif_distance.py / alanine_scan_prep.py (af3_convert는 이번에 사용함 ✓)
- [ ] 기본기: 연구실 서버 접속(ssh/scp)+tmux (docs/03), RFD3 실습 (docs/04)
- [ ] 더 큰 미완: #2 Consensus Docking 본체 / #4 랩 DB / #5 wet / #6 내년 논문

## ❓ 사수에게 확인 (구체화됨)
- [ ] **15 seed × 5 절차 = 서버(15번 제출) vs 로컬 AF3 배치?** (이번에 병목으로 확인됨)
- [ ] 진짜 데이터셋(초안 69 PDB / 5,175 model) 출처·접근법
- [ ] 로컬 AF3 설치 여부 / 서버 ssh 포트 (인증=ID/PW 확인됨)
- [ ] "pseudo" 의미 / Fab vs scFv / CD16a 알로타입(158 V/F) / DockQ용 정답구조 기준
