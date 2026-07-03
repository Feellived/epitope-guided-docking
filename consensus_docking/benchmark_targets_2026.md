# Consensus Docking 벤치마크 세트 — 2026 공개 Ab-Ag 복합체 (leakage-free)

제4회 인공지능 신약개발 해커톤 · 항체-항원 consensus docking.
목적: 같은 Ab-Ag 복합체에 여러 구조예측 모델을 돌려 (신뢰도, DockQ)를 비교.

## 확정 세트 (5개, 전부 RCSB 2026 공개 → 모든 모델에 clean)

| # | PDB | 공개일(RCSB) | 형식 | 항체 사슬 | 항원 (길이, 사슬) | 해상도 | 쓸 사슬(1:1) |
|---|-----|-------------|------|-----------|-------------------|--------|--------------|
| 1 | 10ZO | 2026-03-25 | 나노바디(SD-H) | H=A | Mengla 바이러스 nucleoprotein (~106, B) | 1.65 A | A, B |
| 2 | 9Y0A | 2026-02-11 | Fab (IgE) | H=A, L=B | Bet v 1 자작나무 알레르겐 (~160, C) | 2.12 A | A, B, C |
| 3 | 9OBN | 2026-01-21 | 나노바디(SD-H) | H=E | Pfs48/45 말라리아 항원 (~150, A) | 2.5 A | A, E |
| 4 | 9K6J | 2026-01-07 | Fab | H=A, L=B | SARS-CoV-2 RBD (~271, C) | 2.39 A | A, B, C |
| 5 | 9XF9 | 2026-05-06 | Fab | H=H, L=L | Kasokero 바이러스 엔도뉴클레이스 (~300, A) | 1.85 A | A, H, L |

구성: 3 Fab + 2 나노바디 · 항원 다양(바이러스 NP/알레르겐/기생충 항원/RBD/효소) · 공개일 1~5월 분산 · 전부 X-ray <= 2.5 A.

주의(다중 카피): 10ZO(2카피)·9OBN(4카피)는 결정구조에 복합체가 여러 벌 → 위 "쓸 사슬"의 1:1 한 쌍만 사용.

## 골든/dev 케이스 (벤치마크 아님)
- 9VP0 (Fab + IL-31, X-ray 2.1 A, 공개 2025-10). 사슬 A=항원(IL-31) / H=중쇄 / L=경쇄.
- 파이프라인 개발·디버깅용. 2025-10 공개라 AF3·Boltz-2엔 clean이나 Protenix-v2(2026)엔 노출 가능 → rate 산출은 위 5개(2026)로.

## 방법론 (사수 보고용 핵심)
- leakage 판정 = RCSB `initial_release_date` (2026). SAbDab 날짜 필드는 deposit 미러라 부정확 → 쓰지 말 것.
- 5개 다 2026 공개 → AF3(2021)·Boltz-2(2025중반)·Protenix-v2(2026) 전부 학습 컷오프 이후 = 모델 무관 clean.
- 사슬 역할(H/L/항원, 항원 type=protein)은 SAbDab 표기 기준. 5개 모두 SAbDab 등재.
- 성공 = DockQ >= 0.23. single seed는 노이지 → multi-seed 분포 + ipTM-vs-oracle(Win Rate / Hit-Oracle / Regret).

## 돌릴 모델 (co-folder)
AF3(기준선) · Boltz-2 · Chai-1 · Protenix-v2(최신) · (옵션 tFold-Ag, AF2-Multimer).

## 백업 (검증됨, 교체용)
- 9PIU (Fab+mTOR FRB 92aa, 2026-06, 최소·최신) — 단 SAbDab 미색인 + synthetic Fab
- 9UOI (나노바디 2개+MERS RBD, ternary) / 9SAT, 9MNT (또 다른 CoV/Fab, 약간 중복)
- 큰 항원(참고): 9VDY(338aa), 9TH6(547aa)

## 재현 (RCSB Search API)
POST https://search.rcsb.org/rcsbsearch/v2/query — 필터: `rcsb_accession_info.initial_release_date >= 2026-01-01` + `rcsb_entry_info.resolution_combined <= 2.8` + full_text "nanobody"(또는 "Fab complex"). SAbDab 사슬역할 = GET https://sabdab.opig.stats.ox.ac.uk/api/pdb/pdb_0000<소문자4자리> (`sabdab_chain_type`: H/L/K/N, `type`: SD-H=나노바디/FAB).
