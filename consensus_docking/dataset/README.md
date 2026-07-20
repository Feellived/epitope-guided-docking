# MSA-depth 편향 실험 — 데이터셋 조립 파이프라인

과대표집 항원(RBD·HA·HIV Env) 안에서 **A(지배 에피토프) vs B(비지배)** 항체-항원 복합체를
leakage-clean(post-2023-06)으로 모아, MSA-depth sweep 편향 실험의 데이터셋을 만든다.

설계 근거: 적대검증이 지목한 4구멍(①MSA=fold채널 혼동 ②dominant 정의 순환 ③난이도 confound
④C그룹 교란) 대응 = **같은 과대표집 항원 내 A vs B**(항원·MSA·fold를 상수로 고정) +
**PDB-독립 dominance 라벨**(RBD=ACE2/Barnes class, HA=head/stem, Env=bnAb supersite).

## 파이프라인

```
sabdab2_summary.csv (SAbDab2 전체, date=진짜 구조일)
  └─ count_overrep.py      # 항원별 과대표집 카운트(진단용)
  └─ build_manifest.py     # RBD/HA/Env × post-2023-06 × res≤3.5 → manifest_candidates.csv (375 PDB)
       └─ fetch_structures.py    # RCSB에서 mmCIF 다운로드 → structures/
            └─ classify_epitope.py  # 접촉잔기→에피토프클래스→A/B → manifest_labeled.csv
                 (epitope_defs.py = 잔기정의 + A/B 매핑)
```

## 실행 (서버, boltz 또는 DockQ env = biopython 필요)

```bash
cd ~/projects/bk21-antibody-ml/consensus_docking/dataset
conda activate boltz            # biopython 있는 env
# summary 스냅샷이 없으면 최신 다운로드(라이브 DB):
#   curl -sL https://sabdab.opig.stats.ox.ac.uk/api/download/all-summary -o sabdab2_summary.csv
python fetch_structures.py      # manifest_candidates.csv의 375 PDB 다운로드(~1GB)
python classify_epitope.py      # → manifest_labeled.csv (A/B 라벨 + 요약 출력)
```

## 분류 규칙 (epitope_defs.py)

- **RBD** (Wuhan auth 넘버링 직접): 접촉잔기의 RBM(437-508+417) overlap ≥0.5 → **A**(class1/2);
  class4(369-386,404-408) or class3(N343) 우세 → **B**. RBD 범위(319-541) 밖 접촉만 → `no-RBD-contact`(NTD/S2 결합=배제).
- **HA** (P03437 참조 정렬): 접촉잔기의 head globular 범위(P03437 56-306) 비율 ≥0.6 → **A**(head);
  ≤0.35 → **B**(stem); 사이 → `mixed`(검토).
- **Env** (P04578=HXB2 정렬): bnAb class footprint overlap 최대 클래스 → CD4bs/N332/V2apex = **A**(구조-빈발 supersite),
  나머지 = **B**. `ENV_DOMINANT_MODE`로 dominant 정의 토글(structural=기본 / infection).

## 검증 (2026-07-20, Mac biopython 1.85, 알려진-에피토프 참조 6개)

| 참조 | 항체 | 기대 | 결과 |
|---|---|---|---|
| 6W41 | CR3022 | class4/B | class4/B ✓ |
| 3GBN | CR6261 | stem/B | stem/B ✓ |
| 4FQI | CR9114 | stem/B | stem/B ✓ |
| 2VIR | HC19 | head/A | head/A ✓ |
| 4FP8 | C05 | head/A | head/A ✓ |
| 3NGB | VRC01 | CD4bs/A | CD4bs/A ✓ |

**6/6 정답.** 이후 = A/B 카운트 확인 → 검증 게이트(fold-stability, 난이도 매칭, Neff80 사다리) → lock → depth-sweep 생성.

## 주의 / 한계

- `date` = SAbDab deposition. 엄밀 leakage는 RCSB **release date**로 최종 재확인 권장(deposition≤release라 다소 보수적).
- Env dominance는 미묘(PDB-빈발 bnAb site vs 감염-면역우세 가변루프) → `ENV_DOMINANT_MODE=structural`(모델이 과학습한 것) 기본, 한계 명시.
- HA `mixed`·RBD `ambiguous`·Env `low-classmatch` 플래그 행은 수동 spot-check.
- 다중 Fv/copy 구조는 사슬별로 분리 라벨(한 PDB에 A·B 항체 공존 가능).
