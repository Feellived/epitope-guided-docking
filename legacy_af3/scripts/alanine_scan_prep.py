#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""AF3용 in-silico 알라닌 스캐닝(alanine scanning) 전처리 CLI.

이 스크립트가 하는 일 (한 줄)
-----------------------------
항체 서열에서 "대상 잔기"를 하나씩 알라닌(A)으로 바꾼 **단일치환 변이체**들을 만들고,
각 변이체를 AF3로 다시 예측하기 위한 **입력 파일**(FASTA + 간단한 job JSON)과,
"어떤 입력이 어떤 변이인지" 적은 **매핑표(manifest CSV)**를 한꺼번에 찍어낸다.
야생형(WT, 원본)도 비교 기준으로 함께 만든다.

배경 (왜 필요한가 — 완전 초보용 풀이)
--------------------------------------
- **알라닌 스캐닝:** "이 잔기가 결합에 얼마나 중요한가"를 알아보려고, 그 잔기를
  곁사슬이 거의 없는 알라닌(Ala, 한 글자 코드 ``A``)으로 바꿔보는 고전적 방법.
  바꾼 뒤 결합이 많이 약해지면 그 잔기가 핵심(hotspot)이라고 본다.
- **in silico(전산):** 실제 실험 대신 컴퓨터(여기선 AF3 재예측)로 변이체를 다룬다는 뜻.
- **잔기(residue):** 단백질을 이루는 아미노산 한 개 단위. "100번 잔기" = 100번째 아미노산.
- 변이체가 많아 손으로 서열을 고치다 보면 실수가 잦다. 이 스크립트는 그 **반복·실수
  나기 쉬운 전처리만** 자동화한다. 실제 AF3 제출은 사람이 직접 한다(아래 "범위 밖" 참조).

이 스크립트의 입력
------------------
1) 항체 서열 — 둘 중 하나로 준다:
   - ``--fasta`` : 항체 FASTA 파일. 레코드(>헤더)가 곧 항체의 한 체인이다.
                   여러 레코드면 여러 체인(예: Fab = 중쇄절편 + 경쇄). ``--ab-chains`` 로
                   처리할 체인을 고르거나, 안 주면 모든 레코드를 항체 체인으로 본다.
   - ``--structure`` : .cif/.pdb 구조 파일에서 ``--ab-chains`` 로 지정한 체인의 서열을 뽑는다.
                       (구조 입력일 때만 gemmi 또는 biopython이 필요 — 아래 "의존성" 참조.)
2) 항원 서열 — ``--ag-fasta`` (AF3 job JSON에 함께 넣기 위함). 안 주면 항원 없이
   "항체만" 입력을 만든다(보통은 항원을 함께 넣어 복합체로 예측하므로 지정 권장).
3) 치환 대상 잔기 — 둘 중 하나로 준다:
   - ``--positions`` : "체인:번호" 또는 "번호" 목록. 공백 또는 쉼표로 구분.
                       예) ``B:31 B:100 52`` 또는 ``B:31,B:100,52``.
                       체인을 안 적으면 항체 체인이 하나일 때만 그 체인으로 간주한다.
   - ``--contacts`` : ``cif_distance.py`` 가 만든 계면 접촉 CSV. 이 CSV에서 항체 쪽
                      접촉 잔기를 자동으로 골라 대상으로 삼는다(두 출력 형식 모두 인식).

⚠️ 번호(좌표계) 주의 — 가장 흔하고 위험한 실수
-----------------------------------------------
``--positions`` / ``--contacts`` 의 잔기 번호가 **서열 위치인지 author 잔기번호인지**,
또 **1-기반인지 0-기반인지** 반드시 통일해야 한다. 이 스크립트는 기본적으로
**1-기반 서열 위치**(서열의 N번째 아미노산, 1부터)로 해석한다(``--index-base`` 로 0-기반 전환).

⚠️ **``cif_distance.py`` 의 출력은 author 잔기번호(구조 파일에 적힌 번호)를 쓴다.**
구조의 author 번호와 "서열 위치"는 결손 잔기·삽입코드·번호 오프셋 때문에 **다를 수 있다.**
이 스크립트는 ``--contacts`` CSV의 잔기 라벨에서 **원래 아미노산 이름(예: TYR=Y)** 까지 읽어,
"내가 서열 N번째에서 본 아미노산"과 "CSV가 말하는 아미노산"이 **일치하는지 검사**한다.
   - 일치 → 그 잔기는 안전하게 진행(번호가 맞다는 강한 신호).
   - 불일치 → **건너뛰고(SKIPPED) 시끄럽게 경고**한다(좌표계 어긋남 가능 — 사람이 확인 필요).
``--no-resname-check`` 로 이 안전장치를 끌 수 있으나, 끄면 엉뚱한 잔기를 변이시킬 위험이 커진다.

출력 (``--out`` 폴더 아래)
--------------------------
- ``fasta/`` : 변이체별(그리고 WT) FASTA. 한 파일에 [항원 레코드] + [항체 체인 레코드들].
- ``json/``  : 변이체별(그리고 WT) AF3 job JSON(간단 명세). seed 세트 포함.
- ``manifest.csv`` : 매핑표. 한 행 = 한 변이(또는 WT). 컬럼:
    mutation_id, kind(WT/MUT/SKIPPED), chain, position(1-기반), wt_aa, mut_aa,
    note, fasta_path, json_path
  변이ID 예) ``B_Y100A`` = 체인 B, 100번 타이로신(Y) → 알라닌(A). WT는 ``WT``.
- ``README_notes.txt`` : 사람이 읽는 메모(seed/model 관례, AF3 JSON 형식 가정 등).
  ⚠️ AF3가 먹는 JSON에는 비표준 키를 안 넣는다(서버가 거부할 수 있어서). 메모는 여기 모은다.

⚠️ 범위 밖(이 스크립트가 안 하는 일)
------------------------------------
- **실제 AF3 제출은 수동이다.** 이 스크립트는 입력 파일까지만 만든다. 변이체들을 AF3
  웹서버/실행기에 올려 75모델씩 돌리는 일은 사람이 직접 한다(서버 정책·동시 제출 한도·
  UI 변동 때문에 자동 제출은 깨지기 쉬워 일부러 손에 남긴다).
- **점수 추출·WT 비교는 별도 단계**다(``af3_extract_features.py`` 예정). 진짜 ΔΔG(물리
  역장 기반 결합에너지)도 계산하지 않는다 — 우리는 AF3 점수 기반 "대리지표"만 본다.

AF3 입력 포맷 가정 (버전 따라 다를 수 있음)
-------------------------------------------
- AF3 서버/실행기 입력 포맷은 버전에 따라 다르다. 그래서 **범용 FASTA**와 **간단한 job
  JSON** 양쪽을 만든다.
- 여기서 만드는 JSON은 AlphaFold3 공개 입력 스키마를 따른 형태다(2026-01 기준 공개 문서):
  최상위 ``name``, ``modelSeeds``(정수 seed 배열, 최소 1개), ``sequences``(각 원소가
  ``{"protein": {"id": <체인기호>, "sequence": <서열>}}``), ``dialect``("alphafold3"),
  ``version``(정수 1~4, 문서상 최신 4). 단백질 ``id`` 는 문자열 또는 문자열 배열이 허용된다.
- AF3는 보통 **전체 IgG 포맷을 잘 못 받아** 항체를 Fab/scFv로 넣으므로, 항체 체인이 여러
  개면 sequences에 각각 별도 protein 레코드로 들어간다.
- 체인 순서 관례: **위=항원, 아래=항체** → 결과에서 chain A=항원, chain B=항체. 그래서
  JSON에서도 항원을 먼저(보통 ``A``), 항체 체인을 그다음(``B``, ``C`` ...)에 둔다.
- ⚠️ **실제 서버에 올리기 전에 그 버전의 입력 스펙과 대조해 키 이름·구조·version 값을 맞춰라.**
  특히 ``modelSeeds`` 의 의미(seed 1개 = 모델 몇 개?)는 실행기마다 다를 수 있다(아래 참조).

seed/model 관례 (75모델)
------------------------
우리 파이프라인은 복합체당 **15 seed × 5 model = 75 도킹 모델**을 쓰고, 재현성을 위해
seed를 직접 고정한다(같은 seed = 같은 결과). 이 스크립트는 ``--seeds`` 개수만큼 seed 정수를
JSON의 ``modelSeeds`` 에 채운다(시작값 ``--seed-start``, 1씩 증가). seed당 몇 개 모델을
뽑는지(``--models``)는 실행기 옵션이라 JSON 본문에는 넣지 않고 ``README_notes.txt`` 에만
메모로 남긴다. WT와 모든 변이체에 **동일한 seed 세트**를 써야 비교가 공정하다(이 스크립트는
항상 동일 세트를 쓴다).

의존성
------
- 기본은 **표준 라이브러리만** 사용(argparse, csv, json, sys, re, pathlib, typing 등).
- ``--structure`` 로 구조 파일에서 서열을 뽑을 때만 외부 패키지가 필요하다:
  ``gemmi`` 가 있으면 우선 사용, 없으면 ``biopython``(Bio.PDB)을 시도한다.
  둘 다 없으면 친절한 오류 메시지와 함께 종료한다. (FASTA 입력만 쓰면 외부 의존 없음.)

테스트 상태 (중요)
------------------
- 이 파일은 합성 FASTA/CSV로 기본·엣지 케이스(이미-Ala/글리신/범위밖/다중체인/
  contacts CSV 두 형식/0-기반/잔기이름 불일치 경고)를 확인했고, gemmi가 설치된 환경에서
  합성 PDB로 ``--structure`` 추출 경로도 확인했다. 하지만 **실제 AF3 model.cif·실제
  cif_distance.py 출력으로는 직접 검증하지 못했다.** 특히 실제 AF3 버전의 JSON 입력
  스키마(version 값 등)는 환경/버전에 따라 다를 수 있다. 먼저 작은 예제로 ``manifest.csv`` 와
  생성된 FASTA/JSON 한두 개를 눈으로 검산한 뒤 사용하라.

사용 예시
---------
    # 1) FASTA 항체 + 항원, 대상 잔기 직접 지정
    python alanine_scan_prep.py \
        --fasta data/seqs/항체_WT.fasta \
        --ag-fasta data/seqs/항원.fasta \
        --positions B:31 B:33 B:100 B:101 \
        --seeds 15 --models 5 \
        --out data/alascan/inputs

    # 2) 구조에서 항체 체인 서열 추출 + cif_distance 접촉 CSV로 대상 자동 선정
    #    (cif_distance.py 의 --interface-out 또는 -o 출력 둘 다 인식)
    python alanine_scan_prep.py \
        --structure data/af3_out/WT/seed-1_sample-0/model.cif \
        --ab-chains B \
        --ag-fasta data/seqs/항원.fasta \
        --contacts data/contacts_interface.csv \
        --out data/alascan/inputs

    # 항상 옵션 재확인
    python alanine_scan_prep.py --help
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

# ---------------------------------------------------------------------------
# 단백질 표준 20 아미노산 1-letter 매핑 (3-letter -> 1-letter)
# 구조 파일/접촉 CSV에서 잔기 이름이 보통 3글자(예: ALA, TYR)로 나오므로 변환표가 필요하다.
# ---------------------------------------------------------------------------
THREE_TO_ONE: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
# 표준 1-letter 코드 집합(서열 검증용). 알라닌 = "A".
VALID_AA: frozenset = frozenset(THREE_TO_ONE.values())
ALANINE: str = "A"
GLYCINE: str = "G"

# AF3 job JSON 기본 메타값(버전에 따라 다를 수 있어 옵션으로도 바꿀 수 있게 둠).
# dialect는 공개 문서 기준 'alphafold3' 고정. version은 1~4가 유효(문서상 최신 4)이며
# 환경/실행기에 따라 호환 버전이 다를 수 있어 --json-version 으로 바꾼다.
DEFAULT_DIALECT: str = "alphafold3"
DEFAULT_JSON_VERSION: int = 2  # 보수적 기본값. 최신 실행기는 3/4를 요구할 수 있음(아래 메모/주석).


# ===========================================================================
# 1) 항체/항원 서열 읽기
# ===========================================================================
def read_fasta(path: Path) -> List[Tuple[str, str]]:
    """FASTA 파일을 (헤더, 서열) 목록으로 읽는다(표준 라이브러리만 사용).

    - 헤더는 '>' 다음 첫 공백 전까지를 ID로 본다(나머지 설명은 버린다).
    - 서열은 공백/숫자/'*'(종결자) 등을 제거하고 대문자로 통일한다.
    - 반환 순서는 파일에 나온 순서를 그대로 유지한다(체인 순서가 중요하므로).
    """
    records: List[Tuple[str, str]] = []
    header: Optional[str] = None
    seq_parts: List[str] = []

    def _flush() -> None:
        # 지금까지 모은 한 레코드를 records에 넣는다.
        if header is not None:
            seq = "".join(seq_parts)
            # 알파벳이 아닌 문자(숫자/공백/별표 등) 제거 후 대문자화
            seq = re.sub(r"[^A-Za-z]", "", seq).upper()
            records.append((header, seq))

    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                _flush()
                # '>' 제거 후 첫 토큰을 ID로
                rest = line[1:].strip()
                header = rest.split()[0] if rest else ""
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        _flush()  # 마지막 레코드 처리

    if not records:
        raise ValueError(f"FASTA에서 서열 레코드를 못 찾았습니다: {path}")
    # 빈 서열(헤더만 있고 본문 없음) 방어
    for h, s in records:
        if not s:
            raise ValueError(f"FASTA 레코드 '{h}'의 서열이 비어 있습니다: {path}")
    # 중복 헤더 방어(체인기호로 쓰므로 유일해야 함)
    seen_h: Dict[str, int] = {}
    for h, _ in records:
        seen_h[h] = seen_h.get(h, 0) + 1
    dups = sorted(h for h, c in seen_h.items() if c > 1)
    if dups:
        raise ValueError(
            f"FASTA 헤더 ID가 중복됩니다(체인기호로 써야 해서 유일해야 함): {dups} — {path}"
        )
    return records


def read_chains_from_structure(
    path: Path, chain_ids: Sequence[str]
) -> List[Tuple[str, str]]:
    """구조 파일(.cif/.pdb)에서 지정 체인의 1-letter 서열을 뽑는다.

    gemmi가 있으면 우선 사용하고, 없으면 biopython을 시도한다. 둘 다 없으면
    설치 안내와 함께 예외를 던진다. (이 경로를 쓸 때만 외부 의존이 생긴다.)

    반환: [(체인기호, 서열), ...]  — chain_ids 순서를 따른다.
    표준 20종이 아닌 잔기(HETATM·물·변형아미노산 등)는 건너뛴다(서열에 안 넣는다).
    ⚠️ 이 함수가 만드는 "서열 위치"는 1번부터 차례로 센 값이다(구조의 author 잔기번호와
       다를 수 있음). 이 스크립트 전체가 서열 위치 기준이라 일관된다.
    """
    wanted = list(chain_ids)
    result: List[Tuple[str, str]] = []

    # --- 경로 A: gemmi ---
    try:
        import gemmi  # type: ignore
    except ImportError:
        gemmi = None  # 아래 biopython 경로로 진행

    if gemmi is not None:
        st = gemmi.read_structure(str(path))
        if len(st) == 0:
            raise ValueError(f"구조에 모델이 없습니다: {path}")
        model = st[0]  # 보통 첫 모델만 사용
        # 체인기호 -> 서열
        seqs: Dict[str, List[str]] = {cid: [] for cid in wanted}
        for chain in model:
            if chain.name not in seqs:
                continue
            for residue in chain:
                one = THREE_TO_ONE.get(residue.name.upper())
                if one is not None:
                    seqs[chain.name].append(one)
        for cid in wanted:
            if not seqs[cid]:
                raise ValueError(
                    f"체인 '{cid}'에서 표준 아미노산 잔기를 못 찾았습니다(체인기호 오타? "
                    f"또는 단백질 체인이 아님): {path}"
                )
            result.append((cid, "".join(seqs[cid])))
        return result

    # --- 경로 B: biopython ---
    try:
        from Bio.PDB import MMCIFParser, PDBParser  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "구조 파일에서 서열을 뽑으려면 'gemmi' 또는 'biopython'이 필요합니다. "
            "둘 다 없습니다.\n  pip install gemmi   (권장)\n  pip install biopython\n"
            "FASTA 입력(--fasta)만 쓰면 외부 패키지 없이 동작합니다."
        ) from exc

    suffix = path.suffix.lower()
    if suffix in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", str(path))
    try:
        model = next(iter(structure))  # 첫 모델
    except StopIteration as exc:
        raise ValueError(f"구조에 모델이 없습니다: {path}") from exc
    present = {ch.id: ch for ch in model}
    for cid in wanted:
        if cid not in present:
            raise ValueError(
                f"구조에 체인 '{cid}'가 없습니다. 있는 체인: "
                f"{', '.join(sorted(map(str, present))) or '(없음)'}  파일: {path}"
            )
        ones: List[str] = []
        for residue in present[cid]:
            # 표준 아미노산만(이름이 THREE_TO_ONE 에 있는 잔기)
            resname = residue.get_resname().upper()
            one = THREE_TO_ONE.get(resname)
            if one is not None:
                ones.append(one)
        if not ones:
            raise ValueError(
                f"체인 '{cid}'에서 표준 아미노산 잔기를 못 찾았습니다: {path}"
            )
        result.append((cid, "".join(ones)))
    return result


# ===========================================================================
# 2) 대상 잔기(체인, 위치[, 기대아미노산]) 모으기
#    Target = (chain, position, expected_aa_or_None)
# ===========================================================================
Target = Tuple[str, int, Optional[str]]


def parse_positions(
    tokens: Sequence[str], default_chain: Optional[str]
) -> List[Target]:
    """``--positions`` 토큰들을 (체인, 1-기반위치, None) 목록으로 바꾼다.

    공백과 쉼표 둘 다 구분자로 받는다(문서 예시가 쉼표 구분이라 둘 다 허용).
    허용 형식:
      - "B:31"  -> (체인 B, 31번, None)
      - "31"    -> (default_chain, 31번, None)   # 체인을 안 적으면 기본 체인 사용
    --positions 로 줄 때는 잔기 이름 정보가 없으므로 expected_aa 는 None(검사 안 함).
    default_chain 이 None 인데 체인이 안 적힌 토큰이 있으면 오류(어느 체인인지 모호).
    """
    out: List[Target] = []
    # 토큰 안에 쉼표가 섞여 들어올 수 있으니 한 번 더 분해
    flat: List[str] = []
    for tok in tokens:
        flat.extend(piece for piece in re.split(r"[,\s]+", tok) if piece)
    for tok in flat:
        if ":" in tok:
            chain, num = tok.split(":", 1)
            chain = chain.strip()
            num = num.strip()
        else:
            chain = default_chain or ""
            num = tok
        if not chain:
            raise ValueError(
                f"'{tok}'에 체인이 없습니다. 항체 체인이 여러 개면 'B:31'처럼 "
                f"'체인:번호' 형식으로 주세요."
            )
        try:
            pos = int(num)
        except ValueError as exc:
            raise ValueError(f"위치가 정수가 아닙니다: '{tok}'") from exc
        out.append((chain, pos, None))
    if not out:
        raise ValueError("--positions 에 유효한 잔기가 하나도 없습니다.")
    return out


# cif_distance.py 잔기 라벨 형식: '{체인}/{잔기이름}{author번호}{삽입코드}'  예) 'B/SER100A'
# (cif_distance.py 의 residue_label() 와 동일 형식.)
_RESLABEL_RE = re.compile(
    r"^\s*(?P<chain>[^/]+)/(?P<name>[A-Za-z]{1,3})(?P<num>-?\d+)(?P<icode>[A-Za-z]?)\s*$"
)


def _parse_residue_label(label: str) -> Optional[Tuple[str, str, int]]:
    """cif_distance 라벨('B/SER100A')을 (체인, 잔기이름, author번호)로 푼다.

    형식이 안 맞으면 None. 삽입코드(icode)는 버린다(서열 위치 매핑에 쓸 수 없으므로 —
    삽입코드가 있는 잔기는 author번호↔서열위치 매핑이 특히 위험해 잔기이름 검사로 거른다).
    """
    m = _RESLABEL_RE.match(label)
    if not m:
        return None
    chain = m.group("chain").strip()
    name = m.group("name").upper()
    try:
        num = int(m.group("num"))
    except ValueError:
        return None
    return chain, name, num


def parse_contacts_csv(
    path: Path, antibody_chains: Sequence[str]
) -> List[Target]:
    """``cif_distance.py`` 가 만든 계면 접촉 CSV에서 항체 쪽 접촉 잔기를 뽑는다.

    cif_distance.py 의 **두 가지 출력 형식**을 모두 인식한다(2026-06-29 실측):
      (가) 계면 잔기 목록 CSV (``--interface-out``):
           컬럼 = side, chain_group, residue   / side='B'가 항체 쪽, residue='B/SER100A'
      (나) 잔기쌍 거리 CSV (``-o``/``--out``):
           컬럼 = residue_a, residue_b, distance_A, metric
           관례상 residue_a=항원 쪽(체인 A), residue_b=항체 쪽(체인 B). 양쪽 라벨을 모두
           훑어 antibody_chains 에 속하는 라벨만 항체 잔기로 취한다.

    반환: (체인, author번호, 기대아미노산1글자) 목록.
      - 라벨에 잔기이름이 들어 있어 expected_aa 를 같이 돌려준다(좌표계 검증용).
      - 중복 (체인, 번호)는 한 번만 남긴다(원래 등장 순서 유지).

    ⚠️ cif_distance 의 번호는 **author 잔기번호**다. 이 스크립트는 서열 위치로 변이를
       적용하므로, 호출부에서 expected_aa 로 일치 검사를 한다(불일치 시 건너뜀+경고).
    """
    ab_set = set(antibody_chains)
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames
        if not fieldnames:
            raise ValueError(f"접촉 CSV에 헤더가 없습니다: {path}")
        lower = {fn.lower(): fn for fn in fieldnames}

        # 어떤 형식인지 헤더로 판별하고, 훑을 "라벨 컬럼"들을 정한다.
        label_cols: List[str]
        if "residue" in lower:
            # (가) 계면 잔기 목록 형식
            label_cols = [lower["residue"]]
        elif "residue_a" in lower or "residue_b" in lower:
            # (나) 잔기쌍 형식 — 양쪽 라벨 모두 훑는다.
            label_cols = [lower[c] for c in ("residue_a", "residue_b") if c in lower]
        else:
            raise ValueError(
                "접촉 CSV 형식을 알 수 없습니다. cif_distance.py 출력은 보통\n"
                "  - 계면 목록: side,chain_group,residue (--interface-out)\n"
                "  - 잔기쌍   : residue_a,residue_b,distance_A,metric (-o)\n"
                f"인데, 받은 헤더는 {fieldnames} 입니다.\n"
                "  cif_distance.py 로 다시 뽑거나, --positions 로 직접 (체인:번호)를 주세요."
            )

        seen: set = set()
        out: List[Target] = []
        unparsable = 0
        for row in reader:
            for col in label_cols:
                raw = (row.get(col) or "").strip()
                if not raw:
                    continue
                parsed = _parse_residue_label(raw)
                if parsed is None:
                    unparsable += 1
                    continue
                chain, name, num = parsed
                if chain not in ab_set:
                    continue  # 항원 쪽 라벨은 버림
                expected_aa = THREE_TO_ONE.get(name)  # 비표준 잔기면 None
                key = (chain, num)
                if key in seen:
                    continue
                seen.add(key)
                out.append((chain, num, expected_aa))

    if not out:
        hint = ""
        if unparsable:
            hint = (
                f"\n  (참고: 라벨 {unparsable}개가 'B/SER100A' 형식과 안 맞아 건너뜀 — "
                "cif_distance.py 가 만든 CSV가 맞는지 확인)"
            )
        raise ValueError(
            f"접촉 CSV에서 항체 체인({', '.join(sorted(antibody_chains)) or '?'})의 "
            f"접촉 잔기를 하나도 못 찾았습니다: {path}{hint}\n"
            "  (--ab-chains 가 CSV의 항체 체인기호와 같은지 확인하거나 --positions 사용)"
        )
    return out


# ===========================================================================
# 3) 변이 적용 + 변이ID
# ===========================================================================
def make_mutation_id(chain: str, wt_aa: str, pos: int) -> str:
    """변이ID 문자열을 만든다. 예) (B, 'Y', 100) -> 'B_Y100A' (→A는 알라닌 고정)."""
    return f"{chain}_{wt_aa}{pos}{ALANINE}"


def apply_alanine(
    seq: str, pos_1based: int, index_base: int
) -> Tuple[str, str]:
    """서열의 한 위치를 알라닌으로 치환한다.

    pos_1based: 사용자가 준 위치(기본 1-기반).
    index_base: 사용자가 0-기반으로 줬으면 0, 1-기반이면 1.
    반환: (원래잔기 1-letter, 치환된 새 서열).
    범위를 벗어나거나 비표준 잔기면 예외를 던진다(호출부에서 잡아 건너뜀 처리).
    """
    # 내부 인덱스(0부터)로 변환
    idx = pos_1based - index_base
    if idx < 0 or idx >= len(seq):
        raise IndexError(
            f"위치 {pos_1based}(base={index_base})가 서열 길이 {len(seq)} 범위를 벗어남"
        )
    wt_aa = seq[idx]
    if wt_aa not in VALID_AA:
        raise ValueError(f"위치 {pos_1based}의 잔기 '{wt_aa}'가 표준 20종이 아님")
    new_seq = seq[:idx] + ALANINE + seq[idx + 1:]
    return wt_aa, new_seq


# ===========================================================================
# 4) 출력 파일 만들기 (FASTA / AF3 job JSON)
# ===========================================================================
def _wrap(seq: str, width: int = 60) -> List[str]:
    """서열을 width 글자마다 줄바꿈한 줄 목록으로 반환(FASTA 관례)."""
    return [seq[i:i + width] for i in range(0, len(seq), width)] or [""]


def write_fasta(
    path: Path,
    antigen: Optional[Tuple[str, str]],
    antibody_chains: List[Tuple[str, str]],
) -> None:
    """[항원 레코드] + [항체 체인 레코드들] 을 담은 FASTA 한 개를 쓴다.

    antigen: (체인기호, 서열) 또는 None(항원 없이 항체만).
    antibody_chains: [(체인기호, 서열), ...] — 변이가 이미 적용된 상태로 들어온다.
    체인 순서 관례에 따라 항원을 먼저(위), 항체를 그다음(아래)에 쓴다.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: List[str] = []
    if antigen is not None:
        ag_id, ag_seq = antigen
        lines.append(f">{ag_id}")
        lines.extend(_wrap(ag_seq))
    for cid, seq in antibody_chains:
        lines.append(f">{cid}")
        lines.extend(_wrap(seq))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_job_json(
    name: str,
    antigen: Optional[Tuple[str, str]],
    antibody_chains: List[Tuple[str, str]],
    seeds: List[int],
    dialect: str,
    json_version: int,
) -> dict:
    """AF3 job JSON(공개 스키마 형태)을 dict로 만든다.

    ⚠️ AF3 입력 파서는 보통 **스키마에 없는 최상위 키를 거부**할 수 있어, 사람용 메모는
       여기 넣지 않는다(메모는 README_notes.txt 에 따로 적는다). 실제 서버/실행기 버전의
       스펙과 키 이름·구조·version 값을 반드시 대조해 맞춰라(특히 modelSeeds 의미).
       자세한 형식 가정은 파일 상단 docstring "AF3 입력 포맷 가정" 참조.
    """
    sequences: List[dict] = []
    # 관례상 항원을 먼저(체인 A 자리), 항체 체인들을 그다음에.
    if antigen is not None:
        ag_id, ag_seq = antigen
        sequences.append({"protein": {"id": ag_id, "sequence": ag_seq}})
    for cid, seq in antibody_chains:
        sequences.append({"protein": {"id": cid, "sequence": seq}})

    job = {
        "name": name,
        "modelSeeds": seeds,             # 정수 seed 배열(재현성 위해 고정), 최소 1개
        "sequences": sequences,
        "dialect": dialect,
        "version": json_version,
    }
    return job


def write_json(path: Path, job: dict) -> None:
    """job dict를 들여쓰기된 JSON 파일로 쓴다."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(job, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def write_notes(path: Path, seeds: List[int], models_per_seed: int,
                dialect: str, json_version: int) -> None:
    """사람이 읽는 메모 파일(README_notes.txt)을 쓴다.

    AF3 JSON 본문에는 못 넣는 운영 메모(seed/model 관례, 형식 가정)를 여기 모은다.
    """
    text = (
        "in-silico alanine scan 입력 메모\n"
        "================================\n"
        f"- 고정 seed 세트: {seeds}  (총 {len(seeds)} seed)\n"
        f"- seed당 모델 수(실행기 옵션): {models_per_seed}  "
        f"→ 관례상 {len(seeds)} seed x {models_per_seed} model = "
        f"{len(seeds) * models_per_seed} 모델\n"
        "  (modelSeeds 1개가 모델 몇 개를 뜻하는지는 실행기마다 다를 수 있음. "
        "models_per_seed 는 JSON 본문에 넣지 않았다 — 실행기에서 직접 지정할 것.)\n"
        f"- AF3 job JSON: dialect='{dialect}', version={json_version} "
        "(version은 환경/실행기 버전에 따라 1~4 중 호환값으로 조정 필요).\n"
        "- 체인 순서 관례: 항원 먼저(보통 A), 항체 그다음(B, C ...).\n"
        "- 실제 AF3 제출은 수동. 점수 추출/WT 비교는 af3_extract_features.py(예정).\n"
        "- WT 와 모든 변이체는 동일 seed 세트를 쓴다(공정 비교).\n"
    )
    path.write_text(text, encoding="utf-8")


# ===========================================================================
# 5) 입력 정리 헬퍼
# ===========================================================================
def _resolve_antibody_chains(
    args: argparse.Namespace,
) -> List[Tuple[str, str]]:
    """입력 옵션에서 항체 체인 [(체인기호, 서열), ...] 를 만든다(FASTA 또는 구조)."""
    if args.fasta:
        records = read_fasta(Path(args.fasta))
        if args.ab_chains:
            # 지정한 체인기호만, 지정 순서대로 고른다.
            by_id = {h: s for h, s in records}
            missing = [c for c in args.ab_chains if c not in by_id]
            if missing:
                raise ValueError(
                    f"FASTA에 다음 체인(헤더 ID)이 없습니다: {missing}. "
                    f"있는 ID: {list(by_id)}"
                )
            return [(c, by_id[c]) for c in args.ab_chains]
        # 안 주면 모든 레코드를 항체 체인으로 본다(헤더 ID를 체인기호로 사용).
        return [(h, s) for h, s in records]
    # 구조 입력
    if not args.ab_chains:
        raise ValueError("--structure 를 쓸 때는 --ab-chains 로 항체 체인을 지정해야 합니다.")
    return read_chains_from_structure(Path(args.structure), args.ab_chains)


def _resolve_antigen(args: argparse.Namespace) -> Optional[Tuple[str, str]]:
    """항원 (체인기호, 서열) 또는 None을 만든다."""
    if not args.ag_fasta:
        return None
    records = read_fasta(Path(args.ag_fasta))
    # 항원은 보통 단일 체인. 여러 레코드면 첫 번째만 쓰되 경고.
    if len(records) > 1:
        print(
            f"[경고] 항원 FASTA에 레코드가 {len(records)}개 있습니다. 첫 번째만 "
            f"항원으로 사용합니다(다중 체인 항원이면 코드 수정 필요).",
            file=sys.stderr,
        )
    ag_id_raw, ag_seq = records[0]
    # 항원 체인기호: 관례상 'A'를 권장하되, 사용자가 --ag-chain-id 로 바꿀 수 있게.
    ag_id = args.ag_chain_id or ag_id_raw or "A"
    return (ag_id, ag_seq)


# ===========================================================================
# 6) 메인 로직
# ===========================================================================
def main(argv: Optional[List[str]] = None) -> int:
    """CLI 진입점. 반환값은 프로세스 종료 코드(0=정상)."""
    parser = argparse.ArgumentParser(
        description=(
            "AF3용 in-silico 알라닌 스캐닝 전처리: 대상 잔기를 하나씩 알라닌으로 바꾼 "
            "변이 서열 + AF3 입력(FASTA/JSON) + 매핑표(manifest.csv)를 생성한다. "
            "WT 포함. 실제 AF3 제출은 수동(범위 밖)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # --- 항체 서열 입력(둘 중 하나) ---
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--fasta", type=str, help="항체 FASTA 파일(레코드=체인).")
    src.add_argument(
        "--structure", type=str,
        help="항체 구조 파일(.cif/.pdb). --ab-chains 필요. gemmi/biopython 필요.",
    )
    parser.add_argument(
        "--ab-chains", nargs="+", metavar="CHAIN",
        help=(
            "처리할 항체 체인기호들(예: B 또는 B C). --structure 에는 필수, "
            "--fasta 에서는 일부만 고를 때 사용(안 주면 FASTA의 모든 레코드)."
        ),
    )

    # --- 항원 ---
    parser.add_argument(
        "--ag-fasta", type=str,
        help="항원 FASTA(AF3 입력에 함께 넣음). 안 주면 항체만 입력 생성.",
    )
    parser.add_argument(
        "--ag-chain-id", type=str, default=None,
        help="항원 체인기호(기본: FASTA 헤더 ID, 그래도 없으면 'A'). 관례상 'A' 권장.",
    )

    # --- 치환 대상(둘 중 하나) ---
    tgt = parser.add_mutually_exclusive_group(required=True)
    tgt.add_argument(
        "--positions", nargs="+", metavar="POS",
        help="대상 잔기 목록. '체인:번호'(예: B:31) 또는 '번호'(체인 1개일 때). "
             "공백/쉼표 구분 모두 허용.",
    )
    tgt.add_argument(
        "--contacts", type=str,
        help="cif_distance.py 가 만든 계면 접촉 CSV. 항체 쪽 접촉 잔기를 자동 선정. "
             "(--interface-out 형식·-o 잔기쌍 형식 모두 인식)",
    )

    # --- 번호 좌표계 ---
    parser.add_argument(
        "--index-base", type=int, choices=(0, 1), default=1,
        help="입력 위치 번호의 기준(1=1-기반 서열위치 기본, 0=0-기반).",
    )
    parser.add_argument(
        "--no-resname-check", action="store_true",
        help="--contacts 의 잔기이름(예: TYR)과 서열 위치의 아미노산 일치 검사를 끈다. "
             "(기본은 검사 ON — 불일치면 건너뛰고 경고. 끄면 좌표계 어긋남 위험.)",
    )

    # --- seed/model ---
    parser.add_argument(
        "--seeds", type=int, default=15,
        help="JSON modelSeeds 에 넣을 seed 개수(기본 15, 최소 1). WT/변이체 동일 세트.",
    )
    parser.add_argument(
        "--seed-start", type=int, default=1,
        help="seed 시작 정수(기본 1, 1씩 증가). 재현성 위해 고정값 사용.",
    )
    parser.add_argument(
        "--models", type=int, default=5,
        help="seed당 모델 수(기본 5). JSON엔 안 넣고 README_notes.txt 메모(실행기 옵션). "
             "15x5=75 관례.",
    )

    # --- glycine 정책 ---
    parser.add_argument(
        "--skip-glycine", action="store_true",
        help="원래 잔기가 글리신(G)인 자리도 건너뛴다(곁사슬이 사실상 없어 스캔 의미 적음).",
    )

    # --- JSON 메타(버전 따라 다름) ---
    parser.add_argument(
        "--json-dialect", type=str, default=DEFAULT_DIALECT,
        help=f"AF3 job JSON 'dialect' 값(기본 '{DEFAULT_DIALECT}', 버전 따라 다름).",
    )
    parser.add_argument(
        "--json-version", type=int, default=DEFAULT_JSON_VERSION,
        help=f"AF3 job JSON 'version' 값(기본 {DEFAULT_JSON_VERSION}; 문서상 1~4 유효, "
             "실행기 버전에 맞춰 조정).",
    )

    # --- 출력 ---
    parser.add_argument(
        "--out", type=str, required=True,
        help="출력 폴더(fasta/, json/, manifest.csv, README_notes.txt 가 그 아래 생성됨).",
    )
    parser.add_argument(
        "--no-wt", action="store_true",
        help="WT(야생형) 입력을 만들지 않는다(보통은 만들어 비교 기준으로 씀).",
    )

    args = parser.parse_args(argv)

    # --- 옵션 값 방어 ---
    if args.seeds < 1:
        print("[오류] --seeds 는 1 이상이어야 합니다(AF3 modelSeeds 최소 1).", file=sys.stderr)
        return 2
    if args.models < 1:
        print("[오류] --models 는 1 이상이어야 합니다.", file=sys.stderr)
        return 2

    # --- 입력 정리 ---
    try:
        antibody_chains = _resolve_antibody_chains(args)
        antigen = _resolve_antigen(args)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"[오류] 서열 입력 실패: {exc}", file=sys.stderr)
        return 2

    ab_chain_ids = [cid for cid, _ in antibody_chains]
    ab_seq_by_chain: Dict[str, str] = {cid: seq for cid, seq in antibody_chains}
    default_chain = ab_chain_ids[0] if len(ab_chain_ids) == 1 else None

    # --- 대상 잔기 모으기 ---
    try:
        if args.positions:
            targets = parse_positions(args.positions, default_chain)
        else:
            targets = parse_contacts_csv(Path(args.contacts), ab_chain_ids)
    except (OSError, ValueError) as exc:
        print(f"[오류] 대상 잔기 파싱 실패: {exc}", file=sys.stderr)
        return 2

    # 대상 체인이 항체 체인 목록에 있는지 확인
    bad = sorted({c for c, _, _ in targets if c not in ab_seq_by_chain})
    if bad:
        print(
            f"[오류] 대상 잔기의 체인 {bad}가 항체 체인 {ab_chain_ids}에 없습니다.",
            file=sys.stderr,
        )
        return 2

    # 중복 (체인,위치) 제거(순서 유지). expected_aa 는 처음 것을 유지.
    seen_t: set = set()
    uniq_targets: List[Target] = []
    for chain, pos, exp in targets:
        key = (chain, pos)
        if key not in seen_t:
            seen_t.add(key)
            uniq_targets.append((chain, pos, exp))
    targets = uniq_targets

    out_dir = Path(args.out)
    fasta_dir = out_dir / "fasta"
    json_dir = out_dir / "json"
    out_dir.mkdir(parents=True, exist_ok=True)

    # seed 세트(고정): seed_start 부터 1씩 증가
    seeds = [args.seed_start + i for i in range(args.seeds)]

    # manifest 행 모으기
    manifest_rows: List[Dict[str, str]] = []

    def _write_pair(
        mut_id: str,
        mutated_ab: List[Tuple[str, str]],
    ) -> Tuple[str, str]:
        """주어진 항체 체인 세트로 FASTA+JSON 한 쌍을 쓰고 (fasta경로, json경로) 반환."""
        fasta_path = fasta_dir / f"{mut_id}.fasta"
        json_path = json_dir / f"{mut_id}.json"
        write_fasta(fasta_path, antigen, mutated_ab)
        job = build_job_json(
            name=mut_id,
            antigen=antigen,
            antibody_chains=mutated_ab,
            seeds=seeds,
            dialect=args.json_dialect,
            json_version=args.json_version,
        )
        write_json(json_path, job)
        return str(fasta_path), str(json_path)

    # --- WT(야생형) ---
    if not args.no_wt:
        wt_ab = [(cid, ab_seq_by_chain[cid]) for cid in ab_chain_ids]
        fp, jp = _write_pair("WT", wt_ab)
        manifest_rows.append({
            "mutation_id": "WT",
            "kind": "WT",
            "chain": "",
            "position": "",
            "wt_aa": "",
            "mut_aa": "",
            "note": "야생형 기준(비교용). 변이체와 동일 seed 세트.",
            "fasta_path": fp,
            "json_path": jp,
        })

    # --- 각 변이체 ---
    n_made = 0
    n_skipped = 0
    n_resname_mismatch = 0
    for chain, pos, expected_aa in targets:
        base_seq = ab_seq_by_chain[chain]
        # 먼저 원래 잔기를 확인(범위/표준 검사)
        try:
            wt_aa, mutated_seq = apply_alanine(base_seq, pos, args.index_base)
        except (IndexError, ValueError) as exc:
            # 범위 밖/비표준 -> 건너뛰되 manifest에 SKIPPED로 기록
            n_skipped += 1
            manifest_rows.append({
                "mutation_id": f"{chain}_pos{pos}",
                "kind": "SKIPPED",
                "chain": chain,
                "position": str(pos),
                "wt_aa": "",
                "mut_aa": "",
                "note": f"건너뜀: {exc}",
                "fasta_path": "",
                "json_path": "",
            })
            print(f"[건너뜀] {chain}:{pos} — {exc}", file=sys.stderr)
            continue

        # 좌표계 안전장치: --contacts 라벨의 잔기이름과 서열 위치 아미노산이 다르면
        # author번호↔서열위치 어긋남일 가능성 → 건너뛰고 시끄럽게 경고.
        if (expected_aa is not None and not args.no_resname_check
                and wt_aa != expected_aa):
            n_skipped += 1
            n_resname_mismatch += 1
            manifest_rows.append({
                "mutation_id": f"{chain}_pos{pos}",
                "kind": "SKIPPED",
                "chain": chain,
                "position": str(pos),
                "wt_aa": wt_aa,
                "mut_aa": "",
                "note": (
                    f"건너뜀: 잔기이름 불일치 — CSV는 {expected_aa}, 서열 위치 {pos}은(는) "
                    f"{wt_aa}. author번호↔서열위치 좌표계가 어긋났을 수 있음. "
                    "cif_distance 와 같은 좌표계인지 확인 후 --positions 로 직접 지정하거나 "
                    "--no-resname-check 로 강제 진행."
                ),
                "fasta_path": "",
                "json_path": "",
            })
            print(
                f"[건너뜀/경고] {chain}:{pos} 잔기이름 불일치 "
                f"(CSV={expected_aa}, 서열={wt_aa}) — 좌표계 확인 필요.",
                file=sys.stderr,
            )
            continue

        # 이미 알라닌이면 바꿀 게 없음 -> 건너뜀
        if wt_aa == ALANINE:
            n_skipped += 1
            manifest_rows.append({
                "mutation_id": make_mutation_id(chain, wt_aa, pos),
                "kind": "SKIPPED",
                "chain": chain,
                "position": str(pos),
                "wt_aa": wt_aa,
                "mut_aa": ALANINE,
                "note": "이미 알라닌(A) — 치환할 곁사슬 없음, 건너뜀.",
                "fasta_path": "",
                "json_path": "",
            })
            print(f"[건너뜀] {chain}:{pos} 이미 알라닌(A).", file=sys.stderr)
            continue

        # 글리신 정책(옵션): 건너뛰거나(표시), 아니면 진행하되 manifest에 메모
        gly_note = ""
        if wt_aa == GLYCINE:
            if args.skip_glycine:
                n_skipped += 1
                manifest_rows.append({
                    "mutation_id": make_mutation_id(chain, wt_aa, pos),
                    "kind": "SKIPPED",
                    "chain": chain,
                    "position": str(pos),
                    "wt_aa": wt_aa,
                    "mut_aa": ALANINE,
                    "note": "글리신(G) — 곁사슬 사실상 없음, --skip-glycine 으로 건너뜀.",
                    "fasta_path": "",
                    "json_path": "",
                })
                print(f"[건너뜀] {chain}:{pos} 글리신(G) (--skip-glycine).", file=sys.stderr)
                continue
            gly_note = "원래 글리신(G) — 알라닌과 차이가 작아 스캔 의미가 약할 수 있음. "

        # 변이 적용한 항체 체인 세트 구성(대상 체인만 변이, 나머지는 WT 그대로)
        mut_id = make_mutation_id(chain, wt_aa, pos)
        mutated_ab: List[Tuple[str, str]] = []
        for cid in ab_chain_ids:
            if cid == chain:
                mutated_ab.append((cid, mutated_seq))
            else:
                mutated_ab.append((cid, ab_seq_by_chain[cid]))

        fp, jp = _write_pair(mut_id, mutated_ab)
        n_made += 1
        manifest_rows.append({
            "mutation_id": mut_id,
            "kind": "MUT",
            "chain": chain,
            "position": str(pos),
            "wt_aa": wt_aa,
            "mut_aa": ALANINE,
            "note": gly_note + "단일치환 변이체. WT와 동일 seed 세트.",
            "fasta_path": fp,
            "json_path": jp,
        })

    # --- manifest.csv 쓰기 ---
    manifest_path = out_dir / "manifest.csv"
    fieldnames = [
        "mutation_id", "kind", "chain", "position",
        "wt_aa", "mut_aa", "note", "fasta_path", "json_path",
    ]
    with manifest_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest_rows)

    # --- 사람용 메모 파일 ---
    notes_path = out_dir / "README_notes.txt"
    write_notes(notes_path, seeds, args.models, args.json_dialect, args.json_version)

    # --- 요약 출력 ---
    print()
    print("== 알라닌 스캐닝 전처리 완료 ==")
    print(f"  항체 체인: {', '.join(ab_chain_ids)}")
    if antigen is not None:
        print(f"  항원 체인: {antigen[0]} (길이 {len(antigen[1])})")
    else:
        print("  항원: (없음 — 항체만 입력 생성)")
    print(f"  대상 잔기 수: {len(targets)}")
    print(f"  생성된 변이체: {n_made}개 / 건너뜀: {n_skipped}개"
          f"{' / WT 포함' if not args.no_wt else ''}")
    if n_resname_mismatch:
        print(f"  ⚠️ 잔기이름 불일치로 건너뛴 잔기: {n_resname_mismatch}개 — "
              "좌표계(author번호 vs 서열위치)를 꼭 확인하세요.")
    print(f"  seed 세트(고정): {seeds}  (seed당 {args.models} 모델 관례 → "
          f"{len(seeds) * args.models}모델)")
    print(f"  출력 폴더: {out_dir.resolve()}")
    print(f"    - FASTA: {fasta_dir}/")
    print(f"    - JSON : {json_dir}/")
    print(f"    - 매핑표: {manifest_path}")
    print(f"    - 메모  : {notes_path}")
    print()
    print("다음 단계(수동): 각 JSON/FASTA를 AF3에 제출 → 변이체마다 75모델 예측 → "
          "af3_extract_features.py 로 점수 추출 → WT와 Δscore 비교.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
