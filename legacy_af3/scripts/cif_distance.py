#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""구조 파일(.cif / .pdb)에서 두 체인 사이 잔기 간 거리를 계산하는 CLI.

배경 (왜 이 스크립트가 필요한가)
--------------------------------
AlphaFold3(AF3)로 항체-항원 도킹 모델(model.cif)을 만들면, "어떤 잔기들이
실제로 맞닿아 있는가"(= 계면, interface)를 알고 싶을 때가 많다. 계면 잔기는
- 도킹 포즈가 그럴듯한지 눈으로/수치로 확인할 때,
- 계면 접촉 그래프(메인 프로젝트의 핵심 피처)를 만들 때,
- alanine scanning 같은 후속 분석에서 "어디를 건드릴지" 고를 때
모두 쓰인다.

이 스크립트는 구조 파일 하나를 받아서, 두 체인(기본: A=항원, B=항체) 사이의
모든 잔기쌍에 대해 거리를 재고, 그중 가까운(접촉하는) 쌍만 골라 계면 잔기
목록을 만들어 준다. 결과는 CSV로 저장한다.

⚠️ 체인 배정 관례 (AF3 입력 순서와 맞춰야 함)
--------------------------------------------
우리 파이프라인 관례상 AF3에 "위=항원, 아래=항체" 순서로 넣으면 결과에서
chain A = 항원, chain B = 항체로 매핑된다. 그래서 이 스크립트의 기본값도
--chain-a A (항원), --chain-b B (항체)다.
단, 항체를 Fab/scFv로 넣으면 체인이 여러 개일 수 있다(예: Fab = 중쇄절편 +
경쇄 = 2개 체인). 그런 경우 --chain-b 에 항체 쪽 체인 기호를 정확히 지정해야
한다(여러 체인을 넘기려면 --chain-b H L 처럼 공백으로 나열). 어떤 체인이 있는지
모르면 --list-chains 로 먼저 확인하라.

거리 계산 방식 (--metric)
-------------------------
- ca        : 두 잔기의 Cα(알파 탄소) 원자끼리의 거리. 가장 빠르고 흔한 근사.
              (Cα = 단백질 골격에서 곁사슬이 붙는 중심 탄소. 잔기 1개에 1개)
- centroid  : 두 잔기의 "모든 원자 평균 좌표(무게중심)"끼리의 거리.
              곁사슬 위치까지 어느 정도 반영하는 중간 정도의 방식.
- min-heavy : 두 잔기에 속한 heavy atom(수소가 아닌 무거운 원자) 쌍 중
              "가장 가까운 한 쌍"의 거리. 실제 접촉(원자가 맞닿음)을 가장
              충실히 잡는 방식이지만 가장 느리다. 계면 정의에는 보통 이게 좋다.
              (heavy atom = 수소(H)를 제외한 원자. 결정/예측 구조엔 수소가
              아예 없는 경우가 많아 보통 자동으로 heavy atom만 남는다.)

잔기 번호 (중요)
----------------
출력 CSV의 잔기 번호는 author numbering(저자 번호)을 쓴다. 즉 구조 파일에
원래 적혀 있는 번호(PyMOL 등에서 보이는 그 번호) 그대로다. 내부 일련번호가
아니므로 PyMOL의 `select` 등에 바로 붙여 쓸 수 있다.

의존성 (외부 패키지)
--------------------
- 표준 라이브러리: argparse, csv, sys, math, pathlib, typing.
- 좌표 로드는 외부 패키지가 하나 필요하며, 다음 우선순위로 자동 선택한다:
    1) gemmi  (권장. .cif/.pdb 모두 빠르고 안정적으로 읽음)
         설치: pip install gemmi
    2) Biopython (gemmi가 없을 때 fallback)
         설치: pip install biopython
  둘 다 없으면 친절한 안내 메시지를 내고 종료한다.
  (numpy는 쓰지 않는다 — 좌표 연산은 순수 파이썬 math로 처리해 의존성을 줄였다.)

테스트 상태
-----------
- 로직(거리 metric 3종 ca/centroid/min-heavy, 물·헤테로 필터, 계면 선별, CSV
  출력, author-numbering 라벨, --list-chains, 에러 경로)은 gemmi·Biopython 두
  백엔드 모두에서 합성(가짜) 구조로 end-to-end 실행 검증했다. 두 백엔드가 같은
  결과를 냈고, label_seq 와 auth_seq 가 다른 mmCIF에서도 둘 다 author 번호를
  보고함을 확인했다.
- 단, 실제 AF3 model.cif/실제 PDB 복합체로는 아직 대조하지 못했다. 처음 쓸 때는
  작은/익숙한 구조 하나로 결과를 PyMOL 등과 대조해 본 뒤 본격적으로 사용하라.
  특히 체인 기호·잔기 번호가 의도대로 나오는지 확인할 것.

사용 예시
---------
    # 0) 어떤 체인이 들어있는지 먼저 확인
    python cif_distance.py model.cif --list-chains

    # 1) 기본: A(항원) vs B(항체), min-heavy 5.0 Å 이내 계면 잔기쌍을 CSV로
    python cif_distance.py model.cif --metric min-heavy -o pairs.csv

    # 2) 계면 잔기 유니크 목록도 같이 저장
    python cif_distance.py model.cif -o pairs.csv --interface-out interface.csv

    # 3) 항체가 Fab(중쇄 H + 경쇄 L)일 때
    python cif_distance.py model.cif --chain-a A --chain-b H L -o pairs.csv

    # 4) Cα 거리, 컷오프 8 Å
    python cif_distance.py model.cif --metric ca --cutoff 8.0 -o pairs.csv

    # 5) 컷오프 이내만이 아니라 모든 잔기쌍 거리를 다 저장
    python cif_distance.py model.cif --all-pairs -o all.csv
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

# --------------------------------------------------------------------------
# 자료형 정의 (초보용 설명)
# --------------------------------------------------------------------------
# Atom    : 원자 하나 = (원자이름, 원소기호, x, y, z).
#   - 원자이름(예: 'CA', 'CB', 'N'): Cα 선택에 필요하다.
#   - 원소기호(예: 'C', 'N', 'H'): heavy/수소 구분에 쓴다.
# Residue : 잔기 하나 = (체인기호, author 잔기번호, 삽입코드, 잔기이름, 원자목록).
#   - author 잔기번호: 구조 파일에 적힌 번호(PyMOL에 보이는 그 번호).
#   - 삽입코드(insertion code): 같은 번호를 더 쪼갤 때 붙는 한 글자(보통 빈칸).
#     항체 CDR 같은 곳에서 100A, 100B처럼 나타날 수 있어 함께 보관한다.
Atom = Tuple[str, str, float, float, float]
Residue = Tuple[str, int, str, str, List[Atom]]
Point = Tuple[float, float, float]

# 단백질 표준 아미노산 20종 + 흔한 변형 일부(이름 기준). 분석 대상의 "기본"으로 본다.
# (strip_het 기본 동작에서, 이 집합에 없는 잔기명은 물/이온/리간드로 보고 제거한다.)
STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # 자주 보이는 변형/대체 표기
    "MSE",  # 셀레노메티오닌(메티오닌 대체) — 단백질로 취급
    "SEC", "PYL", "HSD", "HSE", "HSP",
}

# 물 잔기 이름(--strip-water-only 에서 이것만 제거).
WATER_NAMES = {"HOH", "WAT", "DOD"}


# ==========================================================================
# 1) 구조 로드 — gemmi 우선, 없으면 Biopython fallback
# ==========================================================================
def _is_heavy(element: str) -> bool:
    """원소기호를 보고 heavy atom(수소가 아닌 원자)인지 판단한다.

    수소는 'H' 또는 'D'(중수소)로 표기된다. 빈 원소기호는 안전하게 heavy로 본다
    (정보가 없으면 버리지 않는다는 보수적 선택).
    """
    e = (element or "").strip().upper()
    return e not in ("H", "D")


def _load_with_gemmi(path: Path, model_index: int) -> List[Residue]:
    """gemmi로 구조를 읽어 Residue 리스트로 변환한다(원자 이름 포함).

    gemmi는 .cif와 .pdb를 모두 읽을 수 있는 가벼운 구조 라이브러리다.
    model_index: 여러 모델(NMR 등)이 있을 때 몇 번째 모델을 쓸지(0부터).
    """
    import gemmi  # 호출 시점에 import (없으면 caller가 fallback 처리)

    st = gemmi.read_structure(str(path))
    if len(st) == 0:
        raise ValueError("구조에 모델이 하나도 없습니다.")
    if model_index < 0 or model_index >= len(st):
        raise ValueError(
            f"--model {model_index} 가 범위를 벗어났습니다(이 구조의 모델 수={len(st)})."
        )
    model = st[model_index]

    residues: List[Residue] = []
    for chain in model:
        for res in chain:
            # author 잔기번호(seqid.num)와 삽입코드(seqid.icode)를 가져온다.
            seqid = res.seqid
            res_num = int(seqid.num)
            icode = (seqid.icode or " ").strip()
            res_name = res.name.strip().upper()

            atoms: List[Atom] = []
            for atom in res:
                name = atom.name.strip()
                element = atom.element.name if atom.element else ""
                pos = atom.pos
                atoms.append((name, element, float(pos.x), float(pos.y), float(pos.z)))
            residues.append((chain.name, res_num, icode, res_name, atoms))
    return residues


def _load_with_biopython(path: Path, model_index: int) -> List[Residue]:
    """Biopython으로 구조를 읽어 Residue 리스트로 변환한다(fallback, 원자 이름 포함).

    확장자가 .cif/.mmcif면 MMCIFParser, 그 외(.pdb 등)는 PDBParser를 쓴다.
    """
    from Bio.PDB import MMCIFParser, PDBParser  # type: ignore

    suffix = path.suffix.lower()
    if suffix in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    structure = parser.get_structure("s", str(path))
    models = list(structure.get_models())
    if not models:
        raise ValueError("구조에 모델이 하나도 없습니다.")
    if model_index < 0 or model_index >= len(models):
        raise ValueError(
            f"--model {model_index} 가 범위를 벗어났습니다(이 구조의 모델 수={len(models)})."
        )
    model = models[model_index]

    residues: List[Residue] = []
    for chain in model:
        for res in chain:
            # Biopython 잔기 id = (hetflag, resseq, icode).
            #   hetflag: ' '면 표준, 'W'면 물, 'H_...'면 헤테로.
            #   resseq : author 잔기번호.  icode: 삽입코드(없으면 ' ').
            _het_flag, res_num, icode = res.id
            res_name = res.resname.strip().upper()
            icode = (icode or " ").strip()

            atoms: List[Atom] = []
            for atom in res:
                name = atom.get_name().strip()
                # element 속성이 비어 있으면 원자 이름 첫 글자로 보강(완벽하진 않음).
                element = (atom.element or "").strip() or name[:1]
                x, y, z = atom.coord  # numpy 값이지만 float로 캐스팅
                atoms.append((name, element, float(x), float(y), float(z)))
            residues.append((chain.id, int(res_num), icode, res_name, atoms))
    return residues


def load_structure(path: Path, model_index: int = 0) -> Tuple[List[Residue], str]:
    """구조 파일을 읽어 (Residue 리스트, 사용한 백엔드 이름)을 돌려준다.

    gemmi가 설치돼 있으면 gemmi를, 없으면 Biopython을 쓴다. 둘 다 없으면
    설치 안내와 함께 종료한다.
    """
    if not path.exists():
        raise FileNotFoundError(f"입력 파일을 찾을 수 없습니다: {path}")

    # 1순위: gemmi
    try:
        import gemmi  # noqa: F401  (존재 여부만 확인)
        return _load_with_gemmi(path, model_index), "gemmi"
    except ImportError:
        pass

    # 2순위: Biopython
    try:
        import Bio  # noqa: F401  (존재 여부만 확인)
        return _load_with_biopython(path, model_index), "biopython"
    except ImportError:
        pass

    # 둘 다 없음 → 안내 후 종료
    sys.exit(
        "오류: 구조를 읽을 라이브러리가 없습니다.\n"
        "  다음 중 하나를 설치하세요(gemmi 권장):\n"
        "    pip install gemmi\n"
        "    pip install biopython\n"
    )


# ==========================================================================
# 2) 잔기 필터링 / 체인 선택 / 라벨
# ==========================================================================
def filter_residues(
    residues: Sequence[Residue],
    strip_het: bool,
    strip_water_only: bool = False,
) -> List[Residue]:
    """물/헤테로 잔기를 걸러낸다.

    - strip_het=True       : 표준 아미노산이 아닌 잔기를 모두 제거(물·이온·리간드 등).
                             단백질-단백질 계면 분석의 기본값.
    - strip_water_only=True: 물(HOH 등)만 제거하고 다른 헤테로는 남긴다.
                             (리간드가 분석에 필요할 때 쓰는 절충 옵션)
    - 둘 다 False          : 아무것도 거르지 않는다(원본 그대로).
    원자가 0개인 잔기는 좌표가 없으므로 어느 경우든 제거한다.
    """
    out: List[Residue] = []
    for res in residues:
        name = res[3]
        atoms = res[4]
        if not atoms:
            continue  # 좌표 없는 잔기는 의미가 없으니 버린다
        if strip_water_only:
            if name in WATER_NAMES:
                continue
        elif strip_het:
            # 표준 아미노산만 남긴다. (핵산까지 남기려면 이 조건을 넓혀야 함)
            if name not in STANDARD_AA:
                continue
        out.append(res)
    return out


def select_chains(residues: Sequence[Residue], chains: Sequence[str]) -> List[Residue]:
    """주어진 체인 기호들에 속한 잔기만 골라낸다(단순 필터)."""
    wanted = set(chains)
    return [r for r in residues if r[0] in wanted]


def residue_label(res: Residue) -> str:
    """CSV 출력용 잔기 라벨을 만든다. 예: 'A/ASP123' 또는 'B/SER100A'.

    형식 = '{체인}/{잔기이름}{author번호}{삽입코드}'. author numbering 기준이며
    PyMOL select 에 거의 그대로 쓸 수 있다.
    """
    chain, num, icode, name, _atoms = res
    ic = icode if icode and icode != " " else ""
    return f"{chain}/{name}{num}{ic}"


# ==========================================================================
# 3) 거리 계산
# ==========================================================================
def _euclid(a: Point, b: Point) -> float:
    """두 3차원 좌표 사이의 유클리드 거리(Å). 순수 파이썬 계산."""
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def representative_point(atoms: Sequence[Atom], metric: str) -> Optional[Point]:
    """metric이 ca/centroid 일 때, 잔기를 대표하는 점 하나를 계산한다.

    - ca       : 원자 이름이 'CA'인 원자의 좌표(없으면 None → 그 잔기는 건너뜀).
    - centroid : 잔기의 heavy atom 좌표 평균(없으면 전체 원자 평균).
    min-heavy 는 점 하나로 못 줄이므로 여기서 처리하지 않는다(None 반환).
    """
    if metric == "ca":
        for name, _el, x, y, z in atoms:
            if name.strip().upper() == "CA":
                return (x, y, z)
        return None  # Cα가 없는 잔기는 ca 방식에서 제외

    if metric == "centroid":
        heavy = [(x, y, z) for (_n, el, x, y, z) in atoms if _is_heavy(el)]
        pts = heavy if heavy else [(x, y, z) for (_n, _e, x, y, z) in atoms]
        if not pts:
            return None
        n = len(pts)
        sx = sum(p[0] for p in pts)
        sy = sum(p[1] for p in pts)
        sz = sum(p[2] for p in pts)
        return (sx / n, sy / n, sz / n)

    return None  # min-heavy 등은 대표점이 없음


def heavy_coords(atoms: Sequence[Atom]) -> List[Point]:
    """잔기의 heavy atom 좌표 목록을 돌려준다(min-heavy 계산용).

    수소가 모두 빠진 구조(예측 구조 등)에서 heavy가 비는 경우는 거의 없지만,
    혹시 비면 전체 원자로 대체한다.
    """
    heavy = [(x, y, z) for (_n, el, x, y, z) in atoms if _is_heavy(el)]
    if heavy:
        return heavy
    return [(x, y, z) for (_n, _e, x, y, z) in atoms]


def residue_distance(res_a: Residue, res_b: Residue, metric: str) -> Optional[float]:
    """두 잔기 사이 거리를 metric 방식으로 계산한다. 계산 불가 시 None.

    - ca / centroid : 각 잔기의 대표점 사이 유클리드 거리.
    - min-heavy     : 두 잔기 heavy atom 쌍 중 최단 거리.
    """
    atoms_a = res_a[4]
    atoms_b = res_b[4]

    if metric in ("ca", "centroid"):
        pa = representative_point(atoms_a, metric)
        pb = representative_point(atoms_b, metric)
        if pa is None or pb is None:
            return None
        return _euclid(pa, pb)

    if metric == "min-heavy":
        ha = heavy_coords(atoms_a)
        hb = heavy_coords(atoms_b)
        if not ha or not hb:
            return None
        best = float("inf")
        for ca in ha:
            for cb in hb:
                d = _euclid(ca, cb)
                if d < best:
                    best = d
        return best if best != float("inf") else None

    raise ValueError(f"알 수 없는 --metric: {metric}")


# ==========================================================================
# 4) 잔기쌍 거리 계산 + 계면 선별
# ==========================================================================
def compute_pairs(
    res_a_list: Sequence[Residue],
    res_b_list: Sequence[Residue],
    metric: str,
    cutoff: float,
    keep_all: bool,
) -> List[Tuple[Residue, Residue, float]]:
    """체인A 잔기 × 체인B 잔기 모든 조합의 거리를 구한다.

    - keep_all=False(기본): 거리 <= cutoff 인 (접촉) 쌍만 돌려준다(= 계면).
    - keep_all=True       : 모든 쌍을 거리와 함께 돌려준다(필터 없음).
    반환 리스트는 거리 오름차순 정렬.

    안전장치: --chain-a 와 --chain-b 에 같은 체인이 겹쳐 들어오면(예: 둘 다 A)
    "잔기 자기 자신"끼리의 거리 0.0 쌍이 섞여 계면 통계가 왜곡된다. 그래서
    같은 (체인, author번호, 삽입코드, 잔기이름) 인 동일 잔기쌍은 건너뛴다.
    (정상적인 항원 A vs 항체 B처럼 체인이 겹치지 않으면 영향이 전혀 없다.)
    """
    results: List[Tuple[Residue, Residue, float]] = []
    for ra in res_a_list:
        for rb in res_b_list:
            # 같은 잔기끼리(자기 자신) 비교는 거리 0이라 의미가 없으므로 제외.
            # 잔기 식별 = (체인, 번호, 삽입코드, 이름); 좌표(원자목록)는 비교에서 뺀다.
            if ra[0] == rb[0] and ra[1] == rb[1] and ra[2] == rb[2] and ra[3] == rb[3]:
                continue
            d = residue_distance(ra, rb, metric)
            if d is None:
                continue  # 대표점/heavy atom이 없어 계산 불가한 쌍은 건너뜀
            if keep_all or d <= cutoff:
                results.append((ra, rb, d))
    results.sort(key=lambda t: t[2])
    return results


def interface_residue_sets(
    pairs: Sequence[Tuple[Residue, Residue, float]],
) -> Tuple[List[str], List[str]]:
    """접촉 쌍 목록에서 계면 잔기 유니크 라벨을 체인별로 모은다.

    돌려주는 두 리스트 = (체인A 계면 잔기 라벨, 체인B 계면 잔기 라벨).
    중복은 제거하고 보기 좋게 정렬한다.
    """
    set_a: Dict[str, None] = {}
    set_b: Dict[str, None] = {}
    for ra, rb, _d in pairs:
        set_a[residue_label(ra)] = None
        set_b[residue_label(rb)] = None
    return sorted(set_a.keys()), sorted(set_b.keys())


# ==========================================================================
# 5) CSV 저장
# ==========================================================================
def write_pairs_csv(
    path: Path,
    pairs: Sequence[Tuple[Residue, Residue, float]],
    metric: str,
) -> None:
    """잔기쌍 거리 CSV 저장. 컬럼: residue_a, residue_b, distance_A, metric.

    distance_A 단위는 옹스트롬(Å). residue_a/b 는 author numbering 라벨.
    """
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["residue_a", "residue_b", "distance_A", "metric"])
        for ra, rb, d in pairs:
            writer.writerow([residue_label(ra), residue_label(rb), f"{d:.3f}", metric])


def write_interface_csv(
    path: Path,
    set_a: Sequence[str],
    set_b: Sequence[str],
    chain_a_label: str,
    chain_b_label: str,
) -> None:
    """계면 유니크 잔기 목록 CSV 저장. 컬럼: side, chain_group, residue.

    side = 'A'(항원 쪽) / 'B'(항체 쪽). residue 는 author numbering 라벨.
    """
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["side", "chain_group", "residue"])
        for r in set_a:
            writer.writerow(["A", chain_a_label, r])
        for r in set_b:
            writer.writerow(["B", chain_b_label, r])


# ==========================================================================
# 6) CLI
# ==========================================================================
def list_chains(residues: Sequence[Residue]) -> None:
    """구조에 들어 있는 체인 기호와 각 체인의 (필터 전) 잔기 수를 출력한다."""
    counts: Dict[str, int] = {}
    for res in residues:
        counts[res[0]] = counts.get(res[0], 0) + 1
    print("이 구조의 체인 목록 (기호: 잔기수):")
    for name in sorted(counts.keys()):
        print(f"  {name}: {counts[name]}")
    print(
        "\n참고: 우리 관례상 A=항원, B=항체. 항체가 Fab/scFv면 체인이 여러 개일 수 있음.\n"
        "      --chain-b 에 항체 쪽 체인을 정확히 지정하세요(예: --chain-b H L)."
    )


def build_arg_parser() -> argparse.ArgumentParser:
    """argparse 파서 구성. 각 옵션의 한국어 도움말 포함."""
    p = argparse.ArgumentParser(
        prog="cif_distance.py",
        description=(
            "구조 파일(.cif/.pdb)에서 두 체인 사이 잔기 간 거리를 계산하고 "
            "계면(접촉) 잔기를 찾아 CSV로 저장하는 도구. "
            "잔기 번호는 author numbering(구조 파일에 적힌 번호) 기준."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "structure",
        type=Path,
        help="입력 구조 파일 경로(.cif 또는 .pdb). 예: model.cif",
    )
    p.add_argument(
        "--chain-a",
        nargs="+",
        default=["A"],
        metavar="CHAIN",
        help="체인 A(기본 A = 항원) 기호. 여러 개면 공백으로 나열. 기본값: A",
    )
    p.add_argument(
        "--chain-b",
        nargs="+",
        default=["B"],
        metavar="CHAIN",
        help="체인 B(기본 B = 항체) 기호. Fab/scFv면 여러 개일 수 있음(예: H L). 기본값: B",
    )
    p.add_argument(
        "--metric",
        choices=["ca", "centroid", "min-heavy"],
        default="min-heavy",
        help=(
            "거리 계산 방식. ca=Cα-Cα, centroid=잔기 원자중심, "
            "min-heavy=최단 heavy-atom 거리(계면엔 보통 이게 적합). 기본: min-heavy"
        ),
    )
    p.add_argument(
        "--cutoff",
        type=float,
        default=5.0,
        help="접촉(계면)으로 볼 거리 한계값(Å). 이 값 이하인 쌍만 계면으로 본다. 기본: 5.0",
    )
    p.add_argument(
        "-o", "--out",
        type=Path,
        default=None,
        help="잔기쌍 거리 CSV 저장 경로. 생략하면 화면에 요약만 출력.",
    )
    p.add_argument(
        "--interface-out",
        type=Path,
        default=None,
        help="계면 잔기 유니크 목록 CSV 저장 경로(선택).",
    )
    p.add_argument(
        "--all-pairs",
        action="store_true",
        help="cutoff로 거르지 않고 모든 잔기쌍 거리를 저장(파일이 매우 커질 수 있음).",
    )
    p.add_argument(
        "--no-strip-het",
        action="store_true",
        help="물/헤테로(비표준 잔기)를 제거하지 않고 그대로 사용. 기본은 제거함.",
    )
    p.add_argument(
        "--strip-water-only",
        action="store_true",
        help="물(HOH 등)만 제거하고 그 외 헤테로(리간드 등)는 남긴다.",
    )
    p.add_argument(
        "--model",
        type=int,
        default=0,
        help="여러 모델이 든 파일에서 사용할 모델 번호(0부터). 기본: 0",
    )
    p.add_argument(
        "--list-chains",
        action="store_true",
        help="구조의 체인 기호와 잔기 수만 출력하고 종료(분석은 하지 않음).",
    )
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI 진입점. 정상 종료 0, 오류 시 0이 아닌 값."""
    args = build_arg_parser().parse_args(argv)

    # 1) 구조 로드. 실패 시 메시지 출력 후 종료.
    try:
        residues, backend = load_structure(args.structure, args.model)
    except (FileNotFoundError, ValueError) as exc:
        print(f"오류: {exc}", file=sys.stderr)
        return 2
    print(f"[로드] {args.structure}  (백엔드: {backend}, 잔기 {len(residues)}개)")

    # 2) --list-chains 면 체인만 출력하고 끝낸다.
    if args.list_chains:
        list_chains(residues)
        return 0

    # 3) 물/헤테로 필터링.
    if args.no_strip_het:
        filtered = filter_residues(residues, strip_het=False, strip_water_only=False)
    elif args.strip_water_only:
        filtered = filter_residues(residues, strip_het=False, strip_water_only=True)
    else:
        filtered = filter_residues(residues, strip_het=True, strip_water_only=False)

    # 4) 체인 분리.
    res_a = select_chains(filtered, args.chain_a)
    res_b = select_chains(filtered, args.chain_b)

    a_label = "+".join(args.chain_a)
    b_label = "+".join(args.chain_b)
    print(f"[체인] A({a_label}) 잔기 {len(res_a)}개  /  B({b_label}) 잔기 {len(res_b)}개")

    # 어느 한쪽 체인이 비면(기호 오타 등) 안내하고 종료.
    if not res_a or not res_b:
        sys.stdout.flush()  # 위의 [로드]/[체인] stdout 줄을 먼저 비운 뒤 stderr 출력
        print(
            "오류: 한쪽 체인에서 잔기를 찾지 못했습니다. 체인 기호를 확인하세요.\n"
            "      --list-chains 로 이 구조의 실제 체인 기호를 먼저 확인하세요.",
            file=sys.stderr,
        )
        return 3

    # 5) 잔기쌍 거리 계산 + 계면 선별.
    pairs = compute_pairs(
        res_a, res_b, metric=args.metric, cutoff=args.cutoff, keep_all=args.all_pairs
    )

    if args.all_pairs:
        # all-pairs 모드에서도 "계면" 통계는 cutoff 기준으로 따로 센다.
        contact_pairs = [t for t in pairs if t[2] <= args.cutoff]
        print(
            f"[결과] 전체 잔기쌍 {len(pairs)}개 "
            f"(그중 {args.cutoff} Å 이내 {len(contact_pairs)}개, metric={args.metric})"
        )
    else:
        contact_pairs = list(pairs)
        print(f"[결과] {args.cutoff} Å 이내 접촉 잔기쌍 {len(pairs)}개 (metric={args.metric})")

    # 6) 계면 유니크 잔기 집합.
    set_a, set_b = interface_residue_sets(contact_pairs)
    print(f"[계면] 항원 측 {len(set_a)}잔기 / 항체 측 {len(set_b)}잔기 (cutoff={args.cutoff} Å)")

    # 7) CSV 저장 또는 화면 미리보기.
    if args.out is not None:
        write_pairs_csv(args.out, pairs, args.metric)
        print(f"[저장] 잔기쌍 거리 CSV → {args.out}")
    else:
        # 출력 경로가 없으면 가까운 순 상위 몇 개만 미리 보여준다.
        print("\n(미리보기: 가까운 잔기쌍 최대 10개)")
        print(f"  {'residue_a':<20} {'residue_b':<20} {'distance_A':>10}")
        for ra, rb, d in contact_pairs[:10]:
            print(f"  {residue_label(ra):<20} {residue_label(rb):<20} {d:>10.3f}")

    if args.interface_out is not None:
        write_interface_csv(args.interface_out, set_a, set_b, a_label, b_label)
        print(f"[저장] 계면 잔기 목록 CSV → {args.interface_out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
