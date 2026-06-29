#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""AF3 출력 .cif(mmCIF) 구조를 .pdb로 변환하고, 신뢰도 JSON에서 핵심 스칼라만 뽑아 보는 CLI.

배경 (왜 이 스크립트가 필요한가)
--------------------------------
AlphaFold3(이하 AF3)가 항체-항원 도킹을 돌리면 결과 구조가 `.cif`(=mmCIF) 포맷으로 나온다.
최신 뷰어(PyMOL·ChimeraX)는 `.cif`를 잘 읽지만, **오래된 도구나 일부 파이썬 라이브러리·
도킹 프로그램은 `.pdb`만 받거나 `.pdb`에서 더 안정적**이다. 그래서 후속 분석 도구가
`.pdb`를 요구할 때 변환이 필요하다. 이 스크립트는 그 변환을 한 줄 명령으로 해 준다.

추가로, AF3가 함께 내놓는 신뢰도 JSON(요약/상세)에서 `iptm`, `ptm`, `ranking_score`
같은 **값 하나짜리 숫자(스칼라)** 만 빠르게 들여다보는 보조 기능(--summary)도 넣었다.
(상세 행렬 피처 추출은 이 스크립트의 일이 아니다 — 그건 af3_extract_features.py 담당.)

⚠️ .pdb 포맷의 한계 (변환 전에 꼭 알아둘 것)
--------------------------------------------
- `.pdb`는 옛 포맷이라 **체인(chain) 기호가 1글자**다(여기서 체인 = 단백질 사슬 하나,
  예: 항원 사슬, 항체 중쇄, 항체 경쇄). 우리 항원-항체 복합체는 체인이 보통 2~3개라
  A/B/C 한 글자로 충분해 대개 괜찮다. 하지만 체인이 아주 많아 2글자 기호(예: AA, AB)가
  필요한 거대 복합체는 변환 시 깨지거나 잘릴 수 있다.
  - **백엔드별 차이(중요):**
    * gemmi 경로: 2글자 체인이 있어도 보통 에러 없이 변환은 되지만, 1글자로 강제
      재명명되거나 잘릴 수 있어 결과가 원본과 달라질 수 있다(반드시 PyMOL로 확인).
    * Biopython 경로: 2글자 체인이 있으면 저장 단계에서 보통 에러가 난다
      (TypeError: %c requires int or char). 즉 다중문자 체인은 Biopython으로는
      .pdb 저장 자체가 실패할 수 있다.
- 원자/잔기 수가 매우 큰 구조도 컬럼 폭 한계로 문제가 날 수 있다.
- 그래서 변환 후에는 PyMOL로 한 번 열어 정상인지 확인하는 습관을 들이자.

⚠️ 멀티모델 .cif / 큰 파일 주의
-------------------------------
- AF3가 모델당 내놓는 *_model.cif는 보통 모델(구조) 1개만 담는다. 하지만 여러 구조를
  한 파일에 담은 "멀티모델 .cif"가 들어올 수도 있다(예: 외부에서 받은 앙상블).
  - gemmi 경로: gemmi의 write_pdb는 **기본적으로 구조 안의 모든 모델을** PDB로 쓴다
    (각 모델을 MODEL/ENDMDL로 감싼다). 이 스크립트는 기본값에서 **첫 번째 모델만**
    남기고 저장한다(후속 도킹 도구가 단일 모델을 기대하는 경우가 많기 때문).
    --all-models를 주면 모델을 지우지 않고 gemmi가 전부 쓰게 둔다.
  - Biopython 경로(폴백): MMCIFParser는 멀티모델을 읽지만, 이 스크립트는 안전하게
    **항상 첫 모델만** 저장한다(Select로 첫 모델만 골라 저장). --all-models는
    Biopython 경로에서는 무시된다(아래 안내 메시지 출력).
- 파일이 아주 크면 메모리/시간이 든다. --recursive로 폴더를 통째 변환할 때는 개수가
  많을 수 있으니(복합체당 75개 모델!) 먼저 작은 폴더로 시험해 보길 권한다.

의존성 (둘 중 하나만 있으면 됨)
-------------------------------
- 1순위: **gemmi** (가볍고 빠름).  설치:  pip install gemmi
    내부적으로 gemmi.read_structure(...) → structure.write_pdb(...) 를 쓴다.
- 폴백: **Biopython** (gemmi가 없을 때).  설치:  pip install biopython
    내부적으로 Bio.PDB.MMCIFParser → Bio.PDB.PDBIO 를 쓴다.
- --summary(신뢰도 JSON 요약)와 파일 탐색은 **표준 라이브러리만** 쓰므로,
  요약 기능만 쓸 거면 위 외부 패키지 없이도 동작한다.
- 둘 다 없으면 변환 시 친절한 안내 메시지와 함께 종료한다.

테스트 상태
-----------
- 이 파일은 실제 AF3 출력 데이터로 직접 검증하지 못했다(이 환경에 gemmi/Biopython이
  설치돼 있지 않아 변환 결과를 눈으로 확인하지 못함). gemmi/Biopython API 사용법은
  공식 문서를 기준으로 맞췄으나, 버전에 따라 세부 동작이 다를 수 있다. 처음에는
  **샘플 .cif 한두 개로 변환 → PyMOL로 열어 확인**한 뒤 --recursive 일괄 변환에
  쓰는 것을 권한다.
  (참고: --summary/에러처리/파일탐색 등 표준라이브러리 경로는 작성자가 픽스처로
  실행 테스트해 동작을 확인했다. 변환 백엔드 호출만 미실행 검증이다.)

사용 예시
---------
    # 단일 파일 변환: 같은 폴더에 같은 이름의 .pdb 생성
    python scripts/af3_convert.py path/to/<job>_model.cif

    # 출력 경로를 직접 지정
    python scripts/af3_convert.py path/to/<job>_model.cif -o path/to/out.pdb

    # 폴더 통째로(하위까지 재귀) 모든 .cif 변환 (원본 .cif 옆에 .pdb 생성)
    python scripts/af3_convert.py path/to/af3_outputs/ --recursive

    # 재귀 변환하되 결과 .pdb를 한 폴더에 모으기
    python scripts/af3_convert.py path/to/af3_outputs/ --recursive --out-dir path/to/pdbs/

    # 이미 .pdb가 있으면 건너뛰기(다시 안 만듦)
    python scripts/af3_convert.py path/to/af3_outputs/ --recursive --skip-existing

    # 신뢰도 JSON에서 핵심 스칼라만 빠르게 출력 (변환 안 함)
    python scripts/af3_convert.py path/to/<job>_summary_confidences.json --summary
    python scripts/af3_convert.py path/to/af3_outputs/ --recursive --summary
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# --summary로 들여다볼 "핵심 스칼라" 키 후보들.
# AF3 버전/출처(로컬 DeepMind vs AlphaFold Server)에 따라 일부 키는 없을 수 있으므로,
# 존재하는 것만 골라 출력한다(없으면 조용히 건너뜀).
SUMMARY_SCALAR_KEYS: Tuple[str, ...] = (
    "iptm",            # interface pTM — 체인 사이(계면) 정렬 신뢰도 (우리 baseline)
    "ptm",             # predicted TM-score — 구조 전체 정렬 신뢰도
    "ranking_score",   # AF3가 모델을 줄 세울 때 쓰는 종합 점수
    "fraction_disordered",  # 무질서(흐트러진) 영역 비율
    "has_clash",       # 충돌 여부 플래그 (불리언/숫자)
    "num_recycles",    # 리사이클 횟수 (참고용; 버전/출처에 따라 없을 수 있음)
)


# ---------------------------------------------------------------------------
# 1) 변환 백엔드 (gemmi 우선, Biopython 폴백)
# ---------------------------------------------------------------------------

def _convert_with_gemmi(cif_path: Path, pdb_path: Path, all_models: bool) -> None:
    """gemmi로 .cif -> .pdb 변환. gemmi가 없으면 ImportError를 그대로 올린다.

    매개변수
    --------
    cif_path : 입력 mmCIF 경로
    pdb_path : 출력 PDB 경로
    all_models : True면 멀티모델을 그대로 보존(gemmi write_pdb 기본 동작 = 전 모델).
                 False(기본)면 첫 모델만 남기고 저장한다.

    참고: gemmi의 structure.write_pdb는 **기본적으로 구조 안 모든 모델을** 쓴다.
    그래서 "첫 모델만" 원하면 저장 전에 나머지 모델을 지워야 한다.
    """
    import gemmi  # 여기서 import — 없으면 ImportError가 위로 전달되어 폴백을 타게 한다.

    structure = gemmi.read_structure(str(cif_path))

    # 멀티모델 처리: 기본은 첫 모델만 남긴다(.pdb 1글자 체인/원자 폭 한계 + 후속 도구가
    # 단일 모델을 기대하는 경우가 많기 때문). gemmi Structure는 모델들의 시퀀스처럼
    # 동작하며 슬라이스 삭제(del st[1:])를 지원한다.
    if not all_models and len(structure) > 1:
        del structure[1:]  # 첫 모델(index 0)만 남기고 나머지 삭제

    # gemmi 권장: 파일을 읽은 뒤 setup_entities()로 엔티티/서브체인 정보를 정리하면
    # PDB 저장 시 TER 레코드/체인 그룹이 제대로 들어간다("거의 항상 필요"하다고 문서가
    # 안내). 버전에 따라 없거나 시그니처가 다를 수 있어 try로 감싼다.
    try:
        structure.setup_entities()
    except Exception:  # noqa: BLE001 - 정리 단계 실패는 치명적이지 않으니 그대로 저장 시도.
        pass

    structure.write_pdb(str(pdb_path))


def _convert_with_biopython(cif_path: Path, pdb_path: Path) -> None:
    """Biopython으로 .cif -> .pdb 변환(폴백). Biopython이 없으면 ImportError를 올린다.

    멀티모델 안전성: PDBIO의 멀티모델 저장 동작이 버전에 따라 다를 수 있어, 여기서는
    Select로 **첫 번째 모델만** 명시적으로 골라 저장한다(동작을 예측 가능하게).

    다중문자 체인 주의: 입력에 2글자 이상 체인 기호가 있으면 PDBIO 저장 단계에서
    TypeError(%c requires int or char)가 날 수 있다. 이 경우 gemmi 백엔드를 쓰거나
    체인을 1글자로 재명명해야 한다.
    """
    from Bio.PDB import MMCIFParser, PDBIO, Select  # 없으면 ImportError -> 호출부에서 안내

    # QUIET=True: 비표준 레코드 등에 대한 경고를 줄여 초보에게 덜 시끄럽게 한다.
    parser = MMCIFParser(QUIET=True)
    # structure id는 임의 문자열이면 된다(파일 stem 사용).
    structure = parser.get_structure(cif_path.stem, str(cif_path))

    # 구조에 모델이 여러 개일 수 있으니, 첫 모델만 저장하도록 Select를 만든다.
    # (model.id는 보통 0부터 시작. 첫 모델의 실제 id를 읽어 그것만 통과시킨다.)
    try:
        first_model_id = next(iter(structure)).id
    except StopIteration:
        raise ValueError(f"구조에 모델이 하나도 없습니다(빈 구조?): {cif_path}")

    class _FirstModelOnly(Select):
        def accept_model(self, model) -> bool:  # type: ignore[override]
            return model.id == first_model_id

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(pdb_path), select=_FirstModelOnly())


def convert_one(
    cif_path: Path,
    pdb_path: Path,
    backend: str,
    all_models: bool,
) -> str:
    """.cif 한 개를 .pdb로 변환한다.

    매개변수
    --------
    backend : "auto" | "gemmi" | "biopython"
        - "auto": gemmi를 먼저 시도하고, 없으면(ImportError) Biopython으로 폴백.
        - "gemmi"/"biopython": 해당 백엔드만 강제 사용(없으면 에러).
    all_models : gemmi 경로에서만 의미. True면 멀티모델 보존, False면 첫 모델만.
                 (biopython 경로는 항상 첫 모델만 — all_models는 무시됨)

    반환값
    ------
    실제로 사용한 백엔드 이름("gemmi" 또는 "biopython"). 변환 자체가 실패하면 예외를 올린다.
    """
    if backend == "gemmi":
        _convert_with_gemmi(cif_path, pdb_path, all_models)
        return "gemmi"

    if backend == "biopython":
        _convert_with_biopython(cif_path, pdb_path)
        return "biopython"

    # backend == "auto": gemmi 먼저, 없으면 Biopython.
    try:
        _convert_with_gemmi(cif_path, pdb_path, all_models)
        return "gemmi"
    except ImportError:
        # gemmi 미설치 → Biopython 폴백 시도.
        try:
            _convert_with_biopython(cif_path, pdb_path)
            return "biopython"
        except ImportError:
            # 둘 다 없음 → 사용자에게 설치 안내를 담은 명확한 에러를 올린다.
            raise ImportError(
                "변환에 필요한 라이브러리가 없습니다. 다음 중 하나를 설치하세요:\n"
                "    pip install gemmi        # (권장: 가볍고 빠름)\n"
                "    pip install biopython    # (폴백)"
            )


# ---------------------------------------------------------------------------
# 2) 파일 찾기
# ---------------------------------------------------------------------------

def find_cif_files(base: Path, recursive: bool) -> List[Path]:
    """변환 대상 .cif 파일 목록을 모은다.

    - base가 파일이면 그 파일 하나만(확장자가 .cif가 아니어도 사용자가 직접 지정했으니 존중).
    - base가 폴더면: recursive=True면 하위까지(rglob), 아니면 그 폴더 바로 아래만(glob).
    - 대소문자 구분 없이 .cif를 찾는다(.CIF도 포함).
    """
    if base.is_file():
        return [base]

    pattern_glob = base.rglob if recursive else base.glob
    found = set()
    for pat in ("*.cif", "*.CIF"):
        for p in pattern_glob(pat):
            if p.is_file():
                found.add(p)
    # 정렬해 출력 순서를 안정적으로(보기 좋게) 만든다.
    return sorted(found)


def find_summary_jsons(base: Path, recursive: bool) -> List[Path]:
    """--summary 대상 JSON 목록을 모은다(요약/상세 신뢰도 파일).

    파일명 패턴은 환경에 따라 다르므로 넉넉히 잡는다:
    *summary_confidences*.json, *_confidences.json, *full_data*.json.
    (요약 파일이 우선 의미 있지만, 상세 파일에도 같은 스칼라가 들어 있을 수 있다.)
    set으로 모으므로 패턴이 겹쳐도 같은 파일은 한 번만 들어간다.
    """
    if base.is_file():
        return [base]

    pattern_glob = base.rglob if recursive else base.glob
    found = set()
    for pat in (
        "*summary_confidences*.json",
        "*_confidences.json",
        "*full_data*.json",
    ):
        for p in pattern_glob(pat):
            if p.is_file():
                found.add(p)
    return sorted(found)


def default_pdb_path(cif_path: Path, out_dir: Optional[Path]) -> Path:
    """입력 .cif에 대응하는 출력 .pdb 경로를 정한다.

    - out_dir가 없으면: 원본 .cif 옆에 같은 이름 + .pdb.
    - out_dir가 있으면: 그 폴더 안에 같은 파일이름(stem) + .pdb.
      (재귀 변환 시 파일이름이 겹치면 나중 것이 덮어쓸 수 있다 — 아래 main에서 경고.)
    """
    if out_dir is None:
        return cif_path.with_suffix(".pdb")
    return out_dir / (cif_path.stem + ".pdb")


# ---------------------------------------------------------------------------
# 3) --summary: 신뢰도 JSON에서 스칼라만 뽑기
# ---------------------------------------------------------------------------

def extract_summary_scalars(json_path: Path) -> Dict[str, object]:
    """신뢰도 JSON에서 핵심 스칼라(SUMMARY_SCALAR_KEYS)만 골라 dict로 돌려준다.

    - 행렬/배열(예: pae, contact_probs, chain_pair_iptm)은 의도적으로 건너뛴다
      (여기 목적은 "한눈에 보는 요약"이고, 큰 배열은 화면을 더럽히기 때문).
    - 키가 없거나 값이 스칼라가 아니면(딕셔너리/리스트면) 제외한다.
    - JSON이 깨졌거나 읽기 실패면 빈 dict를 돌려주고, 호출부가 알아서 경고한다.
    """
    try:
        with json_path.open("r", encoding="utf-8") as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError, UnicodeDecodeError):
        return {}

    if not isinstance(data, dict):
        # 최상위가 dict가 아니면 우리가 찾는 스칼라 구조가 아니다.
        return {}

    result: Dict[str, object] = {}
    for key in SUMMARY_SCALAR_KEYS:
        if key not in data:
            continue
        value = data[key]
        # 스칼라(불리언/숫자/문자열/None)만 채택. list/dict는 요약에서 제외.
        # bool은 int의 하위형이므로 (bool, int, float) 순서/포함 모두 안전하다.
        if isinstance(value, (bool, int, float, str)) or value is None:
            result[key] = value
    return result


def print_summary(json_paths: List[Path]) -> int:
    """--summary 동작: 각 JSON의 핵심 스칼라를 한 줄씩 보기 좋게 출력한다.

    반환값: 종료 코드(0=정상, 1=하나라도 스칼라 추출 실패).
    """
    if not json_paths:
        print("요약할 신뢰도 JSON을 찾지 못했습니다. (패턴: *summary_confidences*.json 등)")
        return 0

    any_failed = False
    for jp in json_paths:
        scalars = extract_summary_scalars(jp)
        print(f"\n[{jp}]")
        if not scalars:
            print("  (핵심 스칼라를 못 찾음 — 파일이 깨졌거나, JSON이 아니거나, 키가 다를 수 있음)")
            any_failed = True
            continue
        # SUMMARY_SCALAR_KEYS 순서대로 출력(중요한 iptm/ptm/ranking_score가 위로 오도록).
        for key in SUMMARY_SCALAR_KEYS:
            if key in scalars:
                print(f"  {key:>20} : {scalars[key]}")
    print()
    return 1 if any_failed else 0


# ---------------------------------------------------------------------------
# 4) CLI 진입점
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """argparse 파서를 만든다(테스트/도움말에서 재사용하려고 분리)."""
    parser = argparse.ArgumentParser(
        description=(
            "AF3 출력 .cif(mmCIF)를 .pdb로 변환한다. --recursive로 폴더 일괄 변환, "
            "--summary로 신뢰도 JSON의 핵심 스칼라(iptm/ptm/ranking_score 등)만 빠르게 본다. "
            "변환 백엔드는 gemmi(권장) 우선, 없으면 Biopython 폴백."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "path",
        type=str,
        help="입력 .cif 파일, 또는 .cif/JSON이 들어 있는 폴더(--recursive와 함께).",
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        metavar="OUT.pdb",
        help="단일 파일 변환 시 출력 .pdb 경로. (--recursive와 함께 쓰지 말 것)",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default=None,
        metavar="DIR",
        help="변환 결과 .pdb를 모아 둘 폴더. 미지정 시 원본 .cif 옆에 생성.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="path가 폴더일 때 하위까지 재귀적으로 모든 .cif를 일괄 변환.",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="대상 .pdb가 이미 있으면 건너뛴다(다시 만들지 않음).",
    )
    parser.add_argument(
        "--backend",
        choices=("auto", "gemmi", "biopython"),
        default="auto",
        help="변환 백엔드 선택. auto=gemmi 먼저, 없으면 biopython(기본값).",
    )
    parser.add_argument(
        "--all-models",
        action="store_true",
        help=(
            "멀티모델 .cif의 모든 모델을 보존(gemmi 경로 한정; write_pdb 기본 동작). "
            "기본은 첫 모델만 남긴다. biopython 경로에서는 무시되고 항상 첫 모델만 저장."
        ),
    )
    parser.add_argument(
        "--summary",
        action="store_true",
        help=(
            "변환 대신, 신뢰도 JSON(요약/상세)에서 핵심 스칼라만 출력한다. "
            "path가 폴더면 --recursive로 하위까지 훑는다."
        ),
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """CLI 진입점. 반환값은 프로세스 종료 코드(0=정상)."""
    args = build_parser().parse_args(argv)

    base = Path(args.path)

    # --- 공통 입력 검증 ---
    if not base.exists():
        print(f"[오류] 경로가 존재하지 않습니다: {base}", file=sys.stderr)
        return 2

    out_dir: Optional[Path] = Path(args.out_dir) if args.out_dir else None

    # === 모드 1: --summary (변환 안 함, JSON 스칼라만 출력) ===
    if args.summary:
        if base.is_dir():
            json_paths = find_summary_jsons(base, recursive=args.recursive)
        else:
            json_paths = [base]
        return print_summary(json_paths)

    # === 모드 2: 변환 ===
    # 옵션 충돌/주의 안내.
    if args.output and (base.is_dir() or args.recursive):
        print(
            "[오류] -o/--output 은 '단일 파일' 변환에만 씁니다. "
            "폴더/재귀 변환에는 --out-dir 을 쓰세요.",
            file=sys.stderr,
        )
        return 2
    if args.output and args.out_dir:
        print(
            "[오류] -o/--output 과 --out-dir 은 함께 쓸 수 없습니다. 하나만 쓰세요.",
            file=sys.stderr,
        )
        return 2

    # biopython 강제 + --all-models는 의미가 없으므로 미리 알려준다.
    if args.all_models and args.backend == "biopython":
        print(
            "[안내] --all-models 는 gemmi 경로에서만 적용됩니다. "
            "biopython 경로에서는 항상 첫 모델만 저장합니다."
        )

    cif_files = find_cif_files(base, recursive=args.recursive)
    if not cif_files:
        if base.is_dir() and not args.recursive:
            print(
                "변환할 .cif가 이 폴더 바로 아래에 없습니다. "
                "하위 폴더까지 보려면 --recursive 를 붙이세요."
            )
        else:
            print("변환할 .cif 파일을 찾지 못했습니다.")
        return 0

    # 단일 파일을 직접 지정했는데 .cif가 아니면(실수로 .json 등) 가볍게 경고만 한다.
    # (사용자가 일부러 비표준 확장자를 줄 수도 있으니 막지는 않는다.)
    if base.is_file() and base.suffix.lower() not in (".cif",):
        print(
            f"[안내] 확장자가 .cif가 아닙니다({base.suffix}). mmCIF가 맞는지 확인하세요. "
            f"신뢰도 JSON을 보려면 --summary 를 쓰세요."
        )

    # 출력 폴더가 지정됐으면 미리 만들어 둔다.
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)

    # 재귀 + out_dir일 때 파일이름(stem)이 겹치면 덮어쓸 위험이 있어 미리 알린다.
    if out_dir is not None and len(cif_files) > 1:
        stems = [c.stem for c in cif_files]
        if len(set(stems)) != len(stems):
            print(
                "⚠️ 경고: --out-dir 로 모을 때 같은 파일이름(stem)이 여러 개라 "
                "나중 변환이 앞 것을 덮어쓸 수 있습니다. 원본 옆에 두려면 --out-dir 을 빼세요."
            )

    print(f"[정보] 변환 대상 .cif: {len(cif_files)}개")
    print(f"[정보] 백엔드: {args.backend} (auto=gemmi 우선, 없으면 biopython)")
    print()

    converted = 0
    skipped = 0
    failed = 0

    for cif in cif_files:
        # 단일 파일 + -o 지정이면 그 경로, 아니면 규칙대로 정한다.
        if args.output and not args.recursive and base.is_file():
            pdb_path = Path(args.output)
            # -o로 폴더 경로를 줬을 때를 대비: 그 폴더 안에 stem.pdb 로 만든다.
            if pdb_path.is_dir():
                pdb_path = pdb_path / (cif.stem + ".pdb")
        else:
            pdb_path = default_pdb_path(cif, out_dir)

        if args.skip_existing and pdb_path.exists():
            print(f"  [건너뜀] 이미 있음: {pdb_path}")
            skipped += 1
            continue

        # 입력과 출력이 같은 경로면(예: .cif에 -o로 같은 .cif를 준 실수) 덮어쓰지 않는다.
        try:
            same = pdb_path.resolve() == cif.resolve()
        except OSError:
            same = False
        if same:
            print(f"  [실패] 입력과 출력 경로가 같습니다: {cif}", file=sys.stderr)
            failed += 1
            continue

        # 출력 부모 폴더가 없으면 만든다(단일 -o로 새 경로를 준 경우 등).
        pdb_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            used = convert_one(cif, pdb_path, backend=args.backend, all_models=args.all_models)
            print(f"  [변환] {cif}  ->  {pdb_path}   ({used})")
            converted += 1
        except ImportError as exc:
            # 백엔드 자체가 없으면 더 진행해도 전부 실패하므로 즉시 중단한다.
            print(f"\n[오류] {exc}", file=sys.stderr)
            return 3
        except Exception as exc:  # noqa: BLE001 - 변환 실패는 파일마다 다양할 수 있어 폭넓게 잡는다.
            print(f"  [실패] {cif} ({type(exc).__name__}: {exc})", file=sys.stderr)
            # 부분적으로 생성됐을 수 있는 깨진 .pdb는 지워 둔다(혼란 방지).
            try:
                if pdb_path.exists():
                    pdb_path.unlink()
            except OSError:
                pass
            failed += 1

    print()
    print(f"완료: 변환 {converted}개, 건너뜀 {skipped}개, 실패 {failed}개.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
