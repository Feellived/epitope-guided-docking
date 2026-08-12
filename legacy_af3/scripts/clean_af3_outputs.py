#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""AF3 출력 폴더에서 msas / template 등 보조 폴더를 재귀적으로 삭제하는 CLI.

배경 (왜 이 스크립트가 필요한가)
--------------------------------
AlphaFold3(AF3)로 항체-항원 도킹을 돌리면, 각 작업 폴더 안에는 우리가 분석에
쓰는 결과(구조 .cif, *_summary_confidences_#.json, *_full_data_#.json)뿐 아니라
"msas/"(다중서열정렬 입력)나 "template/"(구조 템플릿 입력) 같은 *중간 입력*
폴더도 함께 생긴다. 이 폴더들은 보통 후속 분석에 필요 없는데,
- 파일 수가 매우 많고 용량을 차지하며,
- 결과 파일을 재귀 탐색(os.walk 등)할 때 자꾸 걸려서 스크립트를 느리게 한다.
그래서 분석 전에 이런 보조 폴더만 골라 지우면 편하다.

이 스크립트는 base-dir(기준 폴더) 아래를 재귀적으로 훑어서, 이름이
msas / msa / template / templates 인 디렉터리를 찾아 삭제한다.

⚠️ 위험성 경고 (반드시 읽을 것)
--------------------------------
- 이 스크립트는 폴더를 **통째로(shutil.rmtree)** 지운다. 한 번 지우면 복구할 수
  없다. base-dir를 잘못 지정하면 엉뚱한 폴더가 통째로 날아갈 수 있다.
- 그래서 **기본 동작은 DRY-RUN(예행 연습)** 이다. 즉 "무엇을 지울지"만 화면에
  보여주고 실제로는 아무것도 지우지 않는다.
- 진짜로 지우려면 반드시 `--apply` 플래그를 붙여야 한다. 먼저 --apply 없이
  실행해서 목록을 꼭 눈으로 확인한 다음, 맞으면 --apply로 다시 실행하라.

대상 범위 (무엇을 지우고 무엇을 안 지우나)
------------------------------------------
- 이 도구는 **디렉터리(폴더)만** 대상으로 한다. 이름이 일치하더라도 *파일*이나
  *파일을 가리키는 링크*는 건드리지 않는다(폴더를 정리하는 도구이기 때문).
- 폴더를 가리키는 심볼릭 링크는 대상이 될 수 있으나, 이 경우 링크가 가리키는
  실제 폴더는 절대 건드리지 않고 **링크 자체만 unlink** 한다.

안전장치
--------
- 삭제 대상이 '놓인 위치'가 base-dir 바깥이면 거부한다(경로 탈출 방지).
  단, 링크는 그 위치(링크가 사는 자리)로 판단한다 — 링크가 base-dir 안에 있으면
  바깥을 가리키더라도 링크만 안전하게 제거하고, 가리키는 외부 대상은 손대지 않는다.
- 심볼릭 링크(바로가기)는 따라 들어가지 않는다(os.walk followlinks=False).
  폴더를 가리키는 링크가 대상 이름이면 내용을 지우지 않고 링크 자체만 unlink 한다.
- 실제 삭제 전에 총 개수와 대략적인 용량 합계를 요약해서 보여준다.

의존성
------
- 표준 라이브러리만 사용한다(argparse, os, sys, shutil, pathlib). 외부 패키지 불필요.

테스트 상태
-----------
- 이 파일은 실제 AF3 출력 데이터로 직접 검증하지 못했다(임시 폴더로 동작·안전장치만
  확인). 위험한 작업이므로, 반드시 먼저 DRY-RUN(기본값)으로 목록을 확인하고,
  가능하면 사본·테스트 폴더에서 --apply 동작을 한 번 시험해 본 뒤 실제 데이터에
  사용하라.

사용 예시
---------
    # 1) 먼저 예행 연습: 무엇을 지울지 목록만 본다 (실제 삭제 안 함)
    python clean_af3_outputs.py /path/to/af3_outputs

    # 2) 목록이 맞으면 실제로 삭제
    python clean_af3_outputs.py /path/to/af3_outputs --apply

    # 3) 지울 폴더 이름을 직접 지정 (기본값 대신)
    python clean_af3_outputs.py /path/to/af3_outputs --names msas template
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path
from typing import List, Tuple

# 기본으로 삭제 대상으로 삼을 폴더 이름들.
# (대소문자 구분 없이 비교한다 — 아래 _matches_target 참고)
DEFAULT_TARGET_NAMES: Tuple[str, ...] = ("msas", "msa", "template", "templates")


def _human_size(num_bytes: int) -> str:
    """바이트 수를 사람이 읽기 쉬운 문자열(KB/MB/GB ...)로 바꿔 준다.

    예: 1536 -> "1.5 KB". 화면 출력용일 뿐 계산에는 영향이 없다.
    """
    size = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024.0 or unit == "TB":
            # 바이트 단위는 소수점 없이, 그 외에는 소수점 한 자리로 표시
            if unit == "B":
                return f"{int(size)} {unit}"
            return f"{size:.1f} {unit}"
        size /= 1024.0
    # 위 루프에서 항상 return 되므로 여기까지 오지 않는다(타입체커 안심용)
    return f"{num_bytes} B"


def _dir_size(path: Path) -> int:
    """폴더 안 모든 파일의 용량 합계를 바이트로 계산한다.

    심볼릭 링크는 따라가지 않고(followlinks=False), 링크 파일 자체도 용량 합산에서
    제외한다(엉뚱한 외부 대상의 크기를 더하지 않기 위해서다).
    접근 불가 등 오류가 나는 파일은 조용히 건너뛴다(용량은 어디까지나 참고용).
    """
    total = 0
    for root, _dirs, files in os.walk(path, followlinks=False):
        for name in files:
            fp = Path(root) / name
            try:
                # 링크는 건너뛰고, 실제 파일 크기만 더한다
                if fp.is_symlink():
                    continue
                total += fp.stat().st_size
            except OSError:
                # 권한 문제 등으로 못 읽으면 그냥 무시(용량은 추정치)
                continue
    return total


def _matches_target(name: str, target_names: Tuple[str, ...]) -> bool:
    """폴더 이름이 삭제 대상 이름 목록에 (대소문자 무시하고) 들어있는지 검사."""
    return name.lower() in {t.lower() for t in target_names}


def _is_within(child: Path, parent: Path) -> bool:
    """child가 '놓여 있는 위치'가 parent(base-dir) 안인지 확인한다(경로 탈출 방지).

    중요(심볼릭 링크 처리):
    - child 자체를 resolve()하면 *링크가 가리키는 대상*까지 따라가 버린다. 그러면
      base-dir 안에 있는 링크(예: 이름이 msas 인데 외부 폴더를 가리키는 링크)가
      "밖을 가리킨다"는 이유로 잘못 걸러져, 정작 우리가 지우려는 링크 자체를 못
      지운다. 우리가 하려는 일은 링크가 '놓인 자리'를 보는 것이므로, child는
      resolve하지 않고 그 **부모 폴더만 resolve**한 뒤 이름을 붙여 위치를 잡는다.
    - 이렇게 해도 안전하다: 링크를 지울 때는 unlink로 *링크만* 제거하므로(아래
      delete_target 참고) 링크가 가리키는 외부 대상은 절대 건드리지 않는다.
    - 부모를 resolve하므로, 탐색 경로 중간에 base-dir 밖으로 나가는 링크가 끼어
      있으면(이 경우 followlinks=False라 보통 못 내려가지만) 그래도 걸러진다.
    """
    try:
        parent_resolved = parent.resolve()
        # child 자체는 따라가지 않고(부모만 실제 경로로) 위치를 잡는다.
        child_loc = child.parent.resolve() / child.name
    except OSError:
        return False
    # child의 위치가 parent와 같거나 그 하위인지 검사
    return parent_resolved == child_loc or parent_resolved in child_loc.parents


def find_targets(base_dir: Path, target_names: Tuple[str, ...]) -> List[Path]:
    """base-dir 아래를 재귀 탐색해 삭제 대상 디렉터리 경로 목록을 모은다.

    - os.walk(followlinks=False)로 심볼릭 링크 디렉터리는 따라 들어가지 않는다.
    - 이름이 target_names 중 하나와 일치하는 디렉터리만 모은다.
    - 매칭된 디렉터리는 그 안으로 더 내려가지 않도록 dirs 목록에서 제거한다
      (중첩된 같은 이름 폴더를 중복으로 처리하지 않기 위함).
    - 최종적으로 base-dir 밖에 '놓인' 항목은 안전상 제외한다(아래 _is_within 참고).
    """
    targets: List[Path] = []
    for root, dirs, _files in os.walk(base_dir, topdown=True, followlinks=False):
        # topdown=True 이므로 dirs를 그 자리에서 수정하면 탐색 경로를 제어할 수 있다.
        keep_dirs: List[str] = []
        for d in dirs:
            full = Path(root) / d
            if _matches_target(d, target_names):
                # 대상이면 목록에 추가하고, 그 아래로는 내려가지 않는다.
                targets.append(full)
                # keep_dirs에 넣지 않음 => os.walk가 이 폴더로 들어가지 않음
            else:
                keep_dirs.append(d)
        dirs[:] = keep_dirs  # os.walk가 실제로 방문할 하위 폴더 목록을 갱신

    # 안전장치: 위치가 base-dir 밖인 경로는 제거한다.
    # (링크가 base-dir 안에 살면서 바깥을 가리키는 경우는 여기서 걸러지지 않는다 —
    #  링크 '위치'는 안에 있고, 삭제는 링크만 unlink하므로 안전하다.)
    safe_targets = [t for t in targets if _is_within(t, base_dir)]
    return safe_targets


def summarize(targets: List[Path]) -> Tuple[int, int]:
    """삭제 대상 목록을 화면에 출력하고 (총 개수, 총 용량 바이트)를 돌려준다."""
    total_bytes = 0
    for t in targets:
        if t.is_symlink():
            # 링크는 용량 0으로 보고(실제 대상 크기를 더하지 않음)
            size = 0
            tag = " (심볼릭 링크 -> 링크만 제거)"
        else:
            size = _dir_size(t)
            tag = ""
        total_bytes += size
        print(f"  - {t}{tag}  [{_human_size(size)}]")
    return len(targets), total_bytes


def delete_target(path: Path) -> None:
    """대상 하나를 안전하게 삭제한다.

    - 심볼릭 링크면 폴더 내용을 건드리지 않고 링크 자체만 지운다(unlink).
    - 일반 폴더면 shutil.rmtree로 통째로 지운다.
    - 그 외(예: 일반 파일)는 unlink로 처리한다.
    """
    # 반드시 심볼릭 링크를 가장 먼저 확인한다. (is_dir()는 링크가 폴더를 가리키면
    # True가 되므로, 순서를 바꾸면 링크가 가리키는 외부 폴더를 통째로 지울 위험이
    # 있다. 링크는 unlink로 링크 자체만 제거해야 안전하다.)
    if path.is_symlink():
        # 링크가 가리키는 실제 대상은 건드리지 않고, 링크만 제거
        path.unlink()
    elif path.is_dir():
        # rmtree는 폴더 내용을 따라 들어가지만, 위에서 링크는 이미 걸러냈으므로
        # 여기 도달한 것은 실제 디렉터리다. (방어적으로, 만약을 위한 추가 확인)
        if os.path.islink(path):  # 이론상 도달 불가 — 이중 안전장치
            path.unlink()
        else:
            shutil.rmtree(path)
    else:
        # 디렉터리가 아니면(혹시 같은 이름의 파일) 파일만 제거
        path.unlink()


def main(argv: List[str] | None = None) -> int:
    """CLI 진입점. 반환값은 프로세스 종료 코드(0=정상)."""
    parser = argparse.ArgumentParser(
        description=(
            "AF3 출력 폴더 아래의 msas/template 등 보조 폴더를 재귀적으로 삭제한다. "
            "기본은 DRY-RUN(예행 연습)이며, 실제 삭제는 --apply가 있어야만 한다."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "base_dir",
        type=str,
        help="탐색을 시작할 기준 폴더(이 폴더 아래만 건드린다).",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="이 플래그가 있어야만 실제로 삭제한다. 없으면 목록만 출력(안전 기본값).",
    )
    parser.add_argument(
        "--names",
        nargs="+",
        default=list(DEFAULT_TARGET_NAMES),
        metavar="NAME",
        help=(
            "삭제 대상 폴더 이름들(대소문자 무시). "
            f"기본값: {' '.join(DEFAULT_TARGET_NAMES)}"
        ),
    )
    args = parser.parse_args(argv)

    base_dir = Path(args.base_dir)

    # --- 입력 검증 ---
    if not base_dir.exists():
        print(f"[오류] base-dir가 존재하지 않습니다: {base_dir}", file=sys.stderr)
        return 2
    if not base_dir.is_dir():
        print(f"[오류] base-dir가 폴더가 아닙니다: {base_dir}", file=sys.stderr)
        return 2

    target_names = tuple(args.names)

    print(f"[정보] 기준 폴더(base-dir): {base_dir.resolve()}")
    print(f"[정보] 삭제 대상 폴더 이름: {', '.join(target_names)}")
    print()

    # --- 대상 수집 ---
    targets = find_targets(base_dir, target_names)

    if not targets:
        print("삭제할 대상이 없습니다. (조건에 맞는 폴더를 못 찾음)")
        return 0

    mode = "실제 삭제(--apply)" if args.apply else "DRY-RUN(예행 연습 — 아무것도 안 지움)"
    print(f"== 삭제 대상 목록 [{mode}] ==")
    count, total_bytes = summarize(targets)
    print()
    print(f"총 {count}개 폴더, 약 {_human_size(total_bytes)} (용량은 추정치).")
    print()

    if not args.apply:
        # 기본 안전 동작: 실제로는 지우지 않고 안내만 한다.
        print("⚠️ 지금은 DRY-RUN 모드라 아무것도 삭제하지 않았습니다.")
        print("   위 목록이 맞다면, 같은 명령에 --apply 를 붙여 다시 실행하세요.")
        return 0

    # --- 실제 삭제 단계 ---
    print("⚠️ 실제 삭제를 시작합니다. 이 작업은 되돌릴 수 없습니다.")
    deleted = 0
    failed = 0
    for t in targets:
        # 삭제 직전에 다시 한 번 위치가 base-dir 안인지 확인(이중 안전장치)
        if not _is_within(t, base_dir):
            print(f"  [건너뜀] base-dir 밖이라 거부: {t}", file=sys.stderr)
            failed += 1
            continue
        try:
            delete_target(t)
            print(f"  [삭제됨] {t}")
            deleted += 1
        except OSError as exc:
            print(f"  [실패] {t} ({exc})", file=sys.stderr)
            failed += 1

    print()
    print(f"완료: {deleted}개 삭제, {failed}개 실패.")
    # 실패가 하나라도 있으면 0이 아닌 종료 코드로 알린다.
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
