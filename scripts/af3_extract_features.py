#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""af3_extract_features.py

AlphaFold3(AF3) 도킹 출력 폴더를 재귀 탐색해서, 같은 폴더 안에 있는
`*full_data_*_#.json` 과 `*summary_confidences_*_#.json` 파일 쌍을 인식하고,
항원-항체 '계면(interface)' 피처를 뽑아 CSV로 저장하는 CLI 도구입니다.

이 스크립트가 하는 일을 한 줄로:
  "AF3가 뱉은 수많은 JSON을 돌면서, 항원-항체가 맞닿는 부분의 신뢰도 수치를
   표(CSV) 한 장으로 정리한다."

--------------------------------------------------------------------------
배경 용어 (처음 보는 사람을 위해)
--------------------------------------------------------------------------
- AF3: 단백질 구조를 예측하는 AI. 항원-항체 복합체를 넣으면 둘이 어떻게
  붙는지(도킹) 구조를 만들어 줍니다.
- token / 토큰: AF3 내부에서 잔기(아미노산 1개)를 세는 기본 단위라고
  생각하면 됩니다. 이 스크립트에서는 토큰 = 잔기(residue) 수준 좌표로 다룹니다.
- chain / 체인: 단백질 사슬 1개. AF3에 "위=항원, 아래=항체" 순서로 넣으면
  보통 chain A = 항원, chain B = 항체로 매핑됩니다. 단, 항체를 Fab/scFv로
  넣으면 체인이 여러 개(예: B=중쇄절편, C=경쇄)일 수 있습니다.
- ipTM: AF3가 주는 '복합체 전체' 신뢰도 점수(0~1, 높을수록 좋음). 전역
  점수라서 계면의 국소 디테일은 못 잡는다는 게 이 프로젝트의 출발점입니다.
- PAE (Predicted Aligned Error): 토큰 i 기준으로 토큰 j의 위치 오차 예측
  (Å 단위, 낮을수록 좋음). i->j 와 j->i 가 다른 '비대칭' 행렬입니다.
  그래서 AtoB(항원->항체)와 BtoA(항체->항원)를 둘 다 따로 봅니다.
- contact_probs: 토큰 i와 j가 접촉(contact)할 확률 행렬. 대칭(i,j == j,i)
  이라 한 방향(AtoB)만 봐도 됩니다.
- pLDDT: 원자(또는 잔기)별 국소 구조 신뢰도(0~100, 높을수록 좋음).

--------------------------------------------------------------------------
외부 의존성 (이 두 개만 필요)
--------------------------------------------------------------------------
- numpy
- pandas
표준 라이브러리(os, re, json, argparse, collections)만 추가로 씁니다.
설치:  pip install numpy pandas

--------------------------------------------------------------------------
원본(Colab) 대비 고친 버그 / 추가 기능
--------------------------------------------------------------------------
(1) [핵심] 체인이 정확히 A,B 두 개이고 'A가 먼저 연속, 그 다음 B가 연속'
    이라고 가정하던 부분(`token_chains.index('B')`로 경계 1개를 잡던 방식)을
    완전히 제거했습니다. 이제는 token_chain_ids를 일반적으로 파싱해서
    '항원 체인 집합'과 '항체 체인 집합'에 속하는 토큰 인덱스 마스크를 만들고,
    그 마스크로 계면 교차블록을 뽑습니다. 체인이 B,C,D... 여러 개여도,
    순서가 섞여 있어도(A,B,A,C ...) 동작합니다. 항원 = --antigen-chain,
    항체는 --antibody-chains(미지정 시 항원이 아닌 '나머지 전부').
(2) pLDDT를 '전역 top-N'만 보던 것에 더해, '계면 잔기'에 한정한 pLDDT
    피처도 추가했습니다(--plddt-interface). 계면 잔기 = contact_probs 기준
    접촉 토큰. 전역/계면 둘 다 별도 컬럼으로 저장합니다.
(3) contact/pae 의 top-k 인덱스가 원본에서는 '블록 안 로컬 좌표'였습니다.
    이를 전역(global) 토큰 인덱스로 변환해 저장하고, 컬럼명/주석에 좌표계를
    분명히 적었습니다(컬럼명: *_global_token_idx).
(4) seed 추출이 `_s(\d+)` 정규식 하나에만 의존하던 것을, 파일명 ->
    상위 폴더명 -> summary json 내부 키 순서로 다단계 탐색하도록 바꿨습니다.
    못 찾으면 NaN + 경고. 어디서 찾았는지 seed_source 컬럼으로 추적 가능.
    rank(모델 인덱스) 파싱도 AF3 표준 파일명(예: fold_xxx_full_data_0.json)을
    견고하게 인식합니다.
(5) typ(분류 라벨) 추론을 옵션화했습니다(--typ-level). 기본은 base-dir
    바로 아래 1단계 폴더명(원본 기본 동작). 동작은 determine_typ() docstring 참고.

(6) [검증 단계 신규 수정] contact_probs / pae 행렬의 차원이 토큰 수와 맞는지
    검사합니다. 행렬이 토큰 수보다 작으면 원본/초안은 IndexError로 '실행 전체가
    크래시'했고, 토큰 수보다 크면 '조용히 잘못된 값'을 냈습니다. 이제 둘 다
    경고를 찍고 그 행렬 피처만 건너뜁니다(다른 피처·다른 모델은 정상 처리).
    또한 모델 한 건에서 예상치 못한 예외가 나도 그 모델만 건너뛰고 배치 전체는
    계속 돕니다(수백 개를 무인 실행할 때 한 파일이 전부를 죽이지 않게).

--------------------------------------------------------------------------
⚠️ 경고 (반드시 읽기)
--------------------------------------------------------------------------
- 이 스크립트는 실제 AF3 출력 데이터로 테스트하지 못했습니다(합성 데이터로
  로직·엣지케이스만 검증). JSON 키 이름(token_chain_ids, contact_probs, pae,
  atom_plddts, atom_chain_ids 등)과 파일명 패턴은 일반적인 AF3 관례를 가정한
  것이며, AF3 버전(웹서버 vs 로컬)에 따라 다를 수 있습니다. 처음 한 케이스를
  돌려보고 컬럼이 비어있지 않은지, 값이 말이 되는지 꼭 눈으로 확인하세요.
- contact_probs / pae 는 '토큰(잔기) 단위' 정사각 행렬이라고 가정합니다.
  atom_plddts 는 '원자 단위'라 길이가 다릅니다(체인 매핑은 atom_chain_ids
  사용). 만약 AF3 버전이 다른 단위를 쓰면, (6)의 차원 검사가 그 행렬을 건너뛰니
  관련 컬럼이 비어 보일 수 있습니다 — 그러면 키/단위를 확인하세요.
- --plddt-interface 의 계면 pLDDT는 '원자->토큰 비례 배분'이라는 근사에
  의존합니다(잔기당 원자 수가 달라 정확하지 않음). '대략적 참고치'로만 쓰세요.

--------------------------------------------------------------------------
사용 예시
--------------------------------------------------------------------------
  python af3_extract_features.py \
      --base-dir /path/to/af3_outputs \
      --output-dir /path/to/analysis_out

  # 항체가 Fab(체인 B,C)이고 항원이 A인 경우 (보통 기본값으로 충분):
  python af3_extract_features.py \
      --base-dir ./outputs --output-dir ./analysis \
      --antigen-chain A --antibody-chains B C

  # 계면 pLDDT 피처까지 켜고, 접촉 임계값을 0.2로:
  python af3_extract_features.py \
      --base-dir ./outputs --output-dir ./analysis \
      --plddt-interface --interface-cutoff 0.2
"""

from __future__ import annotations

import argparse
import json
import os
import re
from collections import defaultdict
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd


# =====================================================================
# 1) 작은 수치 유틸리티
# =====================================================================
def topn_mean_1d(a: Any, n: int, largest: bool = True) -> float:
    """1차원으로 펼친 배열에서 상위/하위 n개의 평균을 구한다.

    Args:
        a: 숫자 배열(다차원이면 1차원으로 펼침).
        n: 몇 개를 평균할지.
        largest: True면 '큰 값 n개'(예: contact_probs, pLDDT),
                 False면 '작은 값 n개'(예: PAE는 낮을수록 좋으므로).

    Returns:
        상위/하위 n개의 평균(float). 배열이 비어 있으면 NaN.
    """
    arr = np.asarray(a, dtype=float).ravel()
    if arr.size == 0:
        return float("nan")
    k = min(max(int(n), 1), arr.size)
    if largest:
        part = np.partition(arr, arr.size - k)[-k:]
    else:
        part = np.partition(arr, k - 1)[:k]
    return float(part.mean())


def topk_local_indices(block: np.ndarray, k: int, largest: bool = True) -> List[Tuple[int, int]]:
    """2차원 블록에서 상위/하위 k개 값의 (행, 열) '로컬' 좌표를 정렬해 반환한다.

    여기서 '로컬'이란 입력으로 받은 block 안에서의 행/열 번호라는 뜻입니다.
    바깥에서 이 좌표를 '전역 토큰 인덱스'로 다시 변환해야 합니다
    (local_to_global_indices 참고).

    Args:
        block: 2차원 numpy 배열(교차블록).
        k: 뽑을 개수.
        largest: True면 큰 값, False면 작은 값 기준.

    Returns:
        (row, col) 정수 튜플의 리스트. largest=True면 큰 값부터,
        largest=False면 작은 값부터 정렬.
    """
    arr = np.asarray(block)
    flat = arr.ravel()
    k = min(k, flat.size)
    if k == 0:
        return []
    if largest:
        sel = np.argpartition(flat, flat.size - k)[-k:]
        sel_sorted = sel[np.argsort(flat[sel])[::-1]]
    else:
        sel = np.argpartition(flat, k - 1)[:k]
        sel_sorted = sel[np.argsort(flat[sel])]
    rows, cols = np.unravel_index(sel_sorted, arr.shape)
    return [(int(r), int(c)) for r, c in zip(rows, cols)]


def local_to_global_indices(
    local_pairs: Sequence[Tuple[int, int]],
    row_global: Sequence[int],
    col_global: Sequence[int],
) -> List[List[int]]:
    """블록 안 (행,열) 로컬 좌표를 전역 토큰 인덱스 (행,열) 로 변환한다.

    교차블록은 '항원 토큰들(행) x 항체 토큰들(열)' 식으로 원래 큰 행렬에서
    특정 토큰들만 골라 만든 부분 행렬입니다. 따라서 블록의 행 i는 실제로는
    row_global[i] 번째 토큰, 블록의 열 j는 col_global[j] 번째 토큰입니다.

    Args:
        local_pairs: (row, col) 로컬 좌표 리스트.
        row_global: 블록의 각 행이 가리키는 전역 토큰 인덱스 배열.
        col_global: 블록의 각 열이 가리키는 전역 토큰 인덱스 배열.

    Returns:
        [global_row_token, global_col_token] 리스트의 리스트.
    """
    out: List[List[int]] = []
    for r, c in local_pairs:
        out.append([int(row_global[r]), int(col_global[c])])
    return out


# =====================================================================
# 2) 체인 마스크 / 교차블록 (버그 (1) 핵심 수정)
# =====================================================================
def build_chain_masks(
    token_chain_ids: Sequence[Any],
    antigen_chain: str,
    antibody_chains: Optional[Sequence[str]],
) -> Tuple[np.ndarray, np.ndarray, List[str], List[str]]:
    """token_chain_ids로부터 항원/항체 토큰 마스크를 만든다(연속/정렬 가정 없음).

    원본 코드는 'A가 앞에 쭉, B가 뒤에 쭉'이라고 가정하고 경계 1개
    (b_index)만 잡았습니다. 항체가 Fab/scFv면 체인이 여러 개라 이 가정이
    깨집니다. 여기서는 각 토큰이 어느 체인인지 라벨을 직접 보고, 항원에
    속하면 True, 항체에 속하면 True 하는 식의 불리언 마스크를 만듭니다.
    체인 순서가 섞여 있어도(A,B,A,C ...) 안전합니다.

    Args:
        token_chain_ids: 토큰별 체인 기호 리스트(예: ['A','A',...,'B','B','C',...]).
        antigen_chain: 항원 체인 기호(예: 'A').
        antibody_chains: 항체 체인 기호 리스트. None/빈 값이면 '항원이 아닌 나머지 전부'.

    Returns:
        (antigen_mask, antibody_mask, present_chains, antibody_chains_used)
        - antigen_mask: 길이 = 토큰 수, 항원 토큰이면 True.
        - antibody_mask: 항체 토큰이면 True.
        - present_chains: 실제 등장한 체인 기호(정렬, 중복 제거).
        - antibody_chains_used: 실제로 항체로 사용된 체인 기호들.
    """
    chains = np.asarray([str(c) for c in token_chain_ids])
    present = sorted(set(chains.tolist()))

    antigen_chain = str(antigen_chain)
    if chains.size == 0:
        empty = np.zeros(0, dtype=bool)
        antibody_set = [str(c) for c in antibody_chains] if antibody_chains else []
        return empty, empty, present, antibody_set

    antigen_mask = chains == antigen_chain

    if antibody_chains is None or len(antibody_chains) == 0:
        # 항원이 아닌 모든 체인을 항체로 간주
        antibody_set = [c for c in present if c != antigen_chain]
    else:
        antibody_set = [str(c) for c in antibody_chains]

    # 항원/항체 집합이 겹치면(설정 실수) 항체에서 항원 체인을 빼서 교차블록이
    # 자기 자신을 포함하지 않게 한다.
    antibody_set = [c for c in antibody_set if c != antigen_chain]

    antibody_mask = np.isin(chains, antibody_set)

    return antigen_mask, antibody_mask, present, antibody_set


def is_square_token_matrix(matrix: Any, n_token: int) -> Tuple[bool, np.ndarray]:
    """행렬이 '토큰 수 x 토큰 수' 정사각인지 검사한다 (버그 (6)).

    contact_probs / pae 는 (토큰, 토큰) 정사각이라고 가정합니다. 하지만 AF3
    버전/단위에 따라 차원이 토큰 수와 다를 수 있습니다. 행렬이 토큰 수보다
    작으면 마스크 인덱싱에서 IndexError로 '실행 전체가' 크래시하고, 크면
    '조용히 엉뚱한 토큰'을 가리켜 잘못된 값을 냅니다. 그래서 미리 검사합니다.

    Args:
        matrix: 검사할 행렬(리스트 또는 ndarray).
        n_token: token_chain_ids 길이. 0이면 토큰 길이 검사는 건너뜀(정사각만 확인).

    Returns:
        (ok, mat) — ok가 True면 안전하게 써도 됨. mat는 float ndarray로 변환된 행렬.
    """
    mat = np.asarray(matrix, dtype=float)
    if mat.ndim != 2:
        return False, mat
    if mat.shape[0] != mat.shape[1]:
        return False, mat
    if n_token > 0 and mat.shape[0] != n_token:
        return False, mat
    return True, mat


def cross_block(
    matrix: np.ndarray,
    row_mask: np.ndarray,
    col_mask: np.ndarray,
) -> np.ndarray:
    """정사각 행렬에서 (row_mask 토큰들) x (col_mask 토큰들) 교차블록을 뽑는다.

    예) contact 행렬 M에서 항원 토큰(행) x 항체 토큰(열) 만 잘라낸 부분.
    np.ix_ 를 쓰면 마스크가 가리키는 행/열 조합만 골라낸 2D 부분행렬이 됩니다.

    안전장치: 마스크 인덱스가 행렬 범위를 벗어나면(차원 불일치) 빈 블록을
    돌려줍니다(크래시 방지). 단, 호출 전에 is_square_token_matrix로 거르는
    것을 권장합니다.

    Args:
        matrix: (T, T) 정사각 행렬(T = 토큰 수).
        row_mask: 행으로 쓸 토큰 마스크.
        col_mask: 열로 쓸 토큰 마스크.

    Returns:
        (len(row), len(col)) 부분 행렬. 추출 불가면 (0,0) 빈 배열.
    """
    mat = np.asarray(matrix, dtype=float)
    if mat.ndim != 2:
        return np.empty((0, 0), dtype=float)
    row_idx = np.where(row_mask)[0]
    col_idx = np.where(col_mask)[0]
    if row_idx.size == 0 or col_idx.size == 0:
        return np.empty((0, 0), dtype=float)
    if int(row_idx.max()) >= mat.shape[0] or int(col_idx.max()) >= mat.shape[1]:
        # 차원 불일치 — 크래시 대신 빈 블록
        return np.empty((0, 0), dtype=float)
    return mat[np.ix_(row_idx, col_idx)]


# =====================================================================
# 3) pLDDT (원자 단위) - 전역 + 계면 (버그 (2) 추가)
# =====================================================================
def plddt_topn_global(atom_plddts: Any, ns: Sequence[int]) -> Dict[str, float]:
    """전체 원자 pLDDT에서 상위 N개 평균(전역) 피처를 만든다.

    Returns:
        {'plddt_top1_mean': ..., 'plddt_top3_mean': ..., ...}
    """
    a = np.asarray(atom_plddts, dtype=float).ravel()
    out: Dict[str, float] = {}
    for n in ns:
        out[f"plddt_top{n}_mean"] = topn_mean_1d(a, n, largest=True)
    return out


def plddt_top_indices_with_chain(
    atom_plddts: Any,
    atom_chain_ids: Sequence[Any],
    k: int = 10,
) -> List[List[Any]]:
    """pLDDT 상위 k개 원자를 (체인기호, 그 체인 내 로컬 원자번호)로 반환한다.

    원자 단위라 토큰 인덱스와 직접 비교는 어렵지만, 어느 체인의 원자가
    높은 신뢰도인지 파악하는 용도로 둡니다. (좌표계: 체인 내 로컬 원자 인덱스)

    Returns:
        [[chain_symbol, local_atom_index], ...] (pLDDT 높은 순).
    """
    a = np.asarray(atom_plddts, dtype=float).ravel()
    chains = [str(c) for c in atom_chain_ids]
    if a.size == 0:
        return []

    # 각 원자에 대해 (체인기호, 그 체인 안에서 몇 번째 원자인지) 라벨 부여
    per_chain_counter: Dict[str, int] = defaultdict(int)
    labels: List[Tuple[str, int]] = []
    for c in chains:
        labels.append((c, per_chain_counter[c]))
        per_chain_counter[c] += 1

    # 라벨 수와 pLDDT 수가 어긋나면(데이터 불일치) 짧은 쪽에 맞춘다.
    m = min(a.size, len(labels))
    if m == 0:
        return []
    a = a[:m]
    labels = labels[:m]

    k = min(k, a.size)
    sel = np.argpartition(a, a.size - k)[-k:]
    sel_sorted = sel[np.argsort(a[sel])[::-1]]
    return [[labels[int(i)][0], int(labels[int(i)][1])] for i in sel_sorted]


def map_atoms_to_tokens(
    atom_chain_ids: Sequence[Any],
    token_chain_ids: Sequence[Any],
) -> Optional[np.ndarray]:
    """원자 -> 토큰 매핑을 '체인별 순서'로 추정한다(근사).

    AF3 full_data 에 원자->토큰 직접 매핑이 항상 들어있지는 않습니다.
    여기서는 'atom_chain_ids 와 token_chain_ids 모두 같은 체인 순서로
    나열돼 있다'는 가정 아래, 각 체인 안에서 원자를 토큰에 비례 배분합니다.

    주의: 이는 근사입니다. 잔기당 원자 수가 다르기 때문에 정확하지 않을 수
    있어, 계면 pLDDT는 '대략적 참고치'로만 쓰세요. 매핑이 불가능하면
    None을 반환합니다(계면 pLDDT 계산을 건너뜀).

    Returns:
        길이 = 원자 수, 각 원자가 속한 전역 토큰 인덱스 배열. 또는 None.
    """
    atom_chains = [str(c) for c in atom_chain_ids]
    token_chains = [str(c) for c in token_chain_ids]
    if len(atom_chains) == 0 or len(token_chains) == 0:
        return None

    # 체인별 토큰 전역 인덱스 목록
    token_idx_by_chain: Dict[str, List[int]] = defaultdict(list)
    for ti, c in enumerate(token_chains):
        token_idx_by_chain[c].append(ti)

    # 체인별 원자 개수
    atom_count_by_chain: Dict[str, int] = defaultdict(int)
    for c in atom_chains:
        atom_count_by_chain[c] += 1

    # 모든 원자 체인이 토큰에도 존재해야 매핑 가능
    for c in atom_count_by_chain:
        if c not in token_idx_by_chain or len(token_idx_by_chain[c]) == 0:
            return None

    atom_to_token = np.full(len(atom_chains), -1, dtype=int)
    seen_in_chain: Dict[str, int] = defaultdict(int)
    for ai, c in enumerate(atom_chains):
        toks = token_idx_by_chain[c]
        n_atoms = atom_count_by_chain[c]
        pos = seen_in_chain[c]
        seen_in_chain[c] += 1
        # 체인 안에서 원자 위치를 토큰 개수에 비례시켜 배분(근사)
        t_local = int(pos * len(toks) / max(1, n_atoms))
        t_local = min(t_local, len(toks) - 1)
        atom_to_token[ai] = toks[t_local]
    return atom_to_token


def plddt_interface_features(
    atom_plddts: Any,
    atom_chain_ids: Sequence[Any],
    token_chain_ids: Sequence[Any],
    interface_token_set: Set[int],
    ns: Sequence[int],
) -> Dict[str, float]:
    """계면 토큰에 매핑되는 원자들만으로 pLDDT 통계 피처를 만든다(버그 (2)).

    Args:
        atom_plddts: 원자별 pLDDT.
        atom_chain_ids: 원자별 체인 기호.
        token_chain_ids: 토큰별 체인 기호.
        interface_token_set: '계면 잔기'로 판정된 전역 토큰 인덱스 집합.
        ns: 상위 N 목록.

    Returns:
        {'plddt_iface_mean': ..., 'plddt_iface_top1_mean': ..., 'plddt_iface_n_atoms': ...}
        매핑 불가/계면 없음이면 통계는 NaN(개수는 0).
    """
    out: Dict[str, float] = {"plddt_iface_mean": float("nan"), "plddt_iface_n_atoms": 0.0}
    for n in ns:
        out[f"plddt_iface_top{n}_mean"] = float("nan")

    a = np.asarray(atom_plddts, dtype=float).ravel()
    if a.size == 0 or len(interface_token_set) == 0:
        return out

    atom_to_token = map_atoms_to_tokens(atom_chain_ids, token_chain_ids)
    if atom_to_token is None:
        return out

    # 길이 불일치 방어
    m = min(a.size, atom_to_token.size)
    a = a[:m]
    atom_to_token = atom_to_token[:m]

    sel_mask = np.array([int(t) in interface_token_set for t in atom_to_token], dtype=bool)
    iface_vals = a[sel_mask]
    out["plddt_iface_n_atoms"] = float(iface_vals.size)
    if iface_vals.size == 0:
        return out

    out["plddt_iface_mean"] = float(iface_vals.mean())
    for n in ns:
        out[f"plddt_iface_top{n}_mean"] = topn_mean_1d(iface_vals, n, largest=True)
    return out


# =====================================================================
# 4) 계면 토큰 판정 (버그 (2) 보조)
# =====================================================================
def find_interface_tokens(
    contact_block: np.ndarray,
    antigen_global_idx: np.ndarray,
    antibody_global_idx: np.ndarray,
    contact_threshold: float,
) -> Set[int]:
    """contact 교차블록에서 임계값 이상 접촉하는 토큰들을 계면 토큰으로 모은다.

    교차블록의 (i, j) 가 임계값을 넘으면, 그 행의 항원 토큰과 그 열의 항체
    토큰을 둘 다 '계면 토큰'으로 봅니다. 결과는 전역 토큰 인덱스 집합.

    Args:
        contact_block: 항원(행) x 항체(열) contact 확률 부분행렬.
        antigen_global_idx: 블록 각 행의 전역 토큰 인덱스.
        antibody_global_idx: 블록 각 열의 전역 토큰 인덱스.
        contact_threshold: 접촉 판정 확률 임계값(예: 0.1).

    Returns:
        계면으로 판정된 전역 토큰 인덱스 집합(set of int).
    """
    if contact_block.size == 0:
        return set()
    hit_rows, hit_cols = np.where(contact_block >= contact_threshold)
    iface: Set[int] = set()
    for r in np.unique(hit_rows):
        iface.add(int(antigen_global_idx[r]))
    for c in np.unique(hit_cols):
        iface.add(int(antibody_global_idx[c]))
    return iface


# =====================================================================
# 5) 파일 탐색 / 파싱 (버그 (4) seed·rank 견고화)
# =====================================================================
def collect_rank_pairs(case_dir: str) -> List[Tuple[int, str, str]]:
    """한 폴더 안에서 rank별 (full_data, summary_confidences) 파일 쌍을 모은다.

    AF3 표준 파일명 예: fold_xxx_full_data_0.json / fold_xxx_summary_confidences_0.json
    여기서 끝의 _0, _1 ... 이 rank(모델 인덱스)입니다. 두 종류가 같은 rank로
    모두 있을 때만 분석 대상으로 봅니다.

    Returns:
        [(rank, full_path, summary_path), ...] rank 오름차순.
    """
    full_by_rank: Dict[int, str] = {}
    summ_by_rank: Dict[int, str] = {}

    try:
        names = os.listdir(case_dir)
    except OSError:
        return []

    for fname in names:
        if fname.startswith(".") or fname == "Thumbs.db":
            continue
        # summary_confidences 가 full_data 보다 먼저 와야 한다(이름에 'full_data'가
        # 들어가지 않으므로 충돌은 없지만, 순서를 명확히 하기 위해 summary 먼저 검사).
        m_sm = re.search(r"summary_confidences.*?_(\d+)\.json$", fname)
        if m_sm:
            summ_by_rank[int(m_sm.group(1))] = os.path.join(case_dir, fname)
            continue
        m_fd = re.search(r"full_data.*?_(\d+)\.json$", fname)
        if m_fd:
            full_by_rank[int(m_fd.group(1))] = os.path.join(case_dir, fname)

    ranks = sorted(set(full_by_rank) & set(summ_by_rank))
    return [(r, full_by_rank[r], summ_by_rank[r]) for r in ranks]


def extract_seed(
    full_path: str,
    case_dir: str,
    summary: Optional[Dict[str, Any]],
) -> Tuple[float, str]:
    """seed 값을 여러 곳에서 다단계로 찾는다(버그 (4)).

    탐색 순서:
      1) 파일명 안의 _s<숫자> 또는 seed<숫자> 패턴
      2) 상위 폴더명들 안의 같은 패턴(case_dir 부터 위로)
      3) summary json 내부에 'seed' 비슷한 키가 있으면 그 값
    못 찾으면 (NaN, 'seed<-NOT_FOUND') 반환.

    주의: 2단계(상위 폴더명)는 무관한 조상 폴더의 _s<n>을 잘못 잡을 수 있습니다.
    어디서 찾았는지 seed_source 컬럼으로 항상 추적되니, 의심되면 확인하세요.

    Returns:
        (seed_value_or_nan, note). note는 어디서 찾았는지/실패했는지 설명.
    """
    seed_pat = re.compile(r"(?:_s|seed[_-]?)(\d+)", re.IGNORECASE)

    # 1) 파일명
    m = seed_pat.search(os.path.basename(full_path))
    if m:
        return float(int(m.group(1))), "seed<-filename"

    # 2) 상위 폴더명들 (case_dir 부터 위로)
    path_part = case_dir
    while path_part and path_part not in ("/", ""):
        base = os.path.basename(os.path.normpath(path_part))
        m = seed_pat.search(base)
        if m:
            return float(int(m.group(1))), "seed<-parentdir"
        parent = os.path.dirname(os.path.normpath(path_part))
        if parent == path_part:
            break
        path_part = parent

    # 3) summary json 내부
    if summary:
        for key in summary:
            if "seed" in str(key).lower():
                val = summary[key]
                try:
                    return float(int(val)), f"seed<-summary[{key}]"
                except (TypeError, ValueError):
                    # 값이 리스트/문자열이면 첫 숫자만 시도
                    mm = re.search(r"\d+", str(val))
                    if mm:
                        return float(int(mm.group(0))), f"seed<-summary[{key}]"

    return float("nan"), "seed<-NOT_FOUND"


def determine_typ(case_dir: str, base_dir: str, typ_level: int) -> str:
    """case_dir로부터 typ(분류 라벨, 출력 CSV 파일명에 쓰임)을 정한다(버그 (5)).

    typ_level 의 의미:
      - 0 : base_dir 폴더명 자체를 typ으로 (모든 케이스가 한 CSV로 합쳐짐)
      - 1 : base_dir 바로 아래 1단계 하위 폴더명 (원본 기본 동작)
      - 2 : base_dir 아래 2단계 폴더명, 없으면 가능한 가장 깊은 단계로 fallback
      - n : base_dir 기준 상대경로의 n번째 세그먼트(1부터)

    상대경로가 base_dir 바깥('..')이거나 '.'이면 base_dir 폴더명을 씁니다.

    Returns:
        typ 문자열(파일시스템 안전 문자로 가공하지는 않음 — sanitize_filename 별도).
    """
    base_name = os.path.basename(os.path.normpath(base_dir))
    if typ_level <= 0:
        return base_name

    rel = os.path.relpath(case_dir, base_dir)
    if rel == "." or rel.startswith(".."):
        return base_name

    segments = rel.split(os.sep)
    # typ_level 번째 세그먼트(1-based). 부족하면 마지막 세그먼트로 fallback.
    idx = min(typ_level - 1, len(segments) - 1)
    return segments[idx]


def sanitize_filename(name: str) -> str:
    """typ을 파일명으로 쓸 때 위험한 문자를 _ 로 바꾼다."""
    return re.sub(r"[^0-9A-Za-z가-힣._-]+", "_", str(name)).strip("_") or "case"


# =====================================================================
# 6) 한 케이스(rank 하나) -> 피처 row 만들기
# =====================================================================
def process_pair(
    rank: int,
    full_path: str,
    summary_path: str,
    typ: str,
    case_name: str,
    case_dir: str,
    args: argparse.Namespace,
) -> Optional[Dict[str, Any]]:
    """full_data + summary 한 쌍을 읽어 계면 피처 dict(row) 하나를 만든다.

    Returns:
        피처 컬럼이 담긴 dict. 치명적 오류(full_data 파일 깨짐 등)면 None.
    """
    ns = tuple(args.topn)

    # ----- summary 먼저 로드(seed 탐색에도 쓰임) -----
    summary: Optional[Dict[str, Any]] = None
    try:
        with open(summary_path, "r", encoding="utf-8-sig") as f:
            summary = json.load(f)
    except Exception as exc:  # noqa: BLE001 - 파일 깨짐 등 광범위 방어
        print(f"   - 경고: summary 읽기 실패 (rank {rank}): {exc}")
        summary = None

    # ----- full_data 로드 -----
    try:
        with open(full_path, "r", encoding="utf-8-sig") as f:
            data = json.load(f)
    except Exception as exc:  # noqa: BLE001
        print(f"   - 건너뜀: full_data 읽기 실패 (rank {rank}): {exc}")
        return None

    if not isinstance(data, dict):
        print(f"   - 건너뜀: full_data 최상위가 dict가 아님 (rank {rank})")
        return None

    token_chains = data.get("token_chain_ids", []) or []
    if not token_chains:
        print(f"   - 경고: token_chain_ids 없음 (rank {rank}) -> 체인 마스크 불가, "
              f"계면 피처를 건너뜁니다.")

    # ----- 체인 마스크 (버그 (1) 수정) -----
    antigen_mask, antibody_mask, present_chains, ab_used = build_chain_masks(
        token_chains, args.antigen_chain, args.antibody_chains
    )

    if token_chains:
        if not antigen_mask.any():
            print(f"   - 경고: 항원 체인 '{args.antigen_chain}' 토큰을 찾지 못함. "
                  f"등장 체인={present_chains} (rank {rank})")
        if not antibody_mask.any():
            print(f"   - 경고: 항체 체인 {ab_used} 토큰을 찾지 못함. "
                  f"등장 체인={present_chains} (rank {rank})")

    antigen_global_idx = np.where(antigen_mask)[0]
    antibody_global_idx = np.where(antibody_mask)[0]

    # ----- seed (버그 (4)) -----
    seed_val, seed_note = extract_seed(full_path, case_dir, summary)
    if np.isnan(seed_val):
        print(f"   - 경고: seed를 찾지 못해 NaN 처리 (rank {rank}, {full_path})")

    row: Dict[str, Any] = {
        "typ": typ,
        "name": case_name,
        "case_dir": case_dir,
        "full_data_file": os.path.basename(full_path),
        "model_index": rank,          # = rank, AF3 파일명 끝의 _#
        "seed": seed_val,
        "seed_source": seed_note,
        "antigen_chain": args.antigen_chain,
        "antibody_chains": ",".join(ab_used),
        "present_chains": ",".join(present_chains),
        "n_token": len(token_chains),
        "antigen_n_token": int(antigen_global_idx.size),
        "antibody_n_token": int(antibody_global_idx.size),
    }

    have_chains = antigen_global_idx.size > 0 and antibody_global_idx.size > 0
    n_tok = len(token_chains)

    # 계면 토큰 집합(전역 인덱스) - pLDDT 계면 피처에서 사용
    interface_tokens: Set[int] = set()

    # ----- contact_probs (대칭 -> AtoB 한 방향) -----
    contact_ok = False
    if have_chains and "contact_probs" in data:
        ok_cp, _ = is_square_token_matrix(data["contact_probs"], n_tok)
        if not ok_cp:
            print(f"   - 경고: contact_probs 행렬 차원이 토큰 수({n_tok})와 불일치 "
                  f"-> contact 피처 건너뜀 (rank {rank})")
        else:
            contact_ok = True

    if contact_ok:
        cp_block = cross_block(data["contact_probs"], antigen_mask, antibody_mask)
        cp_flat = cp_block.ravel()
        for n in ns:
            row[f"contact_probs_top{n}_mean"] = topn_mean_1d(cp_flat, n, largest=True)

        # top-10 좌표를 전역 토큰 인덱스로 변환 (버그 (3))
        local_idx = topk_local_indices(cp_block, 10, largest=True)
        global_idx = local_to_global_indices(local_idx, antigen_global_idx, antibody_global_idx)
        # 컬럼명에 좌표계 명시: 각 항목 [antigen_token, antibody_token] = 전역 토큰 인덱스
        row["contact_probs_top10_global_token_idx"] = json.dumps(global_idx)

        # 계면 토큰 판정(버그 (2) 보조)
        interface_tokens = find_interface_tokens(
            cp_block, antigen_global_idx, antibody_global_idx, args.contact_threshold
        )
        row["interface_n_token"] = len(interface_tokens)

    # ----- pae (비대칭 -> AtoB, BtoA 둘 다) -----
    pae_arr: Optional[np.ndarray] = None
    if have_chains and "pae" in data:
        ok_pae, pae_arr = is_square_token_matrix(data["pae"], n_tok)
        if not ok_pae:
            print(f"   - 경고: pae 행렬 차원이 토큰 수({n_tok})와 불일치 "
                  f"-> pae 피처 건너뜀 (rank {rank})")
            pae_arr = None

    if pae_arr is not None:
        # AtoB: 행=항원, 열=항체 / BtoA: 행=항체, 열=항원
        pae_AtoB = cross_block(pae_arr, antigen_mask, antibody_mask)
        pae_BtoA = cross_block(pae_arr, antibody_mask, antigen_mask)
        for n in ns:
            # PAE는 낮을수록 좋으므로 '작은 값 top-n' (largest=False)
            row[f"AtoB_pae_min_top{n}_mean"] = topn_mean_1d(pae_AtoB.ravel(), n, largest=False)
            row[f"BtoA_pae_min_top{n}_mean"] = topn_mean_1d(pae_BtoA.ravel(), n, largest=False)

        # top-10 작은 PAE 좌표 -> 전역 토큰 인덱스(버그 (3))
        atob_local = topk_local_indices(pae_AtoB, 10, largest=False)
        atob_global = local_to_global_indices(atob_local, antigen_global_idx, antibody_global_idx)
        btoa_local = topk_local_indices(pae_BtoA, 10, largest=False)
        btoa_global = local_to_global_indices(btoa_local, antibody_global_idx, antigen_global_idx)
        row["pae_min_top10_global_token_idx"] = json.dumps(
            {
                # AtoB: [antigen_token, antibody_token], BtoA: [antibody_token, antigen_token]
                "AtoB": atob_global,
                "BtoA": btoa_global,
            }
        )

    # ----- pLDDT: 전역(원본 유지) + 계면(버그 (2) 신규) -----
    if "atom_plddts" in data and "atom_chain_ids" in data:
        row.update(plddt_topn_global(data["atom_plddts"], ns))
        row["plddt_top10_chain_local_atom_idx"] = json.dumps(
            plddt_top_indices_with_chain(data["atom_plddts"], data["atom_chain_ids"], 10)
        )

        if args.plddt_interface:
            iface_feats = plddt_interface_features(
                data["atom_plddts"],
                data["atom_chain_ids"],
                token_chains,
                interface_tokens,
                ns,
            )
            row.update(iface_feats)

    # ----- summary 스칼라 병합 (iptm 등) -----
    if isinstance(summary, dict):
        for key, val in summary.items():
            # 리스트/딕트 같은 비스칼라는 JSON 문자열로(컬럼 폭발 방지)
            if isinstance(val, (list, dict)):
                val_to_store: Any = json.dumps(val)
            else:
                val_to_store = val
            if key in row:
                row[f"{key}_summary"] = val_to_store
            else:
                row[key] = val_to_store

    return row


# =====================================================================
# 7) 메인 파이프라인
# =====================================================================
def run(args: argparse.Namespace) -> None:
    """base_dir를 재귀 탐색해 모든 케이스를 처리하고 typ별 CSV로 저장한다."""
    print(f"[*] base-dir   : {args.base_dir}")
    print(f"[*] output-dir : {args.output_dir}")
    print(f"[*] 항원 체인  : {args.antigen_chain}")
    print(f"[*] 항체 체인  : {args.antibody_chains if args.antibody_chains else '(항원 외 전부 자동)'}")
    print(f"[*] top-N      : {args.topn}")
    print(f"[*] typ-level  : {args.typ_level}")
    print(f"[*] 계면 pLDDT : {'ON' if args.plddt_interface else 'OFF'} "
          f"(interface-cutoff={args.contact_threshold})")
    print("-" * 50)

    os.makedirs(args.output_dir, exist_ok=True)

    records_per_typ: Dict[str, List[Dict[str, Any]]] = defaultdict(list)

    case_count = 0
    pair_count = 0
    for case_dir, _subdirs, _files in os.walk(args.base_dir):
        pairs = collect_rank_pairs(case_dir)
        if not pairs:
            continue

        case_count += 1
        typ = determine_typ(case_dir, args.base_dir, args.typ_level)
        case_name = os.path.basename(os.path.normpath(case_dir))
        print(f">> Case: typ='{typ}' | dir='{case_dir}' | ranks={[r for r, _, _ in pairs]}")

        for rank, full_path, summary_path in pairs:
            # 모델 한 건의 예외가 배치 전체를 죽이지 않게 방어 (버그 (6))
            try:
                row = process_pair(
                    rank, full_path, summary_path, typ, case_name, case_dir, args
                )
            except Exception as exc:  # noqa: BLE001
                print(f"   - 건너뜀: 처리 중 예외 (rank {rank}, {full_path}): {exc}")
                row = None
            if row is not None:
                records_per_typ[typ].append(row)
                pair_count += 1

    print("-" * 50)
    print(f"[*] 탐지한 케이스 폴더: {case_count}개, 처리한 모델(rank): {pair_count}개")

    if not records_per_typ:
        print("[!] 분석할 파일 쌍을 찾지 못했습니다. "
              "패턴 확인: '*full_data_*_#.json' + '*summary_confidences_*_#.json'")
        return

    for typ, rows in records_per_typ.items():
        if not rows:
            continue
        df = pd.DataFrame(rows)
        # iptm 있으면 내림차순 정렬(좋은 모델이 위로). 안정 정렬 사용.
        if "iptm" in df.columns:
            df = df.sort_values(by="iptm", ascending=False, kind="mergesort")
        out_path = os.path.join(args.output_dir, f"{sanitize_filename(typ)}_simplified.csv")
        df.to_csv(out_path, encoding="utf-8-sig", index=False)
        print(f"  [저장] {len(df)} rows -> {out_path}")

    print("[*] 모든 처리 완료.")


# =====================================================================
# 8) CLI
# =====================================================================
def build_parser() -> argparse.ArgumentParser:
    """argparse 파서를 구성한다."""
    p = argparse.ArgumentParser(
        description=(
            "AF3 도킹 출력 폴더를 재귀 탐색해 항원-항체 계면 피처를 CSV로 추출합니다. "
            "(full_data + summary_confidences 쌍 인식)"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--base-dir",
        required=True,
        help="재귀 탐색을 시작할 루트 폴더. 이 아래에서 full_data/summary 쌍을 자동 인식.",
    )
    p.add_argument(
        "--output-dir",
        required=True,
        help="결과 CSV([typ]_simplified.csv)를 저장할 폴더. 없으면 생성.",
    )
    p.add_argument(
        "--antigen-chain",
        default="A",
        help="항원 체인 기호. AF3에 '위=항원'으로 넣었으면 보통 A.",
    )
    p.add_argument(
        "--antibody-chains",
        nargs="*",
        default=None,
        help=(
            "항체 체인 기호들(공백 구분, 예: --antibody-chains B C). "
            "미지정 시 항원 체인이 아닌 '나머지 전부'를 항체로 간주. "
            "Fab/scFv는 체인이 여러 개일 수 있으니 이 옵션이 중요."
        ),
    )
    p.add_argument(
        "--topn",
        type=int,
        nargs="+",
        default=[1, 3, 5, 10],
        help="상위 N개 평균 피처를 만들 N 목록.",
    )
    # --interface-cutoff 가 '주' 옵션(기본값 보유). --contact-threshold 는 동일
    # 대상(dest=contact_threshold)을 가리키는 별칭이며 기본값을 따로 두지 않아
    # 한쪽만 줘도, 둘 다 줘도(마지막에 적은 값 우선) 일관되게 동작한다.
    p.add_argument(
        "--interface-cutoff",
        dest="contact_threshold",
        type=float,
        default=0.1,
        help=(
            "계면 토큰 판정용 contact_probs 임계값(0~1). 이 값 이상 접촉하는 "
            "항원/항체 토큰을 계면 잔기로 봄. (별칭: --contact-threshold)"
        ),
    )
    p.add_argument(
        "--contact-threshold",
        dest="contact_threshold",
        type=float,
        help=argparse.SUPPRESS,  # --interface-cutoff 와 동일(별칭)
    )
    p.add_argument(
        "--plddt-interface",
        action="store_true",
        help=(
            "켜면 '계면 잔기'에 한정한 pLDDT 피처(plddt_iface_*)도 추가. "
            "주의: 원자->토큰 매핑은 근사라 참고치로만."
        ),
    )
    p.add_argument(
        "--typ-level",
        type=int,
        default=1,
        help=(
            "CSV 분류 라벨(typ) 추출 단계. 0=base-dir 폴더명(전부 한 파일), "
            "1=base-dir 바로 아래 폴더명(원본 기본), 2 이상=그만큼 깊은 폴더명."
        ),
    )
    return p


def main(argv: Optional[Sequence[str]] = None) -> None:
    """CLI 진입점."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if not os.path.isdir(args.base_dir):
        parser.error(f"--base-dir 가 폴더가 아닙니다: {args.base_dir}")

    run(args)


if __name__ == "__main__":
    main()
