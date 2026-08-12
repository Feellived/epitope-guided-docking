#!/usr/bin/env python3
"""UniProt에서 단백질 서열(또는 특정 도메인 구간)을 받아 FASTA로 저장하는 도구.

매번 코드를 고치지 않고, 명령행 인자로 accession·구간·저장이름만 바꿔 쓴다.
다운로드 + 구간 잘라내기 + 저장을 한 번에 처리한다.

예시:
    # CD16a 세포외 도메인(Extracellular 17-208)만
    python scripts/fetch_uniprot_domain.py --accession P08637 --range 17-208 --name CD16a_FCGR3A_ECD

    # 전체 서열 (구간 생략)
    python scripts/fetch_uniprot_domain.py -a P08637

    # 저장 위치/이름 지정
    python scripts/fetch_uniprot_domain.py -a P01857 -r 1-98 -n IgG1_CH1 -o data/

의존성: 표준 라이브러리만 (urllib). 인터넷 연결 필요.
좌표는 UniProt와 동일하게 1-기반·양끝 포함 (예: 17-208 = 17번부터 208번까지).
"""
from __future__ import annotations

import argparse
import os
import sys
import urllib.error
import urllib.request

UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"


def fetch_fasta(accession: str) -> tuple[str, str]:
    """UniProt에서 canonical FASTA를 받아 (헤더, 서열) 로 반환."""
    url = UNIPROT_FASTA.format(acc=accession)
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            text = resp.read().decode("utf-8")
    except urllib.error.HTTPError as e:
        if e.code == 404:
            sys.exit(f"[오류] accession '{accession}' 을 UniProt에서 못 찾음(404). 철자 확인.")
        sys.exit(f"[오류] UniProt 요청 실패 (HTTP {e.code}).")
    except urllib.error.URLError as e:
        sys.exit(f"[오류] 네트워크 실패: {e.reason}. 인터넷 연결 확인.")

    lines = text.splitlines()
    if not lines or not lines[0].startswith(">"):
        sys.exit("[오류] FASTA 형식이 아님. accession 을 확인하세요.")
    header = lines[0]
    seq = "".join(line.strip() for line in lines[1:] if line.strip())
    return header, seq


def parse_range(range_str: str, seq_len: int) -> tuple[int, int]:
    """'17-208' → (17, 208). 1-기반·양끝 포함. 범위 유효성 검사."""
    try:
        start_s, end_s = range_str.split("-")
        start, end = int(start_s), int(end_s)
    except ValueError:
        sys.exit(f"[오류] --range 형식은 'start-end' (예: 17-208). 받은 값: {range_str!r}")
    if not (1 <= start <= end <= seq_len):
        sys.exit(f"[오류] 구간 {start}-{end} 가 서열 길이(1-{seq_len}) 범위를 벗어남.")
    return start, end


def main() -> None:
    ap = argparse.ArgumentParser(description="UniProt 서열/도메인을 FASTA로 저장")
    ap.add_argument("-a", "--accession", required=True, help="UniProt accession (예: P08637)")
    ap.add_argument("-r", "--range", help="잔기 구간, 1-기반 양끝포함 (예: 17-208). 생략 시 전체")
    ap.add_argument("-n", "--name", help="저장 레코드/파일 이름 (생략 시 accession 기반)")
    ap.add_argument("-o", "--out", default="data",
                    help="저장 폴더, 또는 .fasta 파일 경로 (기본: data/)")
    args = ap.parse_args()

    _, seq = fetch_fasta(args.accession)
    print(f"[*] {args.accession} 전체 길이: {len(seq)}")

    if args.range:
        start, end = parse_range(args.range, len(seq))
        sub = seq[start - 1:end]      # 1-기반 → 0-기반 슬라이스 (양끝 포함)
        tag = f"{start}-{end}"
        print(f"[*] 구간 {tag} 길이: {len(sub)}")
    else:
        sub = seq
        tag = "full"
        print("[*] 전체 서열 사용")

    name = args.name or f"{args.accession}_{tag}"

    # 출력 경로 결정: .fasta/.fa 로 끝나면 파일로, 아니면 폴더로 간주
    if args.out.endswith((".fasta", ".fa")):
        out_path = args.out
        out_dir = os.path.dirname(out_path) or "."
    else:
        out_dir = args.out
        out_path = os.path.join(out_dir, f"{name}.fasta")
    os.makedirs(out_dir, exist_ok=True)

    with open(out_path, "w") as f:
        f.write(f">{name} {args.accession} {tag}\n{sub}\n")
    print(f"[저장] {out_path}")


if __name__ == "__main__":
    main()
