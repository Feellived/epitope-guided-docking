#!/usr/bin/env python3
"""[Fv 파생] runs/<id>/chains.json 의 항체 사슬을 ANARCI로 V-도메인(Fv)만 잘라
runs/<id>_fv/ 세트를 생성. 항원 서열은 그대로. native 구조는 원본으로 symlink.

목적: Fab 대신 VH+VL(Fv)만 예측시켜 constant-domain/elbow 잡음을 제거한 버전을
      나란히 돌리고, "단순화가 실패 케이스를 살리나"를 진단.

파이프라인: prep_chains.py (먼저) -> prep_fv.py -> make_<model>.py (runs/* 전부, _fv 포함)
필요: anarci (pip install anarci; HMMER 의존). 없으면 명확히 에러.
사용: python prep_fv.py --run-dir runs [--only 9W43,9KKJ]
"""
import argparse, glob, json, os, sys


def fv_truncate(seq):
    """ANARCI로 V-도메인 잔기 범위를 찾아 Fv 서열 반환 (실패 시 None)."""
    try:
        from anarci import anarci
    except ImportError:
        sys.exit("[에러] anarci 미설치. `pip install anarci` (HMMER 필요) 후 재실행.")
    numbered, details, _ = anarci([("q", seq)], scheme="imgt", output=False)
    doms = numbered[0]
    if not doms:                      # 항체로 인식 안 됨
        return None
    # 첫 도메인의 원서열 인덱스 범위 = V-domain (constant domain·꼬리 제거)
    _numbering, start, end = doms[0]
    return seq[start:end + 1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", default="runs")
    ap.add_argument("--only", help="쉼표구분 PDB id (대문자)")
    a = ap.parse_args()
    only = set(x.upper() for x in a.only.split(",")) if a.only else None

    made = skipped = 0
    for cj in sorted(glob.glob(os.path.join(a.run_dir, "*", "chains.json"))):
        rd = os.path.dirname(cj)
        if rd.endswith("_fv"):
            continue                                  # 이미 Fv 세트는 건너뜀
        d = json.load(open(cj))
        tid = d["pdb_id"]
        if only and tid.upper() not in only:
            continue
        antibody = set(d["antibody"] if isinstance(d["antibody"], list) else [d["antibody"]])

        new_chains = []
        ok = True
        for c in d["chains"]:
            if c["id"] in antibody:                   # 항체 사슬 -> Fv로 절단
                fv = fv_truncate(c["seq"])
                if not fv:
                    print(f"[경고] {tid} 사슬 {c['id']}: ANARCI 항체 인식 실패 -> Fv 스킵")
                    ok = False; break
                new_chains.append({**c, "seq": fv})
            else:                                     # 항원 -> 그대로
                new_chains.append(dict(c))
        if not ok:
            skipped += 1; continue

        fv_dir = rd + "_fv"
        os.makedirs(fv_dir, exist_ok=True)
        fv_data = {**d, "pdb_id": tid, "chains": new_chains, "fv": True}
        with open(os.path.join(fv_dir, "chains.json"), "w") as f:
            json.dump(fv_data, f, indent=2, ensure_ascii=False)
        with open(os.path.join(fv_dir, f"{tid.lower()}.fasta"), "w") as f:
            for c in new_chains:
                f.write(f">{c['id']}|{c['role']}\n{c['seq']}\n")
        # native 구조 symlink (eval의 DockQ가 찾도록). 계면은 Fv 안이라 full Fab native와 비교 OK.
        for ext in ("pdb", "cif"):
            for cand in (f"{tid}.{ext}", f"{tid.upper()}.{ext}", f"{tid.lower()}.{ext}"):
                src = os.path.join(rd, cand)
                if os.path.exists(src):
                    dst = os.path.join(fv_dir, cand)
                    if not os.path.exists(dst):
                        os.symlink(os.path.relpath(src, fv_dir), dst)
        oldlen = sum(len(c["seq"]) for c in d["chains"] if c["id"] in antibody)
        newlen = sum(len(c["seq"]) for c in new_chains if c["id"] in antibody)
        print(f"[fv] {tid}: 항체 {oldlen}->{newlen} aa -> {fv_dir}")
        made += 1

    print(f"\n완료: Fv 세트 {made}개 생성, {skipped}개 스킵.")


if __name__ == "__main__":
    main()
