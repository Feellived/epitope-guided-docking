#!/usr/bin/env python3
"""[선택 보조] DiscoTope-3.0 으로 항원 표면의 '에피토프다움'을 재서 select_site.py 입력으로 만든다.

왜 이 도구인가 — 우리 설계 철학은 "MSA 말고 **구조·계면**을 보고 고른다"이다. DiscoTope-3.0 은
항원 구조만 보고(역접힘 임베딩 ESM-IF1 기반) 잔기별로 '여기가 항체가 붙을 자리인가'를 낸다.
MSA 를 쓰지 않으므로, MSA 조성이 만든 후보를 고르는 데 써도 **순환 논리가 되지 않는다**.

⚠️ 항원 **단독** 구조를 넣는다(항체를 빼고). 항체가 붙어 있는 복합체를 넣으면 그 계면을 보고
   맞히는 것이 되어 눈먼 선택이 아니다. 기본은 우리 예측물에서 항원 사슬만 뽑아 쓴다
   (--from native 로 실험 구조를 쓸 수도 있으나, 그 경우 '실험 구조를 썼다'를 밝혀야 한다).

═══ 좌표계 주의 (이 파일의 핵심) ═══
DiscoTope 출력의 잔기 번호는 넣어준 PDB 의 번호다. 우리 후보 자리는 **참조서열 0-based 위치**로
적혀 있다. 그대로 맞추면 어긋나므로, 뽑아낸 사슬의 서열을 참조서열에 정렬해 번호를 옮긴다
(msa-depth 의 epitope_recall.posmap 과 같은 방식 — 같은 함수를 그대로 쓴다).

사용(discotope env):
  python scripts/select_score_discotope.py --target 8k3k_D --out results/discotope_8k3k_D.csv
  python scripts/select_score_discotope.py --target 8k3k_D --from native
  python scripts/select_score_discotope.py --target 8k3k_D --pose /경로/model_0.cif
출력: chain_ord,pos,score  (pos = 1-based 참조서열 위치) → select_site.py --method discotope --scores 이 파일
"""
import argparse, csv, glob, json, os, subprocess, sys, tempfile


def load_msad(msad):
    sys.path.insert(0, os.path.join(msad, "scripts"))
    import lib_epitope_recall as ER                      # noqa: E402
    return ER


def pick_pose(msad_data, sites):
    """채점에 쓴 예측물 중 아무거나 하나(항원 사슬 좌표만 필요하다)."""
    base = os.path.join(msad_data, "compreps", "seedrep_cand",
                        sites["model"], sites["target"], sites["depth"])
    ps = sorted(glob.glob(os.path.join(base, "*", "results", "**", "*.cif"), recursive=True))
    return ps[0] if ps else None


def extract_antigen(ER, cj, path, outpdb):
    """구조에서 항원 사슬만 뽑아 PDB 로 저장 → (사슬순번, 실제 잔기 목록, 참조서열) 목록."""
    from Bio.PDB import PDBIO, Select
    m = ER.load(path)
    refs = ER.antigen_refs(cj)
    used, got = set(), []
    for i, ref in enumerate(refs):
        cid, _, rr = ER.best(m, ref, exclude=used)
        if cid is None:
            continue
        used.add(cid); got.append((i, cid, rr, ref))
    if not got:
        sys.exit("!! 항원 사슬을 못 찾았다")

    keep = {cid for _, cid, _, _ in got}

    class S(Select):
        def accept_chain(self, ch):
            return ch.id in keep
        def accept_residue(self, r):
            return "CA" in r

    io = PDBIO(); io.set_structure(m); io.save(outpdb, S())
    return got


def run_discotope(repo, pdb, outdir, struc_type, python_bin):
    script = os.path.join(repo, "src", "predict_webserver.py")
    if not os.path.exists(script):
        sys.exit(f"!! DiscoTope-3.0 을 못 찾았다: {script}\n"
                 f"   설치: bash scripts/run_setup_tools.sh --apply  (또는 --repo 로 경로 지정)")
    cmd = [python_bin, script, "--pdb_or_zip_file", pdb, "--struc_type", struc_type,
           "--out_dir", outdir]
    print("  $ " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2000:]); print(r.stderr[-2000:])
        sys.exit("!! DiscoTope 실행 실패")
    out = glob.glob(os.path.join(outdir, "*.csv"))
    if not out:
        sys.exit(f"!! {outdir} 에 결과 CSV 가 없다")
    return out[0]


def parse_discotope(path):
    """→ {(사슬 id, PDB 잔기번호): 점수}. 열 이름은 판올림마다 조금씩 달라 넉넉히 받는다."""
    rows = list(csv.DictReader(open(path)))
    if not rows:
        sys.exit(f"!! {path} 가 비었다")
    keys = {k.lower(): k for k in rows[0]}
    def find(*cands):
        for c in cands:
            for lk, k in keys.items():
                if c in lk:
                    return k
        return None
    kc, kr = find("chain"), find("res_id", "resid", "residue_number", "position")
    ks = find("discotope-3.0_score", "discotope3_score", "score")
    if not (kr and ks):
        sys.exit(f"!! 열 이름을 못 알아봤다: {list(rows[0])}")
    out = {}
    for r in rows:
        try:
            out[(str(r[kc]).strip() if kc else "", int(float(r[kr])))] = float(r[ks])
        except Exception:
            continue
    print(f"  DiscoTope 잔기 점수 {len(out)}개 (열: chain={kc} res={kr} score={ks})")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--msad", default=os.path.expanduser("~/projects/bk21-msa-depth-bias/pipeline"))
    ap.add_argument("--msad-data", default="/mnt/data/admuser/msadepth")
    ap.add_argument("--repo", default=os.path.expanduser("~/tools/discotope3_web"))
    ap.add_argument("--python", default=sys.executable, help="DiscoTope env 의 python")
    ap.add_argument("--from", dest="src", default="pred", choices=["pred", "native"],
                    help="pred=우리 예측물의 항원(기본·완전한 눈먼 선택) / native=실험 구조")
    ap.add_argument("--pose", default="", help="쓸 구조를 직접 지정")
    ap.add_argument("--struc-type", default="", help="solved|alphafold (생략 시 --from 에 맞춰 자동)")
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    ER = load_msad(a.msad)
    sp = os.path.join(a.msad, "results", f"sites_{a.target}.json")
    if not os.path.exists(sp):
        sys.exit(f"!! {sp} 없음 — msa-depth 에서 먼저: bash analyze_target.sh {a.target}")
    sites = json.load(open(sp))
    cj = json.load(open(os.path.join(a.msad, "targets", a.target, "chains.json")))

    if a.pose:
        src = a.pose
    elif a.src == "native":
        src = os.path.join(a.msad, "targets", a.target, "native.cif")
    else:
        src = pick_pose(a.msad_data, sites)
    if not src or not os.path.exists(src):
        sys.exit("!! 쓸 구조를 못 찾았다 (--pose 로 직접 지정)")
    stype = a.struc_type or ("solved" if a.src == "native" else "alphafold")
    print(f"■ {a.target}   구조 = {src}\n   struc_type = {stype} ({'실험 구조' if a.src=='native' else '우리 예측물'})")
    if a.src == "native":
        print("   ⚠️ 실험 구조를 썼다 — 발표·문서에 그대로 밝힐 것.")

    with tempfile.TemporaryDirectory() as td:
        agpdb = os.path.join(td, f"{a.target}.pdb")
        got = extract_antigen(ER, cj, src, agpdb)
        print(f"   항원 사슬 {len(got)}개 추출 → {os.path.basename(agpdb)}")
        raw = run_discotope(a.repo, agpdb, os.path.join(td, "out"), stype, a.python)
        sc = parse_discotope(raw)

        # PDB 잔기번호 → 참조서열 위치. 사슬 잔기 순서와 서열 정렬로 옮긴다.
        rows = []
        for ordi, cid, rr, ref in got:
            seq = "".join(ER.T3.get(r.get_resname().upper(), "X") for r in rr)
            pm = ER.posmap(seq, ref)            # {사슬 내 index → 참조서열 0-based}
            hit = 0
            for i, r in enumerate(rr):
                if i not in pm:
                    continue
                key = (str(cid), int(r.id[1]))
                v = sc.get(key, sc.get(("", int(r.id[1]))))
                if v is None:
                    continue
                rows.append(dict(chain_ord=ordi, pos=pm[i] + 1, score=round(v, 5)))
                hit += 1
            print(f"   사슬 {cid}(순번 {ordi}): 점수 옮긴 잔기 {hit}/{len(rr)}")

    if not rows:
        sys.exit("!! 옮겨진 점수가 없다 — 사슬 id 나 잔기번호가 안 맞는다")
    out = a.out or f"results/discotope_{a.target}.csv"
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["chain_ord", "pos", "score"])
        w.writeheader(); w.writerows(rows)
    print(f"\n→ {out}  ({len(rows)}잔기)")
    print(f"   다음: python scripts/select_site.py --sites {sp} --method discotope --scores {out}")


if __name__ == "__main__":
    main()
