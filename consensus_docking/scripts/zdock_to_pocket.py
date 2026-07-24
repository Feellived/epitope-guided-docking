#!/usr/bin/env python3
"""ZDOCK 물리 pose → co-folder guided 입력 (diverse 세트). HADDOCK-guided(haddock_to_*)의 ZDOCK·diverse판.

ZDOCK top-rank pose(들)의 항원 계면 잔기(=후보 에피토프)를 co-folder 계면 제약으로 주입:
  - boltz   : pocket YAML  (blind=make_boltz + --use_msa_server 와 동일, constraint 블록만 추가 = 순수 힌트 효과 격리)
  - protenix: pocket JSON  (blind=make_protenix + 서버 MSA, constraint.pocket만 추가)
  - tfold   : ICF .pt(epitope-only) + FASTA(H/L/A)  (blind tFold + --icf --model_version ppi)

⚠️ 힌트 오염: ZDOCK score-pick(=rank1 pose)이 틀릴 수 있음(9kkj oracle=rank66·9y0a=rank3). 그래서
   --topn N(상위 N pose) × [--union(합쳐 넓게) | 기본(pose별 패치)]로 '압축 vs 넓게'를 직접 실험한다.
   top-1이 틀리면 guided도 틀린 자리로 감 → 압축 실험의 핵심 관전 포인트.

입력: zdock/<t>/complex.* (create.pl 산출; chain A=항원 1..N 재넘버, chain B=항체) + runs_diverse/<t>/chains.json
출력: runs_zdock_guided/<t>/{boltz_<t>_patchK.yaml, protenix_<t>_patchK.json, tfold_<t>_patchK.pt, tfold_<t>.fasta}

사용:
  python scripts/zdock_to_pocket.py --targets "9KKJ" --emit boltz                      # 스모크(Boltz만)
  python scripts/zdock_to_pocket.py --targets "9J4C 9KKJ 9Y0A 9W43" --emit boltz,protenix,tfold --topn 1
  python scripts/zdock_to_pocket.py --targets "9J4C 9KKJ 9Y0A 9W43" --topn 10 --union   # 넓게(top10 합집합)
필요: biopython. (tfold emit 시 torch 추가.)
"""
import argparse, glob, json, os, re
import numpy as np
from Bio.PDB import PDBParser
from Bio.Data.IUPACData import protein_letters_3to1

T3 = {k.upper(): v for k, v in protein_letters_3to1.items()}


def load_pdb(p):
    return PDBParser(QUIET=True).get_structure("x", p)[0]


def chain_seq(ch):
    return "".join(T3.get(r.get_resname().upper(), "X") for r in ch if "CA" in r)


def best_ag_chain(model, ag_seq):
    import difflib
    best, bs = None, -1.0
    for ch in model:
        s = chain_seq(ch)
        if not s:
            continue
        r = difflib.SequenceMatcher(None, s, ag_seq).ratio()
        if r > bs:
            bs, best = r, ch.id
    return best


def antigen_epitope(pose, ag_seq, extract_dist):
    """complex.* pose의 항원 계면 잔기 → 항원 서열 1-based 위치 집합.
    항원 사슬은 서열매칭으로 식별, 잔기순서 = 서열순서(zdock_prep가 1..N 재넘버) → index+1 = seq 위치."""
    agc = best_ag_chain(pose, ag_seq)
    if agc is None:
        return set()
    ag_res = [r for r in pose[agc] if "CA" in r]
    ab_ca = [r["CA"].coord for ch in pose if ch.id != agc for r in ch if "CA" in r]
    if not ab_ca or not ag_res:
        return set()
    ab_ca = np.array(ab_ca)
    epi = set()
    for i, r in enumerate(ag_res):
        if np.linalg.norm(ab_ca - r["CA"].coord, axis=1).min() <= extract_dist:
            epi.add(i + 1)                       # 1-based 서열 위치
    return epi


def zdock_poses(zdir, t, topn):
    """zdock/<t>/complex.N[.pdb] 를 순위순(complex.1 = score 1위)으로 top-N. 확장자 유무 모두 대응."""
    od = os.path.join(zdir, t.lower())
    out = []
    for p in glob.glob(os.path.join(od, "complex.*")):
        m = re.search(r"complex\.(\d+)", os.path.basename(p))   # 'complex.10.pdb' → 10
        if m:
            out.append((int(m.group(1)), p))
    out.sort()
    poses = [p for _, p in out]
    return poses[:topn] if topn > 0 else poses


# ── emitters ──
def emit_boltz(outd, t, chains, ag, binder, positions, md, soft, k):
    contacts = "\n".join(f"        - [{ag}, {r}]" for r in sorted(positions))
    L = ["version: 1", "sequences:"]
    for c in chains:
        L.append(f'  - protein: {{ id: {c["id"]}, sequence: "{c["seq"]}" }}')
    L += ["constraints:", "  - pocket:", f"      binder: {binder}", "      contacts:", contacts,
          f"      max_distance: {md}", f"      force: {'false' if soft else 'true'}"]
    open(os.path.join(outd, f"boltz_{t.lower()}_patch{k}.yaml"), "w").write("\n".join(L) + "\n")


def emit_protenix(outd, t, chains, ag, binder, positions, md, k):
    ids = [c["id"] for c in chains]
    ag_ent = ids.index(ag) + 1
    binder_ent = ids.index(binder) + 1
    job = {"name": t.lower(),
           "sequences": [{"proteinChain": {"sequence": c["seq"], "count": 1}} for c in chains],
           "constraint": {"pocket": {
               "binder_chain": {"entity": binder_ent, "copy": 1},
               "contact_residues": [{"entity": ag_ent, "copy": 1, "position": int(r)} for r in sorted(positions)],
               "max_distance": md}}}
    json.dump([job], open(os.path.join(outd, f"protenix_{t.lower()}_patch{k}.json"), "w"), indent=2)


def emit_tfold(outd, t, seqmap, ag, ab, positions, k, wrote_fasta):
    import torch
    L_H = len(seqmap[ab[0]]); L_L = len(seqmap[ab[1]]) if len(ab) >= 2 else 0; L_A = len(seqmap[ag])
    L_ab = L_H + L_L
    if not wrote_fasta[0]:
        recs = [("H", seqmap[ab[0]])] + ([("L", seqmap[ab[1]])] if L_L else []) + [("A", seqmap[ag])]
        with open(os.path.join(outd, f"tfold_{t.lower()}.fasta"), "w") as f:
            for cid, seq in recs:
                f.write(f">{cid}\n{seq}\n")
        wrote_fasta[0] = True
    icf = torch.zeros(L_ab + L_A, dtype=torch.int8)
    for r in sorted(positions):
        i = int(r) - 1
        if 0 <= i < L_A:
            icf[L_ab + i] = 1
    torch.save(icf, os.path.join(outd, f"tfold_{t.lower()}_patch{k}.pt"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", required=True, help="쉼표/공백 구분 PDB id (예: '9J4C 9KKJ')")
    ap.add_argument("--zdir", default="zdock")
    ap.add_argument("--run-dir", default="runs_diverse")
    ap.add_argument("--emit", default="boltz,protenix,tfold", help="쉼표: boltz,protenix,tfold 중 선택")
    ap.add_argument("--topn", type=int, default=1, help="상위 N pose(순위=ZDOCK score). 0=전체")
    ap.add_argument("--union", action="store_true", help="추출 패치들을 하나로 합쳐 단일 제약(넓게)")
    ap.add_argument("--max-patches", type=int, default=1, help="union 아닐 때 pose별 패치 최대 개수")
    ap.add_argument("--extract-dist", type=float, default=10.0, help="에피토프 추출: 항원 CA가 항체 CA와 이 거리(Å) 내")
    ap.add_argument("--max-distance", type=float, default=6.0, help="co-folder pocket 제약 거리(Å)")
    ap.add_argument("--min-contacts", type=int, default=4)
    ap.add_argument("--soft", action="store_true", help="Boltz force:false(소프트). 기본=force:true")
    ap.add_argument("--outdir", default="runs_zdock_guided")
    a = ap.parse_args()
    emit = [e.strip() for e in a.emit.replace(",", " ").split()]
    tgts = [x.upper() for x in a.targets.replace(",", " ").split()]

    for t in tgts:
        cjp = os.path.join(a.run_dir, t.lower(), "chains.json")
        if not os.path.exists(cjp):
            print(f"{t}: chains.json 없음 skip"); continue
        cj = json.load(open(cjp))
        ag = str(cj["antigen"]); ab = cj["antibody"]; ab = [ab] if isinstance(ab, str) else [str(x) for x in ab]
        binder = str(ab[0])
        seqmap = {c["id"]: c["seq"] for c in cj["chains"]}
        ag_seq = seqmap.get(ag, "")
        poses = zdock_poses(a.zdir, t, a.topn)
        if not poses:
            print(f"{t}: zdock/{t.lower()}/complex.* 없음 (ZDOCK run 먼저) skip"); continue

        patches = []
        for p in poses:
            try:
                epi = antigen_epitope(load_pdb(p), ag_seq, a.extract_dist)
            except Exception as e:
                print(f"  {t} {os.path.basename(p)} 파싱 실패: {e}"); continue
            if len(epi) < a.min_contacts:
                continue
            if any(len(epi & q) / max(1, len(epi | q)) > 0.8 for q in patches):
                continue                          # dedup Jaccard>0.8
            patches.append(epi)
            if not a.union and len(patches) >= a.max_patches:
                break
        if not patches:
            print(f"{t}: 유효 패치 없음(min-contacts={a.min_contacts}) skip"); continue
        if a.union:
            patches = [set().union(*patches)]

        outd = os.path.join(a.outdir, t.lower()); os.makedirs(outd, exist_ok=True)
        wrote_fasta = [False]
        for k, pos in enumerate(patches):
            if "boltz" in emit:
                emit_boltz(outd, t, cj["chains"], ag, binder, pos, a.max_distance, a.soft, k)
            if "protenix" in emit:
                emit_protenix(outd, t, cj["chains"], ag, binder, pos, a.max_distance, k)
            if "tfold" in emit:
                emit_tfold(outd, t, seqmap, ag, ab, pos, k, wrote_fasta)
        src = f"top{a.topn} pose" + (" union" if a.union else "")
        sizes = ",".join(str(len(p)) for p in patches)
        print(f"{t}: {src} → 패치 {len(patches)}개(잔기수 {sizes}) emit={emit} "
              f"(binder={binder}, 항원={ag}, force={'false' if a.soft else 'true'}) → {outd}/")


if __name__ == "__main__":
    main()
