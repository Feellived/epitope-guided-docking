#!/usr/bin/env python3
"""[항체-쪽 가설 진단] co-folder 도킹 실패가 '항원 에피토프(자리) 선택' 문제인지
'항체 모델링/배치' 문제인지 분해한다. (2026-07-13, MSA-free 이점 가설 검증용)

각 (타깃, 포맷, 모델)에 대해 pose들을 훑어 두 가지를 따로 잰다:
  - epitope recall : 예측 pose에서 '항체가 닿은 항원 잔기'가 native 에피토프를 얼마나 회복하나.
                     1.0에 가까울수록 = "올바른 자리(RBM 등)는 찾았다".
  - antibody RMSD  : 예측한 항체(Fv/Fab)만 native 항체에 중첩했을 때의 CA-RMSD.
                     작을수록 = "항체 자체(CDR 포함)는 잘 지었다". 배치와 무관하게 fold 품질만 봄.

해석 규칙:
  recall 높음 + DockQ 실패  -> 자리는 맞고 orientation/항체 문제  = 항체-쪽 가설 지지 (Fork B)
  recall 낮음               -> 엉뚱한 자리                       = 에피토프 선택 문제 (Fork A)
  MSA 모델 recall 높은데 실패 + tFold 성공 -> MSA-free 이점(항체-쪽) 강화

사용: conda activate boltz            # Bio.PDB 있는 env
      python scripts/diagnose_antibody_side.py --targets 9K6J
      python scripts/diagnose_antibody_side.py --targets 9K6J,9LY5,9KKJ --formats Fab,Fv
결과: 화면 표 + results/diagnose_<target>.txt (복붙 실패 대비, 파일로도 저장)
"""
import argparse, glob, json, os
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.Align import PairwiseAligner
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.Data.IUPACData import protein_letters_3to1

RUN = os.environ.get("RUN_DIR", "runs_diverse")
POSE_GLOB = {
    "Boltz-2":       "*boltz/**/*_model_*.cif",
    "Chai-1":        "out_chai/pred.model_idx_*.cif",
    "Protenix-base": "out_protenix/**/*_sample_*.cif",
    "AF2-Multimer":  "*colabfold/*_unrelaxed_rank_*.pdb",
    "tFold-Ag":      "out_tfold/*seed*.pdb",
}
MODELS = list(POSE_GLOB)
EPI_THR = 8.0   # CA-CA Å: region 수준 에피토프 정의(측쇄 여유 포함)

T3TO1 = {k.upper(): v for k, v in protein_letters_3to1.items()}
_al = PairwiseAligner()
_al.mode = "global"; _al.match_score = 1; _al.mismatch_score = -1
_al.open_gap_score = -3; _al.extend_gap_score = -0.5


def load(path):
    p = MMCIFParser(QUIET=True) if path.endswith(".cif") else PDBParser(QUIET=True)
    return p.get_structure("x", path)[0]


def seq_ca(chain):
    seq, cas = [], []
    for res in chain:
        if "CA" not in res:
            continue
        aa = T3TO1.get(res.get_resname().upper())
        if aa is None:
            continue
        seq.append(aa); cas.append(res["CA"].get_coord())
    return "".join(seq), (np.array(cas) if cas else np.zeros((0, 3)))


def best_match(model, target_seq, exclude=()):
    """모델 사슬 중 target_seq와 정규화 정체성 최고인 사슬 -> (id, score, seq, ca)."""
    best = (None, -1e9, "", np.zeros((0, 3)))
    for ch in model:
        if ch.id in exclude:
            continue
        s, ca = seq_ca(ch)
        if len(s) < 5:
            continue
        sc = _al.score(s, target_seq) / max(len(s), len(target_seq), 1)
        if sc > best[1]:
            best = (ch.id, sc, s, ca)
    return best


def aln_index_map(s_from, s_to):
    """s_from 잔기 index -> s_to 잔기 index (CA 순서 기준)."""
    if not s_from or not s_to:
        return {}
    aln = _al.align(s_from, s_to)[0]
    m = {}
    for (a0, a1), (b0, b1) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(a1 - a0):
            m[a0 + k] = b0 + k
    return m


def epitope(ag_ca, ab_ca):
    if len(ag_ca) == 0 or len(ab_ca) == 0:
        return set()
    d = np.sqrt(((ag_ca[:, None, :] - ab_ca[None, :, :]) ** 2).sum(-1))
    return set(np.where(d.min(1) <= EPI_THR)[0])


def find_native(rd, T):
    tags = ("boltz", "chai", "protenix", "colabfold", "tfold")
    for f in os.listdir(rd):
        fl = f.lower()
        if fl.endswith((".pdb", ".cif")) and T.lower() in fl and not any(m in fl for m in tags):
            return os.path.join(rd, f)
    for f in os.listdir(rd):
        if f.lower().endswith((".pdb", ".cif")):
            return os.path.join(rd, f)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="9K6J")
    ap.add_argument("--formats", default="Fab,Fv")
    a = ap.parse_args()
    targets = [t.strip().upper() for t in a.targets.split(",")]
    formats = [f.strip() for f in a.formats.split(",")]

    hdr = f"{'target':6} {'fmt':4} {'model':14} {'epi_recall':>10} {'abRMSD@best':>12} {'min_abRMSD':>11} {'poses':>6}"
    print(hdr)
    print("-" * len(hdr))
    for T in targets:
        out_lines = []
        for fmt in formats:
            low = T.lower() + ("_fv" if fmt == "Fv" else "")
            rd = os.path.join(RUN, low)
            cj = os.path.join(rd, "chains.json")
            if not os.path.exists(cj):
                continue
            d = json.load(open(cj))
            ag_id = str(d["antigen"])
            ab = d["antibody"]; ab_ids = [ab] if isinstance(ab, str) else [str(x) for x in ab]
            seqmap = {c["id"]: c["seq"] for c in d.get("chains", [])}
            ag_seq = seqmap.get(ag_id, "")
            ab_seqs = [seqmap.get(i, "") for i in ab_ids]

            natf = find_native(rd, T)
            if not natf:
                continue
            nat = load(natf)
            nat_ids = [c.id for c in nat]
            nag_seq, nag_ca = (seq_ca(nat[ag_id]) if ag_id in nat_ids else best_match(nat, ag_seq)[2:])
            nab = []
            for aid in ab_ids:
                if aid in nat_ids:
                    nab.append(seq_ca(nat[aid]))
            if not nab:
                for s in ab_seqs:
                    r = best_match(nat, s); nab.append((r[2], r[3]))
            nat_ab_ca = np.vstack([c for _, c in nab]) if nab else np.zeros((0, 3))
            epi_nat = epitope(nag_ca, nat_ab_ca)
            if not epi_nat:
                print(f"# {T} {fmt}: native 에피토프 계산 실패(스킵)")
                continue

            for model in MODELS:
                poses = sorted(glob.glob(os.path.join(rd, POSE_GLOB[model]), recursive=True))
                best_recall, rmsd_at_best, min_rmsd, npose = -1.0, float("nan"), float("inf"), 0
                for pose in poses:
                    try:
                        pm = load(pose); npose += 1
                        pid, _, ps, pca = best_match(pm, ag_seq)
                        used = {pid}; pab = []
                        for s in ab_seqs:
                            r = best_match(pm, s, exclude=used)
                            if r[0] is not None:
                                used.add(r[0]); pab.append((r[2], r[3]))
                        pab_ca = np.vstack([c for _, c in pab]) if pab else np.zeros((0, 3))
                        # (1) epitope recall: pred 항원 idx -> native 항원 idx
                        epi_pred = epitope(pca, pab_ca)
                        m = aln_index_map(ps, nag_seq)
                        epi_pred_nat = {m[i] for i in epi_pred if i in m}
                        recall = len(epi_pred_nat & epi_nat) / len(epi_nat)
                        # (2) antibody CA-RMSD: pred 항체 -> native 항체 중첩(배치 무관, fold만)
                        P, Q = [], []
                        for (ps_seq, ps_ca), (ns_seq, ns_ca) in zip(pab, nab):
                            mm = aln_index_map(ps_seq, ns_seq)
                            for pi, ni in mm.items():
                                if pi < len(ps_ca) and ni < len(ns_ca):
                                    P.append(ps_ca[pi]); Q.append(ns_ca[ni])
                        rms = float("nan")
                        if len(P) >= 5:
                            sup = SVDSuperimposer(); sup.set(np.array(Q), np.array(P)); sup.run()
                            rms = sup.get_rms()
                        if recall > best_recall:
                            best_recall, rmsd_at_best = recall, rms
                        if not np.isnan(rms):
                            min_rmsd = min(min_rmsd, rms)
                    except Exception:
                        continue
                if npose == 0:
                    continue
                rab = f"{rmsd_at_best:.2f}" if not np.isnan(rmsd_at_best) else "-"
                mr = f"{min_rmsd:.2f}" if min_rmsd != float("inf") else "-"
                line = f"{T:6} {fmt:4} {model:14} {best_recall:>10.2f} {rab:>12} {mr:>11} {npose:>6}"
                print(line); out_lines.append(line)
        if out_lines:
            os.makedirs("results", exist_ok=True)
            with open(f"results/diagnose_{T}.txt", "w") as fh:
                fh.write(hdr + "\n" + "\n".join(out_lines) + "\n")
            print(f"# saved -> results/diagnose_{T}.txt")
        print()


if __name__ == "__main__":
    main()
