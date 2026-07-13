#!/usr/bin/env python3
"""[어디에 놓았나] 각 모델이 9K6J에서 항체를 항원의 어느 잔기(native 번호)에 놓았는지 덤프.
목적: 실패한 MSA 모델들이 '한 곳(RBD 코어)에 체계적으로 쏠리는가(H1: MSA 오유도)'
      아니면 '흩어지는가(H2: 불안정 도킹)'를 판정.

각 (타깃,포맷,모델)에서 5 pose의 예측 에피토프(항체가 닿은 항원 잔기, native 번호)를
majority-vote(>=절반 pose)로 모아 출력하고, RBD 영역별로 분류한다.

RBD 영역(9K6J 전용 라벨):
  RBM(437-508)  = ACE2/modal 자리(진짜 P5-1C8 에피토프가 여기)
  core(369-386) = CR3022형 보존 코어(공진화 신호 강한 자리 = MSA가 끌 만한 곳)
  class3(337-346,356,440-441,444)
  other         = 그 외

사용: conda activate boltz
      python scripts/where_docked.py --targets 9K6J
결과: 화면 + results/where_docked_<target>.txt
"""
import argparse, glob, json, os
from collections import Counter
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.Align import PairwiseAligner
from Bio.Data.IUPACData import protein_letters_3to1

RUN = "runs_diverse"
POSE_GLOB = {
    "Boltz-2":       "results_boltz/**/*_model_*.cif",
    "Chai-1":        "out_chai/pred.model_idx_*.cif",
    "Protenix-base": "out_protenix/**/*_sample_*.cif",
    "AF2-Multimer":  "results_colabfold/*_unrelaxed_rank_*.pdb",
    "tFold-Ag":      "out_tfold/*seed*.pdb",
}
MODELS = list(POSE_GLOB)
EPI_THR = 8.0

T3TO1 = {k.upper(): v for k, v in protein_letters_3to1.items()}
_al = PairwiseAligner()
_al.mode = "global"; _al.match_score = 1; _al.mismatch_score = -1
_al.open_gap_score = -3; _al.extend_gap_score = -0.5

# RBD 영역 라벨 (9K6J)
def region(rn):
    if 437 <= rn <= 508: return "RBM"
    if (369 <= rn <= 386) or rn in (404, 462, 463, 464): return "core"
    if (337 <= rn <= 346) or rn in (356, 440, 441, 444): return "class3"
    return "other"


def load(path):
    p = MMCIFParser(QUIET=True) if path.endswith(".cif") else PDBParser(QUIET=True)
    return p.get_structure("x", path)[0]


def seq_ca_ids(chain):
    seq, cas, ids = [], [], []
    for res in chain:
        if "CA" not in res: continue
        aa = T3TO1.get(res.get_resname().upper())
        if aa is None: continue
        seq.append(aa); cas.append(res["CA"].get_coord()); ids.append(res.get_id()[1])
    return "".join(seq), (np.array(cas) if cas else np.zeros((0, 3))), ids


def best_match(model, target_seq, exclude=()):
    best = (None, -1e9, "", np.zeros((0, 3)))
    for ch in model:
        if ch.id in exclude: continue
        s, ca, _ = seq_ca_ids(ch)
        if len(s) < 5: continue
        sc = _al.score(s, target_seq) / max(len(s), len(target_seq), 1)
        if sc > best[1]: best = (ch.id, sc, s, ca)
    return best


def aln_map(s_from, s_to):
    if not s_from or not s_to: return {}
    aln = _al.align(s_from, s_to)[0]
    m = {}
    for (a0, a1), (b0, b1) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(a1 - a0): m[a0 + k] = b0 + k
    return m


def contacts_idx(ag_ca, ab_ca):
    if len(ag_ca) == 0 or len(ab_ca) == 0: return set()
    d = np.sqrt(((ag_ca[:, None, :] - ab_ca[None, :, :]) ** 2).sum(-1))
    return set(np.where(d.min(1) <= EPI_THR)[0])


def find_native(rd, T):
    tags = ("boltz", "chai", "protenix", "colabfold", "tfold")
    for f in os.listdir(rd):
        fl = f.lower()
        if fl.endswith((".pdb", ".cif")) and T.lower() in fl and not any(m in fl for m in tags):
            return os.path.join(rd, f)
    for f in os.listdir(rd):
        if f.lower().endswith((".pdb", ".cif")): return os.path.join(rd, f)
    return None


def regbreak(resnums):
    c = Counter(region(r) for r in resnums)
    return f"RBM:{c['RBM']} core:{c['core']} class3:{c['class3']} other:{c['other']}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="9K6J")
    ap.add_argument("--formats", default="Fab,Fv")
    a = ap.parse_args()
    targets = [t.strip().upper() for t in a.targets.split(",")]
    formats = [f.strip() for f in a.formats.split(",")]

    for T in targets:
        out = []
        for fmt in formats:
            low = T.lower() + ("_fv" if fmt == "Fv" else "")
            rd = os.path.join(RUN, low); cj = os.path.join(rd, "chains.json")
            if not os.path.exists(cj): continue
            d = json.load(open(cj))
            ag_id = str(d["antigen"]); ab = d["antibody"]
            ab_ids = [ab] if isinstance(ab, str) else [str(x) for x in ab]
            seqmap = {c["id"]: c["seq"] for c in d.get("chains", [])}
            ag_seq = seqmap.get(ag_id, ""); ab_seqs = [seqmap.get(i, "") for i in ab_ids]

            natf = find_native(rd, T)
            if not natf: continue
            nat = load(natf); nat_ids = [c.id for c in nat]
            if ag_id in nat_ids: nseq, nca, nresnum = seq_ca_ids(nat[ag_id])
            else:
                r = best_match(nat, ag_seq); nseq = r[2]; nca = r[3]; nresnum = list(range(len(nseq)))
            nab_ca = []
            for aid in ab_ids:
                if aid in nat_ids:
                    _, ca, _ = seq_ca_ids(nat[aid]); nab_ca.append(ca)
            nab_ca = np.vstack(nab_ca) if nab_ca else np.zeros((0, 3))
            # native true epitope (auth resnums)
            true_idx = contacts_idx(nca, nab_ca)
            true_res = sorted(nresnum[i] for i in true_idx if i < len(nresnum))
            hdr = f"\n===== {T} {fmt} =====  native 항원 {min(nresnum)}-{max(nresnum)}"
            print(hdr); out.append(hdr)
            tl = f"  TRUE epitope (n={len(true_res)}): {true_res}   [{regbreak(true_res)}]"
            print(tl); out.append(tl)

            for model in MODELS:
                poses = sorted(glob.glob(os.path.join(rd, POSE_GLOB[model]), recursive=True))
                votes = Counter(); npose = 0
                for pose in poses:
                    try:
                        pm = load(pose); npose += 1
                        pid, _, ps, pca = best_match(pm, ag_seq)
                        used = {pid}; pab = []
                        for s in ab_seqs:
                            r = best_match(pm, s, exclude=used)
                            if r[0] is not None: used.add(r[0]); pab.append(r[3])
                        pab_ca = np.vstack(pab) if pab else np.zeros((0, 3))
                        epi_pos = contacts_idx(pca, pab_ca)
                        m = aln_map(ps, nseq)
                        for i in epi_pos:
                            if i in m and m[i] < len(nresnum):
                                votes[nresnum[m[i]]] += 1
                    except Exception:
                        continue
                if npose == 0: continue
                thr = max(1, (npose + 1) // 2)         # majority (>=절반)
                maj = sorted(r for r, v in votes.items() if v >= thr)
                recall = len(set(maj) & set(true_res)) / max(len(true_res), 1)
                line = (f"  {model:14} pred_epitope(>= {thr}/{npose}, n={len(maj)}): {maj}\n"
                        f"                 [{regbreak(maj)}]  true∩recall={recall:.2f}")
                print(line); out.append(line)
        if out:
            os.makedirs("results", exist_ok=True)
            with open(f"results/where_docked_{T}.txt", "w") as fh:
                fh.write("\n".join(out) + "\n")
            print(f"\n# saved -> results/where_docked_{T}.txt")


if __name__ == "__main__":
    main()
