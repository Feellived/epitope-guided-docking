#!/usr/bin/env python3
"""HADDOCK diverse(haddock_diverse/) 포즈 DockQ 채점 — ZDOCK과 동일 렌즈(oracle vs HADDOCK score-pick).
   출력 컬럼 = results/zdock_dockq.csv와 동일:
     target, n_pose, oracle_dockq, oracle_rank, scorepick_dockq, regret, best_score
   → 물리 arm(HADDOCK + ZDOCK)을 results/ 한 폴더에서 나란히 비교.

   방식(자기완결 — lib_structure 의존 없음):
     · 포즈 = haddock_diverse/<t>/run/<stage>/*.pdb.gz (stage 기본 4_emref = 정제 최종 survivor 풀).
     · 항원 사슬 = native 항원 서열(chains.json)과 서열매칭으로 식별(HADDOCK가 A/B 라벨을 바꿔도 안전).
     · DockQ = merged-chain(항원→A, 항체 H+L 병합→B) MAX. ZDOCK/기존 haddock_dockq와 동일 규약.
     · score-pick = HADDOCK score(run/의 tsv, model_name+score, 더 음수=좋음) 1위 pose의 DockQ.
       tsv 없으면 score-pick/regret은 공란(oracle만).

   ⭐ resume/skip: 이미 --out CSV에 있는 타깃은 건너뜀(--force로 재채점). 아직 안 돈 타깃(포즈 없음)도
     건너뛰고 표시만 → HADDOCK 도는 중에 완료된 것만 증분 채점 가능. 타깃마다 CSV flush(중단 안전).

   사용(DockQ 있는 env, 예 boltz):
     python scripts/eval_haddock_diverse.py                       # 전체(configs/targets_benchmark_v2.yaml)
     python scripts/eval_haddock_diverse.py --targets "9K6J"      # 특정 타깃만
     python scripts/eval_haddock_diverse.py --stage 1_rigidbody   # 생성 상한(정제 전) 확인
     python scripts/eval_haddock_diverse.py --force               # 이미 있는 것도 재채점
"""
import argparse, csv, glob, gzip, json, os, re, subprocess, tempfile, difflib
from Bio.PDB import PDBParser, MMCIFParser, PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model as BModel
from Bio.PDB.Chain import Chain as BChain
from Bio.Data.IUPACData import protein_letters_3to1

T3 = {k.upper(): v for k, v in protein_letters_3to1.items()}
DIVERSE = ["9K6J", "9XF9", "9Y0A", "9VXL", "9LY5", "9MQR", "9KKJ", "9W43", "9YC6", "9J4C"]
COLS = ["target", "n_pose", "oracle_dockq", "oracle_rank", "scorepick_dockq", "regret", "best_score"]


def load_any(p):
    h = gzip.open(p, "rt") if p.endswith(".gz") else open(p)
    parser = MMCIFParser(QUIET=True) if p.lower().endswith((".cif", ".cif.gz")) else PDBParser(QUIET=True)
    return parser.get_structure("x", h)[0]


def chain_seq(ch):
    return "".join(T3.get(r.get_resname().upper(), "X") for r in ch if "CA" in r)


def best_ag_chain(model, ag_seq):
    """항원 서열과 가장 잘 맞는 사슬 id (항원 식별; 나머지=항체)."""
    best, bs = None, -1.0
    for ch in model:
        s = chain_seq(ch)
        if not s:
            continue
        r = difflib.SequenceMatcher(None, s, ag_seq).ratio()
        if r > bs:
            bs, best = r, ch.id
    return best


def write_merged(model, ag_ids, ab_ids, out):
    s = Structure("m"); m = BModel(0); s.add(m); A = BChain("A"); Bc = BChain("B"); i = 1; j = 1
    for cid in ag_ids:
        if cid not in model:
            continue
        for r in model[cid]:
            if "CA" not in r:
                continue
            rr = r.copy(); rr.id = (" ", i, " "); A.add(rr); i += 1
    for cid in ab_ids:
        if cid not in model:
            continue
        for r in model[cid]:
            if "CA" not in r:
                continue
            rr = r.copy(); rr.id = (" ", j, " "); Bc.add(rr); j += 1
    if i == 1 or j == 1:
        return False
    m.add(A); m.add(Bc); io = PDBIO(); io.set_structure(s); io.save(out); return True


def dockq(model_pdb, native_pdb):
    try:
        r = subprocess.run(["DockQ", model_pdb, native_pdb], capture_output=True, text=True, timeout=180)
    except Exception:
        return None
    out = (r.stdout or "") + (r.stderr or "")
    vals = re.findall(r"\bDockQ\b[:\s]+([0-9]*\.[0-9]+)", out)
    return max(float(v) for v in vals) if vals else None


def target_info(t, run_dir):
    rd = os.path.join(run_dir, t.lower())
    cjp = os.path.join(rd, "chains.json")
    if not os.path.exists(cjp):
        return None
    cj = json.load(open(cjp))
    ag = str(cj["antigen"]); ab = cj["antibody"]; ab = [ab] if isinstance(ab, str) else [str(x) for x in ab]
    seqmap = {c["id"]: c["seq"] for c in cj["chains"]}
    nat = None
    for f in sorted(os.listdir(rd)):
        low = f.lower()
        if (t.lower() in low and low.endswith((".pdb", ".cif"))
                and not any(m in low for m in ("boltz", "chai", "protenix", "colabfold", "tfold", "antigen", "antibody", "_fv"))):
            nat = os.path.join(rd, f); break
    return dict(rd=rd, ag=ag, ab=ab, ag_seq=seqmap.get(ag, ""), native=nat)


def haddock_scores(run_dir):
    """run/의 tsv(model_name+score)에서 {pose 파일 basename(.gz 제외): score}. 더 음수=좋음."""
    tsv = None
    cands = glob.glob(os.path.join(run_dir, "*clustfcc*", "*.tsv")) + glob.glob(os.path.join(run_dir, "**", "*.tsv"), recursive=True)
    for f in cands:
        try:
            hdr = next(csv.reader(open(f), delimiter="\t"))
            if {"model_name", "score"} <= set(hdr):
                tsv = f; break
        except Exception:
            pass
    if not tsv:
        return {}
    sc = {}
    for row in csv.DictReader(open(tsv), delimiter="\t"):
        try:
            sc[os.path.basename(row["model_name"])] = float(row["score"])
        except (ValueError, KeyError):
            pass
    return sc


def score_target(t, haddock_dir, run_dir, stage, maxn):
    run = os.path.join(haddock_dir, t.lower(), "run")
    poses = sorted(glob.glob(os.path.join(run, stage, "*.pdb.gz")) + glob.glob(os.path.join(run, stage, "*.pdb")))
    if not poses:
        return None, "포즈없음(미실행)"
    ni = target_info(t, run_dir)
    if ni is None or ni["native"] is None:
        return None, "native/chains.json 없음"
    scores = haddock_scores(run)
    # HADDOCK score 순으로 정렬(있으면), 없으면 파일순. score-pick = 정렬 1위.
    def key(p):
        return scores.get(os.path.basename(p).replace(".gz", ""), float("inf"))
    ordered = sorted(poses, key=key)
    if maxn > 0:
        ordered = ordered[:maxn]
    with tempfile.TemporaryDirectory() as td:
        natm = os.path.join(td, "nat.pdb")
        if not write_merged(load_any(ni["native"]), [ni["ag"]], ni["ab"], natm):
            return None, "native merge 실패"
        best_q, best_i, sp_q, n = -1.0, -1, None, 0
        for k, pose in enumerate(ordered):
            try:
                pm = load_any(pose)
                agc = best_ag_chain(pm, ni["ag_seq"])
                abc = [ch.id for ch in pm if ch.id != agc]
                mp = os.path.join(td, "p.pdb")
                if not write_merged(pm, [agc], abc, mp):
                    continue
                q = dockq(mp, natm)
                if q is None:
                    continue
                n += 1
                if q > best_q:
                    best_q, best_i = q, k + 1
                if k == 0:
                    sp_q = q
            except Exception:
                continue
    if best_q < 0:
        return None, "채점된 pose 없음"
    bs = min(scores.values()) if scores else ""
    reg = round(best_q - sp_q, 3) if sp_q is not None else ""
    row = dict(target=t, n_pose=n, oracle_dockq=round(best_q, 3), oracle_rank=best_i,
               scorepick_dockq=round(sp_q, 3) if sp_q is not None else "", regret=reg,
               best_score=round(bs, 2) if bs != "" else "")
    return row, None


def load_done(out):
    if not os.path.exists(out):
        return {}
    d = {}
    for r in csv.DictReader(open(out)):
        if r.get("oracle_dockq", "") != "":
            d[r["target"].upper()] = r
    return d


def flush(out, done):
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS); w.writeheader()
        for t in sorted(done, key=lambda x: -float(done[x]["oracle_dockq"])):
            w.writerow(done[t])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--haddock-dir", default="haddock_diverse")
    ap.add_argument("--run-dir", default="runs_diverse")
    ap.add_argument("--stage", default="4_emref")
    ap.add_argument("--targets", default=" ".join(DIVERSE))
    ap.add_argument("--max", type=int, default=0, help=">0이면 score 상위 N개만 채점(oracle 상한 낮아질 수 있음). 0=전체.")
    ap.add_argument("--out", default="results/haddock_dockq_diverse.csv")
    ap.add_argument("--force", action="store_true", help="이미 CSV에 있는 타깃도 재채점")
    a = ap.parse_args()
    if subprocess.run(["which", "DockQ"], capture_output=True).returncode != 0:
        raise SystemExit("!! DockQ 없음 — DockQ 있는 env(예: boltz)에서 실행")

    done = {} if a.force else load_done(a.out)
    if done:
        print(f"[resume] 이미 채점된 {len(done)}타깃 skip: {', '.join(sorted(done))}")
    print(f"[HADDOCK diverse DockQ] stage={a.stage}  dir={a.haddock_dir}  (oracle vs HADDOCK score-pick)")
    print(f"{'target':8}{'n_pose':>7}{'oracle(rank)':>16}{'score-pick':>12}{'regret':>8}")
    print("-" * 55)
    for t in a.targets.split():
        t = t.upper()
        if not a.force and t in done:
            continue
        row, why = score_target(t, a.haddock_dir, a.run_dir, a.stage, a.max)
        if row is None:
            print(f"{t:8} skip — {why}"); continue
        done[t] = {k: str(v) for k, v in row.items()}
        flush(a.out, done)  # 타깃마다 저장(중단 안전)
        print(f"{t:8}{row['n_pose']:>7}{row['oracle_dockq']:>10}(r{row['oracle_rank']}){str(row['scorepick_dockq']):>12}{str(row['regret']):>8}")
    print(f"\n→ {a.out} ({len(done)}행). oracle≫score-pick = HADDOCK score도 pose 선택 실패(=ZDOCK·ipTM과 같은 '생성 O·선택 X').")


if __name__ == "__main__":
    main()
