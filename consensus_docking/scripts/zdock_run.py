#!/usr/bin/env python3
"""[ZDOCK 2/2] 실행 + DockQ 채점 (diverse Fab, native-split bound).
  zdock/<id>/ (zdock_prep.py 산출) 별로:
    mark_sur(항체·항원) → block.pl(비-CDR, 항체) → zdock -R 항체 -L 항원 → create.pl → 상위 pose DockQ vs native.
  DockQ = merged-chain(항원->A, 항체->B, 재넘버), MAX. HADDOCK(haddock_dockq.py)과 동일 규약.
  출력:
    results/zdock_dockq.csv        — 타깃별 요약 (oracle vs score-pick = ipTM 대체 렌즈)
    results/zdock_poses_<id>.csv   — 포즈별 (rank, zdock_score, dockq)  ← Phase 1 재랭커 학습용

환경변수:
  ZDOCK_HOME (기본 ~/tools/zdock3.0.2_linux_x64) · LIBG2C (기본 ~/tools/lib, libg2c.so.0 위치)
  ZDOCK_FLAGS (기본 "-F -N 100": -F=수용체 고정(블로킹시), -N=출력 포즈수)
실행(DockQ 있는 env, 예 boltz):  python scripts/zdock_run.py [--only 9K6J] [--runs runs_diverse] [--zdir zdock]
필요: DockQ CLI, biopython. mark_sur/zdock/create.pl은 서브프로세스(LD_LIBRARY_PATH 자동 주입).
"""
import argparse, csv, glob, json, os, re, shlex, subprocess, tempfile
from Bio.PDB import PDBParser, MMCIFParser, PDBIO
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model as BModel
from Bio.PDB.Chain import Chain as BChain

ZHOME = os.environ.get("ZDOCK_HOME", os.path.expanduser("~/tools/zdock3.0.2_linux_x64"))
LIBG2C = os.environ.get("LIBG2C", os.path.expanduser("~/tools/lib"))
ZFLAGS = shlex.split(os.environ.get("ZDOCK_FLAGS", "-F -N 100"))
ENV = {**os.environ, "LD_LIBRARY_PATH": LIBG2C + ":" + os.environ.get("LD_LIBRARY_PATH", "")}

def load(p):
    return (MMCIFParser(QUIET=True) if p.endswith(".cif") else PDBParser(QUIET=True)).get_structure("x", p)[0]

def write_merged(model, ag_ids, ab_ids, out):
    s = Structure("m"); m = BModel(0); s.add(m); A = BChain("A"); Bc = BChain("B"); i = 1; j = 1
    for cid in ag_ids:
        if cid not in model: continue
        for r in model[cid]:
            if "CA" not in r: continue
            rr = r.copy(); rr.id = (" ", i, " "); A.add(rr); i += 1
    for cid in ab_ids:
        if cid not in model: continue
        for r in model[cid]:
            if "CA" not in r: continue
            rr = r.copy(); rr.id = (" ", j, " "); Bc.add(rr); j += 1
    if i == 1 or j == 1: return False
    m.add(A); m.add(Bc); io = PDBIO(); io.set_structure(s); io.save(out); return True

def dockq(model_pdb, native_pdb):
    try:
        r = subprocess.run(["DockQ", model_pdb, native_pdb], capture_output=True, text=True, timeout=180)
    except Exception:
        return None
    out = (r.stdout or "") + (r.stderr or "")
    vals = re.findall(r"\bDockQ\b[:\s]+([0-9]*\.[0-9]+)", out)
    return max(float(v) for v in vals) if vals else None

def sh(cmd, **kw):
    return subprocess.run(cmd, env=ENV, capture_output=True, text=True, **kw)

def native_merged(rd, tid, td):
    cj = json.load(open(os.path.join(rd, "chains.json")))
    ag = str(cj["antigen"]); ab = cj["antibody"]; ab = [ab] if isinstance(ab, str) else [str(x) for x in ab]
    nat = None
    for f in sorted(os.listdir(rd)):
        low = f.lower()
        if (tid.lower() in low and low.endswith((".pdb", ".cif"))
                and not any(m in low for m in ("boltz", "chai", "protenix", "colabfold", "tfold", "antigen", "antibody", "_fv"))):
            nat = os.path.join(rd, f); break
    if not nat: return None
    out = os.path.join(td, "nat.pdb")
    return out if write_merged(load(nat), [ag], ab, out) else None

def parse_scores(zout):
    """zdock.out: 헤더 뒤 각 예측줄 마지막 컬럼 = ZDOCK score, 순서=순위(최고 먼저)."""
    scores = []
    for line in open(zout):
        p = line.split()
        if len(p) >= 7:
            try:
                scores.append(float(p[-1]))
            except ValueError:
                pass
    return scores

def run_target(tid, zdir, rd, out_dir):
    od = os.path.join(zdir, tid)
    ab = os.path.join(od, "antibody.pdb"); ag = os.path.join(od, "antigen.pdb"); blk = os.path.join(od, "block_noncdr.txt")
    if not all(os.path.exists(x) for x in (ab, ag, blk)):
        print(f"  {tid}: prep 산출 없음 (zdock_prep.py 먼저)"); return None
    ab_m = os.path.join(od, "antibody_m.pdb"); ag_m = os.path.join(od, "antigen_m.pdb")
    sh([os.path.join(ZHOME, "mark_sur"), ab, ab_m])
    sh([os.path.join(ZHOME, "mark_sur"), ag, ag_m])
    if not (os.path.exists(ab_m) and os.path.exists(ag_m)):
        print(f"  {tid}: mark_sur 실패 (libg2c LD_LIBRARY_PATH 확인)"); return None
    ab_b = os.path.join(od, "antibody_block.pdb")
    r = sh([os.path.join(ZHOME, "block.pl"), ab_m, blk]); open(ab_b, "w").write(r.stdout)
    if os.path.getsize(ab_b) == 0:
        print(f"  {tid}: block.pl 출력 비어있음\n{r.stderr[:200]}"); return None
    zout = os.path.join(od, "zdock.out")
    zr = sh([os.path.join(ZHOME, "zdock"), "-R", ab_b, "-L", ag_m, *ZFLAGS, "-o", zout])
    if not os.path.exists(zout) or os.path.getsize(zout) == 0:
        print(f"  {tid}: zdock 실패\n  cmd flags={ZFLAGS}\n{zr.stdout[-300:]}\n{zr.stderr[-300:]}"); return None
    for f in ("create_lig", "uniCHARMM", "create.pl"):     # create.pl은 cwd에 이 파일들을 요구
        lk = os.path.join(od, f)
        if not os.path.exists(lk):
            os.symlink(os.path.join(ZHOME, f), lk)
    sh(["./create.pl", os.path.basename(zout)], cwd=od)
    poses = sorted(glob.glob(os.path.join(od, "complex.*")),   # create.pl 산출 complex.1, complex.2, ...
                   key=lambda p: int((re.findall(r"\d+", os.path.basename(p)) or [0])[0]))
    scores = parse_scores(zout)
    if not poses:
        print(f"  {tid}: create.pl pose 없음"); return None

    rows = []
    with tempfile.TemporaryDirectory() as td:
        natm = native_merged(rd, tid, td)
        if natm is None:
            print(f"  {tid}: native merge 실패"); return None
        best_q = -1.0; best_i = -1; sp_q = None
        for k, pose in enumerate(poses):
            try:
                pm = os.path.join(td, "p.pdb")
                if not write_merged(load(pose), ["A"], ["B"], pm):
                    continue
                q = dockq(pm, natm)
                if q is None:
                    continue
                sc = scores[k] if k < len(scores) else ""
                rows.append({"rank": k + 1, "zdock_score": sc, "dockq": round(q, 3)})
                if q > best_q:
                    best_q = q; best_i = k + 1
                if k == 0:
                    sp_q = q
            except Exception:
                continue
    if best_q < 0:
        print(f"  {tid}: 채점된 pose 없음"); return None

    with open(os.path.join(out_dir, f"zdock_poses_{tid}.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["rank", "zdock_score", "dockq"]); w.writeheader(); w.writerows(rows)
    reg = round(best_q - sp_q, 3) if sp_q is not None else ""
    print(f"  {tid:8} pose {len(rows):3}  oracle {best_q:.3f}(rank{best_i})  score-pick {sp_q if sp_q is None else round(sp_q,3)}  regret {reg}")
    return dict(target=tid, n_pose=len(rows), oracle_dockq=round(best_q, 3), oracle_rank=best_i,
                scorepick_dockq=round(sp_q, 3) if sp_q is not None else "", regret=reg,
                best_score=round(scores[0], 2) if scores else "")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default="runs_diverse")
    ap.add_argument("--zdir", default="zdock")
    ap.add_argument("--out", default="results/zdock_dockq.csv")
    ap.add_argument("--only", help="쉼표구분 PDB id")
    a = ap.parse_args()
    if subprocess.run(["which", "DockQ"], capture_output=True).returncode != 0:
        raise SystemExit("!! DockQ 없음 — DockQ 있는 env(예: boltz)에서 실행")
    only = set(x.upper() for x in a.only.split(",")) if a.only else None
    out_dir = os.path.dirname(a.out) or "."; os.makedirs(out_dir, exist_ok=True)
    print(f"ZDOCK_HOME={ZHOME}  flags={ZFLAGS}")
    print(f"{'target':8} {'pose':>4}  {'oracle(rank)':>14}  {'score-pick':>10}  regret")
    print("-" * 60)
    rows = []
    for od in sorted(glob.glob(os.path.join(a.zdir, "*"))):
        tid = os.path.basename(od)
        if not os.path.isdir(od):
            continue
        if only and tid.upper() not in only:
            continue
        rd = os.path.join(a.runs, tid)
        if not os.path.isdir(rd):
            print(f"  {tid}: runs_diverse/{tid} 없음"); continue
        r = run_target(tid, a.zdir, rd, out_dir)
        if r:
            rows.append(r)
    with open(a.out, "w", newline="") as f:
        cols = ["target", "n_pose", "oracle_dockq", "oracle_rank", "scorepick_dockq", "regret", "best_score"]
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
    print(f"\n→ {a.out} ({len(rows)}행). 성공기준 DockQ 0.23/0.49/0.80. oracle≫score-pick = ZDOCK score가 pose 선택 실패(재랭커 여지).")

if __name__ == "__main__":
    main()
