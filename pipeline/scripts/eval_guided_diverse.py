#!/usr/bin/env python3
"""guided co-folder(diverse) DockQ vs blind — 물리 에피토프 힌트가 co-folder를 native로 정교화하나.
   ZDOCK-guided(run_zdock_guided.sh) · HADDOCK-guided(추후) 공용. 모델별(boltz/protenix/tfold).
   guided_best = 패치×pose 중 최고 DockQ(= 클러스터/패치 oracle 상한). 실배포 '어느 패치인가'는 재랭커 몫.
   blind = runs_diverse/<t>/out_<model> (기존 5모델 예측). 항체 H+L 병합, 항원=서열매칭.

   출력: results/guided_<model>_<tag>.csv  (target, blind_best, guided_best, diff, n_pose_blind, n_pose_guided)
   resume: 이미 CSV에 있으면 skip(--force). 타깃마다 flush.

   사용(DockQ env, 예 boltz):
     python scripts/eval_guided_diverse.py --model boltz --guided-dir runs_zdock_guided --targets "9J4C 9KKJ 9Y0A 9W43"
     python scripts/eval_guided_diverse.py --model boltz --guided-dir runs_zdock_guided_union --tag union ...   # 넓게 버전
     python scripts/eval_guided_diverse.py --model protenix ...   /   --model tfold ...
"""
import argparse, csv, glob, os, tempfile
from haddock_dockq_diverse import load_any, write_merged, dockq, best_ag_chain, target_info, DIVERSE

# 모델별 (blind 하위폴더, pose glob, guided 하위폴더 패턴)
# 모델별 (blind 하위폴더, pose glob, guided 하위폴더 패턴). blind 폴더는 out_<model>로 통일.
SPEC = {
    "boltz":    ("out_boltz",    "*_model_*.cif",  "results_patch*"),
    "protenix": ("out_protenix", "*_sample_*.cif", "results_patch*"),
    "tfold":    ("out_tfold",    "*.pdb",          "out_tfold_patch*"),
}
COLS = ["target", "model", "blind_best", "guided_best", "diff", "n_blind", "n_guided"]


def best_over(poses, ni, natm, td):
    bq = None; n = 0
    for p in poses:
        try:
            pm = load_any(p)
            agc = best_ag_chain(pm, ni["ag_seq"])
            abc = [ch.id for ch in pm if ch.id != agc]
            mp = os.path.join(td, "m.pdb")
            if not write_merged(pm, [agc], abc, mp):
                continue
            q = dockq(mp, natm); n += 1
            if q is not None and (bq is None or q > bq):
                bq = q
        except Exception:
            continue
    return bq, n


def load_done(out):
    if not os.path.exists(out):
        return {}
    return {r["target"].upper(): r for r in csv.DictReader(open(out)) if r.get("guided_best", "") not in ("", "-")}


def flush(out, done):
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS); w.writeheader()
        for t in done:
            w.writerow(done[t])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True, choices=list(SPEC))
    ap.add_argument("--guided-dir", default="runs_zdock_guided")
    ap.add_argument("--run-dir", default="runs_diverse")
    ap.add_argument("--targets", default=" ".join(DIVERSE))
    ap.add_argument("--tag", default="", help="출력 파일 태그(예: union). 기본=guided-dir 이름")
    ap.add_argument("--out", default="")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    blind_sub, pose_glob, guided_sub = SPEC[a.model]
    tag = a.tag or os.path.basename(a.guided_dir.rstrip("/"))
    out = a.out or f"results/guided_{a.model}_{tag}.csv"

    done = {} if a.force else load_done(out)
    if done:
        print(f"[resume] skip {len(done)}: {', '.join(done)}")
    print(f"[guided {a.model} | blind={a.run_dir}/*/{blind_sub} vs guided={a.guided_dir}/*/{guided_sub}]")
    print(f"{'target':8}{'blind':>9}{'guided':>9}{'diff':>9}  (n_blind/n_guided)")
    print("-" * 52)
    for t in a.targets.split():
        t = t.upper()
        if not a.force and t in done:
            continue
        ni = target_info(t, a.run_dir)
        if ni is None or ni["native"] is None:
            print(f"{t:8} native 없음 skip"); continue
        blind = glob.glob(os.path.join(a.run_dir, t.lower(), blind_sub, "**", pose_glob), recursive=True)
        guided = glob.glob(os.path.join(a.guided_dir, t.lower(), guided_sub, "**", pose_glob), recursive=True)
        if not guided:
            print(f"{t:8} guided pose 없음(미실행) skip"); continue
        with tempfile.TemporaryDirectory() as td:
            natm = os.path.join(td, "nat.pdb")
            if not write_merged(load_any(ni["native"]), [ni["ag"]], ni["ab"], natm):
                print(f"{t:8} native merge 실패 skip"); continue
            bq, nb = best_over(blind, ni, natm, td)
            gq, ng = best_over(guided, ni, natm, td)
        bs = f"{bq:.3f}" if bq is not None else "-"
        gs = f"{gq:.3f}" if gq is not None else "-"
        ds = f"{gq - bq:+.3f}" if (gq is not None and bq is not None) else "-"
        print(f"{t:8}{bs:>9}{gs:>9}{ds:>9}  ({nb}/{ng})")
        done[t] = dict(target=t, model=a.model, blind_best=bs, guided_best=gs, diff=ds, n_blind=nb, n_guided=ng)
        flush(out, done)
    print(f"\n→ {out} ({len(done)}행). guided_best > blind = 물리 에피토프 힌트가 co-folder를 native로 끌어올림.")
    print("  (⚠️ guided_best=패치 oracle 상한 — 실배포는 '어느 패치' 선택문제=재랭커 몫)")


if __name__ == "__main__":
    main()
