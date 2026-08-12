#!/usr/bin/env python3
"""[자리 그림] 팔마다 최고 자세를 골라 plot_epitope.py 로 넘긴다 — 한 세션·같은 각도.

plot_demo.py 가 '항체를 얹은 그림'이라면 이쪽은 '항원 표면에 자리를 칠한 그림'이다.
자세(DockQ)가 낮아도 **자리가 옮겨간 것**은 그대로 보이므로 이 연구의 주장에 맞는다.

사용 (DockQ env):
  python scripts/plot_sites.py --targets "8k3k_D 8sis_HL 8ume_HL" \\
      --outpng fig --pymol ~/miniconda3/envs/pymol/bin/pymol
"""
import argparse, csv, glob, json, os, subprocess, sys, tempfile

ARM_EN = {"blind": "original-MSA runs (other model, larger pose budget)",
          "noconstraint": "no constraint",
          "fullmsa": "original-MSA site", "sizematch": "area-matched decoy",
          "ours": "composition site"}


def load_pf(msad):
    sys.path.insert(0, msad)
    import pose_features as PF                      # noqa: E402
    return PF


def poses_of(d):
    return sorted(glob.glob(os.path.join(d, "**", "*.cif"), recursive=True))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", required=True)
    ap.add_argument("--arms", default="noconstraint fullmsa sizematch ours",
                    help="blind 는 다른 모델·다른 자세 예산이라 기본 비교에서 제외한다")
    ap.add_argument("--outdir", default="runs_sites_guided")
    ap.add_argument("--outpng", default="fig")
    ap.add_argument("--msad", default=os.path.expanduser("~/projects/bk21-msa-depth-bias/pipeline"))
    ap.add_argument("--msad-data", default="/mnt/data/admuser/msadepth")
    ap.add_argument("--models", default="boltz protenix")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--pymol", default="pymol")
    ap.add_argument("--no-fp", action="store_true",
                    help="잘못 짚은 부분을 칠하지 않는다(그림이 깨끗해진다)")
    ap.add_argument("--hit-only", action="store_true", help="맞힌 부분만 칠한다")
    ap.add_argument("--main-patch", action="store_true", help="흩어진 잔기를 버리고 큰 덩어리만")
    ap.add_argument("--repr", default="surface", choices=["surface", "cartoon"])
    ap.add_argument("--show-ab", action="store_true", help="예측 항체를 반투명으로 얹는다")
    ap.add_argument("--show-native-ab", action="store_true", help="정답 항체도 흐릿하게")
    ap.add_argument("--suffix", default="", help="출력 이름에 붙일 꼬리표 (두 판을 함께 보관할 때)")
    ap.add_argument("--csv", default="results/site_figures.csv")
    a = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    script = os.path.join(here, "plot_epitope.py")
    if not os.path.exists(script):
        sys.exit(f"!! plot_epitope.py 가 없다: {script}")
    PF = load_pf(a.msad)
    models = a.models.replace(",", " ").split()
    arms = a.arms.replace(",", " ").split()
    os.makedirs(a.outpng, exist_ok=True)
    rows = []

    for t in a.targets.replace(",", " ").split():
        sp = os.path.join(a.msad, "results", f"sites_{t}.json")
        cjp = os.path.join(a.msad, "targets", t, "chains.json")
        native = os.path.join(a.msad, "targets", t, "native.cif")
        if not all(os.path.exists(p) for p in (sp, cjp, native)):
            print(f"{t}: 입력 파일 없음 — 건너뜀"); continue
        sj, cj = json.load(open(sp)), json.load(open(cjp))
        tr = PF.native_true(cj, native, a.cutoff)
        if tr is None:
            print(f"{t}: 정답 결합자리 계산 실패 — 건너뜀"); continue
        true = tr[0]
        try:
            pop = PF.ES.popular_refset(cj, "")
        except Exception:
            pop = None

        preds = []
        with tempfile.TemporaryDirectory() as td:
            for arm in arms:
                if arm == "blind":
                    base = os.path.join(a.msad_data, "compreps", "seedrep_cand",
                                        sj["model"], t, sj["depth"])
                    ps = [p for d in sorted(glob.glob(os.path.join(base, "seedfull_r*")))
                          for p in poses_of(d)]
                else:
                    ps = [p for m in models
                          for p in poses_of(os.path.join(a.outdir, t.lower(), arm, f"out_{m}"))]
                best = None
                for p in ps:
                    try:
                        met = PF.pose_all_metrics(cj, p, a.cutoff, true, pop)
                    except Exception:
                        continue
                    rc = (met or {}).get("recall")
                    if rc is not None and (best is None or rc > best[0]):
                        best = (rc, p)
                if best is None:
                    print(f"{t} / {arm}: 채점된 자세 없음 — 건너뜀"); continue
                preds.append((ARM_EN.get(arm, arm), best[1]))
                rows.append(dict(target=t, arm=arm, recall=f"{best[0]:.4f}",
                                 n_pose=len(ps), pose=best[1]))
                print(f"{t:10}{arm:12}recall={best[0]:.3f}  n={len(ps):>3}  {best[1]}")

        if not preds:
            continue
        cl = [a.pymol, "-cq", script, "--", "--native", native, "--chains", cjp,
              "--out-prefix", os.path.join(a.outpng, f"site_{t}{a.suffix}"),
              "--cutoff", str(a.cutoff)]
        for flag, on in (("--no-fp", a.no_fp), ("--hit-only", a.hit_only),
                         ("--main-patch", a.main_patch), ("--show-ab", a.show_ab),
                         ("--show-native-ab", a.show_native_ab)):
            if on:
                cl.append(flag)
        cl += ["--repr", a.repr]
        for lab, p in preds:
            cl += ["--pred", f"{lab}={p}"]
        r = subprocess.run(cl, capture_output=True, text=True)
        print(r.stdout.strip() or r.stderr[-500:])

    if rows:
        os.makedirs(os.path.dirname(a.csv) or ".", exist_ok=True)
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
        print(f"\n→ {a.csv}")


if __name__ == "__main__":
    main()
