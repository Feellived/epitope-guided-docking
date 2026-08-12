#!/usr/bin/env python3
"""[데모 채점] 세 팔의 DockQ·결합자리 겹침을 한 표로 비교한다.

  (a) blind    제약 없음 — msa-depth 본 검정의 '원래 MSA(seedfull)' 실행
  (b) fullmsa  원래 MSA 가 간 자리로 유도
  (c) ours     우리가 고른 후보로 유도

판정은 **(c) > (b)**. (c) > (a) 만 보이면 '제약을 주면 좋아진다'와 구별되지 않는다.

채점 자체는 bk21-msa-depth-bias 의 pose_features 를 그대로 쓴다(같은 DockQ·같은 결합자리
정의라야 팔끼리 비교가 성립한다). 그래서 --msad 로 그 pipeline 경로를 받아 import 한다.

사용(DockQ env):
  python scripts/eval_demo.py --targets "8k3k_D 8sit_HL" --outdir runs_sites_guided
  python scripts/eval_demo.py --targets 8k3k_D --csv results/demo_dockq.csv
"""
import argparse, csv, glob, json, os, statistics as st, sys, tempfile

ARMS = ["blind", "noconstraint", "fullmsa", "sizematch", "ours"]


def load_pf(msad):
    sys.path.insert(0, msad)
    import pose_features as PF                      # noqa: E402
    return PF


def poses_of(d):
    return sorted(glob.glob(os.path.join(d, "**", "*.cif"), recursive=True))


def blind_dir(msad_data, target, model, depth):
    """msa-depth 의 원래 MSA 실행 = compreps/seedrep_cand/<model>/<타깃>/<깊이>/seedfull_r*/"""
    base = os.path.join(msad_data, "compreps", "seedrep_cand", model, target, depth)
    return sorted(glob.glob(os.path.join(base, "seedfull_r*")))


def score(PF, cj, native, poses, cutoff, true, pop, td):
    natm = PF.native_merged(cj, native, td)
    out = []
    for p in poses:
        try:
            met = PF.pose_all_metrics(cj, p, cutoff, true, pop)
            pm = PF.pose_merged(cj, p, td)
            q = PF.dockq(pm, natm) if pm else None
        except Exception:
            continue
        out.append((q, (met or {}).get("recall")))
    return out


def agg(vals, i):
    v = [x[i] for x in vals if x[i] is not None]
    return (max(v), st.median(v), len(v)) if v else (None, None, 0)


def fmt(x):
    return f"{x:.3f}" if isinstance(x, float) else "  -  "


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", required=True)
    ap.add_argument("--outdir", default="runs_sites_guided")
    ap.add_argument("--msad", default=os.path.expanduser("~/projects/bk21-msa-depth-bias/pipeline"))
    ap.add_argument("--msad-data", default="/mnt/data/admuser/msadepth")
    ap.add_argument("--models", default="boltz protenix")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--csv", default="results/demo_dockq.csv")
    ap.add_argument("--arms", default=" ".join(ARMS),
                    help="채점할 팔. 기본은 다섯 팔 전부(예전과 같다). "
                         "후보별 재도킹처럼 'ours' 만 있는 폴더는 --arms ours 로 "
                         "주면 blind(자세 40개)를 폴더마다 다시 채점하지 않는다.")
    a = ap.parse_args()
    PF = load_pf(a.msad)
    models = a.models.replace(",", " ").split()
    arms = [x for x in a.arms.replace(",", " ").split() if x in ARMS]
    bad = [x for x in a.arms.replace(",", " ").split() if x not in ARMS]
    if bad:
        raise SystemExit(f"!! 모르는 팔: {bad}  (있는 것 = {ARMS})")
    if not arms:
        raise SystemExit("!! 채점할 팔이 없다")
    rows = []

    print(f"{'타깃':11}{'모델':9}{'팔':9}{'자세':>5}{'DockQ최고':>10}{'DockQ중앙':>10}"
          f"{'자리겹침최고':>12}{'자리겹침중앙':>12}")
    print("-" * 80)
    for t in a.targets.replace(",", " ").split():
        sites_p = os.path.join(a.msad, "results", f"sites_{t}.json")
        if not os.path.exists(sites_p):
            print(f"{t}: {sites_p} 없음 — 건너뜀"); continue
        sj = json.load(open(sites_p))
        cjp = os.path.join(a.msad, "targets", t, "chains.json")
        native = os.path.join(a.msad, "targets", t, "native.cif")
        if not (os.path.exists(cjp) and os.path.exists(native)):
            print(f"{t}: chains.json/native.cif 없음 — 건너뜀"); continue
        cj = json.load(open(cjp))
        tr = PF.native_true(cj, native, a.cutoff)
        if tr is None:
            print(f"{t}: native 결합자리 계산 실패 — 건너뜀"); continue
        true = tr[0]
        try:
            pop = PF.ES.popular_refset(cj, "")
        except Exception:
            pop = None

        with tempfile.TemporaryDirectory() as td:
            for m in models:
                for arm in arms:
                    if arm == "blind":
                        ds = blind_dir(a.msad_data, t, sj["model"], sj["depth"])
                        ps = [p for d in ds for p in poses_of(d)]
                        note = f"(msa-depth {sj['model']} {sj['depth']} seedfull)"
                        if m != models[0]:
                            continue          # blind 는 모델당 중복이 아니라 한 번만
                    else:
                        ps = poses_of(os.path.join(a.outdir, t.lower(), arm, f"out_{m}"))
                        note = ""
                    if not ps:
                        print(f"{t:11}{m:9}{arm:9}    0   (자세 없음) {note}")
                        continue
                    vals = score(PF, cj, native, ps, a.cutoff, true, pop, td)
                    dq_max, dq_med, n = agg(vals, 0)
                    rc_max, rc_med, _ = agg(vals, 1)
                    label = m if arm != "blind" else f"{sj['model']}*"
                    print(f"{t:11}{label:9}{arm:9}{n:>5}{fmt(dq_max):>10}{fmt(dq_med):>10}"
                          f"{fmt(rc_max):>12}{fmt(rc_med):>12}  {note}")
                    rows.append(dict(target=t, model=label, arm=arm, n_pose=n,
                                     dockq_max=fmt(dq_max), dockq_med=fmt(dq_med),
                                     recall_max=fmt(rc_max), recall_med=fmt(rc_med)))

    if rows:
        os.makedirs(os.path.dirname(a.csv) or ".", exist_ok=True)
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        print(f"\n→ {a.csv}")
    print("\n판정: (c) ours > (b) fullmsa 여야 '조성 다양화가 옳은 자리를 올렸다'가 증명된다.")
    print("      (c) > (a) blind 만으로는 '제약을 주면 좋아진다'와 구별되지 않는다.")


if __name__ == "__main__":
    main()
