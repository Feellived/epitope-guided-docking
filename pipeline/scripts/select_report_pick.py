#!/usr/bin/env python3
"""[대칭 비교] 정답을 보지 않고 고른 자세의 성적을 낸다.

왜 필요한가 — 데모 표의 값은 자세 5개 중 **최고**이고, 그건 정답을 봐야 고를 수 있다.
AlphaFold Server 는 자기 순위로 1순위를 내놓으므로, 그것과 맞대려면 우리 쪽도 정답 없이
고른 자세여야 한다. 여기서는 AbEpiScore 가 가장 높은 자세 하나를 골라 그 자세의
결합 자리 겹침과 DockQ 를 보고한다. 천장(최고값)도 같이 찍어 얼마나 잃는지 보인다.

입력은 eval_abepiscore.py score 가 낸 CSV(path,target,run,pose,score) 다.
팔 이름은 경로에서 읽는다: runs_.../<타깃>/<팔>/out_boltz/...

⚠️ DockQ CLI 가 있는 env 에서 돌린다(= conda activate boltz).

사용:
  python -u scripts/select_report_pick.py --scores results/abepi_*.csv
"""
import argparse, csv, glob, json, os, re, sys, tempfile


def arm_of(path):
    m = re.search(r"/([^/]+)/([^/]+)/out_[a-z]+/", path.replace("\\", "/"))
    return (m.group(1), m.group(2)) if m else (None, None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", nargs="+", required=True, help="eval_abepiscore.py score 결과 CSV")
    ap.add_argument("--msad", default=os.path.expanduser("~/projects/bk21-msa-depth-bias/pipeline"))
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--out", default="results/picked_by_abepi.csv")
    a = ap.parse_args()

    sys.path.insert(0, a.msad)
    import pose_features as PF                                   # noqa: E402

    canon = {d.lower(): d for d in os.listdir(os.path.join(a.msad, "targets"))}
    rows = []
    for pat in a.scores:
        for f in glob.glob(pat):
            rows += list(csv.DictReader(open(f)))
    if not rows:
        raise SystemExit("!! 점수 CSV 가 비었다")

    grp = {}
    for r in rows:
        p = r.get("path") or ""
        tdir, arm = arm_of(p)
        if not tdir:
            continue
        try:
            s = float(r["score"])
        except (TypeError, ValueError):
            continue
        grp.setdefault((tdir, arm), []).append((s, p))
    print(f"점수 {len(rows)}행 · (타깃, 팔) 묶음 {len(grp)}개\n")

    print(f"{'타깃':<10}{'팔':<14}{'자세':>4}{'고른 겹침':>10}{'천장 겹침':>10}"
          f"{'고른 DockQ':>11}{'천장 DockQ':>11}")
    print("-" * 72)
    out = []
    for (tdir, arm) in sorted(grp):
        tgt = canon.get(tdir.lower())
        if not tgt:
            print(f"{tdir:<10}{arm:<14}  targets/ 에 없는 이름 — 건너뜀"); continue
        cj = json.load(open(os.path.join(a.msad, "targets", tgt, "chains.json")))
        nat = os.path.join(a.msad, "targets", tgt, "native.cif")
        tr = PF.native_true(cj, nat, a.cutoff)
        if tr is None:
            print(f"{tgt:<10}{arm:<14}  native 결합자리 계산 실패 — 건너뜀"); continue
        true = tr[0]
        try:
            pop = PF.ES.popular_refset(cj, "")
        except Exception:
            pop = None
        poses = sorted(grp[(tdir, arm)], reverse=True)            # 점수 높은 순
        best_path = poses[0][1]
        with tempfile.TemporaryDirectory() as td:
            natm = PF.native_merged(cj, nat, td)
            def sc(p):
                met = PF.pose_all_metrics(cj, p, a.cutoff, true, pop) or {}
                pm = PF.pose_merged(cj, p, td)
                return met.get("recall"), (PF.dockq(pm, natm) if pm else None)
            rc_pick, dq_pick = sc(best_path)
            allm = [sc(p) for _, p in poses]
        rc_top = max((x for x, _ in allm if x is not None), default=None)
        dq_top = max((y for _, y in allm if y is not None), default=None)
        f = lambda v: f"{v:.3f}" if isinstance(v, float) else "  -  "
        print(f"{tgt:<10}{arm:<14}{len(poses):>4}{f(rc_pick):>10}{f(rc_top):>10}"
              f"{f(dq_pick):>11}{f(dq_top):>11}")
        out.append(dict(target=tgt, arm=arm, n_pose=len(poses),
                        recall_pick=rc_pick, recall_top=rc_top,
                        dockq_pick=dq_pick, dockq_top=dq_top,
                        picked=os.path.basename(best_path)))

    if out:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
        print(f"\n→ {a.out}  ({len(out)}행)")
    print("\n'고른' 값이 정답을 보지 않고 얻는 성적이고, '천장'은 정답을 봐야 얻는 값이다.")


if __name__ == "__main__":
    main()
