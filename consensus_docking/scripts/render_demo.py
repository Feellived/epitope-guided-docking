#!/usr/bin/env python3
"""[데모 그림] 유도 도킹 전/후를 PyMOL 그림으로 — 팔마다 최고 자세를 골라 나란히 붙인다.

무엇을 보여주나:
  native 항원(회색) 위에 **진짜 결합자리를 진홍**으로 칠하고, 모델이 놓은 항체를 주황으로 얹는다.
  주황이 진홍 위에 앉으면 성공, 딴 데 가 있으면 실패 — 숫자를 안 봐도 보인다.

왜 세 팔인가 (기본값):
  blind    힌트 없음        = 돌리기 전
  fullmsa  원래 MSA가 간 자리로 유도  = **대조**(제약을 주면 아무 자리나 좋아지는지)
  ours     우리가 고른 후보로 유도    = 돌린 후
  판정은 ours > fullmsa 다. blind 와만 비교하면 "제약을 주면 좋아진다"와 구별되지 않으므로
  가운데 칸을 빼지 말 것.

자세 선택 = 그 팔 안에서 DockQ 최고(--by recall 로 자리겹침 기준으로 바꿀 수 있다).
**정답을 보고 고른 값이므로 그림은 '천장'이다** — 캡션에 그대로 쓸 것.
고른 파일 경로를 CSV로 남긴다(같은 그림을 다시 만들 수 있게).

사용 (DockQ env, pymol-open-source 필요):
  # 좋아진 두 종
  python scripts/render_demo.py --targets "8k3k_D 8sis_HL" --outpng fig
  # 팔·기준 바꾸기
  python scripts/render_demo.py --targets 8k3k_D --arms "blind sizematch ours" --by recall
  # 자세만 고르고 렌더링은 건너뛰기(경로 확인용)
  python scripts/render_demo.py --targets 8k3k_D --dry-run
"""
import argparse, csv, glob, json, os, subprocess, sys, tempfile

# 캡션은 DejaVu 폰트로 얹히는데 한글 글리프가 없어 네모로 나온다 → 영문 표기
ARM_EN = {"blind": "no hint", "noconstraint": "no constraint",
          "fullmsa": "original-MSA site", "sizematch": "size-matched decoy",
          "ours": "OURS"}


def load_pf(msad):
    sys.path.insert(0, msad)
    import pose_features as PF                      # noqa: E402
    return PF


def poses_of(d):
    return sorted(glob.glob(os.path.join(d, "**", "*.cif"), recursive=True))


def arm_poses(a, target, arm, sj, models):
    """(경로, 모델라벨) 목록. blind 는 msa-depth 의 원래 MSA 실행에서 가져온다."""
    if arm == "blind":
        base = os.path.join(a.msad_data, "compreps", "seedrep_cand",
                            sj["model"], target, sj["depth"])
        return [(p, sj["model"]) for d in sorted(glob.glob(os.path.join(base, "seedfull_r*")))
                for p in poses_of(d)]
    out = []
    for m in models:
        out += [(p, m) for p in poses_of(os.path.join(a.outdir, target.lower(), arm, f"out_{m}"))]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", required=True)
    ap.add_argument("--arms", default="blind fullmsa ours")
    ap.add_argument("--outdir", default="runs_sites_guided")
    ap.add_argument("--outpng", default="fig")
    ap.add_argument("--msad", default=os.path.expanduser("~/projects/bk21-msa-depth-bias/pipeline"))
    ap.add_argument("--msad-data", default="/mnt/data/admuser/msadepth")
    ap.add_argument("--models", default="boltz protenix")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--by", choices=["dockq", "recall"], default="dockq")
    ap.add_argument("--surface", action="store_true", help="항원을 표면으로 (에피토프가 더 잘 보인다)")
    ap.add_argument("--width", type=int, default=1400)
    ap.add_argument("--height", type=int, default=1100)
    ap.add_argument("--pymol", default="pymol")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--csv", default="results/demo_figures.csv")
    a = ap.parse_args()

    here = os.path.dirname(os.path.abspath(__file__))
    pose_script = os.path.join(here, "pymol_pose.py")
    if not os.path.exists(pose_script):
        sys.exit(f"!! pymol_pose.py 가 없다: {pose_script}")
    PF = load_pf(a.msad)
    models = a.models.replace(",", " ").split()
    arms = a.arms.replace(",", " ").split()
    os.makedirs(a.outpng, exist_ok=True)
    picked = []

    for t in a.targets.replace(",", " ").split():
        sp = os.path.join(a.msad, "results", f"sites_{t}.json")
        cjp = os.path.join(a.msad, "targets", t, "chains.json")
        native = os.path.join(a.msad, "targets", t, "native.cif")
        for p in (sp, cjp, native):
            if not os.path.exists(p):
                sys.exit(f"!! {t}: 없는 파일 {p}")
        sj, cj = json.load(open(sp)), json.load(open(cjp))
        tr = PF.native_true(cj, native, a.cutoff)
        if tr is None:
            print(f"{t}: native 결합자리 계산 실패 — 건너뜀"); continue
        true = tr[0]
        try:
            pop = PF.ES.popular_refset(cj, "")
        except Exception:
            pop = None

        with tempfile.TemporaryDirectory() as td:
            natm = PF.native_merged(cj, native, td)
            for arm in arms:
                ps = arm_poses(a, t, arm, sj, models)
                if not ps:
                    print(f"{t} / {arm}: 자세 없음 — 건너뜀"); continue
                best = None
                for path, mlab in ps:
                    try:
                        met = PF.pose_all_metrics(cj, path, a.cutoff, true, pop)
                        pm = PF.pose_merged(cj, path, td)
                        dq = PF.dockq(pm, natm) if pm else None
                    except Exception:
                        continue
                    rc = (met or {}).get("recall")
                    key = dq if a.by == "dockq" else rc
                    if key is None:
                        continue
                    if best is None or key > best[0]:
                        best = (key, dq, rc, path, mlab)
                if best is None:
                    print(f"{t} / {arm}: 채점된 자세 없음 — 건너뜀"); continue
                _, dq, rc, path, mlab = best
                png = os.path.join(a.outpng, f"demo_{t}_{arm}.png")
                cap = (f"{t}  |  {ARM_EN.get(arm, arm)}  |  {mlab}  |  "
                       f"DockQ={dq:.3f}  epitope-recall={rc:.3f}  (ceiling)")
                print(f"{t:10}{arm:14}{mlab:9}DockQ={dq:.3f}  recall={rc:.3f}  {path}")
                picked.append(dict(target=t, arm=arm, model=mlab,
                                   dockq=f"{dq:.4f}", recall=f"{rc:.4f}",
                                   pose=path, png=png))
                if a.dry_run:
                    continue
                cl = [a.pymol, "-cq", pose_script, "--", "--pred", path,
                      "--native", native, "--chains", cjp, "--out", png,
                      "--label", cap, "--width", str(a.width), "--height", str(a.height)]
                if a.surface:
                    cl.append("--surface")
                r = subprocess.run(cl, capture_output=True, text=True)
                if not os.path.exists(png):
                    print(r.stdout[-600:]); print(r.stderr[-600:])
                    sys.exit(f"!! 렌더링 실패 {t}/{arm}. pymol 설치 확인: {a.pymol} -cq -d 'print(1)'")

    if not picked:
        sys.exit("!! 고른 자세가 없다. --outdir / --msad-data 경로를 볼 것.")
    os.makedirs(os.path.dirname(a.csv) or ".", exist_ok=True)
    with open(a.csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(picked[0])); w.writeheader(); w.writerows(picked)
    print(f"\n→ 고른 자세 목록: {a.csv}")

    # 타깃별로 팔들을 가로로 붙여 한 장으로 (Pillow 있을 때만)
    if not a.dry_run:
        try:
            from PIL import Image
        except Exception:
            print("[정보] Pillow 없음 → 칸별 PNG 만 생성.")
            return
        for t in sorted({x["target"] for x in picked}):
            fs = [x["png"] for x in picked if x["target"] == t and os.path.exists(x["png"])]
            if len(fs) < 2:
                continue
            ims = [Image.open(f).convert("RGB") for f in fs]
            h = max(i.height for i in ims)
            sh = [i.resize((int(i.width * h / i.height), h)) for i in ims]
            canvas = Image.new("RGB", (sum(i.width for i in sh), h), (20, 20, 20))
            x = 0
            for i in sh:
                canvas.paste(i, (x, 0)); x += i.width
            out = os.path.join(a.outpng, f"demo_{t}_panels.png")
            canvas.save(out)
            print(f"→ {out}  ({len(fs)}칸)")


if __name__ == "__main__":
    main()
