#!/usr/bin/env python3
"""summary.csv -> (A) ipTM vs DockQ 산점도(모델별 facet + Spearman rho + 0.23선)
                  (B) 모델 x 타깃 best-DockQ 히트맵.
필요: matplotlib, numpy (없으면 pip install matplotlib). scipy 있으면 사용, 없으면 rho 수동계산.
사용(consensus_docking 루트):
  python scripts/plot_iptm_dockq.py --csv results/summary.csv --phase Benchmark
결과: results/iptm_dockq_scatter.png, results/dockq_heatmap.png
"""
import argparse, csv, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DEF_MODELS = "Boltz-2,Chai-1,Protenix,AF2-Multimer,tFold-Ag"
THR = 0.23


def spearman(x, y):
    if len(x) < 3:
        return float("nan")
    try:
        from scipy.stats import spearmanr
        r = spearmanr(x, y).correlation
        return float(r) if r == r else float("nan")
    except Exception:
        x = np.asarray(x, float); y = np.asarray(y, float)
        rx = x.argsort().argsort().astype(float)
        ry = y.argsort().argsort().astype(float)
        rx -= rx.mean(); ry -= ry.mean()
        den = math.sqrt((rx**2).sum()) * math.sqrt((ry**2).sum())
        return float((rx*ry).sum()/den) if den > 0 else float("nan")


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def load(csv_path, phase):
    rows = []
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            if phase and r.get("phase") != phase:
                continue
            dq = fnum(r.get("dockq_abag"))
            if dq is None:
                continue
            rows.append({"model": r["model"], "target": r["target"],
                         "dockq": dq, "gipTM": fnum(r.get("global_iptm"))})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="results/summary.csv")
    ap.add_argument("--phase", default="Benchmark")
    ap.add_argument("--models", default=DEF_MODELS, help="쉼표구분; CSV에 있는 것만 사용")
    ap.add_argument("--out-dir", default="results")
    a = ap.parse_args()

    rows = load(a.csv, a.phase)
    if not rows:
        raise SystemExit(f"[에러] {a.csv} (phase={a.phase}) 유효행 없음")
    want = [m.strip() for m in a.models.split(",")]
    present = [m for m in want if any(r["model"] == m for r in rows)]
    missing = [m for m in want if m not in present]
    if missing:
        print(f"[경고] CSV에 없는 모델(스킵): {missing}")
    models = present
    targets = sorted({r["target"] for r in rows if r["model"] in models})
    os.makedirs(a.out_dir, exist_ok=True)

    # ---- A) ipTM vs DockQ 산점도 (모델별 facet) ----
    n = len(models)
    ncol = min(5, n); nrow = math.ceil(n / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6*ncol, 3.4*nrow),
                             sharex=True, sharey=True, squeeze=False)
    axl = [axes[i//ncol][i % ncol] for i in range(nrow*ncol)]
    for k, m in enumerate(models):
        ax = axl[k]
        sub = [r for r in rows if r["model"] == m and r["gipTM"] is not None]
        xs = [r["gipTM"] for r in sub]; ys = [r["dockq"] for r in sub]
        cols = ["#2ca02c" if y >= THR else "#999999" for y in ys]  # 성공=초록
        ax.scatter(xs, ys, s=28, c=cols, alpha=0.8, edgecolors="none")
        ax.axhline(THR, ls="--", color="#d62728", lw=1)
        rho = spearman(xs, ys)
        ax.set_title(f"{m}   rho={rho:.2f}", fontsize=10)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel("global ipTM (confidence)")
        if k % ncol == 0:
            ax.set_ylabel("DockQ (Ab-Ag accuracy)")
    for k in range(len(models), nrow*ncol):
        axl[k].axis("off")
    fig.suptitle("ipTM vs DockQ — each point = one pose  (green = success DockQ>=0.23; dashed = 0.23)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    p1 = os.path.join(a.out_dir, "iptm_dockq_scatter.png")
    fig.savefig(p1, dpi=150, bbox_inches="tight"); plt.close(fig)

    # ---- B) 모델 x 타깃 best-DockQ 히트맵 ----
    mat = np.full((len(models), len(targets)), np.nan)
    for i, m in enumerate(models):
        for j, t in enumerate(targets):
            vals = [r["dockq"] for r in rows if r["model"] == m and r["target"] == t]
            if vals:
                mat[i, j] = max(vals)
    fig, ax = plt.subplots(figsize=(1.3*len(targets)+2.5, 0.72*len(models)+2))
    im = ax.imshow(mat, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(targets))); ax.set_xticklabels(targets)
    ax.set_yticks(range(len(models))); ax.set_yticklabels(models)
    for i in range(len(models)):
        for j in range(len(targets)):
            v = mat[i, j]
            if not np.isnan(v):
                ok = v >= THR
                ax.text(j, i, f"{v:.2f}" + ("*" if ok else ""),
                        ha="center", va="center", fontsize=9,
                        fontweight="bold" if ok else "normal",
                        color="black")
    plt.colorbar(im, ax=ax, label="best DockQ", fraction=0.046, pad=0.04)
    ax.set_title("Best DockQ per model x target  (* = success >= 0.23)", fontsize=11)
    fig.tight_layout()
    p2 = os.path.join(a.out_dir, "dockq_heatmap.png")
    fig.savefig(p2, dpi=150, bbox_inches="tight"); plt.close(fig)

    print(f"[저장] A 산점도 -> {p1}")
    print(f"[저장] B 히트맵 -> {p2}")
    print(f"모델 {models} | 타깃 {targets}")


if __name__ == "__main__":
    main()
