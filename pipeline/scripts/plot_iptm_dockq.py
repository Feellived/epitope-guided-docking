#!/usr/bin/env python3
"""summary.csv -> (A) ipTM vs DockQ 산점도(모델별 facet + Spearman rho + 0.23선)
                  (B) 모델 x 타깃 best-DockQ 히트맵 (그룹 브래킷: 행=MSA/no-MSA, 열=nanobody/Fab).
co-folder(서열->복합체) 모델만. 필요: matplotlib, numpy (scipy 있으면 rho에 사용).
사용(pipeline 루트):
  python scripts/plot_iptm_dockq.py --csv results/summary.csv --phase Benchmark
결과: results/iptm_dockq_scatter.png, results/dockq_heatmap.png
"""
import argparse, csv, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# co-folder 모델만 (Chai-1-noMSA 포함해야 no-MSA 그룹이 채워짐)
MODEL_ORDER = ["Boltz-2", "Chai-1", "Protenix", "AF2-Multimer", "tFold-Ag", "Chai-1-noMSA"]
MODEL_MSA = {"Boltz-2": "MSA", "Chai-1": "MSA", "Protenix": "MSA", "AF2-Multimer": "MSA",
             "tFold-Ag": "no-MSA", "Chai-1-noMSA": "no-MSA"}
TARGET_TYPE = {"10ZO": "nanobody", "9OBN": "nanobody", "9K6J": "Fab", "9XF9": "Fab", "9Y0A": "Fab"}
TYPE_ORDER = ["nanobody", "Fab"]
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
        rx = x.argsort().argsort().astype(float); ry = y.argsort().argsort().astype(float)
        rx -= rx.mean(); ry -= ry.mean()
        den = math.sqrt((rx**2).sum()) * math.sqrt((ry**2).sum())
        return float((rx*ry).sum()/den) if den > 0 else float("nan")


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def spans(labels):
    """연속 그룹 구간 [(start_idx, end_idx, label), ...]."""
    out, i = [], 0
    while i < len(labels):
        j = i
        while j + 1 < len(labels) and labels[j + 1] == labels[i]:
            j += 1
        out.append((i, j, labels[i])); i = j + 1
    return out


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
    ap.add_argument("--models", default=",".join(MODEL_ORDER))
    ap.add_argument("--out-dir", default="results")
    a = ap.parse_args()

    rows = load(a.csv, a.phase)
    if not rows:
        raise SystemExit(f"[에러] {a.csv} (phase={a.phase}) 유효행 없음")
    want = [m.strip() for m in a.models.split(",")]
    models = [m for m in want if any(r["model"] == m for r in rows)]
    missing = [m for m in want if m not in models]
    if missing:
        print(f"[경고] CSV에 없는 모델(스킵): {missing}")
    # 모델 정렬: MSA 그룹 먼저(MODEL_ORDER 순)
    models.sort(key=lambda m: (MODEL_MSA.get(m, "zz") != "MSA",
                               MODEL_ORDER.index(m) if m in MODEL_ORDER else 99))
    # 타깃 정렬: nanobody 먼저
    tset = sorted({r["target"] for r in rows if r["model"] in models})
    targets = sorted(tset, key=lambda t: (TYPE_ORDER.index(TARGET_TYPE.get(t, "Fab"))
                                          if TARGET_TYPE.get(t, "Fab") in TYPE_ORDER else 9, t))
    os.makedirs(a.out_dir, exist_ok=True)

    # ---- A) ipTM vs DockQ 산점도 (모델별 facet) ----
    n = len(models); ncol = min(6, n); nrow = math.ceil(n / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4*ncol, 3.3*nrow),
                             sharex=True, sharey=True, squeeze=False)
    axl = [axes[i//ncol][i % ncol] for i in range(nrow*ncol)]
    for k, m in enumerate(models):
        ax = axl[k]
        sub = [r for r in rows if r["model"] == m and r["gipTM"] is not None]
        xs = [r["gipTM"] for r in sub]; ys = [r["dockq"] for r in sub]
        cols = ["#2ca02c" if y >= THR else "#999999" for y in ys]
        ax.scatter(xs, ys, s=28, c=cols, alpha=0.8, edgecolors="none")
        ax.axhline(THR, ls="--", color="#d62728", lw=1)
        ax.set_title(f"{m} [{MODEL_MSA.get(m,'?')}]   rho={spearman(xs, ys):.2f}", fontsize=9)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xlabel("global ipTM")
        if k % ncol == 0:
            ax.set_ylabel("DockQ (Ab-Ag)")
    for k in range(len(models), nrow*ncol):
        axl[k].axis("off")
    fig.suptitle("ipTM vs DockQ per co-folder — each point = one pose "
                 "(green = success DockQ>=0.23; dashed line = 0.23)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    p1 = os.path.join(a.out_dir, "iptm_dockq_scatter.png")
    fig.savefig(p1, dpi=150, bbox_inches="tight"); plt.close(fig)

    # ---- B) best-DockQ 히트맵 + 그룹 브래킷 ----
    mat = np.full((len(models), len(targets)), np.nan)
    for i, m in enumerate(models):
        for j, t in enumerate(targets):
            vals = [r["dockq"] for r in rows if r["model"] == m and r["target"] == t]
            if vals:
                mat[i, j] = max(vals)
    fig, ax = plt.subplots(figsize=(1.35*len(targets)+3.8, 0.8*len(models)+3.2))
    fig.subplots_adjust(left=0.26, top=0.76, right=0.92, bottom=0.16)
    im = ax.imshow(mat, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(targets))); ax.set_xticklabels(targets)
    ax.set_yticks(range(len(models))); ax.set_yticklabels(models)
    for i in range(len(models)):
        for j in range(len(targets)):
            v = mat[i, j]
            if not np.isnan(v):
                ok = v >= THR
                ax.text(j, i, f"{v:.2f}" + ("*" if ok else ""), ha="center", va="center",
                        fontsize=9, fontweight="bold" if ok else "normal", color="black")
    cb = plt.colorbar(im, ax=ax, label="best DockQ", shrink=0.6, aspect=28, pad=0.02)
    cb.ax.tick_params(labelsize=9)

    # 열(타깃) 그룹 브래킷 — 위 (x=data, y=axes-frac)
    txf = ax.get_xaxis_transform()
    for c0, c1, lab in spans([TARGET_TYPE.get(t, "?") for t in targets]):
        yb, ycap, yt = 1.075, 1.025, 1.13
        ax.plot([c0-0.4, c1+0.4], [yb, yb], transform=txf, color="k", lw=2, clip_on=False)
        ax.plot([c0-0.4, c0-0.4], [yb, ycap], transform=txf, color="k", lw=2, clip_on=False)
        ax.plot([c1+0.4, c1+0.4], [yb, ycap], transform=txf, color="k", lw=2, clip_on=False)
        ax.text((c0+c1)/2, yt, lab, transform=txf, ha="center", va="bottom",
                fontsize=12, fontweight="bold", clip_on=False)
    # 행(모델) 그룹 브래킷 — 왼쪽 (x=axes-frac, y=data)
    tyf = ax.get_yaxis_transform()
    for r0, r1, lab in spans([MODEL_MSA.get(m, "?") for m in models]):
        xb, xcap, xt = -0.24, -0.21, -0.33
        ax.plot([xb, xb], [r0-0.4, r1+0.4], transform=tyf, color="k", lw=2, clip_on=False)
        ax.plot([xb, xcap], [r0-0.4, r0-0.4], transform=tyf, color="k", lw=2, clip_on=False)
        ax.plot([xb, xcap], [r1+0.4, r1+0.4], transform=tyf, color="k", lw=2, clip_on=False)
        ax.text(xt, (r0+r1)/2, lab, transform=tyf, ha="center", va="center",
                rotation=90, fontsize=12, fontweight="bold", clip_on=False)

    fig.suptitle("Co-folder (sequence->complex) models only  —  Best DockQ per model x complex "
                 "(* = success >= 0.23)", fontsize=12, y=0.97)
    fig.text(0.5, 0.03, "rows: MSA vs no-MSA   |   cols: nanobody vs Fab   |   "
             "tFold-Ag = antibody MSA-free / antigen MSA ; Chai-1-noMSA = full no-MSA",
             ha="center", fontsize=8, color="#555555")
    p2 = os.path.join(a.out_dir, "dockq_heatmap.png")
    fig.savefig(p2, dpi=150, bbox_inches="tight"); plt.close(fig)

    print(f"[저장] A 산점도 -> {p1}")
    print(f"[저장] B 히트맵 -> {p2}")
    print(f"모델 {models} | 타깃 {targets}")


if __name__ == "__main__":
    main()
