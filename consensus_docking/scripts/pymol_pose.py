#!/usr/bin/env python3
"""[PyMOL] 예측 pose를 native에 정렬해 도킹 성패를 그림 한 장으로.

핵심 아이디어: 항원(antigen)을 기준으로 예측 구조를 native 위에 구조정렬(cealign)로 겹친 뒤,
'모델이 놓은 항체'가 '진짜 에피토프'에 앉았는지 본다.
  - native 항원   = 회색 (공통 프레임)
  - native 에피토프 = 진홍 (항원 중 native 항체와 5Å 이내 접촉 잔기)
  - native 항체   = 연하늘 반투명 (정답 pose = '유령')
  - 예측 항체     = 주황 (모델의 pose)
겹치면 성공, 빨강 패치를 벗어나면 실패가 한눈에 보인다.

예측 항원 사슬은 자동 식별(각 예측 사슬을 native 항원에 cealign → 가장 잘 맞는 것 = 항원).
따라서 AF2-M(A/B/C 재라벨)·tFold(H/L/A)처럼 사슬 문자가 native와 달라도 동작.

헤드리스 실행 (pymol-open-source):
  # (a) summary.csv에서 best-DockQ seed 자동 선택 — 권장
  pymol -cq scripts/pymol_pose.py -- --run-dir runs/9k6j --target 9K6J \
        --model tFold-Ag --from-csv results/summary.csv --out fig/tfold_9k6j.png
  # (b) 예측 파일 직접 지정
  pymol -cq scripts/pymol_pose.py -- --run-dir runs/9k6j --pred runs/9k6j/out_tfold/xxx_seed3.pdb \
        --label "tFold-Ag 9K6J DockQ=0.49" --out fig/tfold_9k6j.png

필요: pymol-open-source (conda -c conda-forge). 캡션은 Pillow 있으면 자동으로 그림 상단에 얹음(없으면 생략).
"""
import argparse
import csv
import glob
import json
import os
import sys

from pymol import cmd

# 모델별 예측 pose glob (run-dir 기준). {seed}=summary.csv의 seed 열, {low}=소문자 PDB id
PRED_GLOB = {
    "Boltz-2":       "results_boltz/**/*_model_{seed}.cif",
    "Chai-1":        "out_chai/pred.model_idx_{seed}.cif",
    "Chai-1-noMSA":  "out_chai_nomsa/pred.model_idx_{seed}.cif",
    "Protenix":      "out_protenix/**/{low}_sample_{seed}.cif",
    "Protenix-base": "out_protenix/**/{low}_sample_{seed}.cif",
    "AF2-Multimer":  "results_colabfold/*_model_{seed}_*.pdb",
    "tFold-Ag":      "out_tfold/*seed{seed}*.pdb",
}


def argv_after_dashdash():
    a = sys.argv
    return a[a.index("--") + 1:] if "--" in a else a[1:]


def pick(patterns):
    for p in patterns:
        h = sorted(glob.glob(p, recursive=True))
        if h:
            return h[0]
    return None


def best_from_csv(csv_path, phase, model, target):
    """summary.csv에서 (phase,model,target)의 best-DockQ 행 → (seed, dockq)."""
    best = None
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            if r.get("model") != model or r.get("target") != target:
                continue
            if phase and r.get("phase") != phase:
                continue
            try:
                dq = float(r.get("dockq_abag"))
            except (TypeError, ValueError):
                continue
            seed = r.get("seed", "0")
            if best is None or dq > best[1]:
                best = (seed, dq)
    return best


def chains_of(sel):
    return [c for c in cmd.get_chains(sel) if c and c.strip()]


def detect_antigen(native_ag_sel):
    """예측(prd)의 각 사슬을 native 항원에 cealign → 가장 잘 맞는 사슬 = 예측 항원."""
    best = None
    for c in chains_of("prd"):
        try:
            r = cmd.cealign(native_ag_sel, f"prd and chain {c}")
        except Exception:
            continue
        score = (r.get("alignment_length", 0), -r.get("RMSD", 99.0))
        if best is None or score > best[0]:
            best = (score, c)
    return best[1] if best else None


def add_caption(png_path, text):
    """Pillow 있으면 그림 상단에 검은 띠 + 흰 글씨 캡션."""
    if not text:
        return
    try:
        from PIL import Image, ImageDraw, ImageFont
    except Exception:
        print("[정보] Pillow 없음 → 캡션 생략 (파일명으로 구분).")
        return
    try:
        img = Image.open(png_path).convert("RGB")
        W, _ = img.size
        bar = max(34, W // 26)
        canvas = Image.new("RGB", (W, img.height + bar), (20, 20, 20))
        canvas.paste(img, (0, bar))
        d = ImageDraw.Draw(canvas)
        try:
            font = ImageFont.truetype("DejaVuSans-Bold.ttf", int(bar * 0.62))
        except Exception:
            font = ImageFont.load_default()
        d.text((12, bar * 0.16), text, fill=(255, 255, 255), font=font)
        canvas.save(png_path)
    except Exception as e:
        print(f"[정보] 캡션 실패({e}) → 원본 유지.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", default=".")
    ap.add_argument("--pred", help="예측 구조 파일(.cif/.pdb) 또는 glob (미지정 시 --model+--from-csv로 자동)")
    ap.add_argument("--native", help="정답 구조 (기본: run-dir 최상위 *.pdb/*.cif)")
    ap.add_argument("--chains", help="chains.json (기본: run-dir/chains.json)")
    ap.add_argument("--model", help="PRED_GLOB 키 (Boltz-2/Chai-1/Chai-1-noMSA/Protenix/AF2-Multimer/tFold-Ag)")
    ap.add_argument("--target", help="PDB id (예: 9K6J) — from-csv/glob에 사용")
    ap.add_argument("--seed", help="pose seed (미지정+from-csv면 best-DockQ 자동)")
    ap.add_argument("--from-csv", help="summary.csv — best seed·DockQ 자동 선택 + 캡션")
    ap.add_argument("--phase", default="Benchmark")
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", default="", help="캡션 (미지정+from-csv면 자동 생성)")
    ap.add_argument("--surface", action="store_true", help="항원 표면 표시(에피토프 강조)")
    ap.add_argument("--width", type=int, default=1500)
    ap.add_argument("--height", type=int, default=1100)
    a = ap.parse_args(argv_after_dashdash())
    rd = a.run_dir
    low = (a.target or os.path.basename(rd.rstrip("/"))).lower()

    # --- 파일 경로 확정 ---
    native = a.native or pick([f"{rd}/*.pdb", f"{rd}/*.cif"])
    chains_path = a.chains or os.path.join(rd, "chains.json")
    if not native or not os.path.exists(native):
        sys.exit(f"[에러] native 구조 못 찾음 (run-dir={rd}) → --native 로 지정")
    if not os.path.exists(chains_path):
        sys.exit(f"[에러] chains.json 없음: {chains_path}")

    seed, auto_dq = a.seed, None
    if a.from_csv and a.model and a.target:
        b = best_from_csv(a.from_csv, a.phase, a.model, a.target)
        if b:
            if seed is None:
                seed = b[0]
            auto_dq = b[1]
        else:
            print(f"[경고] summary.csv에 {a.model}/{a.target} 행 없음 → seed={seed}")
    if seed is None:
        seed = "0"

    pred = a.pred
    if not pred and a.model:
        if a.model not in PRED_GLOB:
            sys.exit(f"[에러] --model '{a.model}' 미지원. 지원: {list(PRED_GLOB)}")
        g = PRED_GLOB[a.model].format(seed=seed, low=low)
        pred = pick([os.path.join(rd, g)])
        if not pred:
            sys.exit(f"[에러] 예측 pose 못 찾음: {rd}/{g}\n    → --pred 로 직접 지정하세요.")
    elif pred:
        pred = pick([pred]) or pred
    if not pred or not os.path.exists(pred):
        sys.exit(f"[에러] 예측 파일 없음: {pred} (--pred 또는 --model+--from-csv 필요)")

    cj = json.load(open(chains_path))
    ag = str(cj["antigen"])
    ab = cj["antibody"]
    ab = [ab] if isinstance(ab, str) else [str(x) for x in ab]

    # --- 로드 + 항원 기준 정렬 ---
    cmd.load(native, "nat")
    cmd.load(pred, "prd")
    nat_ag_sel = f"nat and chain {ag}"
    if cmd.count_atoms(nat_ag_sel) == 0:
        sys.exit(f"[에러] native에 항원 사슬 {ag} 없음. chains.json 확인.")
    prd_ag = detect_antigen(nat_ag_sel)
    if not prd_ag:
        sys.exit("[에러] 예측 항원 사슬 식별(cealign) 실패.")
    cmd.cealign(nat_ag_sel, f"prd and chain {prd_ag}")  # 최종 정렬(항원 겹침)
    prd_ab = [c for c in chains_of("prd") if c != prd_ag]

    # --- 스타일 ---
    ab_nat = "nat and (" + " or ".join(f"chain {c}" for c in ab) + ")"
    ab_prd = "prd and (" + " or ".join(f"chain {c}" for c in prd_ab) + ")"
    cmd.hide("everything")
    cmd.bg_color("white")
    cmd.set("cartoon_transparency", 0.0)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.set("cartoon_fancy_helices", 1)

    cmd.show("cartoon", nat_ag_sel)
    cmd.color("grey80", nat_ag_sel)
    cmd.select("epi", f"({nat_ag_sel}) within 5 of ({ab_nat})")
    if cmd.count_atoms("epi"):
        cmd.color("firebrick", "epi")
    if a.surface:
        cmd.set("transparency", 0.45)
        cmd.show("surface", nat_ag_sel)

    cmd.show("cartoon", ab_nat)
    cmd.color("skyblue", ab_nat)
    cmd.set("cartoon_transparency", 0.55, ab_nat)

    cmd.show("cartoon", ab_prd)
    cmd.color("orange", ab_prd)

    cmd.hide("everything", f"prd and chain {prd_ag}")   # 예측 항원은 native와 겹치므로 숨김
    cmd.orient(nat_ag_sel + f" or {ab_nat}")

    # --- 렌더 + 저장 ---
    outdir = os.path.dirname(os.path.abspath(a.out))
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    cmd.ray(a.width, a.height)
    cmd.png(a.out, dpi=150)
    pse = a.out.rsplit(".", 1)[0] + ".pse"
    cmd.save(pse)

    label = a.label
    if not label and a.model and a.target:
        tag = ""
        if auto_dq is not None:
            tag = f"  DockQ={auto_dq:.2f}  ({'SUCCESS' if auto_dq >= 0.23 else 'FAIL'})"
        label = f"{a.model}  {a.target}{tag}"
    add_caption(a.out, label)

    print(f"[저장] {a.out}  |  .pse {pse}")
    print(f"  native 항원={ag} → 예측 항원={prd_ag} | 예측 항체={prd_ab} | seed={seed}"
          + (f" | DockQ={auto_dq:.3f}" if auto_dq is not None else ""))


# PyMOL이 스크립트를 실행하면 __name__ == "pymol" (직접 python 실행이면 "__main__")
if __name__ in ("__main__", "pymol"):
    main()
