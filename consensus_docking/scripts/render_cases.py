#!/usr/bin/env python3
"""[PyMOL] 한 복합체를 여러 조건으로 나란히 그린다 — 카메라 각도를 고정해서.

각 패널은 같은 항원 위에 조건만 다른 예측 항체를 얹은 것이다. 각도가 다르면 비교가
성립하지 않으므로 **정답 구조에서 각도를 한 번 잡고 모든 패널에 그대로 쓴다.**

  회색      항원 표면
  진홍      정답 결합 자리 (정답 구조에서 항체와 5 Å 이내로 닿는 항원 잔기)
  연하늘    정답 항체 (반투명 — '있어야 할 자리')
  주황      예측된 항체

보는 법은 초록 면적이 아니라 **주황이 진홍 위에 앉았는가**이다.

⚠️ boltz 환경이 아니라 pymol 환경에서 돈다. 우리 코드에 의존하지 않는다.

사용:
  ~/miniconda3/envs/pymol/bin/pymol -cq render_cases.py -- \\
      --target 9zdu_HL --arms "noconstraint ours" --out fig/9zdu

  # 네 조건 전부
  ... -- --target 8q7s_H --arms "noconstraint fullmsa sizematch ours" --out fig/8q7s_h
"""
import argparse, glob, json, os, sys

from pymol import cmd

MSAD = os.path.expanduser("~/projects/bk21-msa-depth-bias/pipeline")
RUNS = os.path.expanduser("~/projects/bk21-antibody-ml/consensus_docking/runs_sites_guided_honest")
LABEL = {"noconstraint": "제약 없음", "fullmsa": "원래 MSA 자리",
         "sizematch": "같은 크기 다른 자리", "ours": "조성이 고른 자리"}


def argv_after_dashdash():
    a = sys.argv
    return a[a.index("--") + 1:] if "--" in a else a[1:]


def best_pose(tgt, arm, pick):
    """그 팔의 자세 파일 하나. pick='first' 면 model_0, 숫자면 그 번호."""
    pat = os.path.join(RUNS, tgt.lower(), arm, "out_boltz", "results",
                       "**", "*model*.cif")
    hits = sorted(glob.glob(pat, recursive=True))
    if not hits:
        return None
    if pick.isdigit():
        want = [h for h in hits if f"model_{pick}." in h]
        if want:
            return want[0]
    return hits[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--arms", default="noconstraint ours")
    ap.add_argument("--pick", default="first", help="자세 번호. 기본은 첫 번째")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--out", required=True, help="접두사. <접두사>_<팔>.png 로 저장")
    ap.add_argument("--width", type=int, default=1200)
    ap.add_argument("--height", type=int, default=1000)
    a = ap.parse_args(argv_after_dashdash())

    tgt = a.target
    cj = json.load(open(os.path.join(MSAD, "targets", tgt, "chains.json")))
    ag = cj["antigen"]; ag = [ag] if isinstance(ag, str) else [str(x) for x in ag]
    ab = cj["antibody"]; ab = [ab] if isinstance(ab, str) else [str(x) for x in ab]
    native = os.path.join(MSAD, "targets", tgt, "native.cif")
    if not os.path.exists(native):
        sys.exit(f"!! native 가 없다: {native}")

    cmd.set("ray_opaque_background", 1)
    cmd.bg_color("white")
    cmd.set("surface_quality", 1)
    cmd.set("transparency", 0.55)
    cmd.set("cartoon_transparency", 0.55, "nat_ab")

    # ── 정답 구조로 각도를 한 번 잡는다 ────────────────────────────────────
    cmd.load(native, "nat")
    ag_sel = " or ".join(f"chain {c}" for c in ag)
    ab_sel = " or ".join(f"chain {c}" for c in ab)
    cmd.create("nat_ag", f"nat and polymer and ({ag_sel})")
    cmd.create("nat_ab", f"nat and polymer and ({ab_sel})")
    cmd.delete("nat")
    cmd.select("epi", f"nat_ag within {a.cutoff} of nat_ab")
    n_epi = cmd.count_atoms("epi and name CA")
    print(f"[{tgt}] 항원={ag} 항체={ab} · 정답 결합자리 {n_epi}잔기")
    if n_epi == 0:
        sys.exit("!! 정답 결합자리가 비었다 — chains.json 의 사슬 기호를 확인할 것")

    cmd.hide("everything")
    cmd.show("surface", "nat_ag")
    cmd.color("grey80", "nat_ag")
    cmd.color("deeppink", "epi")
    cmd.show("cartoon", "nat_ab")
    cmd.color("lightblue", "nat_ab")

    cmd.orient("nat_ag or nat_ab")        # ★ 여기서 잡은 각도를 끝까지 쓴다
    cmd.zoom("nat_ag or nat_ab", 3)
    VIEW = cmd.get_view()

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    made = []
    for arm in a.arms.split():
        p = best_pose(tgt, arm, a.pick)
        if not p:
            print(f"  ! {arm}: 자세 파일이 없다 — 건너뜀"); continue
        cmd.delete("prd*")
        cmd.load(p, "prd")
        # 예측 항원을 정답 항원에 겹쳐 같은 좌표계로 옮긴다
        cmd.create("prd_ag", f"prd and polymer and ({ag_sel})")
        cmd.create("prd_ab", f"prd and polymer and ({ab_sel})")
        if cmd.count_atoms("prd_ag") == 0 or cmd.count_atoms("prd_ab") == 0:
            print(f"  ! {arm}: 예측에서 항원/항체 사슬을 못 찾았다 "
                  f"(사슬 = {sorted(set(cmd.get_chains('prd')))}) — 건너뜀")
            continue
        try:
            rms = cmd.cealign("nat_ag", "prd_ag")["RMSD"]
        except Exception:
            rms = cmd.align("prd_ag", "nat_ag")[0]
        cmd.matrix_copy("prd_ag", "prd_ab")

        cmd.hide("everything", "prd*")
        cmd.show("cartoon", "prd_ab")
        cmd.color("orange", "prd_ab")
        cmd.set("cartoon_transparency", 0, "prd_ab")

        cmd.set_view(VIEW)                 # ★ 각도 복원 — 패널마다 동일
        png = f"{a.out}_{arm}.png"
        cmd.png(png, width=a.width, height=a.height, dpi=200, ray=1)
        print(f"  {arm:<14} {LABEL.get(arm, arm):<16} 정렬 RMSD {rms:.2f} Å  → {png}")
        made.append(png)

    print(f"\n완료 {len(made)}장. 가로로 붙이려면:")
    print(f"  montage {' '.join(made)} -tile {len(made)}x1 -geometry +6+6 {a.out}_all.png")


main()
