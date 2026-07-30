#!/usr/bin/env python3
"""[PyMOL] 항원 표면에 '정답 결합자리'와 '예측 결합자리'를 직접 칠해 자리 이동을 보인다.

왜 이 그림인가:
  이 연구의 주장은 자세(DockQ)가 아니라 **자리**(epitope recall)다. 항체를 통째로 그리면
  자세가 틀린 경우 그림이 설득력을 잃지만, **항원 표면만 놓고 두 패치를 칠하면**
  "옳은 자리로 옮겨갔다"가 그대로 보인다. 색칠된 넓이가 곧 recall 이다.

색 (파스텔 · 그림 안에 글자를 넣지 않으므로 캡션에 옮겨 적을 것):
  회색   = 항원 표면
  분홍   = 정답 결합자리인데 예측이 놓친 부분      (false negative)
  파랑   = 예측이 짚었으나 정답이 아닌 부분        (false positive)
  초록   = 둘이 겹친 부분 = 맞힌 부분              (true positive)
  → 초록이 많고 분홍이 적을수록 좋다.

⚠️ 여러 조건을 **한 세션에서** 그려 카메라 각도를 동일하게 고정한다. 각도가 다르면
   전/후 비교가 성립하지 않는다.

헤드리스 실행:
  pymol -cq pymol_epitope.py -- --native targets/8k3k_D/native.cif \\
        --chains targets/8k3k_D/chains.json --out-prefix fig/site_8k3k_D \\
        --pred "no constraint=runs/.../a.cif" --pred "composition site=runs/.../b.cif"
"""
import argparse
import json
import os
import sys

from pymol import cmd

GREY, MISS, FALSE, HIT = "grey80", "0xEDA8C4", "0x9EC5E8", "0xA6D8B4"


def argv_after_dashdash():
    a = sys.argv
    return a[a.index("--") + 1:] if "--" in a else a[1:]


def chains_of(sel):
    return [c for c in cmd.get_chains(sel) if c and c.strip()]


def load_chains(cj):
    """native의 실제 사슬 ID. msa-depth는 src_chains, consensus_docking은 antigen/antibody."""
    def lst(v):
        return [v] if isinstance(v, str) else [str(x) for x in v]
    src = cj.get("src_chains") or {}
    ag = [str(x) for x in src.get("antigen", [])] or lst(cj["antigen"])
    ab = [str(x) for x in (list(src.get("H", [])) + list(src.get("L", [])))] or lst(cj["antibody"])
    return ag, ab


def detect_antigen(obj, nat_ag_sel):
    """예측 구조에서 항원 사슬을 고른다 — native 항원과 CA 수가 가장 가까운 사슬."""
    n = cmd.count_atoms(nat_ag_sel + " and name CA")
    best = None
    for c in chains_of(obj):
        k = cmd.count_atoms(f"{obj} and chain {c} and name CA")
        if not k:
            continue
        try:
            rms = cmd.cealign(nat_ag_sel, f"{obj} and chain {c}").get("RMSD", 99.0)
        except Exception:
            rms = 99.0
        cand = (abs(k - n), rms, c)
        if best is None or cand < best:
            best = cand
    return best[2] if best else None


def resset(sel):
    out = set()
    cmd.iterate(sel + " and name CA", "out.add((chain, resi))", space={"out": out})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--native", required=True)
    ap.add_argument("--chains", required=True)
    ap.add_argument("--pred", action="append", default=[],
                    help='"라벨=경로" 형식. 여러 번 줄 수 있고 모두 같은 각도로 그린다')
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--width", type=int, default=1300)
    ap.add_argument("--height", type=int, default=1100)
    a = ap.parse_args(argv_after_dashdash())

    cj = json.load(open(a.chains))
    agl, abl = load_chains(cj)
    cmd.load(a.native, "nat")
    ag_sel = "nat and (" + " or ".join(f"chain {c}" for c in agl) + ")"
    ab_sel = "nat and (" + " or ".join(f"chain {c}" for c in abl) + ")"
    if cmd.count_atoms(ag_sel) == 0 or cmd.count_atoms(ab_sel) == 0:
        sys.exit(f"[에러] 사슬 매핑 실패. native 사슬 = {chains_of('nat')} · 기대 = {agl}/{abl}")

    # 정답 결합자리
    cmd.select("true_epi", f"byres (({ag_sel}) within {a.cutoff} of ({ab_sel}))")
    true_res = resset("true_epi")
    if not true_res:
        sys.exit("[에러] 정답 결합자리를 찾지 못했다.")

    # 표시 기본값 — 항원 표면만 남긴다
    cmd.hide("everything")
    cmd.bg_color("white")
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.set("surface_quality", 1)
    cmd.set("ray_opaque_background", 1)
    cmd.show("surface", ag_sel)
    cmd.color(GREY, ag_sel)

    # 카메라는 정답 자리가 정면에 오도록 한 번만 정하고 모든 칸에 같은 값을 쓴다
    cmd.orient("true_epi")
    cmd.zoom(ag_sel, 3.0)
    view = cmd.get_view()

    print(f"정답 결합자리 {len(true_res)}잔기 · 조건 {len(a.pred)}개")
    for spec in a.pred:
        label, path = spec.split("=", 1)
        tag = "".join(ch if ch.isalnum() else "_" for ch in label)
        cmd.delete("prd")
        cmd.load(path, "prd")
        pag = detect_antigen("prd", ag_sel)
        if not pag:
            print(f"  !! {label}: 예측 항원 사슬 식별 실패 — 건너뜀"); continue
        cmd.cealign(ag_sel, f"prd and chain {pag}")
        pab = [c for c in chains_of("prd") if c != pag]
        if not pab:
            print(f"  !! {label}: 예측 항체 사슬 없음 — 건너뜀"); continue
        pab_sel = "prd and (" + " or ".join(f"chain {c}" for c in pab) + ")"

        # 예측 결합자리 = native 항원 표면 위에서 예측 항체와 가까운 잔기
        cmd.select("pred_epi", f"byres (({ag_sel}) within {a.cutoff} of ({pab_sel}))")
        pred_res = resset("pred_epi")
        hit = true_res & pred_res
        rec = len(hit) / len(true_res)
        prec = len(hit) / len(pred_res) if pred_res else 0.0

        cmd.color(GREY, ag_sel)
        cmd.color(MISS, "true_epi")
        cmd.color(FALSE, "pred_epi and not true_epi")
        cmd.color(HIT, "pred_epi and true_epi")
        cmd.set_view(view)
        out = f"{a.out_prefix}_{tag}.png"
        os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
        cmd.png(out, width=a.width, height=a.height, dpi=300, ray=1)
        print(f"  {label:24} 예측 {len(pred_res):>3}잔기 · 맞힌 {len(hit):>3} · "
              f"recall {rec:.3f} · precision {prec:.3f} → {out}")

    # 정답만 칠한 참조 그림 한 장
    cmd.color(GREY, ag_sel)
    cmd.color(MISS, "true_epi")
    cmd.set_view(view)
    cmd.png(f"{a.out_prefix}_reference.png", width=a.width, height=a.height, dpi=300, ray=1)
    print(f"  {'true epitope only':24} → {a.out_prefix}_reference.png")


main()
