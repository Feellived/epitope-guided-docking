#!/usr/bin/env python3
"""패치별 HADDOCK 제약(ambig) + run.cfg (기존 run.cfg 템플릿 재사용: sampling·run_dir·ambig·select만 교체).
   CDR active(antibody.actpass) ↔ 패치잔기 active(chain B). 사용(haddock3 env):
   python scripts/make_haddock_patch_restraints.py --dir haddock/8sit --samples 100 --refine 20
   (2026-07-16: --tolerance 추가 — flexref/emref CNS 정제 5% abort 방지, 기본 50%)"""
import argparse, glob, os, re, shutil, subprocess

def restraint_cmd(ab_ap, ag_ap):
    """active_passive_to_ambig 호출 커맨드 (standalone or haddock3-restraints 서브커맨드)."""
    tail = [ab_ap, ag_ap, "--segid-one", "A", "--segid-two", "B"]
    if shutil.which("active_passive_to_ambig"):
        return ["active_passive_to_ambig"] + tail
    if shutil.which("haddock3-restraints"):
        return ["haddock3-restraints", "active_passive_to_ambig"] + tail
    raise SystemExit("active_passive_to_ambig / haddock3-restraints 둘 다 없음 — HADDOCK3 env 확인(HAD_ENV).")

def _has_param_in_block(text, mod, param):
    """cfg의 [mod] 블록 안에 param이 이미 있나 (중복 삽입 방지)."""
    m = re.search(r'\[' + re.escape(mod) + r'\](.*?)(\n\[|\Z)', text, re.S)
    return bool(m and re.search(r'\b' + re.escape(param) + r'\s*=', m.group(1)))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True); ap.add_argument("--samples", type=int, default=100)
    ap.add_argument("--refine", type=int, default=20)
    ap.add_argument("--tolerance", type=int, default=50,
                    help="CNS 정제(flexref/emref/rigidbody) 실패 허용 %% (기본 5%%면 소수 실패에 abort)")
    a = ap.parse_args()
    d = a.dir
    tmpl = os.path.join(d, "run.cfg")
    if not os.path.exists(tmpl):
        cands = [c for c in glob.glob(os.path.join(d, "*.cfg")) if "patch" not in c]
        tmpl = cands[0] if cands else None
    if not tmpl: raise SystemExit(f"{d}/run.cfg (템플릿) 없음 — make_haddock_restraints 먼저")
    cfg_text = open(tmpl).read()
    ab_actpass = os.path.join(d, "antibody.actpass")
    if not os.path.exists(ab_actpass): raise SystemExit(f"{ab_actpass} 없음")
    patches = sorted(glob.glob(os.path.join(d, "patches", "patch_*.txt")))
    if not patches: raise SystemExit(f"{d}/patches/patch_*.txt 없음 (partition_surface 먼저)")
    for pf in patches:
        pi = os.path.basename(pf).replace("patch_", "").replace(".txt", "")
        resid = open(pf).read().split()
        if not resid: continue
        ag_ap = os.path.join(d, f"antigen_patch_{pi}.actpass")
        open(ag_ap, "w").write(" ".join(resid) + "\n\n")   # active=패치잔기, passive 빈줄
        ambig = os.path.join(d, f"ambig_patch_{pi}.tbl")
        with open(ambig, "w") as f:
            subprocess.run(restraint_cmd(ab_actpass, ag_ap), stdout=f, check=True)
        txt = cfg_text
        txt = re.sub(r'run_dir\s*=\s*"[^"]*"', 'run_dir = "%s"' % os.path.abspath(os.path.join(d, "patch_"+pi, "run")), txt)
        txt = re.sub(r'\bsampling\s*=\s*\d+', 'sampling = %d' % a.samples, txt)
        txt = re.sub(r'ambig_fname\s*=\s*"[^"]*"', 'ambig_fname = "%s"' % os.path.abspath(ambig), txt)
        if re.search(r'\bselect\s*=\s*\d+', txt):
            txt = re.sub(r'\bselect\s*=\s*\d+', 'select = %d' % a.refine, txt)
        elif '[seletop]' in txt:
            txt = txt.replace('[seletop]', '[seletop]\nselect = %d' % a.refine, 1)
        # CNS 정제 실패 허용치 상향 (flexref/emref/rigidbody 의 '5% 초과 실패 → abort' 방지, 멱등)
        for mod in ("rigidbody", "flexref", "emref"):
            if ('[%s]' % mod) in txt and not _has_param_in_block(txt, mod, "tolerance"):
                txt = txt.replace('[%s]' % mod, '[%s]\ntolerance = %d' % (mod, a.tolerance), 1)
        outcfg = os.path.join(d, f"run_patch_{pi}.cfg")
        open(outcfg, "w").write(txt)
        print(f"patch_{pi}: {len(resid)}잔기 → {outcfg} (sampling {a.samples}, refine {a.refine}, tolerance {a.tolerance}%)")

if __name__ == "__main__": main()
