#!/usr/bin/env python3
"""Boltz 예측 1건 평가 -> DockQ + 신뢰도 파싱 -> summary.csv 자동 기록.

서버(conda dock)에서 실행. 토큰 불필요.
필요: DockQ (CLI), pyyaml(선택; targets.yaml 읽을 때). run-dir 안에
  <TARGET>.pdb (정답),  boltz_<target>.yaml (입력, 사슬 순서용),
  *boltz/**/*_model_<seed>.cif,  confidence_*_model_<seed>.json

예:
  python eval_boltz.py --target 9VP0 --run-dir . --targets ../targets.yaml \
      --seed 0 --phase Preliminary --model Boltz-2 --csv ../results/summary.csv
targets.yaml 없거나 못 읽으면:
  ... --antigen A --antibody H,L
"""
import argparse, glob, json, os, re, subprocess, sys
import lib_log


def find_one(patterns):
    for p in patterns:
        hits = sorted(glob.glob(p, recursive=True))
        if hits:
            return hits[0]
    return None


def load_targets(path, tid):
    try:
        import yaml
    except ImportError:
        return None, None
    with open(path) as f:
        data = yaml.safe_load(f)
    for grp in (data or {}).values():
        for t in grp:
            if str(t.get("id", "")).upper() == tid.upper():
                antibody = list(t["antibody"].values())
                antigen = t["antigen"]["chain"]
                return antigen, antibody
    return None, None


def yaml_chain_order(path):
    """boltz yaml의 id: 등장 순서 = pair_chains_iptm 인덱스 순서."""
    order = []
    with open(path) as f:
        for line in f:
            m = re.search(r"id:\s*([A-Za-z0-9])", line)
            if m:
                order.append(m.group(1))
    return order


def _interfaces(out):
    inter, cur = {}, None
    for line in out.splitlines():
        m = re.search(r"Native chains:\s*([A-Za-z0-9]),\s*([A-Za-z0-9])", line)
        if m:
            cur = frozenset([m.group(1), m.group(2)])
            continue
        m = re.search(r"^\s*DockQ:\s*([\d.]+)", line)
        if m and cur is not None:
            inter[cur] = float(m.group(1))
            cur = None
    return inter


def parse_dockq(model_cif, native_pdb, mapping=None):
    def run(mp):
        cmd = ["DockQ", model_cif, native_pdb]
        if mp:
            cmd += ["--mapping", mp]
        return subprocess.run(cmd, capture_output=True, text=True).stdout
    out = run(mapping)                       # 예측 사슬 -> native 같은 사슬 고정 (단일카피·cognate 다중카피)
    inter = _interfaces(out)
    if not inter and mapping:                # 강제매핑이 native 계면 못 찾음(다중카피 비-cognate) -> auto 재시도
        out2 = run(None)
        i2 = _interfaces(out2)
        if i2:
            return i2, out2
    return inter, out


def parse_conf(conf_path, order, antigen, antibody):
    with open(conf_path) as f:
        c = json.load(f)
    giptm = c.get("iptm")
    pci = c.get("pair_chains_iptm", {}) or {}
    idx = {ch: i for i, ch in enumerate(order)}
    ai = idx.get(antigen)
    vals = []
    for ab in antibody:
        bi = idx.get(ab)
        if ai is None or bi is None:
            continue
        for a, b in [(ai, bi), (bi, ai)]:
            try:
                vals.append(float(pci[str(a)][str(b)]))
            except (KeyError, TypeError, ValueError):
                pass
    return giptm, (min(vals) if vals else None), vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)
    ap.add_argument("--run-dir", default=".")
    ap.add_argument("--targets", default="targets.yaml")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--phase", default="Benchmark")
    ap.add_argument("--model", default="Boltz-2")
    ap.add_argument("--csv", default="results/summary.csv")
    ap.add_argument("--antigen")
    ap.add_argument("--antibody")
    a = ap.parse_args()
    rd, tid = a.run_dir, a.target

    antigen = antibody = None
    if os.path.exists(a.targets):
        antigen, antibody = load_targets(a.targets, tid)
    if a.antigen:
        antigen = a.antigen
    if a.antibody:
        antibody = a.antibody.split(",")
    if not antigen or not antibody:
        sys.exit("[에러] antigen/antibody 사슬 미상 -> --antigen A --antibody H,L 로 지정")

    native = find_one([f"{rd}/{tid}.pdb", f"{rd}/{tid}.cif",
                       f"{rd}/{tid.upper()}.pdb", f"{rd}/{tid.upper()}.cif",
                       f"{rd}/{tid.lower()}.pdb", f"{rd}/{tid.lower()}.cif"])
    yml = find_one([f"{rd}/boltz_{tid.lower()}.yaml", f"{rd}/boltz_*.yaml"])
    cif = find_one([f"{rd}/*boltz/**/*_model_{a.seed}.cif"])
    conf = find_one([f"{rd}/*boltz/**/confidence_*_model_{a.seed}.json"])
    miss = [n for n, v in [("native", native), ("yaml", yml), ("cif", cif), ("conf", conf)] if not v]
    if miss:
        sys.exit(f"[에러] 파일 못 찾음: {miss} (run-dir={rd})")

    order = yaml_chain_order(yml)
    _m = "".join(order)
    inter, dockq_out = parse_dockq(cif, native, mapping=f"{_m}:{_m}")
    giptm, pair, pair_vals = parse_conf(conf, order, antigen, antibody)

    abag = {ab: inter.get(frozenset([antigen, ab])) for ab in antibody}
    have = [v for v in abag.values() if v is not None]
    if not have and len(antibody) == 1 and len(inter) == 1:
        # 나노바디 + 다중카피: auto-map이 다른 카피 사슬로 맵핑 -> 유일 계면이 곧 Ab-Ag
        only = next(iter(inter.values()))
        abag = {antibody[0]: only}
        have = [only]
    if not have:
        sys.exit(f"[에러] Ab-Ag interface DockQ 파싱 실패. antigen={antigen} antibody={antibody}\n"
                 f"DockQ가 본 interface: {[tuple(k) for k in inter]}\n--- DockQ 원문 ---\n{dockq_out}")
    dockq_min = min(have)
    success = dockq_min >= 0.23

    detail = " / ".join(f"{antigen}-{ab} {abag[ab]:.3f}" if abag[ab] is not None else f"{antigen}-{ab} NA"
                        for ab in antibody)
    pr = f"{min(pair_vals):.2f}-{max(pair_vals):.2f}" if pair_vals else "NA"
    notes = f"DockQ(Ab-Ag) {detail}; pair ipTM {pr}; interfaces={len(inter)}"

    row = {"phase": a.phase, "model": a.model, "target": tid, "seed": a.seed,
           "dockq_abag": round(dockq_min, 3),
           "global_iptm": round(giptm, 3) if giptm is not None else "",
           "pair_iptm_abag": round(pair, 3) if pair is not None else "",
           "success": success, "notes": notes}
    saved = lib_log.append_row(a.csv, row)
    print("=== 기록됨 ->", a.csv, "===")
    for k in lib_log.FIELDS:
        print(f"  {k}: {saved[k]}")


if __name__ == "__main__":
    main()
