#!/usr/bin/env python3
"""analyze_haddock_score_rank.py의 'DockQ 전부 실패' 원인 진단. 예외 안 삼키고 각 단계 출력.
   사용: python scripts/analyze_score_rank_debug.py [8SDF] [--stage 4_emref]"""
import argparse, csv, glob, os, subprocess, sys, tempfile, traceback
import lib_structure as B
import haddock_dockq as HD


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target", nargs="?", default="8SDF")
    ap.add_argument("--stage", default="4_emref")
    ap.add_argument("--n", type=int, default=3, help="검사할 pose 수")
    a = ap.parse_args()
    t = a.target.upper()
    rundir = os.path.join("haddock", t.lower(), "run"); emdir = os.path.join(rundir, a.stage)
    print(f"===== 진단: {t} (stage={a.stage}) =====")

    # 1) 실제 emref 파일 vs TSV model_name 대조
    real = sorted(os.path.basename(p) for p in glob.glob(os.path.join(emdir, "*.pdb.gz")))
    print(f"[1] 실제 {a.stage}/*.pdb.gz 파일: {len(real)}개  예: {real[:3]}")
    tsvs = glob.glob(os.path.join(rundir, "*clustfcc*", "*.tsv")) + glob.glob(os.path.join(rundir, "**", "*.tsv"), recursive=True)
    tsvs = [f for f in dict.fromkeys(tsvs)]
    tsv = None
    for f in tsvs:
        try:
            hdr = next(csv.reader(open(f), delimiter="\t"))
            if {"model_name", "score", "cluster_id"} <= set(hdr): tsv = f; break
        except Exception: pass
    print(f"[2] clustfcc TSV: {tsv}")
    if tsv is None: print("    ✗ TSV 없음 → 종료"); return
    rows = list(csv.DictReader(open(tsv), delimiter="\t"))
    mns = [r["model_name"] for r in rows]
    print(f"    TSV 행 {len(rows)}개  model_name 예: {mns[:3]}")
    # 경로 구성 3가지 방식 존재 여부
    for m in mns[:3]:
        c_gz = os.path.join(emdir, m + ".gz"); c_bare = os.path.join(emdir, m)
        c_base = os.path.join(emdir, os.path.basename(m) + (".gz" if not m.endswith(".gz") else ""))
        print(f"    '{m}': +'.gz' 존재={os.path.exists(c_gz)} | bare 존재={os.path.exists(c_bare)} | basename+gz 존재={os.path.exists(c_base)}")
    nres = sum(1 for m in mns if os.path.exists(os.path.join(emdir, m + ".gz")) or os.path.exists(os.path.join(emdir, m)))
    print(f"    → TSV model_name 중 emdir 파일로 해석되는 수: {nres}/{len(mns)}")

    # 3) native
    ni = B.native_info(t); print(f"[3] native_info: {'OK' if ni else 'None'}")
    if ni is None: return
    with tempfile.TemporaryDirectory() as td:
        try:
            natm = HD.native_merged(t, td)
        except Exception:
            print("    ✗ native_merged 예외:"); traceback.print_exc(); return
        print(f"    native_merged → {natm}  존재={os.path.exists(natm) if natm else None}  "
              f"크기={os.path.getsize(natm) if natm and os.path.exists(natm) else '-'}")
        if natm is None: return

        # 4) pose N개를 단계별로 (예외 그대로 노출)
        for m in mns[:a.n]:
            print(f"\n[4] pose model_name='{m}'")
            pose = os.path.join(emdir, m + ".gz")
            if not os.path.exists(pose): pose = os.path.join(emdir, m)
            if not os.path.exists(pose):
                print(f"    ✗ 파일 없음: {pose}"); continue
            try:
                pm = HD.load_gz(pose); print(f"    load_gz OK  chains={[c.id for c in pm]}")
                rec, pid, abids = HD.id_and_recall(pm, ni)
                print(f"    id_and_recall OK  pid={pid}  abids={abids}  recall={rec:.2f}")
                if not abids: print("    ✗ abids 비어있음 → None 반환됨"); continue
                mp = os.path.join(td, "m.pdb"); HD.write_merged(pm, pid, abids, mp)
                print(f"    write_merged OK  크기={os.path.getsize(mp)}")
                r = subprocess.run(["DockQ", mp, natm], capture_output=True, text=True, timeout=180)
                print(f"    DockQ returncode={r.returncode}")
                print(f"    STDOUT(앞600):\n{(r.stdout or '')[:600]}")
                if r.stderr: print(f"    STDERR(앞600):\n{r.stderr[:600]}")
            except Exception:
                print("    ✗ 예외 발생:"); traceback.print_exc()

    # 5) DockQ 버전
    try:
        v = subprocess.run(["DockQ", "--version"], capture_output=True, text=True, timeout=30)
        print(f"\n[5] DockQ --version: {(v.stdout or v.stderr).strip()[:200]}")
    except Exception as e:
        print(f"[5] DockQ --version 실패: {e}")


if __name__ == "__main__": main()
