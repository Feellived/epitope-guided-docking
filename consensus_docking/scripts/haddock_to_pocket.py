#!/usr/bin/env python3
"""HADDOCK 물리 pose → co-folder guided 입력 (diverse). zdock_to_pocket.py의 HADDOCK판.

════════ 왜 이렇게 만들었나 (나중에 이 파일만 봐도 알도록 근거를 남김) ════════

■ 왜 클러스터(seletopclusts) 대표가 아니라 '서로 다른 자리 전부'인가
  옛 RBD HADDOCK-guided(haddock_to_boltz_pocket.py)는 HADDOCK의 FCC 클러스터 대표를 패치로 썼다.
  그런데 클러스터링은 원래 '배포용 선택 장치'(여러 pose가 같은 자리에 모임 = 재현성 있는 결합모드)이고,
  힌트 생성의 목적(= 진짜 자리를 후보 목록에 넣기)과 맞지 않는다. 어느 게 정답인지는 뒤에서
  best-over-patches + 재랭커가 정하면 되므로, 힌트 단계에 필요한 건 '신뢰도'가 아니라 '커버리지'다.

  ⭐ 우리 Exp4 자체 데이터가 이걸 증명함:
    (1) singleton 배제 — "best pose 9/10이 singleton이라 cluster 모드가 배제".
        near-native는 희귀해서 대개 클러스터를 못 이루는데, 클러스터링이 바로 그걸 탈락시킨다
        = 가장 필요한 것을 정확히 버리는 구조.
    (2) 점수 컷 손실 — rigidbody raw엔 near-native가 있는데 seletop(상위 200) 통과 후 학살:
        9ML9 0.99→emref 0.41 · 8SIQ 0.33→0.06 · 9ML8 0.22→0.06.
  즉 클러스터 기반 힌트는 (점수 컷 + 클러스터 컷) 두 번 걸러진 잔여물에서 나온다 → 힌트로 부적합.

■ 왜 stage 기본이 1_rigidbody(정제 전)인가
  위 (2)와 같은 이유. near-native가 살아있는 건 점수 컷 이전의 rigidbody raw다.
  힌트에는 pose 품질이 필요 없고 '접촉 자리(어느 항원 잔기냐)'만 필요하므로,
  정제 안 된 rigidbody 좌표로도 에피토프 위치 추출엔 충분하다 → 커버리지 최대.
  (정제본으로 보고 싶으면 --stage 4_emref)

■ 파일명 숫자는 순위가 아님
  haddock3의 <stage>/xxx_N.pdb.gz 의 N은 생성 순서일 뿐 점수 순위가 아니다(Exp4 확인).
  힌트 목적상 순서는 중요치 않다(어차피 겹치는 자리는 병합하고 서로 다른 자리를 모두 취함).
  다만 --max-patches 상한에 걸릴 때 어떤 자리가 살아남는지에만 영향 → 걱정되면 상한을 올릴 것.

■ 클러스터 정보를 버리라는 뜻은 아님
  '그 자리에 몇 개 pose가 모였나'(클러스터 크기)는 힌트 필터로는 나쁘지만
  나중에 재랭커의 피처로는 유용하다. 여기서 안 쓸 뿐 폐기하지 말 것.

════════ 동작 ════════
  입력: haddock_diverse/<t>/run/<stage>/*.pdb.gz  +  runs_diverse/<t>/chains.json
  처리: pose마다 항원 계면 잔기(=후보 에피토프) 추출 → 겹침 80%↑ 만 병합 → 서로 다른 자리 전부
  출력: runs_haddock_guided/<t>/{boltz_<t>_patchK.yaml, protenix_<t>_patchK.json, tfold_<t>_patchK.pt}
        (emit 포맷·pocket 규약은 zdock_to_pocket.py와 동일 — 함수 재사용)
  이후: 패치마다 co-folder 1회 실행 → best-over-patches. (어느 패치가 맞는지 고르는 건 재랭커 몫)

사용:
  python scripts/haddock_to_pocket.py --targets "9K6J 9Y0A" --emit boltz
  python scripts/haddock_to_pocket.py --targets 9K6J --emit boltz,protenix,tfold --max-patches 20
  python scripts/haddock_to_pocket.py --targets 9K6J --stage 4_emref     # 정제본에서 뽑고 싶을 때
필요: biopython (tfold emit 시 torch 추가).
"""
import argparse, glob, gzip, json, os, re
from Bio.PDB import PDBParser

# 에피토프 추출·emit 로직은 ZDOCK판과 완전히 동일 → 재사용(중복 구현 금지)
from zdock_to_pocket import antigen_epitope, emit_boltz, emit_protenix, emit_tfold

DIVERSE = ["9K6J", "9XF9", "9Y0A", "9VXL", "9LY5", "9MQR", "9KKJ", "9W43", "9YC6", "9J4C"]


def load_gz(p):
    """HADDOCK pose는 .pdb.gz(또는 .pdb)."""
    h = gzip.open(p, "rt") if p.endswith(".gz") else open(p)
    return PDBParser(QUIET=True).get_structure("x", h)[0]


def haddock_poses(haddock_dir, t, stage, pool):
    """haddock_diverse/<t>/run/<stage>/*.pdb[.gz]. 파일명 숫자순(=생성순, 순위 아님). pool>0이면 앞에서 N개만."""
    d = os.path.join(haddock_dir, t.lower(), "run", stage)
    ps = glob.glob(os.path.join(d, "*.pdb.gz")) + glob.glob(os.path.join(d, "*.pdb"))
    def num(p):
        m = re.findall(r"(\d+)", os.path.basename(p))
        return int(m[-1]) if m else 0
    ps.sort(key=num)
    return ps[:pool] if pool > 0 else ps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default=" ".join(DIVERSE), help="쉼표/공백 구분 (기본 diverse 10 전체)")
    ap.add_argument("--haddock-dir", default="haddock_diverse")
    ap.add_argument("--run-dir", default="runs_diverse")
    ap.add_argument("--stage", default="1_rigidbody", help="pose 단계. 기본=정제 전(near-native 생존, 커버리지 최대)")
    ap.add_argument("--emit", default="boltz,protenix,tfold")
    ap.add_argument("--pool", type=int, default=0, help="스캔할 pose 수 상한(0=전부). 느리면 300 정도로 제한")
    ap.add_argument("--max-patches", type=int, default=20, help="서로 다른 자리 최대 개수(연산량 상한)")
    ap.add_argument("--extract-dist", type=float, default=10.0, help="항원 CA가 항체 CA와 이 거리(Å) 내면 계면")
    ap.add_argument("--max-distance", type=float, default=6.0, help="co-folder pocket 제약 거리(Å)")
    ap.add_argument("--min-contacts", type=int, default=4)
    ap.add_argument("--soft", action="store_true", help="Boltz force:false(소프트). 기본=force:true")
    ap.add_argument("--outdir", default="runs_haddock_guided")
    a = ap.parse_args()
    emit = [e.strip() for e in a.emit.replace(",", " ").split()]
    tgts = [x.upper() for x in a.targets.replace(",", " ").split()]

    for t in tgts:
        cjp = os.path.join(a.run_dir, t.lower(), "chains.json")
        if not os.path.exists(cjp):
            print(f"{t}: chains.json 없음 skip"); continue
        cj = json.load(open(cjp))
        ag = str(cj["antigen"]); ab = cj["antibody"]; ab = [ab] if isinstance(ab, str) else [str(x) for x in ab]
        binder = str(ab[0])
        seqmap = {c["id"]: c["seq"] for c in cj["chains"]}
        ag_seq = seqmap.get(ag, "")

        poses = haddock_poses(a.haddock_dir, t, a.stage, a.pool)
        if not poses:
            print(f"{t}: {a.haddock_dir}/{t.lower()}/run/{a.stage}/ pose 없음 (HADDOCK 먼저) skip"); continue

        # 서로 다른 자리 전부 — 겹침 80%↑ 만 병합(중복 제거일 뿐 클러스터링 아님)
        patches = []
        for p in poses:
            try:
                epi = antigen_epitope(load_gz(p), ag_seq, a.extract_dist)
            except Exception:
                continue
            if len(epi) < a.min_contacts:
                continue
            if any(len(epi & q) / max(1, len(epi | q)) > 0.8 for q in patches):
                continue
            patches.append(epi)
            if len(patches) >= a.max_patches:
                break
        if not patches:
            print(f"{t}: 유효 패치 없음(min-contacts={a.min_contacts}) skip"); continue

        outd = os.path.join(a.outdir, t.lower()); os.makedirs(outd, exist_ok=True)
        wrote_fasta = [False]
        for k, pos in enumerate(patches):
            if "boltz" in emit:
                emit_boltz(outd, t, cj["chains"], ag, binder, pos, a.max_distance, a.soft, k)
            if "protenix" in emit:
                emit_protenix(outd, t, cj["chains"], ag, binder, pos, a.max_distance, k)
            if "tfold" in emit:
                emit_tfold(outd, t, seqmap, ag, ab, pos, k, wrote_fasta)
        sizes = ",".join(str(len(p)) for p in patches)
        print(f"{t}: {a.stage} {len(poses)} pose에서 서로 다른 자리 {len(patches)}개(잔기수 {sizes}) "
              f"emit={emit} (binder={binder}, 항원={ag}, force={'false' if a.soft else 'true'}) → {outd}/")


if __name__ == "__main__":
    main()
