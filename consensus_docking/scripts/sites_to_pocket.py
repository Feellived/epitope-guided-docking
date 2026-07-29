#!/usr/bin/env python3
"""[유도 입력] MSA 조성이 만든 결합자리 후보 → co-folder pocket 제약 입력.

haddock_to_pocket.py / zdock_to_pocket.py 의 형제. 힌트의 **출처만 다르다**:
  · zdock/haddock 판 = 물리 도커 pose 에서 뽑은 표면 패치
  · 이 파일        = bk21-msa-depth-bias 본 검정에서 'MSA 조성을 바꿔' 얻은 결합자리 후보
emit 포맷·pocket 규약은 zdock_to_pocket.py 의 함수를 그대로 재사용한다(중복 구현 금지).

═══ 세 팔(arm) — 왜 셋인가 ═══
"제약을 주면 좋아진다"는 이미 알려진 사실이라 그것만 보이면 우리 기여가 증명되지 않는다.
증명해야 하는 것은 **조성 다양화가 옳은 자리를 후보에 올렸다**이므로 대조가 하나 더 필요하다.
  (a) blind    제약 없음. 기존 예측 그대로(여기서 만들지 않는다)
  (b) fullmsa  **원래 MSA 가 간 자리**로 유도 = 우리 파이프라인 없이 얻는 최선의 힌트
  (c) ours     우리가 고른 후보로 유도
판정은 (c) > (b). (c) > (a) 만으로는 '제약 효과'와 구별되지 않는다.

═══ 잔기 번호 변환 ═══
msa-depth 쪽 잔기 키 = (항원 사슬 순번, 그 사슬 참조서열의 0-based 위치)
co-folder pocket 제약 = (사슬 id, 서열 1-based 위치)
→ 사슬 순번을 chains.json 의 role=="antigen" 순서로 id 에 대응시키고, 위치에 +1.

사용(consensus_docking/ 에서):
  python scripts/sites_to_pocket.py --sites ~/projects/bk21-msa-depth-bias/pipeline/results/sites_8k3k_D.json \\
      --chains ~/projects/bk21-msa-depth-bias/pipeline/targets/8k3k_D/chains.json \\
      --arm ours --cand 1 --emit boltz,protenix
  python scripts/sites_to_pocket.py --sites ... --chains ... --arm fullmsa --emit boltz,protenix
필요: biopython 불필요(좌표를 안 읽는다). tfold emit 시 torch.
"""
import argparse, json, os, random, sys

from zdock_to_pocket import emit_boltz, emit_protenix, emit_tfold


def load_chains(path):
    """msa-depth targets/<t>/chains.json → (co-folder emit 이 쓰는 형식, 항원 사슬 id 목록).

    msa-depth 형식: {"chains":[{"id","seq","role":"antigen"|"heavy"|"light"}, ...]}
    emit_* 가 기대하는 형식: chains=[{"id","seq"}...] + 항원 id + binder(항체 첫 사슬) id
    """
    cj = json.load(open(path))
    chains = [{"id": str(c["id"]), "seq": c["seq"]} for c in cj["chains"]]
    ag_ids = [str(c["id"]) for c in cj["chains"] if c.get("role") == "antigen"]
    ab_ids = [str(c["id"]) for c in cj["chains"] if c.get("role") in ("heavy", "light")]
    if not ag_ids or not ab_ids:
        sys.exit(f"!! {path} 에서 항원/항체 사슬을 못 찾았다 (role 키 확인)")
    return chains, ag_ids, ab_ids


def emit_one(a, d, t, c, chains, ag_ids, ab_ids, k):
    """후보 하나 → co-folder 입력. 파일명 접미사 patch<k>. 성공 1 / 건너뜀 0."""
    emit = [e.strip() for e in a.emit.replace(",", " ").split()]

    # (사슬 순번, 0-based) → 사슬 id 별 1-based 위치
    bych = {}
    for ordi, pos0 in c["residues"]:
        if ordi >= len(ag_ids):
            print(f"  ! 후보 {c['cand']}: 항원 사슬 순번 {ordi} 가 chains.json 의 항원 수"
                  f"({len(ag_ids)})를 넘는다 — 건너뜀")
            return 0
        bych.setdefault(ag_ids[ordi], []).append(int(pos0) + 1)
    if not bych:
        print(f"  ! 후보 {c['cand']}: 잔기가 없다 — 건너뜀")
        return 0

    if len(bych) > 1:
        # emit_* 는 pocket contacts 를 사슬 하나로 받는다. 항원이 여러 사슬이면 잔기가 가장 많은
        # 사슬만 쓰고 그 사실을 알린다(데모 타깃은 단일 항원 사슬로 고르는 것을 권장).
        keep = max(bych, key=lambda x: len(bych[x]))
        print(f"  ⚠️ 항원 사슬이 여러 개다({', '.join(bych)}) — 잔기가 가장 많은 {keep} 만 제약으로 쓴다 "
              f"({len(bych[keep])}/{sum(len(v) for v in bych.values())}잔기).")
        bych = {keep: bych[keep]}

    ag = list(bych)[0]
    pos = sorted(bych[ag])
    binder = ab_ids[0]
    outd = os.path.join(a.outdir, t.lower(), a.arm)
    os.makedirs(outd, exist_ok=True)
    seqmap = {ch["id"]: ch["seq"] for ch in chains}
    wrote = [False]
    if "boltz" in emit:
        emit_boltz(outd, t, chains, ag, binder, pos, a.max_distance, a.soft, k)
    if "protenix" in emit:
        emit_protenix(outd, t, chains, ag, binder, pos, a.max_distance, k)
    if "tfold" in emit:
        emit_tfold(outd, t, seqmap, ag, ab_ids, pos, k, wrote)

    print(f"{t} [{a.arm}] 후보 {c['cand']} (patch{k}) · 항원 {ag} {len(pos)}잔기 · binder {binder} "
          f"· force={'false' if a.soft else 'true'} · emit={emit}")
    print(f"  잔기(1-based): {','.join(map(str, pos[:40]))}{' ...' if len(pos) > 40 else ''}")
    print(f"  → {outd}/")
    return 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sites", required=True, help="msa-depth results/sites_<타깃>.json")
    ap.add_argument("--chains", required=True, help="msa-depth targets/<타깃>/chains.json")
    ap.add_argument("--arm", required=True,
                    choices=["ours", "fullmsa", "all", "sizematch", "noconstraint"],
                    help="all = 후보 전부 방출(AbEpiTope 사전 스크리닝용). 파일명은 patch<후보번호>. "
                         "sizematch = 크기 맞춘 대조군. noconstraint = 제약 없는 같은 모델 기준선")
    ap.add_argument("--seed", type=int, default=0, help="sizematch 의 무작위 씨앗(재현용)")
    ap.add_argument("--cand", type=int, default=0,
                    help="--arm ours 일 때 쓸 후보 번호. rank_sites.py 의 PICK 값")
    ap.add_argument("--emit", default="boltz,protenix")
    ap.add_argument("--max-distance", type=float, default=6.0)
    ap.add_argument("--soft", action="store_true", help="Boltz force:false(소프트). 기본=force:true")
    ap.add_argument("--outdir", default="runs_sites_guided")
    a = ap.parse_args()
    emit = [e.strip() for e in a.emit.replace(",", " ").split()]

    d = json.load(open(a.sites))
    t = d["target"]
    chains, ag_ids, ab_ids = load_chains(a.chains)

    if a.arm == "all":
        # 후보 전부 — AbEpiTope 사전 스크리닝(후보마다 값싸게 한 번 모델링해 계면을 채점)용.
        n = 0
        for c in d["candidates"]:
            n += emit_one(a, d, t, c, chains, ag_ids, ab_ids, k=c["cand"])
        print(f"{t} [all] 후보 {n}개 방출 → {os.path.join(a.outdir, t.lower(), 'all')}/")
        return

    if a.arm == "noconstraint":
        # ── 제약 없는 기준선 (같은 모델·같은 MSA 경로·같은 예산) ──────────────
        # 왜 필요한가: msa-depth 의 blind 실행은 protenix 이고 MSA 도 자체 a3m 이다. 그래서
        # boltz 로 유도한 팔과 견주면 '모델이 달라서'와 'MSA 출처가 달라서'가 섞인다.
        # 여기서 같은 모델·같은 경로로 제약만 뺀 실행을 만들면 네 팔이 완전히 맞춰진다.
        outd = os.path.join(a.outdir, t.lower(), a.arm)
        os.makedirs(outd, exist_ok=True)
        emit_l = [e.strip() for e in a.emit.replace(",", " ").split()]
        if "boltz" in emit_l:
            L = ["version: 1", "sequences:"]
            for ch in chains:
                L.append(f'  - protein: {{ id: {ch["id"]}, sequence: "{ch["seq"]}" }}')
            open(os.path.join(outd, f"boltz_{t.lower()}_patch0.yaml"), "w").write("\n".join(L) + "\n")
        if "protenix" in emit_l:
            job = {"name": t.lower(),
                   "sequences": [{"proteinChain": {"sequence": ch["seq"], "count": 1}} for ch in chains]}
            json.dump([job], open(os.path.join(outd, f"protenix_{t.lower()}_patch0.json"), "w"), indent=2)
        print(f"{t} [noconstraint] 제약 없음 · 사슬 {len(chains)}개 · emit={emit_l}")
        print(f"  → {outd}/")
        return

    if a.arm == "sizematch":
        # ── 크기 맞춘 대조군 ──────────────────────────────────────────────────
        # 왜 필요한가: 우리 후보와 원래 MSA 자리는 잔기 수가 다를 수 있다(8k3k_D 는 62 대 43).
        # 넓은 패치는 우연히 정답과 겹칠 확률이 높으므로, "자리가 옳아서 이겼나 그냥 넓어서 이겼나"가
        # 갈리지 않는다. 그래서 원래 MSA 자리를 무작위 잔기로 **우리 후보와 같은 크기까지** 늘린
        # 팔을 만든다. 이것까지 이기면 크기 설명이 배제된다.
        # ⚠️ 채우는 잔기는 항원 서열에서 무작위로 고른다 — 정답 구조를 보지 않는다.
        base = next((x for x in d["candidates"] if x["from_full_msa"]), None)
        tgtc = next((x for x in d["candidates"] if x["cand"] == a.cand), None)
        if base is None or tgtc is None:
            sys.exit("!! sizematch 에는 원래 MSA 후보와 --cand <우리 후보> 가 둘 다 필요하다")
        need = len(tgtc["residues"]) - len(base["residues"])
        cur = {tuple(x) for x in base["residues"]}
        if need <= 0:
            print(f"  원래 MSA 자리({len(cur)})가 이미 우리 후보({len(tgtc['residues'])}) 이상 — 채우지 않는다")
        else:
            rng = random.Random(a.seed)
            lens = {}
            for i, cid in enumerate(ag_ids):
                lens[i] = len(next(ch["seq"] for ch in chains if ch["id"] == cid))
            pool = [(i, p) for i in lens for p in range(lens[i]) if (i, p) not in cur]
            rng.shuffle(pool)
            add = pool[:need]
            cur |= set(add)
            print(f"  크기 맞춤: 원래 MSA 자리 {len(base['residues'])}잔기 + 무작위 {len(add)}잔기 "
                  f"= {len(cur)} (우리 후보 {len(tgtc['residues'])}와 동일)")
        c = dict(base, cand=f"{base['cand']}+무작위", residues=[list(x) for x in sorted(cur)])
    elif a.arm == "fullmsa":
        c = next((x for x in d["candidates"] if x["from_full_msa"]), None)
        if c is None:
            sys.exit(f"!! {t}: 원래 MSA 가 속한 후보가 없다 — 이 타깃은 (b) 대조군을 만들 수 없다.")
    else:
        if not a.cand:
            sys.exit("!! --arm ours 에는 --cand <후보번호> 가 필요하다 (rank_sites.py 의 PICK)")
        c = next((x for x in d["candidates"] if x["cand"] == a.cand), None)
        if c is None:
            sys.exit(f"!! {t}: 후보 {a.cand} 가 없다")
        if c["from_full_msa"]:
            print(f"  ⚠️ 고른 후보 {a.cand} 에 원래 MSA 가 포함돼 있다 — (b)와 같은 자리라 "
                  f"이 타깃은 데모로 부적합하다(대조와 구별이 안 됨).")

    emit_one(a, d, t, c, chains, ag_ids, ab_ids, k=0)


if __name__ == "__main__":
    main()
