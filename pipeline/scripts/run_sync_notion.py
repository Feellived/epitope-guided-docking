#!/usr/bin/env python3
"""summary.csv -> Notion 실험 로그 DB 동기화 (새 행만 push).

토큰은 환경변수 NOTION_API_TOKEN 에서만 읽음 (하드코딩·공용서버 저장 금지).
개인 Mac 등 토큰이 안전한 곳에서 실행:
  export NOTION_API_TOKEN=ntn_xxx
  python run_sync_notion.py --csv results/summary.csv
이미 push한 행은 results/.notion_synced.txt 로 추적해 중복 방지.
"""
import argparse, csv, json, os, sys, urllib.request, urllib.error

DATA_SOURCE = "afec6f93-c8eb-4497-ae05-4f1bb24c7df9"   # 기본 = Experiment 1 실험 로그 DB
DS_EXP2 = "68dc0f8e-3a4c-4bd3-bc52-7e1c3e1a52ec"        # Experiment 2 (Fv/다양-항원) 실험 로그 DB
V = "2025-09-03"
TRUE = ("1", "true", "yes", "o", "y", "✅")


def key(r):
    return "|".join([r["date"], r["model"], r["target"], str(r["seed"]), r["phase"]])


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def infer_format(r):
    """CSV의 format 열 또는 phase(Diverse-Fab/Fv)에서 Fab/Fv 추론. 없으면 None."""
    f = (r.get("format") or "").strip()
    if f in ("Fab", "Fv"):
        return f
    p = (r.get("phase") or "").lower()
    if "fv" in p:
        return "Fv"
    if "fab" in p:
        return "Fab"
    return None


def props(r):
    def t(s):
        return [{"text": {"content": s or ""}}]
    p = {
        "실험": {"title": t(f'{r["model"]} · {r["target"]} · seed{r["seed"]}')},
        "구분": {"select": {"name": r["phase"]}},
        "날짜": {"date": {"start": r["date"]}},
        "모델": {"select": {"name": r["model"]}},
        "타깃": {"select": {"name": r["target"]}},
        "seed": {"number": int(float(r["seed"]))},
        "성공": {"checkbox": str(r.get("success", "")).strip().lower() in TRUE},
        "비고": {"rich_text": t(r.get("notes", ""))},
    }
    fmt = infer_format(r)
    if fmt:                       # 포맷 열이 있는 DB(Exp2)에서만 채워짐; 없으면 Notion이 무시
        p["포맷"] = {"select": {"name": fmt}}
    for ck, nk in [("dockq_abag", "DockQ(Ab-Ag)"), ("global_iptm", "global ipTM"),
                   ("pair_iptm_abag", "pair ipTM(Ab-Ag)")]:
        v = num(r.get(ck))
        if v is not None:
            p[nk] = {"number": v}
    return p


def post(token, p, data_source):
    body = {"parent": {"type": "data_source_id", "data_source_id": data_source}, "properties": p}
    req = urllib.request.Request(
        "https://api.notion.com/v1/pages", data=json.dumps(body).encode("utf-8"), method="POST",
        headers={"Authorization": f"Bearer {token}", "Notion-Version": V, "Content-Type": "application/json"})
    try:
        urllib.request.urlopen(req)
        return True, ""
    except urllib.error.HTTPError as e:
        return False, e.read().decode()[:200]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="results/summary.csv")
    ap.add_argument("--state", default="results/.notion_synced.txt")
    ap.add_argument("--data-source", default=DATA_SOURCE,
                    help="대상 DB data_source. 'exp2' = Experiment 2 DB. 기본 = Experiment 1.")
    a = ap.parse_args()
    ds = DS_EXP2 if a.data_source == "exp2" else a.data_source
    token = os.environ.get("NOTION_API_TOKEN")
    if not token:
        sys.exit("[에러] 환경변수 NOTION_API_TOKEN 없음. 개인 Mac에서 `export NOTION_API_TOKEN=...` 후 실행.")
    synced = set()
    if os.path.exists(a.state):
        synced = set(x for x in open(a.state).read().splitlines() if x)
    with open(a.csv) as f:
        rows = list(csv.DictReader(f))
    new = 0
    with open(a.state, "a") as st:
        for r in rows:
            k = key(r)
            if k in synced:
                continue
            ok, err = post(token, props(r), ds)
            if ok:
                st.write(k + "\n")
                st.flush()
                new += 1
                print("pushed:", k)
            else:
                print("FAIL :", k, "->", err)
    print(f"done. {new} new rows pushed to Notion.")


if __name__ == "__main__":
    main()
