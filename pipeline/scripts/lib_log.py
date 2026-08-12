"""summary.csv 한 줄 append (다른 스크립트가 import).
CSV = 숫자 기록의 source of truth. Notion은 이걸 읽어 동기화."""
import csv, os
from datetime import date as _date

FIELDS = ["date", "phase", "model", "target", "seed",
          "dockq_abag", "global_iptm", "pair_iptm_abag", "success", "notes"]

# 모델 표기 정규화.
#   기록하는 쪽(eval_*.py)과 읽는 쪽(analyze_*.py)이 서로 다른 이름을 써서
#   집계가 조용히 비는 사고가 있었다. 정본 이름은 오른쪽 값이고,
#   왼쪽의 옛 표기로 적힌 과거 CSV 행도 그대로 읽히도록 둘 다 매핑한다.
MODEL_ALIASES = {
    "protenix":       "Protenix-base",
    "protenix-base":  "Protenix-base",
    "tfold":          "tFold-Ag",
    "tfold-ag":       "tFold-Ag",
    "boltz":          "Boltz-2",
    "boltz-2":        "Boltz-2",
    "chai":           "Chai-1",
    "chai-1":         "Chai-1",
    "chai-1-nomsa":   "Chai-1-noMSA",
    "af2-multimer":   "AF2-Multimer",
    "colabfold":      "AF2-Multimer",
}


def canon_model(name):
    """모델 이름을 정본 표기로 바꾼다. 모르는 이름은 그대로 돌려준다."""
    if not name:
        return name
    return MODEL_ALIASES.get(str(name).strip().lower(), str(name).strip())


def append_row(csv_path, row):
    out = {k: row.get(k, "") for k in FIELDS}
    if not out["date"]:
        out["date"] = _date.today().isoformat()
    out["model"] = canon_model(out["model"])
    new = not os.path.exists(csv_path)
    d = os.path.dirname(csv_path)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(csv_path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        if new:
            w.writeheader()
        w.writerow(out)
    return out
