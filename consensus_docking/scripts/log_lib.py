"""summary.csv 한 줄 append (다른 스크립트가 import).
CSV = 숫자 기록의 source of truth. Notion은 이걸 읽어 동기화."""
import csv, os
from datetime import date as _date

FIELDS = ["date", "phase", "model", "target", "seed",
          "dockq_abag", "global_iptm", "pair_iptm_abag", "success", "notes"]


def append_row(csv_path, row):
    out = {k: row.get(k, "") for k in FIELDS}
    if not out["date"]:
        out["date"] = _date.today().isoformat()
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
