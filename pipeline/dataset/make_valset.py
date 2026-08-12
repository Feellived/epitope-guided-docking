#!/usr/bin/env python3
"""검증용 미니 매니페스트: 알려진-에피토프 참조 + 실제 후보 샘플. (Mac에서 분류기 QA)"""
import pandas as pd, re, csv
# 참조(정답 아는 것): 6W41=CR3022 class4(B), 3GBN=CR6261 stem(B), 4FQI=CR9114 stem(B), 3NGB=VRC01 CD4bs(A)
REF = {"6W41":("RBD","class4/B"),"3GBN":("HA","stem/B"),"4FQI":("HA","stem/B"),"3NGB":("Env","CD4bs/A")}
df = pd.read_csv("sabdab2_summary.csv", low_memory=False)
def pdb4(x):
    m = re.match(r"pdb_0000([a-z0-9]{4})$", str(x)); return m.group(1).upper() if m else str(x)
df["p4"] = df["PDB"].map(pdb4)
rows = []
for p4,(antg,exp) in REF.items():
    sub = df[df["p4"]==p4]
    if len(sub)==0: print(f"  참조 {p4} summary에 없음"); continue
    r = sub.iloc[0]
    rows.append(dict(pdb=p4, pdb_full=r["PDB"], antigen=antg, Hchain=str(r["Hchain"]).strip(),
                     Lchain=str(r["Lchain"]).strip(), antigen_chain=str(r["antigen_chain"]).strip(),
                     nb=False, date=r["date"], res=r["resolution"], antigen_name=r["antigen_name"],
                     organism=r["organism"], expect=exp))
# 실제 후보 샘플(항원별 8개)
mani = pd.read_csv("manifest_candidates.csv")
for antg in ["RBD","HA","Env"]:
    s = mani[mani.antigen==antg].head(8).copy(); s["expect"]="(sample)"
    rows += s.to_dict("records")
out = pd.DataFrame(rows)
out.to_csv("valset.csv", index=False)
print(f"valset {len(out)} rows → valset.csv (참조 {len(REF)} + 샘플 24)")
print(out[["pdb","antigen","Hchain","Lchain","antigen_chain","expect"]].to_string(index=False))
