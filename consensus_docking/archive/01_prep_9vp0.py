#!/usr/bin/env python3
"""9VP0 (Fab-IL31) 표준 입력 준비.
   PDB 다운로드 -> 사슬 서열 추출 -> 모든 모델 공용 FASTA + Boltz YAML 생성.
   서버 ~/dock/9vp0 에서 실행. 필요: pip install biopython
"""
import os
import urllib.request
from Bio import SeqIO

PDB = "9VP0"
pdb_file = f"{PDB}.pdb"

# 1) 정답 구조 다운로드 (DockQ 정답으로도 쓰임)
if not os.path.exists(pdb_file):
    urllib.request.urlretrieve(
        f"https://files.rcsb.org/download/{PDB}.pdb", pdb_file)
    print(f"downloaded {pdb_file}")

# 2) 사슬별 서열 추출 (ATOM 기준)
recs = {r.annotations["chain"]: str(r.seq)
        for r in SeqIO.parse(pdb_file, "pdb-atom")}
print("=== chains ===")
for c, s in recs.items():
    print(f">{c}  (len {len(s)})")

# 3) 표준 FASTA (사슬 ID 보존 -> DockQ 매칭 자동)
with open(f"{PDB.lower()}.fasta", "w") as f:
    for c, s in recs.items():
        f.write(f">{c}\n{s}\n")

# 4) Boltz용 YAML
with open(f"boltz_{PDB.lower()}.yaml", "w") as f:
    f.write("version: 1\nsequences:\n")
    for c, s in recs.items():
        f.write(f'  - protein: {{ id: {c}, sequence: "{s}" }}\n')

print(f"\nsaved: {PDB.lower()}.fasta, boltz_{PDB.lower()}.yaml")
print("A=IL-31 항원 / H=Fab 중쇄 / L=Fab 경쇄  (항체-항원 = A-H, A-L)")
