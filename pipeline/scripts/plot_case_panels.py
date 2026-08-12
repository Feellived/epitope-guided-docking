"""[구조 3패널] 정답 구조로 카메라를 한 번 잡고 세 예측을 같은 각도로 그린다.

⚠️ 반드시 **한 세션 안에서** 세 장을 렌더한다. PyMOL 의 orient 는 주축 부호가 기기·빌드마다
   달라, 다른 기기에서 따로 그린 패널은 각도가 맞지 않는다(실제로 겪었다).

환경변수:
  NAT   정답 구조 · AGN/ABN  정답의 항원/항체 사슬
  P1 P2 P3   예측 세 개 (AlphaFold 3 · 제약 없음 · 선택 자리 순)
  AG1 AB1 ... 각 예측의 항원/항체 사슬 (비우면 AGN/ABN 과 같게)
  OUT   저장 경로 (3600x1000 로 이어 붙인다)

사용:
  NAT=native.cif AGN=G ABN=H P1=af3.cif AG1=B AB1=A \
  P2=noconstraint.cif P3=ours.cif AG2=G AB2=H AG3=G AB3=H \
  OUT=fig6.png pymol -cq plot_case_panels.py
"""
import os
from pymol import cmd
E=os.environ
sel=lambda cs: " or ".join(f"chain {c}" for c in cs.split())
cmd.set("auto_zoom",0); cmd.set("ray_opaque_background",1); cmd.bg_color("white")
cmd.set("surface_quality",1); cmd.set("transparency",0.55)
cmd.set_color("c_epi",[0.86,0.12,0.45]); cmd.set_color("c_ag",[0.80,0.80,0.80])
cmd.set_color("c_true",[0.62,0.78,0.92]); cmd.set_color("c_pred",[0.95,0.55,0.10])
AGN,ABN=E["AGN"],E["ABN"]
cmd.load(E["NAT"],"nat")
cmd.create("nat_ag",f"nat and polymer and ({sel(AGN)})")
cmd.create("nat_ab",f"nat and polymer and ({sel(ABN)})")
cmd.delete("nat"); cmd.select("epi",f"nat_ag within 5.0 of nat_ab")
print("정답 결합자리", cmd.count_atoms("epi and name CA"), "잔기")
cmd.hide("everything"); cmd.show("surface","nat_ag")
cmd.color("c_ag","nat_ag"); cmd.color("c_epi","epi")
cmd.show("cartoon","nat_ab"); cmd.color("c_true","nat_ab")
cmd.set("cartoon_transparency",0.55,"nat_ab")
cmd.orient("nat_ag or nat_ab"); cmd.zoom("nat_ag or nat_ab",3)
V=cmd.get_view()                      # ★ 이 각도를 세 패널에 그대로 쓴다
outs=[]
for i in (1,2,3):
    f=E.get(f"P{i}")
    if not f: print(f"  ! P{i} 없음 — 건너뜀"); continue
    ag=E.get(f"AG{i}") or AGN; ab=E.get(f"AB{i}") or ABN
    if not os.path.exists(f):
        print(f"  ! P{i}: 파일이 없다 — 건너뜀  {f}"); continue
    cmd.delete("prd*"); cmd.load(f,"prd")
    cmd.create("prd_ag",f"prd and polymer and ({sel(ag)})")
    cmd.create("prd_ab",f"prd and polymer and ({sel(ab)})")
    if not cmd.count_atoms("prd_ag") or not cmd.count_atoms("prd_ab"):
        print(f"  ! P{i}: 사슬을 못 찾았다 (있는 것 {sorted(set(cmd.get_chains('prd')))})"); continue
    try: rms=cmd.cealign("nat_ag","prd_ag")["RMSD"]
    except Exception: rms=cmd.align("prd_ag","nat_ag")[0]
    cmd.matrix_copy("prd_ag","prd_ab")
    cmd.hide("everything","prd*"); cmd.show("cartoon","prd_ab")
    cmd.color("c_pred","prd_ab"); cmd.set("cartoon_transparency",0.25,"prd_ab")
    cmd.set_view(V)
    p=f"{E['OUT']}.p{i}.png"; cmd.png(p,width=1200,height=1000,dpi=200,ray=1)
    outs.append(p); print(f"  패널 {i}  정렬 RMSD {rms:.2f} Å → {p}")
try:
    from PIL import Image
    ims=[Image.open(p).convert("RGB") for p in outs]
    W=sum(i.width for i in ims); H=max(i.height for i in ims)
    m=Image.new("RGB",(W,H),"white"); x=0
    for i in ims: m.paste(i,(x,0)); x+=i.width
    m.save(E["OUT"]); print("→", E["OUT"], m.size)
except Exception as e:
    print("이어 붙이기 실패 — 패널은 저장됨:", e)
