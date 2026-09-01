#!/usr/bin/env bash
# Runs on apo receptors (crystallographic covalent inhibitor stripped, free catalytic Cys):
#  (A) pH 5.5 vs 7.0 sensitivity for FP-2/FP-3 (16-lead XP redock)
#  (B) OPT1 full covalent cascade: CovDock into FP-2, FP-3, cathepsin K, cathepsin L
set -euo pipefail
: "${SCHRODINGER:?}"
cd "${WORKDIR:-.}"

declare -A APO CEN CYS
APO[FP2]=3BPF_apoA.pdb; CEN[FP2]="-57.249588584905666, -0.92186258490566, -15.734413641509432"; CYS[FP2]=42
APO[FP3]=3BPM_apoA.pdb; CEN[FP3]="4.8004247, 16.389353699999997, -23.692555199999997"; CYS[FP3]=51

# ---------- (A) apo prep + grid + XP redock at both pH ----------
for T in FP2 FP3; do
  for PH in 55 70; do
    ph=$( [ "$PH" = 55 ] && echo 5.5 || echo 7.0 )
    tag=${T}_apo_pH${PH}
    "$SCHRODINGER/utilities/prepwizard" -fillsidechains -disulfides -rehtreat \
      -epik_pH $ph -epik_pHt 0.5 -propka_pH $ph -minimize_adj_h -f OPLS_2005 \
      -WAIT "${APO[$T]}" ${tag}_prep.maegz
    cat > ${tag}_grid.in <<EOF
FORCEFIELD   OPLS_2005
GRID_CENTER   ${CEN[$T]}
GRIDFILE   ${tag}_grid.zip
INNERBOX   10, 10, 10
OUTERBOX   30, 30, 30
RECEP_FILE   ${tag}_prep.maegz
EOF
    "$SCHRODINGER/glide" ${tag}_grid.in -OVERWRITE -HOST "localhost:2" -WAIT
    cat > ${tag}_xp.in <<EOF
GRIDFILE ${tag}_grid.zip
LIGANDFILE leads16.maegz
PRECISION XP
POSES_PER_LIG 1
EOF
    "$SCHRODINGER/glide" ${tag}_xp.in -OVERWRITE -HOST "localhost:2" -NJOBS 2 -WAIT
  done
done

"$SCHRODINGER/run" python3 - <<'PY'
from schrodinger import structure
import re, csv, collections, math
def sid(t):
    m=re.match(r"(SANC\d+)",t); return m.group(1) if m else t
data=collections.defaultdict(dict)
for T in ("FP2","FP3"):
    for PH,lab in (("55","pH5.5"),("70","pH7.0")):
        best={}
        try:
            for st in structure.StructureReader(f"{T}_apo_pH{PH}_xp_pv.maegz"):
                g=st.property.get("r_i_glide_gscore")
                if g is None: continue
                s=sid(st.title); g=float(g)
                if s not in best or g<best[s]: best[s]=g
        except Exception as e: print("miss",T,PH,e); continue
        for s,g in best.items(): data[(T,s)][lab]=g
with open("pH_sensitivity_apo.csv","w",newline="") as o:
    w=csv.writer(o); w.writerow(["target","SANC_id","XP_pH7.0","XP_pH5.5","delta(5.5-7.0)"])
    for (T,s),d in sorted(data.items()):
        a,b=d.get("pH7.0"),d.get("pH5.5")
        w.writerow([T,s,a,b,(b-a) if a is not None and b is not None else None])
def spearman(xs,ys):
    n=len(xs)
    if n<3: return None
    rx={v:i for i,v in enumerate(sorted(range(n),key=lambda k:xs[k]))}
    ry={v:i for i,v in enumerate(sorted(range(n),key=lambda k:ys[k]))}
    return 1-6*sum((rx[i]-ry[i])**2 for i in range(n))/(n*(n*n-1))
for T in ("FP2","FP3"):
    ps=[(d["pH7.0"],d["pH5.5"]) for (t,s),d in data.items() if t==T and "pH7.0" in d and "pH5.5" in d]
    if len(ps)>=3:
        rho=spearman([p[0] for p in ps],[p[1] for p in ps])
        md=sum(abs(p[1]-p[0]) for p in ps)/len(ps)
        print(f"[APO] {T}: n={len(ps)} Spearman={rho:.3f} mean|delta|={md:.2f}")
print("wrote pH_sensitivity_apo.csv")
PY

# ---------- (B) OPT1 covalent cascade on apo receptors ----------
declare -A REC RCYS
REC[FP2]=FP2_apo_pH70_prep.maegz; RCYS[FP2]=42
REC[FP3]=FP3_apo_pH70_prep.maegz; RCYS[FP3]=51
REC[CatK]=1ATK_prep.maegz;        RCYS[CatK]=25
REC[CatL]=5MQY_prep.maegz;        RCYS[CatL]=25
for T in FP2 FP3 CatK CatL; do
  c=${RCYS[$T]}
  "$SCHRODINGER/covalent_docking" -jobname OPT1_${T}_cov \
    -receptor "${REC[$T]}" \
    -grid_center_asl "res.num $c and (atom.ptype \" SG \")" \
    -reactive_residue_asl "res.num $c and (atom.ptype \" SG \")" \
    -reaction_type michael_addition -mode VSGB_ENRICH \
    opt1_lig.maegz -HOST "localhost:2" -WAIT || echo "WARN $T covdock nonzero"
done
"$SCHRODINGER/run" python3 - <<'PY'
from schrodinger import structure
import glob, csv
rows=[]
for T in ("FP2","FP3","CatK","CatL"):
    best=None
    for f in glob.glob(f"OPT1_{T}_cov*-out.maegz"):
        for st in structure.StructureReader(f):
            for k in ("r_i_docking_score","r_i_glide_gscore"):
                if k in st.property:
                    v=float(st.property[k]); best=v if best is None or v<best else best
    rows.append([T,best])
with open("OPT1_covalent_cascade.csv","w",newline="") as o:
    w=csv.writer(o); w.writerow(["target","OPT1_CovDock"]); w.writerows(rows)
print("OPT1 covalent cascade:",rows)
PY
echo "DONE"
