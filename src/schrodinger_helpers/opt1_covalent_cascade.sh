#!/usr/bin/env bash
# OPT1 covalent cascade. covalent_docking is called positionally: rec lig chain:resnum.
# Cathepsins re-prepped apo (free Cys25). FP-2/FP-3 use the apo pH-7 preps (free Cys).
set -euo pipefail
: "${SCHRODINGER:?}"
cd "${WORKDIR:-.}"

# prep apo cathepsins (free Cys25) at pH 7
for pdb in 1ATK 5MQY; do
  "$SCHRODINGER/utilities/prepwizard" -fillsidechains -disulfides -rehtreat \
    -epik_pH 7.0 -epik_pHt 0.5 -propka_pH 7.0 -minimize_adj_h -f OPLS_2005 \
    -WAIT ${pdb}_apoA.pdb ${pdb}_apo_prep.maegz
done

declare -A REC RES
REC[FP2]=FP2_apo_pH70_prep.maegz; RES[FP2]="A:42"
REC[FP3]=FP3_apo_pH70_prep.maegz; RES[FP3]="A:51"
REC[CatK]=1ATK_apo_prep.maegz;    RES[CatK]="A:25"
REC[CatL]=5MQY_apo_prep.maegz;    RES[CatL]="A:25"

for T in FP2 FP3 CatK CatL; do
  "$SCHRODINGER/covalent_docking" \
    "${REC[$T]}" opt1_lig.maegz "${RES[$T]}" \
    -mode enrichment -rxn_type "Michael Addition" \
    -jobname OPT1_${T}_cov -HOST "localhost:2" -WAIT || echo "WARN $T covdock nonzero"
done

"$SCHRODINGER/run" python3 - <<'PY'
from schrodinger import structure
import glob, csv
rows=[]
for T in ("FP2","FP3","CatK","CatL"):
    best=None
    for f in glob.glob(f"OPT1_{T}_cov*-out.maegz"):
        for st in structure.StructureReader(f):
            for k in ("r_i_docking_score","r_i_glide_gscore","r_i_covdock_Prime_Energy"):
                if k in st.property:
                    v=float(st.property[k]); best=v if best is None or v<best else best
                    break
    rows.append([T,best])
with open("OPT1_covalent_cascade.csv","w",newline="") as o:
    w=csv.writer(o); w.writerow(["target","OPT1_CovDock"]); w.writerows(rows)
print("OPT1 covalent cascade:",rows)
PY
echo "OPT1_COVDOCK_DONE"
