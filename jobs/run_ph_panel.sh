#!/usr/bin/env bash
# pH sensitivity of the reported panel, run under the same standardized
# catalytic-dyad convention as the primary screen.
#
# Falcipains act in the acidic food vacuole, so the panel is repeated with both
# the receptor and the ligands prepared at pH 5.5 and compared with the pH 7.0
# protocol used for the reported screen. The catalytic dyad is held at the
# standardized assignment (neutral cysteine thiol, neutral HIE histidine) at both
# pH values, so the comparison reports the effect of pH on the remainder of the
# site and on ligand ionization rather than on the dyad model.
set -euo pipefail
SCHRO="/c/Program Files/Schrodinger2021-2"
cd "$(dirname "$0")"

echo "=== receptor preparation at pH 5.5 ==="
for T in FP2 FP3; do
  "$SCHRO/utilities/prepwizard" "${T}_chainA_in.mae" "${T}_pH55_prep.mae" \
      -propka_pH 5.5 -fillsidechains -rmsd 0.30 -f 2005 -noepik -NOJOBID
done

echo "=== enforce the standardized catalytic dyad at pH 5.5 ==="
"$SCHRO/run.exe" python3 enforce_dyad.py

echo "=== grids at the published centres ==="
cat > FP2_pH55_grid.in <<'EOF'
FORCEFIELD   OPLS_2005
GRID_CENTER   -53.134676400000004, -6.734088999999999, -18.7730746
GRIDFILE   FP2_pH55_grid.zip
INNERBOX   10, 10, 10
OUTERBOX   30, 30, 30
RECEP_FILE   FP2_pH55_standard.mae
EOF
cat > FP3_pH55_grid.in <<'EOF'
FORCEFIELD   OPLS_2005
GRID_CENTER   7.065, 13.502, -16.81
GRIDFILE   FP3_pH55_grid.zip
INNERBOX   10, 10, 10
OUTERBOX   30, 30, 30
RECEP_FILE   FP3_pH55_standard.mae
EOF
for T in FP2 FP3; do
  "$SCHRO/glide.exe" "${T}_pH55_grid.in" -OVERWRITE -HOST localhost:2 -WAIT
done

echo "=== ligand preparation at pH 5.5 ==="
"$SCHRO/ligprep" -ismi panel.smi -omae panel_pH55.maegz \
    -i 2 -ph 5.5 -pht 2.0 -s 8 -bff 14 -NJOBS 4 -HOST localhost:4 -WAIT

echo "=== dock the panel at both targets, pH 5.5 ==="
for T in FP2 FP3; do
  cat > "${T}_pH55_xp.in" <<EOF
FORCEFIELD   OPLS_2005
GRIDFILE   ${T}_pH55_grid.zip
LIGANDFILE   panel_pH55.maegz
POSTDOCK_XP_DELE   0.5
PRECISION   XP
WRITE_XP_DESC   False
EOF
  "$SCHRO/glide.exe" "${T}_pH55_xp.in" -OVERWRITE -HOST localhost:4 -WAIT
done
echo "PH_PANEL_DONE"
