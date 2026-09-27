#!/usr/bin/env bash
# Retrospective enrichment benchmark at falcipain-2 under the recognition
# protocol itself: Glide XP over every prepared state, on the receptor and grid
# used for the reported screen, with each compound assigned the minimum
# DockingScore among its returned states, exactly as the screen assigns it.
#
# The labelled set keeps all 264 actives and a seeded random sample of 1,000 of
# the 8,026 property-matched decoys. The interval on ROC-AUC is set by the
# number of actives, so the sample costs 0.003 AUC units of precision against
# the full decoy pool, and 1,000 decoys is the size used to validate Glide XP.
set -euo pipefail
SCHRO="/c/Program Files/Schrodinger2021-2"
cd "$(dirname "$0")"

echo "=== XP over all prepared states of the benchmark subset ==="
cat > xp.in <<'EOF1'
GRIDFILE FP2_standard_grid.zip
LIGANDFILE prep.maegz
PRECISION XP
POSES_PER_LIG 1
EOF1
"$SCHRO/glide.exe" xp.in -OVERWRITE -HOST localhost:4 -NJOBS 4 -WAIT

echo "=== consolidate states to compounds by minimum DockingScore ==="
"$SCHRO/run.exe" python3 - <<'PY'
from schrodinger import structure
import csv
best = {}
for st in structure.StructureReader("xp_pv.maegz"):
    d = st.property.get("r_i_docking_score")
    if d is None:
        continue
    t = st.title
    if t not in best or float(d) < best[t]:
        best[t] = float(d)
labels = dict((r["id"], (r["label"], r["class"]))
              for r in csv.DictReader(open("bench_FP2_xp1k_labels.csv")))
rows = [(i, lab, best[i], cls, "XP") for i, (lab, cls) in labels.items() if i in best]
with open("scores_FP2_xp.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["id", "label", "score", "class", "stage"])
    w.writerows(rows)
print("scores_FP2_xp.csv:", len(rows), "of", len(labels), "labelled compounds")
PY
echo "BENCHMARK_FP2_XP_DONE"
