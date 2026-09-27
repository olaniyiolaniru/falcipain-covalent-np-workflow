#!/usr/bin/env bash
# Retrospective enrichment benchmark at falcipain-3, run on the same receptor
# and grid as the reported screen and through the same funnel as the
# falcipain-2 benchmark: Glide HTVS over all prepared states, then Glide SP on
# the top 40% by HTVS GlideScore, with the best score per compound retained.
set -euo pipefail
SCHRO="/c/Program Files/Schrodinger2021-2"
cd "$(dirname "$0")"

echo "=== HTVS on all prepared benchmark states ==="
cat > htvs.in <<'EOF2'
GRIDFILE FP3_apo_grid.zip
LIGANDFILE prep.maegz
PRECISION HTVS
POSES_PER_LIG 1
EOF2
"$SCHRO/glide.exe" htvs.in -OVERWRITE -HOST localhost:4 -NJOBS 4 -WAIT

echo "=== select top 40% by HTVS GlideScore ==="
"$SCHRO/run.exe" python3 - <<'PY'
from schrodinger import structure
import csv, math
best = {}
for st in structure.StructureReader("htvs_pv.maegz"):
    g = st.property.get("r_i_glide_gscore")
    if g is None:
        continue
    t = st.title
    if t not in best or float(g) < best[t]:
        best[t] = float(g)
with open("htvs_all_scores.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["id", "htvs_gscore"])
    for k, v in best.items():
        w.writerow([k, v])
keep = set(k for k, _ in sorted(best.items(), key=lambda kv: kv[1])[:max(1, int(round(0.40 * len(best))))])
open("sp_survivors.txt", "w").write("\n".join(sorted(keep)) + "\n")
n = 0
with structure.StructureWriter("sp_input.maegz") as out:
    for st in structure.StructureReader("prep.maegz"):
        if st.title in keep:
            out.append(st); n += 1
print("HTVS %d -> SP %d" % (len(best), len(keep)))
print("SP input: %d" % n)
PY

echo "=== SP on survivors ==="
cat > sp.in <<'EOF3'
GRIDFILE FP3_apo_grid.zip
LIGANDFILE sp_input.maegz
PRECISION SP
POSES_PER_LIG 1
EOF3
"$SCHRO/glide.exe" sp.in -OVERWRITE -HOST localhost:4 -NJOBS 4 -WAIT

echo "=== combine funnel ranking and labels ==="
"$SCHRO/run.exe" python3 - <<'PY'
from schrodinger import structure
import csv
sp = {}
for st in structure.StructureReader("sp_pv.maegz"):
    g = st.property.get("r_i_glide_gscore")
    if g is None:
        continue
    t = st.title
    if t not in sp or float(g) < sp[t]:
        sp[t] = float(g)
htvs = dict((r["id"], float(r["htvs_gscore"])) for r in csv.DictReader(open("htvs_all_scores.csv")))
labels = dict((r["id"], (r["label"], r["class"])) for r in csv.DictReader(open("bench_FP3_labels.csv")))
rows = []
for i, (lab, cls) in labels.items():
    s = sp.get(i, htvs.get(i))
    if s is None:
        continue
    rows.append((i, lab, s, cls, "SP" if i in sp else "HTVS"))
with open("scores_FP3_standard.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["id", "label", "score", "class", "stage"]); w.writerows(rows)
print("scores_FP3_standard.csv:", len(rows))
PY
echo "BENCHMARK_FP3_DONE"
