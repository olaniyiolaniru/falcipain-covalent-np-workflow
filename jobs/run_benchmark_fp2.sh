#!/usr/bin/env bash
# Retrospective enrichment benchmark of the FINAL standardized FP-2 protocol.
#
# Ligand sets, decoy construction and the HTVS -> SP(top 40%) funnel are those of
# the original benchmark; the only change is the receptor/grid, which is now the
# standardized catalytic-dyad receptor used for the reported screen. Metrics are
# ROC-AUC, PR-AUC, EF1%, EF5% and BEDROC(alpha=20) with bootstrap 95% CIs.
set -euo pipefail
SCHRO="/c/Program Files/Schrodinger2021-2"
cd "$(dirname "$0")"

echo "=== HTVS on all prepared benchmark states ==="
cat > htvs.in <<'EOF'
GRIDFILE FP2_standard_grid.zip
LIGANDFILE prep.maegz
PRECISION HTVS
POSES_PER_LIG 1
EOF
"$SCHRO/glide.exe" htvs.in -OVERWRITE -HOST localhost:4 -NJOBS 4 -WAIT

echo "=== select top 40% by HTVS GlideScore ==="
"$SCHRO/run.exe" python3 - <<'PY'
from schrodinger import structure
best = {}
for st in structure.StructureReader("htvs_pv.maegz"):
    g = st.property.get("r_i_glide_gscore")
    if g is None:
        continue
    t = st.title
    g = float(g)
    if t not in best or g < best[t]:
        best[t] = g
items = sorted(best.items(), key=lambda x: x[1])
n = max(1, int(0.40 * len(items)))
keep = set(t for t, _ in items[:n])
open("htvs_all_scores.csv", "w").write("id,htvs_gscore\n" + "\n".join("%s,%s" % (t, g) for t, g in items))
open("sp_survivors.txt", "w").write("\n".join(keep))
print("HTVS %d -> SP %d" % (len(items), len(keep)))
PY

"$SCHRO/run.exe" python3 - <<'PY'
from schrodinger import structure
keep = set(open("sp_survivors.txt").read().split())
w = structure.StructureWriter("sp_input.maegz")
n = 0
for st in structure.StructureReader("prep.maegz"):
    if st.title in keep:
        w.append(st)
        n += 1
w.close()
print("SP input:", n)
PY

echo "=== SP on survivors ==="
cat > sp.in <<'EOF'
GRIDFILE FP2_standard_grid.zip
LIGANDFILE sp_input.maegz
PRECISION SP
POSES_PER_LIG 1
EOF
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
labels = dict((r["id"], (r["label"], r["class"])) for r in csv.DictReader(open("bench_FP2_labels.csv")))
rows = []
for i, (lab, cls) in labels.items():
    s = sp.get(i, htvs.get(i))
    if s is None:
        continue
    rows.append((i, lab, s, cls, "SP" if i in sp else "HTVS"))
with open("scores_FP2_standard.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["id", "label", "score", "class", "stage"])
    w.writerows(rows)
print("scores_FP2_standard.csv:", len(rows))
with open("scores_FP2_standard_cov.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["id", "label", "score"])
    for i, lab, s, cls, stg in rows:
        if lab == "1" and cls != "covalent":
            continue
        w.writerow([i, 1 if lab == "1" else 0, s])
PY

echo "=== enrichment metrics with bootstrap 95% CIs ==="
{
  echo "=== Standardized FP-2 receptor, all actives (covalent + non-covalent) ==="
  "$SCHRO/run.exe" python3 enrichment_metrics.py scores_FP2_standard.csv
  echo "=== Standardized FP-2 receptor, covalent actives only ==="
  "$SCHRO/run.exe" python3 enrichment_metrics.py scores_FP2_standard_cov.csv
} | tee benchmark_FP2_standard_results.txt
echo "BENCHMARK_DONE"
