#!/usr/bin/env bash
# Full-scale FP-2 retrospective benchmark: 264 actives (<=10 uM) + 8026 property-matched decoys.
# LigPrep (NJOBS 8) -> Glide HTVS (all) -> SP (top 40%) -> combined ranking -> enrichment + CIs.
set -euo pipefail
: "${SCHRODINGER:?}"
cd "${WORKDIR:-.}"

# --- LigPrep (Ionizer, pH 7, 1 stereoisomer), full concurrency ---
"$SCHRODINGER/ligprep" -ismi bench_FP2.smi -omae prep.maegz \
  -i 2 -ph 7.0 -pht 2.0 -s 1 -bff 16 -NJOBS 8 -HOST "localhost:8" -WAIT
echo "LIGPREP_DONE $(("$SCHRODINGER/run" python3 -c "from schrodinger import structure;print(sum(1 for _ in structure.StructureReader('prep.maegz')))"))"

# --- HTVS on all ---
cat > htvs.in <<'EOF'
GRIDFILE FP2_grid.zip
LIGANDFILE prep.maegz
PRECISION HTVS
POSES_PER_LIG 1
EOF
"$SCHRODINGER/glide" htvs.in -OVERWRITE -HOST "localhost:8" -NJOBS 8 -WAIT

# --- top 40% by HTVS gscore -> SP ---
"$SCHRODINGER/run" python3 - <<'PY'
from schrodinger import structure
best={}
for st in structure.StructureReader("htvs_pv.maegz"):
    g=st.property.get("r_i_glide_gscore")
    if g is None: continue
    t=st.title; g=float(g)
    if t not in best or g<best[t]: best[t]=g
items=sorted(best.items(),key=lambda x:x[1])
n=max(1,int(0.40*len(items)))
keep=set(t for t,_ in items[:n])
open("htvs_all_scores.csv","w").write("id,htvs_gscore\n"+"\n".join(f"{t},{g}" for t,g in items))
open("sp_survivors.txt","w").write("\n".join(keep))
print(f"HTVS {len(items)} -> SP {len(keep)}")
PY
"$SCHRODINGER/run" python3 - <<'PY'
from schrodinger import structure
keep=set(open("sp_survivors.txt").read().split())
w=structure.StructureWriter("sp_input.maegz"); n=0
for st in structure.StructureReader("prep.maegz"):
    if st.title in keep: w.append(st); n+=1
w.close(); print("SP input:",n)
PY

# --- SP on survivors ---
cat > sp.in <<'EOF'
GRIDFILE FP2_grid.zip
LIGANDFILE sp_input.maegz
PRECISION SP
POSES_PER_LIG 1
EOF
"$SCHRODINGER/glide" sp.in -OVERWRITE -HOST "localhost:8" -NJOBS 8 -WAIT

# --- combined funnel ranking + labels ---
"$SCHRODINGER/run" python3 - <<'PY'
from schrodinger import structure
import csv
sp={}
for st in structure.StructureReader("sp_pv.maegz"):
    g=st.property.get("r_i_glide_gscore")
    if g is None: continue
    t=st.title
    if t not in sp or float(g)<sp[t]: sp[t]=float(g)
htvs=dict((r["id"],float(r["htvs_gscore"])) for r in csv.DictReader(open("htvs_all_scores.csv")))
labels=dict((r["id"],(r["label"],r["class"])) for r in csv.DictReader(open("bench_FP2_labels.csv")))
rows=[]
for i,(lab,cls) in labels.items():
    s=sp.get(i, htvs.get(i))
    if s is None: continue
    rows.append((i,lab,s,cls,"SP" if i in sp else "HTVS"))
with open("scores_FP2_full.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["id","label","score","class","stage"]); w.writerows(rows)
print("scores_FP2_full.csv:",len(rows))
# covalent-only subset
with open("scores_FP2_full_cov.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["id","label","score"])
    for i,lab,s,cls,stg in rows:
        if lab=="1" and cls!="covalent": continue
        w.writerow([i,1 if lab=="1" else 0,s])
PY

# --- enrichment metrics with bootstrap CIs ---
echo "=== FULL benchmark (all actives) ==="
"$SCHRODINGER/run" python3 enrichment_metrics.py scores_FP2_full.csv | tee benchmark_FP2_full_results.txt
echo "=== FULL benchmark (covalent actives only) ===" | tee -a benchmark_FP2_full_results.txt
"$SCHRODINGER/run" python3 enrichment_metrics.py scores_FP2_full_cov.csv | tee -a benchmark_FP2_full_results.txt
echo "FULL_PIPELINE_DONE"
