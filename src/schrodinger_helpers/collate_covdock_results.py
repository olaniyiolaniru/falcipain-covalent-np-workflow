#!/usr/bin/env python3
# Final collation: same-mode CovDock ranking + covalent selectivity. Run with $SCHRODINGER/run python3
from schrodinger import structure
import glob, os, csv

def best(f):
    b = None
    for st in structure.StructureReader(f):
        s = st.property.get("r_i_docking_score")
        if s is not None and (b is None or s < b):
            b = float(s)
    return round(b, 2) if b is not None else None

data = {}
for tgt in ["FP2", "FP3", "CatK", "CatL"]:
    for f in glob.glob(f"{tgt}_SANC*-out.maegz"):
        sid = os.path.basename(f).split("_")[1].split("-")[0]
        data.setdefault(sid, {})[tgt] = best(f)

# ---- same-mode CovDock ranking ----
rows = []
for sid, d in data.items():
    f2, f3 = d.get("FP2"), d.get("FP3")
    m = round((f2 + f3) / 2, 2) if f2 is not None and f3 is not None else None
    rows.append((sid, f2, f3, m))
rows.sort(key=lambda r: (r[3] if r[3] is not None else 99))
with open("covdock_samemode_final.csv", "w", newline="") as o:
    w = csv.writer(o); w.writerow(["rank", "SANC_id", "CovDock_FP2", "CovDock_FP3", "mean", "mode"])
    for i, (sid, f2, f3, m) in enumerate(rows, 1):
        w.writerow([i, sid, f2, f3, m, "enrichment"])
print("=== SAME-MODE CovDock (uniform enrichment) ranking ===")
for i, (sid, f2, f3, m) in enumerate(rows, 1):
    print(f"  {i:2d}. {sid}  FP2 {f2}  FP3 {f3}  mean {m}")

# ---- covalent selectivity ----
with open("covalent_selectivity_final.csv", "w", newline="") as o:
    w = csv.writer(o); w.writerow(["SANC_id", "cov_FP2", "cov_FP3", "cov_CatK", "cov_CatL",
                                   "bestFP", "bestCat", "cov_gap"])
    print("\n=== COVALENT selectivity (bestFP - bestCat; negative = falcipain-preferring) ===")
    sel = []
    for sid, d in sorted(data.items()):
        f2, f3, ck, cl = d.get("FP2"), d.get("FP3"), d.get("CatK"), d.get("CatL")
        bF = min([x for x in (f2, f3) if x is not None], default=None)
        bC = min([x for x in (ck, cl) if x is not None], default=None)
        gap = round(bF - bC, 2) if bF is not None and bC is not None else None
        w.writerow([sid, f2, f3, ck, cl, bF, bC, gap])
        sel.append((sid, gap))
        print(f"  {sid}: bestFP {bF}  CatK {ck} CatL {cl}  gap {gap}")
    ng = [g for _, g in sel if g is not None and g < 0]
    print(f"\n  falcipain-preferring (covalent gap<0): {len(ng)}/{len([g for _,g in sel if g is not None])} with complete data")
    miss = [sid for sid, d in data.items() if len(d) < 4]
    if miss:
        print("  compounds with incomplete covalent data:", ", ".join(sorted(miss)))
