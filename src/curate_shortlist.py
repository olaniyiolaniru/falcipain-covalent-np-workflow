"""
Curate the CovDock shortlist from the dual-active compounds.

Goal = best OUTPUT, not most compounds:
  - dual-target-weighted (mean of FP-2/FP-3 XP GlideScore)
  - chemotype diversity (cap 2 compounds per generic Murcko scaffold, so the
    list isn't 6 near-identical caffeoyl glycosides)
  - force-include the tier-2 high-reactivity warhead subgroup
  - PAINS + MW kept visible (annotated, never silently dropped)
"""
import csv
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

CORE_N = 12   # diverse score-ranked core; tier-2 appended on top

rows = list(csv.DictReader(open("dual_target_shortlist.csv")))
for r in rows:
    m = Chem.MolFromSmiles(r["SMILES"])
    r["scaffold"] = (Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(
        MurckoScaffold.GetScaffoldForMol(m))) if m else r["SANC_id"])
    r["mean"] = float(r["mean"])

rows.sort(key=lambda r: r["mean"])   # most negative first

picked, scaf_count = [], {}
for r in rows:
    s = r["scaffold"]
    if scaf_count.get(s, 0) >= 2:      # diversity cap
        continue
    picked.append(r); scaf_count[s] = scaf_count.get(s, 0) + 1
    if len(picked) >= CORE_N:
        break

# force-include tier-2 high-reactivity warheads that dock both targets
# (a positive GlideScore = non-binding; such a ligand yields no Glide poses and
#  crashes the covalent_docking subjob, so it must be excluded from the set)
have = {r["SANC_id"] for r in picked}
tier2 = [r for r in rows if r["max_tier"] == "2" and r["SANC_id"] not in have
         and float(r["GlideScore_FP2"]) < 0 and float(r["GlideScore_FP3"]) < 0]
picked += tier2

picked.sort(key=lambda r: r["mean"])
cols = ["SANC_id", "GlideScore_FP2", "GlideScore_FP3", "mean", "max_tier",
        "pains_flag", "MW", "cLogP", "warheads", "SMILES"]
with open("covdock_shortlist.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
    w.writeheader()
    for r in picked:
        w.writerow(r)

print(f"CovDock shortlist: {len(picked)} compounds "
      f"({CORE_N} diverse score-ranked core + {len(tier2)} added tier-2)")
print(f"  PAINS-flagged: {sum(1 for r in picked if r['pains_flag'])} "
      f"(kept, flagged) | MW>550: {sum(1 for r in picked if float(r['MW'])>550)}")
print("\n  SANC        FP2    FP3    mean   tier pains        MW    warhead")
for r in picked:
    print("  %-10s %5.1f %5.1f %6.2f   %s   %-11s %5.0f  %s" % (
        r["SANC_id"], float(r["GlideScore_FP2"]), float(r["GlideScore_FP3"]),
        r["mean"], r["max_tier"], r["pains_flag"] or "-", float(r["MW"]),
        r["warheads"]))
