"""Merge FP-2 and FP-3 XP rankings into a dual-target-weighted shortlist."""
import csv

def load(path, key):
    d = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            d[r["SANC_id"]] = r
    return d

fp2 = load("FP2_XP_ranked.csv", "FP2")
fp3 = load("FP3_XP_ranked.csv", "FP3")

both = sorted(set(fp2) & set(fp3))
rows = []
for sid in both:
    s2 = float(fp2[sid]["GlideScore_FP2_XP"])
    s3 = float(fp3[sid]["GlideScore_FP3_XP"])
    m = fp2[sid]
    rows.append({
        "SANC_id": sid,
        "GlideScore_FP2": round(s2, 3),
        "GlideScore_FP3": round(s3, 3),
        "mean": round((s2 + s3) / 2, 3),
        "worse": round(max(s2, s3), 3),   # least-negative = weakest of the two
        "warheads": m["warheads"], "max_tier": m["max_tier"],
        "pains_flag": m["pains_flag"], "MW": m["MW"], "cLogP": m["cLogP"],
        "HBD": m["HBD"], "HBA": m["HBA"], "TPSA": m["TPSA"], "SMILES": m["SMILES"],
    })

# dual-target rank: by mean score (both must be present => already dual-active)
rows.sort(key=lambda r: r["mean"])
cols = ["dual_rank", "SANC_id", "GlideScore_FP2", "GlideScore_FP3", "mean",
        "worse", "warheads", "max_tier", "pains_flag", "MW", "cLogP",
        "HBD", "HBA", "TPSA", "SMILES"]
with open("dual_target_shortlist.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    for i, r in enumerate(rows, 1):
        r["dual_rank"] = i
        w.writerow(r)

print(f"dual-active (docked to BOTH targets): {len(rows)} compounds")
print(f"  FP-2 only: {len(set(fp2)-set(fp3))} | FP-3 only: {len(set(fp3)-set(fp2))}")
tier2 = [r for r in rows if r["max_tier"] == "2"]
clean = [r for r in rows if not r["pains_flag"]]
print(f"  tier-2 warheads in dual set: {len(tier2)} | PAINS-free in dual set: {len(clean)}")
print("\n  rank SANC        FP2     FP3    mean   tier pains       warhead")
for r in rows[:20]:
    print("  %2d  %-10s %6.2f %6.2f %6.2f   %s   %-11s %s" % (
        r["dual_rank"], r["SANC_id"], r["GlideScore_FP2"], r["GlideScore_FP3"],
        r["mean"], r["max_tier"], r["pains_flag"] or "-", r["warheads"]))
