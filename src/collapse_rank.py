"""
collapse_rank.py  -- Funnel stage-1 post-processing for the covalent NP screen.

Collapses LigPrep-expanded states (tautomers/stereoisomers/protomers) back to
ONE row per parent SANCDB compound, keeping each compound's single best
(most negative) GlideScore, and merges warhead/property metadata.

Reusable for FP-2 and FP-3, and for SP or XP pose files -- same pipeline.

Run with the Schrodinger Python:
  & "$env:SCHRODINGER\run.exe" python3 collapse_rank.py <pv.maegz> <label> <out.csv>

Parent key = the SANCxxxxx id parsed from each pose title, so it is robust to
whatever suffixes Epik/LigPrep append to state titles.
"""
import sys, re, csv, os
from schrodinger import structure

META_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "sancdb_hits_accessible.csv")

def load_meta(path):
    meta = {}
    if not os.path.exists(path):
        return meta
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            m = re.search(r'(SANC\d+)', row["id"])
            if m:
                meta[m.group(1)] = row
    return meta

def parent_id(title):
    m = re.search(r'(SANC\d+)', title)
    return m.group(1) if m else title

def main():
    pv, label, out = sys.argv[1], sys.argv[2], sys.argv[3]
    meta = load_meta(META_CSV)

    best = {}   # parent -> best score
    nstates = {}
    for i, st in enumerate(structure.StructureReader(pv)):
        ds = st.property.get('r_i_docking_score')
        if ds is None:            # receptor entry (no docking score) -> skip
            continue
        pid = parent_id(st.title)
        nstates[pid] = nstates.get(pid, 0) + 1
        if pid not in best or ds < best[pid]:
            best[pid] = ds

    ranked = sorted(best.items(), key=lambda kv: kv[1])

    cols = ["rank", "SANC_id", f"GlideScore_{label}", "n_states",
            "warheads", "max_tier", "pains_flag", "MW", "cLogP",
            "HBD", "HBA", "TPSA", "SMILES"]
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for rank, (pid, score) in enumerate(ranked, 1):
            m = meta.get(pid, {})
            w.writerow([rank, pid, round(score, 3), nstates[pid],
                        m.get("warheads", ""), m.get("max_tier", ""),
                        m.get("pains_flag", ""), m.get("MW", ""),
                        m.get("cLogP", ""), m.get("HBD", ""), m.get("HBA", ""),
                        m.get("TPSA", ""), m.get("SMILES", "")])

    print(f"[{label}] {len(ranked)} unique parent compounds "
          f"from {sum(nstates.values())} docked states")
    print(f"  best GlideScore {ranked[0][1]:.3f} ({ranked[0][0]}), "
          f"median {ranked[len(ranked)//2][1]:.3f}, "
          f"worst {ranked[-1][1]:.3f}")
    print(f"  written -> {out}")

if __name__ == "__main__":
    main()
