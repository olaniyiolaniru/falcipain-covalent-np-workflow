"""Export every matched-analogue pose from the wider pose-sampling runs.

The manuscript comparison retains one pose per prepared state. To test whether
the direction of each matched contrast survives broader conformational
sampling, the same twelve prepared states were re-docked at both grids while
retaining up to ten poses per state (XPn10). Enhanced Glide sampling is not
available in XP mode, so the independent sampling axis is a wider SP run with
four-fold enhanced sampling (SPes4). SP and XP scores are not comparable to one
another; what the SP run tests is whether the SIGN of each matched contrast
survives a different sampling protocol and scoring function.
Every returned pose of every run is exported here.

  run.exe python3 extract_matched_poses.py --runs ../jobs/matched_pose_sensitivity \
      --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import csv
import hashlib
import json

from schrodinger import structure

RUNS = [("FP2", "XPn10", "FP2_matched_n10_pv.maegz"),
        ("FP3", "XPn10", "FP3_matched_n10_pv.maegz"),
        ("FP2", "SPes4", "FP2_matched_sp4_pv.maegz"),
        ("FP3", "SPes4", "FP3_matched_sp4_pv.maegz")]

TITLES = {"Parent_SANC00867", "AnalogueA", "AnalogueB", "OPT1"}


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--runs", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    rows, missing = [], []
    for target, run, fname in RUNS:
        path = a.runs / fname
        if not path.is_file():
            missing.append(fname)
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        rank = {}
        for n, st in enumerate(structure.StructureReader(str(path)), 1):
            if n == 1:
                continue
            title = st.title
            rank[title] = rank.get(title, 0) + 1
            rows.append(dict(
                target=target, run=run, source=fname, source_sha256=digest,
                pose_index=n, title=title,
                variant=st.property.get("s_lp_Variant", title),
                pose_rank_within_compound=rank[title],
                formal_charge=st.formal_charge,
                DockingScore=st.property.get("r_i_docking_score"),
                GlideScore=st.property.get("r_i_glide_gscore"),
                glide_emodel=st.property.get("r_i_glide_emodel"),
                glide_energy=st.property.get("r_i_glide_energy"),
                in_matched_set=title in TITLES))

    if not rows:
        raise SystemExit("no pose viewers found in %s" % a.runs)
    with (a.out / "Matched_pose_sensitivity.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    summary = dict(
        runs=sorted({(r["target"], r["run"]) for r in rows}),
        poses=len(rows),
        missing_pose_viewers=missing,
        poses_per_run={"%s_%s" % (t, r): sum(1 for x in rows if x["target"] == t and x["run"] == r)
                       for t, r, _ in RUNS if (t, r) in {(x["target"], x["run"]) for x in rows}},
        compounds=sorted({r["title"] for r in rows}))
    summary["runs"] = ["%s_%s" % (t, r) for t, r in summary["runs"]]
    (a.out / "Matched_pose_sensitivity.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
