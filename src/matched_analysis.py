"""Matched-analogue effects, state sensitivity and pose-robustness summary.

Consumes only supplied tables:
  Matched_state_records.csv       every returned matched pose record
  Matched_pose_contacts.csv       contacts of every returned matched pose
  Matched_core_displacement.csv   atom-mapped core displacement per rule
  Matched_halogen_bond_check.csv  every Cl...N/O/S contact with its geometry
  Matched_pose_sensitivity.csv    optional wider pose sampling, if supplied

Reports, for each target and architecture, the Cl-minus-OH change under every
selection rule available, so that the direction of each contrast can be judged
against the rule used rather than against one chosen number.

  python matched_analysis.py --data ../exports --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import json
import pandas as pd

ARCH = {"SANC00867": ("glycosylated", "OH"), "Parent_SANC00867": ("glycosylated", "OH"),
        "AnalogueA": ("glycosylated", "Cl"),
        "AnalogueB": ("aglycone_ester", "OH"), "OPT1": ("aglycone_ester", "Cl")}
PAIRS = [("glycosylated", "Parent_SANC00867", "AnalogueA"),
         ("aglycone_ester", "AnalogueB", "OPT1")]


def contrasts(frame, field, label):
    out = []
    for target in sorted(frame.target.unique()):
        block = frame[frame.target.eq(target)]
        picked = {}
        for name in block.title.unique():
            sub = block[block.title.eq(name)]
            picked[name] = sub.loc[sub[field].idxmin()]
        for arch, oh, cl in PAIRS:
            if oh not in picked or cl not in picked:
                continue
            out.append(dict(target=target, architecture=arch, rule=label,
                            OH_member=oh, Cl_member=cl,
                            OH_value=float(picked[oh][field]), Cl_value=float(picked[cl][field]),
                            Cl_minus_OH=float(picked[cl][field] - picked[oh][field]),
                            OH_pose_record=int(picked[oh].pose_record),
                            Cl_pose_record=int(picked[cl].pose_record),
                            OH_states_returned=int(len(block[block.title.eq(oh)])),
                            Cl_states_returned=int(len(block[block.title.eq(cl)]))))
    return out


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--data", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    poses = pd.read_csv(a.data / "Matched_pose_contacts.csv")
    poses = poses[poses.title.isin(ARCH)]
    rows = contrasts(poses, "DockingScore", "minimum_DockingScore")
    rows += contrasts(poses, "GlideScore", "minimum_GlideScore")

    wide = a.data / "Matched_pose_sensitivity.csv"
    if wide.is_file():
        w = pd.read_csv(wide)
        w = w[w.title.isin(ARCH)]
        for run in sorted(w.run.unique()):
            block = w[w.run.eq(run)].rename(columns={"pose_index": "pose_record"})
            rows += contrasts(block, "DockingScore", "minimum_DockingScore_" + run)
            rows += contrasts(block, "GlideScore", "minimum_GlideScore_" + run)

    eff = pd.DataFrame(rows)
    eff.to_csv(a.out / "Matched_effects_by_rule.csv", index=False)

    arch_rows = []
    for rule in eff.rule.unique():
        block = eff[eff.rule.eq(rule)].set_index(["target", "architecture"])
        for target in sorted({t for t, _ in block.index}):
            try:
                g = block.loc[(target, "glycosylated")]
                ag = block.loc[(target, "aglycone_ester")]
            except KeyError:
                continue
            arch_rows.append(dict(target=target, rule=rule,
                                  glycosylated_change=g.Cl_minus_OH,
                                  aglycone_change=ag.Cl_minus_OH,
                                  difference_of_effects=ag.Cl_minus_OH - g.Cl_minus_OH,
                                  architecture_OH=ag.OH_value - g.OH_value,
                                  architecture_Cl=ag.Cl_value - g.Cl_value))
    pd.DataFrame(arch_rows).to_csv(a.out / "Matched_architecture_effects.csv", index=False)

    core = pd.read_csv(a.data / "Matched_core_displacement.csv")
    hal = pd.read_csv(a.data / "Matched_halogen_bond_check.csv")

    signs = {}
    for (target, arch), block in eff.groupby(["target", "architecture"]):
        vals = block.Cl_minus_OH
        signs["%s_%s" % (target, arch)] = dict(
            rules=int(len(vals)),
            min=float(vals.min()), max=float(vals.max()),
            sign_consistent=bool((vals > 0).all() or (vals < 0).all()),
            magnitude_min=float(vals.abs().min()))

    summary = dict(
        rules=sorted(eff.rule.unique()),
        contrast_sign_by_target_architecture=signs,
        core_displacement_rmsd_A={
            "%s_%s_%s" % (r.target, r.architecture, r.selection_rule): r.core_displacement_rmsd_A
            for r in core.itertuples()},
        halogen_bond_contacts_screened=int(len(hal)),
        halogen_bonds_found=int(hal.halogen_bond.sum()),
        halogen_max_angle_deg=float(hal.C_Cl_acceptor_angle_deg.max()),
        halogen_min_distance_A=float(hal.Cl_acceptor_distance_A.min()),
        interpretation=("Score differences accompany a change of binding mode: the "
                        "mapped core relocates between the two members."))
    (a.out / "Matched_analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
