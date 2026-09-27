"""Assemble the property profile of the reported panel.

The panel is defined by the stated recognition rule - parents in the top decile
at both receptors - so this script reads that panel from the receptor-normalized
ranking rather than from a fixed list, joins the QikProp record for every member,
and reports each endpoint as a range over the prepared-state records of that
molecule. Members whose QikProp record was generated separately, because they
entered the panel under the receptor-normalized rule, are merged from their own
QikProp output file and flagged with their source.

  python panel_properties.py --paired ../exports/Recognition_objectives_paired.csv \
      --workbook ../data/Prior_Supplementary_Dataset.xlsx \
      --extra ../jobs/qikprop_gap/gap_prep.CSV --names ../exports/Library_compound_names.csv \
      --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import json

import pandas as pd

ENDPOINTS = ["mol_MW", "QPlogPo/w", "QPlogS", "QPPCaco", "QPPMDCK", "QPlogBB",
             "PercentHumanOralAbsorption", "QPlogHERG", "QPlogKhsa",
             "RuleOfFive", "RuleOfThree", "#stars"]
PERCENTILE_CUT = 90.0


def clean(series):
    return series.astype(str).str.replace("_minRM1.pdb", "", regex=False)


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--paired", type=Path, required=True)
    p.add_argument("--workbook", type=Path, required=True)
    p.add_argument("--extra", type=Path, nargs="*", default=[])
    p.add_argument("--names", type=Path)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    paired = pd.read_csv(a.paired)
    panel = paired[paired.worse_percentile >= PERCENTILE_CUT].copy()
    panel = panel.sort_values("mean_percentile", ascending=False)
    members = list(panel.compound)

    qp = pd.read_excel(a.workbook, sheet_name="QikProp_shortlist")
    qp["compound"] = clean(qp.molecule)
    qp["qikprop_source"] = "shortlist QikProp run"
    frames = [qp]
    for extra in a.extra:
        if not Path(extra).is_file():
            continue
        e = pd.read_csv(extra)
        e["compound"] = clean(e.molecule)
        e["qikprop_source"] = Path(extra).name
        frames.append(e)
    qp = pd.concat(frames, ignore_index=True, sort=False)

    names = {}
    if a.names and a.names.is_file():
        nf = pd.read_csv(a.names)
        names = dict(zip(nf.compound, nf.name))

    rows, missing = [], []
    for c in members:
        block = qp[qp.compound.eq(c)]
        if block.empty:
            missing.append(c)
            continue
        pr = panel[panel.compound.eq(c)].iloc[0]
        row = dict(compound=c, name=names.get(c, ""),
                   mean_percentile=round(float(pr.mean_percentile), 2),
                   worse_percentile=round(float(pr.worse_percentile), 2),
                   heavy_atoms=int(pr.heavy_atoms),
                   qikprop_records=int(len(block)),
                   qikprop_source=";".join(sorted(set(block.qikprop_source))))
        for e in ENDPOINTS:
            if e not in block.columns:
                continue
            lo, hi = float(block[e].min()), float(block[e].max())
            key = e.replace("/", "_").replace("#", "n_")
            if abs(hi - lo) < 5e-4:
                row[key] = round(lo, 3)
            else:
                row[key + "_min"] = round(lo, 3)
                row[key + "_max"] = round(hi, 3)
        rows.append(row)

    frame = pd.DataFrame(rows)
    frame.to_csv(a.out / "Panel_property_profile.csv", index=False)

    summary = dict(
        panel_rule=("parents in the top decile at both receptors "
                    "(worse within-target percentile >= %.0f)" % PERCENTILE_CUT),
        panel_size=len(members),
        panel_members=members,
        members_with_qikprop=int(len(frame)),
        members_without_qikprop=missing,
        endpoints=ENDPOINTS,
        reporting_rule="range over every prepared-state record for the molecule")
    (a.out / "Panel_property_summary.json").write_text(json.dumps(summary, indent=2),
                                                        encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(frame.to_string(index=False))


if __name__ == "__main__":
    main()
