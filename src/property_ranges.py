"""Export property predictions with an explicit record-selection rule.

Several molecules have more than one prepared-state QikProp record. Rather than
quoting one record without saying which, every reported endpoint is given as
the range over all records for that molecule, together with the record count.

  python property_ranges.py --workbook <xlsx> --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import json
import pandas as pd

ENDPOINTS = ["mol_MW", "QPlogPo/w", "QPlogS", "QPPCaco", "QPlogBB",
             "PercentHumanOralAbsorption", "QPlogHERG", "RuleOfFive", "RuleOfThree", "#stars"]
EXCLUDE = {"OPT3"}


def block(frame, label):
    rows = []
    for name, g in frame.groupby("molecule", sort=True):
        key = str(name).replace("_minRM1.pdb", "")
        if key in EXCLUDE:
            continue
        row = {"molecule": key, "set": label, "records": len(g)}
        for e in ENDPOINTS:
            if e not in g.columns:
                continue
            row[e + "_min"] = float(g[e].min())
            row[e + "_max"] = float(g[e].max())
        rows.append(row)
    return rows


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--workbook", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    short = pd.read_excel(a.workbook, sheet_name="QikProp_shortlist")
    anal = pd.read_excel(a.workbook, sheet_name="QikProp_analogues")
    rows = block(short, "shortlist") + block(anal, "analogue")
    frame = pd.DataFrame(rows)
    frame.to_csv(a.out / "Property_ranges.csv", index=False)
    excluded = sorted({str(x).replace("_minRM1.pdb", "") for x in anal.molecule} & EXCLUDE)
    summary = dict(molecules=len(frame), endpoints=ENDPOINTS,
                   selection_rule="range over every prepared-state record for the molecule",
                   excluded=excluded,
                   excluded_reason="never docked and outside the four-member matched design")
    (a.out / "Property_ranges.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    cols = ["molecule", "records", "mol_MW_min", "QPPCaco_min", "QPPCaco_max",
            "PercentHumanOralAbsorption_min", "PercentHumanOralAbsorption_max",
            "RuleOfFive_min", "RuleOfFive_max"]
    print(frame[frame.molecule.isin(["SANC00867", "OPT1"])][cols].to_string(index=False))


if __name__ == "__main__":
    main()
