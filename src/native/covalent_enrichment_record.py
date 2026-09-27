"""Record the covalent enrichment outcome for a single compound.

The covalent stage and the reported recognition panel are defined by different
rules, so they select overlapping but not identical compound sets. This script
runs the covalent record for one compound at every target under the identical
enrichment protocol used for the rest of the stage, and reports what that
protocol returned at each, including the targets where it returned no product
pose, so the coverage of the covalent stage can be read exactly.

The record stops at the enrichment stage. It is deliberately not merged into the
covalent stage table, whose complexes were all carried through Prime REAL_MIN
optimisation and post-reaction scoring.

  run.exe python3 covalent_enrichment_record.py --runs <covalent-run-directory> \
      --compound SANC01001 --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import csv
import json
import re

from schrodinger import structure
from schrodinger.structutils import measure

TARGETS = [("FP2", 42, "falcipain-2"), ("FP3", 51, "falcipain-3"),
           ("CatK", 25, "cathepsin K"), ("CatL", 25, "cathepsin L")]
EXIT = re.compile(r"backend died with exit code (\d+)", re.I)


def read_best(path, cys_resnum):
    best, records = None, 0
    for index, st in enumerate(structure.StructureReader(str(path)), 1):
        score = st.property.get("r_i_docking_score")
        lig_atom = st.property.get("i_cdock_lig_atom")
        if score is None or lig_atom is None:
            continue
        records += 1
        if best is not None and score >= best["enrichment_DockingScore"]:
            continue
        sg = [a for a in st.atom
              if a.pdbname.strip() == "SG" and a.resnum == cys_resnum]
        if not sg:
            continue
        cb = st.atom[int(lig_atom)]
        best = dict(enrichment_pose_record=index,
                    enrichment_DockingScore=round(float(score), 4),
                    enrichment_GlideScore=round(
                        float(st.property.get("r_i_sample_glide_gscore",
                                              st.property.get("r_i_glide_gscore",
                                                              float("nan")))), 4),
                    enrichment_minimisation_converged=int(
                        st.property.get("b_ff_Minimization_Converged", 0)),
                    enrichment_S_Cbeta_A=round(
                        float(measure.measure_distance(sg[0], cb)), 4),
                    Cbeta_atom_index=int(cb.index), SG_atom_index=int(sg[0].index))
    if best is not None:
        best["records_returned"] = records
    return best


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--runs", type=Path, required=True)
    p.add_argument("--compound", required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    rows = []
    for target, cys, label in TARGETS:
        out_file = a.runs / ("%s_%s-out.maegz" % (target, a.compound))
        log_file = a.runs / ("%s_%s.log" % (target, a.compound))
        row = dict(target=target, target_name=label, compound=a.compound,
                   reaction="Michael Addition", reactive_residue=cys,
                   covdock_mode="enrichment")
        if out_file.is_file():
            best = read_best(out_file, cys)
            if best is None:
                row.update(status="no record carried both a score and an attachment")
            else:
                row.update(status="completed")
                row.update(best)
        else:
            reason = ""
            if log_file.is_file():
                text = log_file.read_text(errors="replace")
                m = EXIT.search(text)
                if m:
                    reason = " (backend exit code %s)" % m.group(1)
            row.update(status="no product pose returned" + reason)
        rows.append(row)

    fields = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with (a.out / "Covalent_panel_extension.csv").open("w", newline="",
                                                       encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})

    completed = [r["target"] for r in rows if r["status"] == "completed"]
    summary = dict(
        compound=a.compound,
        protocol=("CovDock enrichment mode, Michael addition at the catalytic cysteine, "
                  "identical settings to the covalent panel"),
        stage_reported="enrichment only",
        not_carried_further=("this compound was not carried through Prime REAL_MIN "
                             "optimisation or post-reaction scoring, so it is not part "
                             "of the covalent stage table"),
        targets_attempted=[r["target"] for r in rows],
        targets_returning_a_product_pose=completed,
        targets_returning_nothing=[r["target"] for r in rows
                                   if r["status"] != "completed"])
    (a.out / "Covalent_panel_extension.json").write_text(json.dumps(summary, indent=2),
                                                          encoding="utf-8")
    for r in rows:
        print("%-5s %-14s %s%s" % (r["target"], r["compound"], r["status"],
                                   ("  DS %.4f" % r["enrichment_DockingScore"])
                                   if r.get("enrichment_DockingScore") is not None else ""))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
