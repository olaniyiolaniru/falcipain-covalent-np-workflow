"""Build the covalent pose manifest for one or more CovDock enrichment outputs.

The manifest names, for each target and compound, the pose record that the
optimisation and post-reaction scoring stages act on. The selection rule is the
one used for the rest of the panel: the record with the minimum DockingScore
among every record the enrichment run returned.

For each selected record the script locates the catalytic cysteine S-gamma in
the complex by residue number and atom name, takes the attacked ligand carbon
from the CovDock `i_cdock_lig_atom` property, and measures the S-C distance, so
that every field is read from the structure rather than assumed.

  run.exe python3 covalent_pose_manifest_builder.py --runs <covalent-run-directory> \
      --compound SANC01001 --out <covalent-run-directory>/covalent_pose_manifest.csv
"""
from pathlib import Path
from argparse import ArgumentParser
import csv

from schrodinger import structure
from schrodinger.structutils import measure

# target, catalytic cysteine residue number, receptor file used by the panel
TARGETS = [("FP2", 42, "FP2_rec.mae"),
           ("FP3", 51, "FP3_rec.maegz"),
           ("CatK", 25, "CatK_rec.maegz"),
           ("CatL", 25, "CatL_rec.maegz")]

COMMANDLINE = ('covalent_docking.exe <receptor> <ligand> A:%d '
               '-rxn_type "Michael Addition" -mode enrichment')


def best_record(path, cys_resnum):
    best = None
    for index, st in enumerate(structure.StructureReader(str(path)), 1):
        score = st.property.get("r_i_docking_score")
        lig_atom = st.property.get("i_cdock_lig_atom")
        if score is None or lig_atom is None:
            continue
        if best is not None and score >= best["score"]:
            continue
        sg = [a for a in st.atom
              if a.pdbname.strip() == "SG" and a.resnum == cys_resnum]
        if not sg:
            continue
        cb = st.atom[int(lig_atom)]
        best = dict(record=index, score=float(score),
                    sample_gscore=float(st.property.get("r_i_sample_glide_gscore",
                                                        st.property.get("r_i_glide_gscore"))),
                    minimization_converged=int(st.property.get("b_ff_Minimization_Converged", 0)),
                    rms_derivative=float(st.property.get("r_ff_RMS_Derivative", float("nan"))),
                    SG_index=int(sg[0].index), Cbeta_index=int(cb.index),
                    SG_Cbeta_A=float(measure.measure_distance(sg[0], cb)))
    return best


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--runs", type=Path, required=True)
    p.add_argument("--compound", required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()

    rows, missing = [], []
    for target, cys, _rec in TARGETS:
        source = "%s_%s-out.maegz" % (target, a.compound)
        path = a.runs / source
        if not path.is_file():
            missing.append(source)
            continue
        best = best_record(path, cys)
        if best is None:
            missing.append(source + " (no record carried both a score and an attachment)")
            continue
        rows.append(dict(target=target, compound=a.compound, source=source,
                         record=best["record"], score=best["score"],
                         sample_gscore=best["sample_gscore"],
                         minimization_converged=best["minimization_converged"],
                         rms_derivative=best["rms_derivative"],
                         SG_index=best["SG_index"], SG_residue=cys,
                         Cbeta_index=best["Cbeta_index"],
                         SG_Cbeta_A=best["SG_Cbeta_A"],
                         commandline=COMMANDLINE % cys))

    a.out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["target", "compound", "source", "record", "score", "sample_gscore",
              "minimization_converged", "rms_derivative", "SG_index", "SG_residue",
              "Cbeta_index", "SG_Cbeta_A", "commandline"]
    with a.out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    for r in rows:
        print("%-5s %s  record %2d  DS %8.4f  S-Cbeta %.3f A"
              % (r["target"], r["compound"], r["record"], r["score"], r["SG_Cbeta_A"]))
    print("wrote %d manifest rows to %s" % (len(rows), a.out))
    if missing:
        print("not present:", ", ".join(missing))


if __name__ == "__main__":
    main()
