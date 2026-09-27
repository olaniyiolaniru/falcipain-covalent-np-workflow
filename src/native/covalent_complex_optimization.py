"""Run the covalent-complex Prime REAL_MIN optimisation stage.

This script performs the optimisation stage. Post-reaction XP scoring of the
optimised complexes is a separate stage, implemented in
covalent_postreaction_scoring.py, and the two are released and documented
separately.

Requires the licensed Schrodinger suite and a directory containing the prepared
pose files named in data/covalent_pose_manifest.csv.

  run.exe python3 covalent_complex_optimization.py --input-dir <covalent-pose-directory>
"""
from pathlib import Path
from argparse import ArgumentParser
import csv
import json
import time
import traceback

from schrodinger import structure
from schrodinger.structutils import analyze
from schrodinger.application.prime.packages import Prime
from schrodinger.application.prime.packages import utilities as psp_util

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data" / "covalent_pose_manifest.csv"
OUTPUT = ROOT / "data" / "covalent_complex_optimization"


def graph(st):
    """Atoms and bonds of a structure, for checking the graph is unchanged."""
    atoms = [(a.index, a.element, a.formal_charge) for a in st.atom]
    bonds = sorted((min(b.atom1.index, b.atom2.index),
                    max(b.atom1.index, b.atom2.index), b.order) for b in st.bond)
    return atoms, bonds


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--manifest", type=Path, default=MANIFEST)
    p.add_argument("--input-dir", type=Path, required=True,
                   help="directory holding the prepared pose files named in the manifest")
    p.add_argument("--output", type=Path, default=OUTPUT)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)

    rows = list(csv.DictReader(a.manifest.open()))
    if not rows:
        raise ValueError("the pose manifest is empty")
    print("optimising %d covalent complexes" % len(rows))

    for row in rows:
        key = row["target"] + "_" + row["compound"]
        work = a.output / key
        work.mkdir(exist_ok=True)
        source = a.input_dir / row["source"]
        start = time.time()
        report = dict(target=row["target"], compound=row["compound"],
                      source_file=row["source"], success=False)
        try:
            ct = list(structure.StructureReader(str(source)))[int(row["record"]) - 1]
            sg, cb = int(row["SG_index"]), int(row["Cbeta_index"])
            atom = ct.atom[sg]

            # minimise the reacting residue and everything within 3 A of it
            asl = "chain.name " + atom.chain + " and res.num " + str(atom.resnum)
            select = ("asl=fillres within 3.0 (fillres (atom.i_cdock_lig_attach 1-1000"
                      " OR (" + asl + ")))")
            selected = analyze.evaluate_asl(ct, select[4:])
            assert sg in selected and cb in selected

            params = {"PRIME_TYPE": "REAL_MIN", "MINIM_NITER": "5",
                      "ADD_MISSING_SIDE_CHAINS": "no", "OPLS_VERSION": "OPLS4",
                      "SELECT": select}
            ct.write(str(work / "input.maegz"))
            initial = graph(ct)
            with psp_util.IOToFile(str(work / "Prime_optimization.log"),
                                   str(work / "Prime_optimization.err")):
                output = Prime.run_prime(key, "REAL_MIN", ct, Prime.dict_to_args(params))

            # the optimisation must not change the molecular graph
            result = output[0]
            assert graph(result) == initial
            result.write(str(work / "optimized.maegz"))
            result.write(str(work / "optimized.pdb"))
            report.update(success=True, parameters=params,
                          elapsed_seconds=time.time() - start)
        except Exception as exc:
            report.update(error=str(exc), traceback=traceback.format_exc(),
                          elapsed_seconds=time.time() - start)
        (work / "result.json").write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
