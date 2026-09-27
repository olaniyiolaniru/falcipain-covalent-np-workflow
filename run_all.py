"""Regenerate every open analysis and figure from the supplied inputs.

One documented sequence, from `data/` into `outputs/`. Nothing in `data/` is
modified or overwritten. Stages that require a licensed Schrodinger installation
are not run here: they live in `src/native/`, they produced the files in `data/`,
and each one names its inputs and outputs in MANIFEST.csv.

  python run_all.py                 # run everything
  python run_all.py --list          # show the sequence without running it
  python run_all.py --only quantum  # run one stage

Exit status is non-zero if any stage fails.
"""
from pathlib import Path
from argparse import ArgumentParser
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
OUT = ROOT / "outputs"
SRC = ROOT / "src"

STAGES = [
    ("provenance", "Regenerate the selection trail from the source SDF",
     [sys.executable, str(SRC / "library_provenance.py"),
      "--sdf", str(DATA / "SANCDB_1012_records.sdf"), "--out", str(OUT)]),
    ("descriptors", "Open descriptors for the retained parents",
     [sys.executable, str(SRC / "library_descriptors.py"),
      "--sdf", str(DATA / "Retained_185.sdf"),
      "--out", str(OUT / "Library_descriptors.csv")]),
    ("recognition", "Receptor-normalized ranking and its sensitivity",
     [sys.executable, str(SRC / "receptor_normalised_ranking.py"),
      "--states", str(DATA / "Primary_state_records.csv"),
      "--failures", str(DATA / "Primary_state_failures.csv"),
      "--descriptors", str(DATA / "Library_descriptors.csv"),
      "--out", str(OUT)]),
    ("quantum", "Electronic descriptors and fragment sums from the atomic charges",
     [sys.executable, str(SRC / "quantum_analysis.py"),
      "--data", str(DATA), "--out", str(OUT)]),
    ("matched", "Matched-analogue contrasts under every selection rule",
     [sys.executable, str(SRC / "matched_analysis.py"),
      "--data", str(DATA), "--out", str(OUT)]),
    ("properties", "Panel property profile and two-platform concordance",
     [sys.executable, str(SRC / "panel_properties.py"),
      "--paired", str(DATA / "Recognition_objectives_paired.csv"),
      "--workbook", str(DATA / "Panel_property_predictions.xlsx"),
      "--names", str(DATA / "Library_compound_names.csv"),
      "--extra", str(DATA / "Panel_property_extra_qikprop.csv"),
      "--out", str(OUT)]),
    ("benchmark_xp", "Retrospective enrichment under the recognition protocol, falcipain-2",
     [sys.executable, str(SRC / "enrichment_metrics.py"),
      "--scores", str(DATA / "Benchmark_scores_FP2_XP.csv"),
      "--labels", str(DATA / "Benchmark_labels_FP2_XP.csv"),
      "--receptor", "FP-2, standardized neutral catalytic dyad, scored by Glide XP "
                    "over every prepared state as the recognition screen does",
      "--tag", "FP2_XP", "--out", str(OUT)]),
    ("benchmark_xp_fp3", "Retrospective enrichment under the recognition protocol, falcipain-3",
     [sys.executable, str(SRC / "enrichment_metrics.py"),
      "--scores", str(DATA / "Benchmark_scores_FP3_XP.csv"),
      "--labels", str(DATA / "Benchmark_labels_FP3.csv"),
      "--receptor", "FP-3, standardized neutral catalytic dyad, scored by Glide XP "
                    "over every prepared state as the recognition screen does",
      "--tag", "FP3_XP", "--out", str(OUT)]),
    ("benchmark", "Retrospective enrichment at falcipain-2, both attrition treatments",
     [sys.executable, str(SRC / "enrichment_metrics.py"),
      "--scores", str(DATA / "Benchmark_scores_FP2_standard.csv"),
      "--labels", str(DATA / "Benchmark_labels_FP2.csv"),
      "--receptor", "FP-2, standardized neutral catalytic dyad, the receptor and "
                    "grid used for the reported screen",
      "--out", str(OUT)]),
    ("benchmark_fp3", "Retrospective enrichment at falcipain-3, both attrition treatments",
     [sys.executable, str(SRC / "enrichment_metrics.py"),
      "--scores", str(DATA / "Benchmark_scores_FP3_standard.csv"),
      "--labels", str(DATA / "Benchmark_labels_FP3.csv"),
      "--receptor", "FP-3, standardized neutral catalytic dyad, the receptor and "
                    "grid used for the reported screen",
      "--tag", "FP3", "--out", str(OUT)]),
    ("precision", "Docking precision and decoy-pool size, separated",
     [sys.executable, str(SRC / "benchmark_precision_comparison.py"),
      "--xp", str(DATA / "Benchmark_scores_FP2_XP.csv"),
      "--funnel", str(DATA / "Benchmark_scores_FP2_standard.csv"),
      "--labels", str(DATA / "Benchmark_labels_FP2_XP.csv"),
      "--out", str(OUT)]),
    ("dyad", "Paired comparison of the two catalytic-dyad receptors",
     [sys.executable, str(SRC / "benchmark_dyad_comparison.py"),
      "--standard", str(DATA / "Benchmark_scores_FP2_standard.csv"),
      "--alternative", str(DATA / "Benchmark_scores_FP2_ion_pair.csv"),
      "--out", str(OUT)]),
    ("figures", "Figures 1 to 4, Figure S1 and the table-of-contents graphic",
     [sys.executable, str(SRC / "build_figures.py"),
      "--data", str(OUT), "--structures", str(DATA / "structures"),
      "--out", str(OUT / "figures")]),
]

# tables the later stages read that are supplied rather than regenerated here
SUPPLIED = [
    "Matched_pose_contacts.csv", "Matched_core_displacement.csv",
    "Matched_halogen_bond_check.csv", "Matched_pose_sensitivity.csv",
    "Matched_pose_clusters.csv", "Matched_pose_family_summary.csv",
    "Quantum_job_manifest.csv", "Quantum_atom_charges.csv",
    "Covalent_stage_records.csv", "Covalent_geometry_audit.csv",
    "Dyad_sensitivity.json",
    "Recognition_objectives_paired.csv", "Recognition_ranking_summary.json",
    "Library_compound_names.csv",
]


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--list", action="store_true")
    p.add_argument("--only", nargs="*")
    a = p.parse_args()

    if a.list:
        for name, note, cmd in STAGES:
            print("%-12s %s" % (name, note))
            print("             %s" % " ".join(cmd[1:]))
        return 0

    OUT.mkdir(exist_ok=True)
    (OUT / "figures").mkdir(exist_ok=True)
    for name in SUPPLIED:
        src = DATA / name
        if src.is_file():
            (OUT / name).write_bytes(src.read_bytes())

    failures = []
    run = 0
    for name, note, cmd in STAGES:
        if a.only and name not in a.only:
            continue
        run += 1
        print("\n=== %s : %s" % (name, note), flush=True)
        start = time.time()
        rc = subprocess.run(cmd).returncode
        print("--- %s %s in %.1fs" % (name, "ok" if rc == 0 else "FAILED",
                                      time.time() - start), flush=True)
        if rc != 0:
            failures.append(name)

    print("\n%d stage(s) run, %d failed%s"
          % (run, len(failures), (": " + ", ".join(failures)) if failures else ""))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
