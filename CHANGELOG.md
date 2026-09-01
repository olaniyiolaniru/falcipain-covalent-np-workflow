# Release notes

Reproducibility package accompanying the paper *Reactivity Over Recognition: Covalent Screening of
Natural Products against Falcipain-2 and Falcipain-3*.

## Contents

- **Warhead-filtering cascade** (`src/michael_filter.py`, `src/triage_beta.py`): SMARTS-based
  Michael-acceptor classification (tier 1 / tier 2) and the β-carbon accessibility rule that
  reduce the South African Natural Compounds Database to the 185-compound working set.
- **Analysis and collation scripts** (`src/`, `src/schrodinger_helpers/`): dual-target rank
  consolidation, shortlist curation, figure generation, and parsers for the Glide/CovDock/Prime
  and Jaguar outputs.
- **Retrospective enrichment benchmark** (`src/benchmark/`, `data/validation/`): ChEMBL
  falcipain-2 actives and property-matched decoys, an HTVS→SP funnel docking, and enrichment
  metrics (EF, ROC-AUC, PR-AUC, BEDROC) with bootstrap 95 % confidence intervals for falcipain-2
  and falcipain-3; a warhead-matched decoy control; a docking-box size control; a
  recognition-ablation robustness check; and a pH 5.5 / 7.0 sensitivity analysis.
- **DFT reactivity data** (`data/validation/`, `data/leadopt_warhead_tuning.csv`): global
  electrophilicity ω, condensed β-carbon Fukui indices, methanethiolate reaction energies,
  population-scheme and basis-set robustness, and a ΔSCF electrophilicity cross-check.
- **Covalent pose coordinates** (`data/poses/`): CovDock covalent complexes of SANC00867 at
  FP-2 (Cys42) and FP-3 (Cys51), with a render script.
- **Regression tests and checksums** (`tests/`, `SHA256SUMS.txt`).

## Reproducibility

`enrichment_metrics.py` and the RDKit filtering cascade are pure open-source Python. The
docking, MM-GBSA, and quantum-chemistry steps require a licensed Schrödinger suite and Jaguar;
the helper scripts that drive and parse those jobs are under `src/schrodinger_helpers/`.

## Notes

Transition-state searches for the model thia-Michael reaction did not yield a one-step saddle
in implicit solvent, so no activation-barrier values are reported; warhead reactivity is
characterized thermodynamically and by conceptual-DFT descriptors (see the paper's Methods and
Limitations).
