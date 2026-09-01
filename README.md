# Reactivity-weighted covalent natural-product screening workflow
**Authors:** Olaniyi Victor Olaniru, Clement Odunayo Ajiboye. **Corresponding author:** co.ajiboye@ui.edu.ng

A reproducible, reactivity-weighted *in silico* workflow for discovering **covalent
natural-product inhibitors**, demonstrated against the *Plasmodium falciparum*
haemoglobinases **falcipain-2 and falcipain-3**. The pipeline combines a covalent-warhead
filtering cascade, dual recognition and covalent docking, Prime MM-GBSA rescoring, a
multi-descriptor DFT treatment of warhead reactivity, ADMET prediction, and a covalent
counter-screen against human cathepsins K and L.

This repository holds the **open-source, tool-agnostic components** (the RDKit filtering
cascade and the analysis/collation scripts). The docking, MM-GBSA and DFT steps use the
commercial Schrödinger suite (Glide, CovDock, Prime) and Jaguar; the helper scripts that
drive and parse those jobs are provided under `src/schrodinger_helpers/` for users with a
licence.

## What the workflow does

```
library (SMILES/SDF)
   └─ michael_filter.py        SMARTS-based Michael-acceptor / warhead classifier (tier 1 / tier 2)
        └─ triage_beta.py      β-carbon accessibility filter (accessible vs sterically buried)
             └─ [Glide XP]     dual FP-2 / FP-3 recognition docking
                  └─ [CovDock] covalent docking at the catalytic cysteine (Michael addition)
                       └─ [Prime MM-GBSA] relative rescoring
                            └─ [Jaguar DFT] ω, condensed Fukui f⁺, methanethiolate ΔE  (rank + tune)
                                 └─ [Glide/CovDock] human cathepsin K/L counter-screen
```

The two RDKit stages (`michael_filter.py`, `triage_beta.py`) reproduce the reduction of
the South African Natural Compounds Database from 1,012 entries to the 185-compound
β-accessible working set. The DFT reactivity engine and its use for prospective warhead
tuning are the methodological core of the accompanying paper.

## Installation

```bash
python -m pip install -r requirements.txt
```

Open-source dependencies only (RDKit, NumPy, pandas, openpyxl, matplotlib). Python ≥ 3.9.

## Usage

```bash
# 1. Isolate Michael-acceptor warheads (tier 1/2), annotate PAINS, keep everything
python src/michael_filter.py -i library.sdf -o hits.csv --sdf-out hits.sdf

# 2. Keep only β-accessible warheads (β-carbon bears ≥1 H)
python src/triage_beta.py hits.csv        # writes hits_accessible.csv

# 3. (Schrödinger) dock, covalently dock, MM-GBSA, DFT (see src/schrodinger_helpers/)
#    then collate:
"$SCHRODINGER/run" python3 src/schrodinger_helpers/collate_covdock_results.py
"$SCHRODINGER/run" python3 src/schrodinger_helpers/compute_reaction_energies.py --dft-dir <jaguar_out_dir>
```

`data/` contains the processed datasets for the falcipain demonstration (warhead hits,
185-compound working set, Hammett DFT series, uniform-mode covalent-docking ranking,
covalent selectivity, and the matched-pair scores). The complete 18-sheet workbook is
distributed with the paper as `Supplementary_Dataset.xlsx`.

## Citation

If you use this workflow, please cite the accompanying paper (details on acceptance) and
this repository. A `CITATION.cff` is provided; an archival release is deposited on Zenodo
(DOI assigned at release).

## Licence

MIT (see `LICENSE`). The Schrödinger-dependent helper scripts require a valid Schrödinger
licence to run; the scripts themselves are released under the same terms.
