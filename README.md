# Structure-Based Prioritization and Electronic Tuning of Cinnamate-Bearing Natural Products at Falcipain-2 and Falcipain-3

Analysis code and data for **"Structure-Based Prioritization and Electronic
Tuning of Cinnamate-Bearing Natural Products at Falcipain-2 and Falcipain-3"**
(*ACS Omega*), release **1.2.1**.

## What is here

| Path | Contents |
| --- | --- |
| `run_all.py` | One documented sequence that regenerates every open analysis and figure from `data/` into `outputs/` |
| `src/` | Open analysis scripts (RDKit, pandas, numpy, scipy, matplotlib only) |
| `src/native/` | Stages that require a licensed Schrodinger installation; these produced the tables in `data/` |
| `jobs/` | The docking-job drivers, so the screening, benchmark and pH stages can be repeated exactly |
| `data/` | Inputs and the extracted records every reported number is computed from |
| `outputs/` | Written by `run_all.py`; diff it against `data/` to confirm a reproduction |
| `figures/` | The article figures, the Supporting Information figure and the table-of-contents graphic |
| `tests/` | Checks of the scientific outputs and of release provenance |
| `MANIFEST.csv` | Article section, reported result, input file, output record, script and checksum |

Nothing in `data/` is modified by the analyses. All generated output goes to
`outputs/`.

## Reproduce the analyses

```bash
python -m pip install -r requirements.txt
python run_all.py            # regenerate every open analysis and figure
python run_all.py --list     # show the sequence without running it
python -m pytest tests -q    # check the outputs against the published values
```

## The catalytic-dyad convention

Both falcipain receptors used for the reported screen and the reported matched
comparison carry one catalytic-dyad assignment: a neutral cysteine thiol with an
explicit S-gamma hydrogen, and a neutral HIE histidine.
`data/Receptor_preparation_record.csv` states that assignment for every docking
stage in the study, including the alternative ion-pair receptor used for the
sensitivity experiment and the covalent-stage receptors, which carry their own
convention because covalent docking assigns the reacting cysteine itself.

The sensitivity experiment is released as a script. `src/native/dyad_sensitivity.py`
takes the two pose-viewer files and reports the contrast under each receptor, the
residue-level difference between them, and the resulting shift.
`src/native/dyad_library_effect.py` runs the same comparison across every parent
the screen docked, which is how the size of the term is established rather than
illustrated. `src/native/receptor_diff.py` performs the residue-level comparison
on any two prepared receptors, so the same control can be run on a different pair.

## Stages

| Stage | Script | Regenerates |
| --- | --- | --- |
| `provenance` | `src/library_provenance.py` | The disposition of all 1,012 source records, the duplicate-to-parent mapping, matched atom indices and the distinct-site count |
| `descriptors` | `src/library_descriptors.py` | Open descriptors and heavy-atom counts for the 185 retained parents |
| `recognition` | `src/receptor_normalised_ranking.py` | Within-target percentiles, z-scores and ordinal ranks, the six dual-target composites, their Spearman agreement, the Pareto front and the size diagnostics |
| `benchmark_xp` | `src/enrichment_metrics.py` | **The benchmark behind the reported ROC-AUC.** Glide XP over every prepared state at falcipain-2, on the receptor and grid of the screen, each compound taking its minimum DockingScore, with bootstrap intervals and both treatments of the compounds that returned no pose |
| `benchmark_xp_fp3` | `src/enrichment_metrics.py` | The same protocol-matched benchmark at falcipain-3 |
| `precision` | `src/benchmark_precision_comparison.py` | Docking precision and decoy-pool size separated, each measured with the other held fixed |
| `benchmark` | `src/enrichment_metrics.py` | Lower-precision sensitivity benchmark at falcipain-2: an HTVS-then-SP funnel over all 8,026 decoys, reported alongside the protocol-matched run |
| `benchmark_fp3` | `src/enrichment_metrics.py` | The same lower-precision funnel at falcipain-3 |
| `dyad` | `src/benchmark_dyad_comparison.py` | The two catalytic-dyad receptors restricted to the compounds both scored, with the interval on the difference from a bootstrap that resamples compounds jointly |
| `quantum` | `src/quantum_analysis.py` | omega, f+beta, omega*f+beta and the fragment partitions, from the atomic charge table |
| `matched` | `src/matched_analysis.py` | The Cl minus OH contrasts under every selection rule, the core displacements and the halogen-bond screen |
| `properties` | `src/panel_properties.py`, `src/property_ranges.py`, `src/property_concordance.py` | The panel property profile, the matched-set ranges and the two-platform concordance |
| `dyad` | `src/native/dyad_sensitivity.py`, `src/native/dyad_library_effect.py` | The catalytic-dyad term on the matched pair and across all 175 screened parents |
| `figures` | `src/build_figures.py` | Figures 1 to 4, Figure S1 and the table-of-contents graphic |

`src/build_figures.py` is the single source of truth for the panel letters used
in the article captions: the captions are written against what this script draws.

## Stages that need a licensed Schrodinger installation

The scripts in `src/native/` produced the tables in `data/`. They are released so
that the extraction from the native output files can be audited, and so that the
same extraction can be run against a fresh calculation. They require Schrodinger
2021-2 and are invoked through `$SCHRODINGER/run.exe python3`.

The docking-job drivers in `jobs/` carry the exact Glide and Protein Preparation
Wizard settings used for the screen, the retrospective benchmarks and the pH 5.5
panel.

## Data and provenance

Every released file is listed in `MANIFEST.csv` with the article section it
supports, the input it was produced from, the script that produced it and its
SHA-256. `SHA256SUMS.txt` carries the checksum of every file in the tree.

The benchmark is reported at two docking precisions. The protocol-matched runs,
`Benchmark_summary_FP2_XP.json` and `Benchmark_summary_FP3_XP.json`, score every
prepared state with Glide XP exactly as the recognition screen does, and they are
the runs behind the article's reported enrichment. The files without the `_XP`
suffix are the lower-precision HTVS-then-SP sensitivity benchmark over the full
decoy pool, kept so the protocol-matched result can be read against seven times
the decoys.

The benchmark files in `data/` are of two kinds. The inputs are the label tables
`Benchmark_labels_FP2.csv` and `Benchmark_labels_FP3.csv`, which list every
labelled compound including those that returned no pose, and the raw funnel
output for each receptor, `Benchmark_scores_FP2_standard.csv`,
`Benchmark_scores_FP2_ion_pair.csv` and `Benchmark_scores_FP3_standard.csv`. The
reference outputs those inputs produce are `Benchmark_summary*.json`,
`Benchmark_curves*.csv`, `Benchmark_scores.csv`, `Benchmark_scores_FP3.csv` and
`Benchmark_dyad_comparison.*`, kept so a reproduction can be checked against them.

- Repository: https://github.com/olaniyiolaniru/falcipain-covalent-np-workflow
- Archived release: https://doi.org/10.5281/zenodo.22132513
- Preprint: https://doi.org/10.26434/chemrxiv.15008280

## Licence

Code is released under the MIT licence (see `LICENSE`). The SANCDB source
records are redistributed under the terms of the South African Natural Compounds
Database; cite the SANCDB papers if you use them.
