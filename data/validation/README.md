# Validation data

Computational validation analyses for the workflow. All values are reproducible from the
scripts in `../../src/benchmark/` (open-source) and `../../src/schrodinger_helpers/` (licensed
Schrödinger/Jaguar).

| file | description |
|---|---|
| `scores_FP2_full.csv` | Funnel (HTVS→SP) docking scores + active/decoy labels for the full FP-2 benchmark |
| `scores_FP3.csv` | Funnel docking scores + labels for the FP-3 benchmark |
| `benchmark_FP2_full_results.txt` | EF / ROC-AUC / PR-AUC / BEDROC with bootstrap 95% CIs |
| `benchmark_FP3_results.txt` | FP-3 enrichment metrics (docked ligands and full-set treatments) |
| `actives_FP2_classified.csv` | ChEMBL FP-2 actives (≤10 µM), covalent/noncovalent classified |
| `decoys_FP2.csv` | Property-matched decoys and their matched active |
| `boxsize_control.csv` | Enrichment recomputed with oversized (>20 Å) actives removed |
| `recognition_ablation.csv` | Consensus ranking with the recognition component removed |
| `dscf_electrophilicity.csv` | ΔSCF cross-check of the frontier-orbital electrophilicity ordering |
| `pH_sensitivity_apo.csv` | Per-compound XP scores at pH 5.5 and 7.0 (apo FP-2/FP-3), + deltas |
| `OPT1_covalent_cascade.csv` | OPT1 CovDock scores at FP-2, FP-3, cathepsin K, cathepsin L |

## Headline results
- **Retrospective enrichment (FP-2):** ROC-AUC 0.533 (95% CI 0.498–0.569); recognition docking gives
  limited enrichment, so the workflow ranks by reactivity-weighted consensus rather than GlideScore.
- **Retrospective enrichment (FP-3):** ROC-AUC 0.58 (95% CI 0.44–0.71) over docked ligands (22 actives,
  350 decoys), matching the limited FP-2 enrichment.
- **Box-size control:** ROC-AUC is unchanged (0.538) when oversized actives are removed, so the limited
  enrichment is intrinsic to recognition docking, not a docking-box artifact.
- **pH robustness:** Spearman ρ = 0.77 (FP-2) and 0.97 (FP-3) between pH 7.0 and pH 5.5 rankings.
- **Reactivity descriptors:** the frontier-orbital and ΔSCF electrophilicity orderings agree (r = 0.999).

## Reproduce
```
# open-source: rebuild actives/decoys + metrics
$SCHRODINGER/run python3 ../../src/benchmark/fetch_props.py
$SCHRODINGER/run python3 ../../src/benchmark/match_decoys.py
python3 ../../src/benchmark/enrichment_metrics.py scores_FP2_full.csv
# licensed: docking funnel / pH / OPT1 (see src/schrodinger_helpers/)
```
