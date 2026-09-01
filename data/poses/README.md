# Covalent pose coordinates (SANC00867)

Best-scoring CovDock covalent complexes of SANC00867 at the catalytic cysteines of
falcipain-2 (FP-2, 3BPF) and falcipain-3 (FP-3, 3BPM), generated with Schrödinger
CovDock ("Michael Addition" reaction type, uniform mode). These are the coordinates
underlying Figure 3; `render_panels.py` reproduces both 3D panels from them
(`pymol -cq render_panels.py`, open-source PyMOL).

| file | target | catalytic Cys | CovDock score (kcal/mol) | C–S bond (Å) |
|---|---|---|---|---|
| `SANC00867_FP2_covalent.pdb` | falcipain-2 | Cys42 | −7.40 | 2.17 |
| `SANC00867_FP3_covalent.pdb` | falcipain-3 | Cys51 | −6.39 | 1.90 |

Key measured contacts (best pose):

- **FP-2:** covalent bond Cys42:SG–ligand C3 (2.17 Å); oxyanion-hole H-bond to Gln36:NE2
  (2.82 Å); backbone contact Cys39:N (3.02 Å); catalytic His174:ND1 (3.15 Å).
- **FP-3:** covalent bond Cys51:SG–ligand C3 (1.90 Å); oxyanion-hole H-bond to Gln45:NE2
  (3.00 Å); Asp44 side-chain and Gly92/Trp215 backbone subsite contacts (2.9–3.4 Å).

The covalent ligand is written as residue `UNK`. `render_panels.py` renders both targets
in one pass (`Figure3_FP2.png`, `Figure3_FP3.png`); `Figure3_covalent_poses.png` is the
assembled two-panel figure.
