#!/usr/bin/env python3
"""
triage_beta.py - split Michael-acceptor hits into ACCESSIBLE vs BURIED.

An accessible warhead requires the beta-carbon (the conjugate-addition site,
i.e. the alkene carbon NOT bonded to the carbonyl/EWG) to carry >=1 hydrogen.
Fully-substituted beta-carbons (ring-fused steroidal/triterpenoid enones) are
sterically blocked and flagged BURIED.

Reads the hits CSV (needs a SMILES column), writes:
  *_accessible.csv  - realistic working set
and prints the split.
"""
import sys, csv
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

infile  = sys.argv[1] if len(sys.argv) > 1 else "sancdb_hits.csv"
outfile = infile.rsplit(".",1)[0] + "_accessible.csv"

# Accessible-warhead SMARTS: beta-carbon must be H1 or H2.
ACCESSIBLE = {
    "acc_ab_unsat_carbonyl": "[CX3;H1,H2]=[CX3][CX3]=[OX1]",
    "acc_exo_methylene":     "[CH2]=[CX3][CX3]=[OX1]",
    "acc_vinyl_sulfone":     "[CX3;H1,H2]=[CX3][SX4](=[OX1])(=[OX1])",
    "acc_vinyl_nitrile":     "[CX3;H1,H2]=[CX3][CX2]#[NX1]",
}
# Aggressive tier-2 rings kept regardless (intrinsically reactive)
KEEP_RINGS = {
    "maleimide":    "O=C1C=CC(=O)N1",
    "para_quinone": "O=C1C=CC(=O)C=C1",
    "ortho_quinone":"O=C1C=CC=CC1=O",
}
acc_q  = {k: Chem.MolFromSmarts(v) for k,v in ACCESSIBLE.items()}
ring_q = {k: Chem.MolFromSmarts(v) for k,v in KEEP_RINGS.items()}

rows = list(csv.DictReader(open(infile)))
smi_col = "SMILES" if "SMILES" in rows[0] else list(rows[0])[1]

acc_rows, n_buried = [], 0
for r in rows:
    mol = Chem.MolFromSmiles(r[smi_col])
    if mol is None:
        continue
    acc_hits = [k for k,q in acc_q.items() if mol.HasSubstructMatch(q)]
    ring_hits= [k for k,q in ring_q.items() if mol.HasSubstructMatch(q)]
    if acc_hits or ring_hits:
        r["accessible_warheads"] = ";".join(acc_hits + ring_hits)
        acc_rows.append(r)
    else:
        n_buried += 1

fn = list(rows[0].keys())
if "accessible_warheads" not in fn:
    fn = fn + ["accessible_warheads"]
with open(outfile, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fn); w.writeheader()
    for r in acc_rows: w.writerow(r)

print(f"input hits        : {len(rows)}")
print(f"ACCESSIBLE (kept) : {len(acc_rows)}")
print(f"BURIED   (dropped): {n_buried}")
print(f"-> {outfile}")
