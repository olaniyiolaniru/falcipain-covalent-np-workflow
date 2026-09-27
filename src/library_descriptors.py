"""Calculate open molecular descriptors for the retained SANCDB structures.

Identifiers are taken from the SDF in the order the records actually carry
them: the lowercase `id` property written by the filtering step, then the
molecule title. If neither is present the script fails with the record index
rather than writing a blank identifier, because a blank identifier silently
breaks every downstream join.

The heavy-atom count required by the size-scaling analysis in Section 3.1 is
included.

  python library_descriptors.py --sdf data/retained_185.sdf --out outputs/Library_descriptors.csv
"""
from pathlib import Path
from argparse import ArgumentParser
import re
import sys

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")
ID_PROPS = ("id", "compound", "ID", "SANC_id", "_Name")


def identifier(mol, index):
    """Prefer a value that carries a SANCDB accession; fall back to any non-empty one."""
    values = []
    for prop in ID_PROPS:
        if prop == "_Name":
            value = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        else:
            value = mol.GetProp(prop) if mol.HasProp(prop) else ""
        value = (value or "").strip()
        if value:
            values.append(value)
    for value in values:
        m = re.search(r"(SANC\d+)", value)
        if m:
            return m.group(1)
    if values:
        return values[0]
    raise SystemExit(
        "record %d carries none of the identifier properties %s; refusing to write a "
        "blank identifier" % (index, ", ".join(ID_PROPS)))


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--sdf", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    if not a.sdf.is_file():
        raise SystemExit("input SDF not found: %s" % a.sdf)
    a.out.parent.mkdir(parents=True, exist_ok=True)

    rows, skipped = [], 0
    for index, mol in enumerate(Chem.SDMolSupplier(str(a.sdf), removeHs=False), 1):
        if mol is None:
            skipped += 1
            continue
        rows.append({
            "compound": identifier(mol, index),
            "heavy_atoms": mol.GetNumHeavyAtoms(),
            "MW": round(Descriptors.MolWt(mol), 3),
            "cLogP": round(Crippen.MolLogP(mol), 3),
            "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 2),
            "HBD": Lipinski.NumHDonors(mol),
            "HBA": Lipinski.NumHAcceptors(mol),
            "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "formal_charge": Chem.GetFormalCharge(mol),
            "SMILES": Chem.MolToSmiles(mol),
        })
    if not rows:
        raise SystemExit("no readable records in %s" % a.sdf)
    frame = pd.DataFrame(rows)
    duplicates = frame.compound[frame.compound.duplicated()].tolist()
    if duplicates:
        raise SystemExit("duplicate identifiers in %s: %s" % (a.sdf, sorted(set(duplicates))))
    frame.to_csv(a.out, index=False)
    print("wrote %d descriptor records to %s (%d unreadable records skipped)"
          % (len(frame), a.out, skipped), file=sys.stderr)


if __name__ == "__main__":
    main()
