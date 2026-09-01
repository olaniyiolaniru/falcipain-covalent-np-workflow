#!/usr/bin/env python3
"""
michael_filter.py
-----------------
Isolate genuine Michael-acceptor (and optionally epoxide) warheads from a
natural-product library, ahead of covalent docking against falcipain-2/-3.

Design notes
------------
* All alkene atoms in the warhead SMARTS are written as ALIPHATIC sp2 carbon
  ([CX3], uppercase), so AROMATIC C=C (benzene rings, coumarin pyranones, etc.)
  do NOT match. This is deliberate: ring-aromatic "enones" are poor electrophiles.
* The general alpha,beta-unsaturated carbonyl pattern intentionally leaves the
  non-carbonyl substituent unconstrained, so it captures enones, enals,
  acrylates/acrylamides, unsaturated acids AND unsaturated lactones in one rule.
* Hits are annotated, never silently dropped, and a coarse REACTIVITY TIER is
  assigned. Do NOT run a blanket PAINS removal on this set: most of your
  intended warheads are PAINS-flagged by design. PAINS is annotated only.
* The actual reactivity ranking is the job of your DFT step (electrophilicity
  index omega, Fukui f+, thiol DeltaE). This filter is a topological pre-screen.

Usage
-----
    python michael_filter.py -i library.smi  -o hits.csv
    python michael_filter.py -i library.sdf  -o hits.csv --sdf-out hits.sdf
    python michael_filter.py -i library.csv  -o hits.csv --smiles-col canonical_smiles --id-col id
    python michael_filter.py -i library.smi  -o hits.csv --include-epoxide

Input formats: .smi / .txt (SMILES [whitespace ID] per line), .csv, .sdf
"""

import argparse
import os
import sys
import csv

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")  # silence parse warnings; we handle failures ourselves

# ---------------------------------------------------------------------------
# Warhead definitions
# ---------------------------------------------------------------------------
# tier 1 = tunable / medicinal-chemistry-friendly electrophiles
# tier 2 = intrinsically aggressive / promiscuous (keep, but flag for scrutiny)
MICHAEL_ACCEPTORS = {
    # --- tier 1 -----------------------------------------------------------
    # General alpha,beta-unsaturated carbonyl: enone, enal, acrylate,
    # acrylamide, unsaturated acid, unsaturated lactone (chalcones, curcumin,
    # cinnamaldehyde, sesquiterpene lactones, etc.)
    "ab_unsat_carbonyl":      ("[CX3]=[CX3][CX3]=[OX1]", 1),
    # Terminal exocyclic methylene conjugated to C=O
    # (alpha-methylene-gamma-butyrolactone: parthenolide, helenalin, etc.)
    "exo_methylene_carbonyl": ("[CH2]=[CX3][CX3]=[OX1]", 1),
    # Acrylonitrile / vinyl nitrile
    "vinyl_nitrile":          ("[CX3]=[CX3][CX2]#[NX1]", 1),
    # --- tier 2 -----------------------------------------------------------
    # Vinyl sulfone (classic falcipain warhead; reactive but well precedented)
    "vinyl_sulfone":          ("[CX3]=[CX3][SX4](=[OX1])(=[OX1])", 2),
    # Maleimide (very reactive)
    "maleimide":              ("O=C1C=CC(=O)N1", 2),
    # para-quinone
    "para_quinone":           ("O=C1C=CC(=O)C=C1", 2),
    # ortho-quinone
    "ortho_quinone":          ("O=C1C=CC=CC1=O", 2),
}

# Optional secondary covalent class (E64 is an epoxysuccinate). Off by default.
EPOXIDE = {
    "epoxide":                ("[OX2r3]1[#6r3][#6r3]1", 2),
}

# PAINS annotation only (never auto-removed here)
from rdkit.Chem import FilterCatalog
_pains_params = FilterCatalog.FilterCatalogParams()
_pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
_PAINS = FilterCatalog.FilterCatalog(_pains_params)


def compile_patterns(include_epoxide: bool):
    src = dict(MICHAEL_ACCEPTORS)
    if include_epoxide:
        src.update(EPOXIDE)
    compiled = {}
    for name, (smarts, tier) in src.items():
        q = Chem.MolFromSmarts(smarts)
        if q is None:
            sys.exit(f"FATAL: bad SMARTS for {name}: {smarts}")
        compiled[name] = (q, tier)
    return compiled


_normalizer = rdMolStandardize.Normalizer()
_lfc = rdMolStandardize.LargestFragmentChooser()


def clean_mol(mol):
    """Strip salts/solvents, normalize. Returns cleaned mol or None."""
    if mol is None:
        return None
    try:
        mol = _lfc.choose(mol)        # keep largest fragment (desalt)
        mol = _normalizer.normalize(mol)
        Chem.SanitizeMol(mol)
        return mol
    except Exception:
        return None


def analyze(mol, compiled):
    """Return dict of warhead annotations, or None if no warhead present."""
    hits, tiers, total = [], set(), 0
    for name, (q, tier) in compiled.items():
        n = len(mol.GetSubstructMatches(q, uniquify=True))
        if n:
            hits.append(f"{name}:{n}")
            tiers.add(tier)
            total += n
    if not hits:
        return None
    pains = _PAINS.GetFirstMatch(mol)
    return {
        "warheads": ";".join(hits),
        "n_warheads": total,
        "max_tier": max(tiers),                       # 1 tunable, 2 aggressive
        "pains_flag": pains.GetDescription() if pains else "",
        "MW": round(Descriptors.MolWt(mol), 1),
        "cLogP": round(Crippen.MolLogP(mol), 2),
        "HBD": rdMolDescriptors.CalcNumHBD(mol),
        "HBA": rdMolDescriptors.CalcNumHBA(mol),
        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 1),
        "ArRings": rdMolDescriptors.CalcNumAromaticRings(mol),
    }


def read_input(path, smiles_col, id_col):
    """Yield (id, mol) pairs from .smi/.txt, .csv, or .sdf."""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".smi", ".txt"):
        with open(path) as fh:
            for i, line in enumerate(fh):
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                smi = parts[0]
                cid = parts[1] if len(parts) > 1 else f"mol_{i}"
                yield cid, Chem.MolFromSmiles(smi)
    elif ext == ".csv":
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            if smiles_col not in reader.fieldnames:
                sys.exit(f"FATAL: --smiles-col '{smiles_col}' not in CSV header {reader.fieldnames}")
            for i, row in enumerate(reader):
                cid = row.get(id_col, f"mol_{i}") if id_col else f"mol_{i}"
                yield cid, Chem.MolFromSmiles(row[smiles_col])
    elif ext == ".sdf":
        supplier = Chem.SDMolSupplier(path)
        for i, mol in enumerate(supplier):
            cid = mol.GetProp("_Name") if (mol and mol.HasProp("_Name") and mol.GetProp("_Name")) else f"mol_{i}"
            yield cid, mol
    else:
        sys.exit(f"FATAL: unsupported input extension '{ext}' (use .smi/.txt/.csv/.sdf)")


def main():
    ap = argparse.ArgumentParser(description="Michael-acceptor warhead filter for NP libraries")
    ap.add_argument("-i", "--input", required=True, help="library file (.smi/.txt/.csv/.sdf)")
    ap.add_argument("-o", "--output", required=True, help="output CSV of hits")
    ap.add_argument("--sdf-out", help="optional SDF of hits (carries warhead tags)")
    ap.add_argument("--smiles-col", default="SMILES", help="SMILES column name for CSV input")
    ap.add_argument("--id-col", default="", help="ID column name for CSV input")
    ap.add_argument("--include-epoxide", action="store_true",
                    help="also flag epoxides (E64-type secondary warhead)")
    ap.add_argument("--tier1-only", action="store_true",
                    help="keep only tunable (tier-1) warheads; drop quinones/maleimides/sulfones")
    ap.add_argument("--max-mw", type=float, default=None,
                    help="optional MW ceiling (NPs run large; default: no limit)")
    args = ap.parse_args()

    compiled = compile_patterns(args.include_epoxide)

    seen_inchikey = set()
    n_read = n_valid = n_hit = n_written = 0
    rows = []
    sdf_writer = Chem.SDWriter(args.sdf_out) if args.sdf_out else None

    for cid, mol in read_input(args.input, args.smiles_col, args.id_col):
        n_read += 1
        mol = clean_mol(mol)
        if mol is None:
            continue
        n_valid += 1

        ann = analyze(mol, compiled)
        if ann is None:
            continue
        n_hit += 1

        if args.tier1_only and ann["max_tier"] != 1:
            continue
        if args.max_mw is not None and ann["MW"] > args.max_mw:
            continue

        # de-duplicate by InChIKey
        try:
            ik = Chem.MolToInchiKey(mol)
        except Exception:
            ik = Chem.MolToSmiles(mol)
        if ik in seen_inchikey:
            continue
        seen_inchikey.add(ik)

        smi = Chem.MolToSmiles(mol)
        row = {"id": cid, "SMILES": smi, **ann}
        rows.append(row)
        n_written += 1

        if sdf_writer:
            mol.SetProp("_Name", str(cid))
            for k, v in ann.items():
                mol.SetProp(k, str(v))
            sdf_writer.write(mol)

    if sdf_writer:
        sdf_writer.close()

    fieldnames = ["id", "SMILES", "warheads", "n_warheads", "max_tier",
                  "pains_flag", "MW", "cLogP", "HBD", "HBA", "TPSA", "ArRings"]
    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"read              : {n_read}")
    print(f"valid (parsed)    : {n_valid}")
    print(f"warhead-bearing   : {n_hit}")
    print(f"written (unique)  : {n_written}")
    print(f"-> {args.output}" + (f"  +  {args.sdf_out}" if args.sdf_out else ""))


if __name__ == "__main__":
    main()