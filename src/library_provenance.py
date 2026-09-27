"""Regenerate the complete chemical selection trail from the source SDF.

One row per source record, with an explicit disposition and the reason for it,
the InChIKey used for deduplication, the parent a duplicate was merged into,
the matched atom indices of every motif query, and a count of distinct
electrophilic sites rather than a count of overlapping SMARTS matches.

Distinct site definition: one site per unique alkene carbon-carbon pair that a
motif query matches. A p-quinone therefore contributes two sites, a maleimide
one, and a molecule matched by both the general enone query and a ring query at
the same alkene contributes one.

The beta-hydrogen rule is a topological substitution filter applied to the
two-dimensional structure; three-dimensional accessibility is resolved by the
docking calculations that follow it.

  python library_provenance.py --sdf ../../data/SANCDB_1012_records.sdf --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import csv
import hashlib
import json
import re

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, FilterCatalog
from rdkit.Chem.MolStandardize import rdMolStandardize

RDLogger.DisableLog("rdApp.*")

# name -> (SMARTS, tier, alkene atom-index pairs within the match)
MOTIFS = {
    "ab_unsat_carbonyl":      ("[CX3]=[CX3][CX3]=[OX1]", 1, ((0, 1),)),
    "exo_methylene_carbonyl": ("[CH2]=[CX3][CX3]=[OX1]", 1, ((0, 1),)),
    "vinyl_nitrile":          ("[CX3]=[CX3][CX2]#[NX1]", 1, ((0, 1),)),
    "vinyl_sulfone":          ("[CX3]=[CX3][SX4](=[OX1])(=[OX1])", 2, ((0, 1),)),
    "maleimide":              ("O=C1C=CC(=O)N1", 2, ((2, 3),)),
    "para_quinone":           ("O=C1C=CC(=O)C=C1", 2, ((2, 3), (6, 7))),
    "ortho_quinone":          ("O=C1C=CC=CC1=O", 2, ((2, 3), (4, 5))),
}

# beta-carbon must carry at least one hydrogen for the acyclic/vinyl motifs
BETA_H = {
    "ab_unsat_carbonyl":      "[CX3;H1,H2]=[CX3][CX3]=[OX1]",
    "exo_methylene_carbonyl": "[CH2]=[CX3][CX3]=[OX1]",
    "vinyl_nitrile":          "[CX3;H1,H2]=[CX3][CX2]#[NX1]",
    "vinyl_sulfone":          "[CX3;H1,H2]=[CX3][SX4](=[OX1])(=[OX1])",
}
RING_MOTIFS = {"maleimide", "para_quinone", "ortho_quinone"}

_norm = rdMolStandardize.Normalizer()
_largest = rdMolStandardize.LargestFragmentChooser()
_params = FilterCatalog.FilterCatalogParams()
_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
_PAINS = FilterCatalog.FilterCatalog(_params)


def clean(mol):
    mol = _largest.choose(mol)
    mol = _norm.normalize(mol)
    Chem.SanitizeMol(mol)
    return mol


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--sdf", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    queries = {k: (Chem.MolFromSmarts(v[0]), v[1], v[2]) for k, v in MOTIFS.items()}
    beta_q = {k: Chem.MolFromSmarts(v) for k, v in BETA_H.items()}
    ring_q = {k: Chem.MolFromSmarts(MOTIFS[k][0]) for k in RING_MOTIFS}

    rows, sites_rows = [], []
    seen = {}
    counts = dict(records=0, parse_failure=0, no_motif=0, duplicate=0,
                  failed_selection=0, retained=0)

    supplier = Chem.ForwardSDMolSupplier(str(a.sdf), removeHs=False, sanitize=False)
    for index, mol in enumerate(supplier, start=1):
        counts["records"] += 1
        row = dict(record=index, source_id="", title="", inchikey="",
                   status="", reason="", duplicate_of="",
                   motif_matches="", matched_atom_indices="",
                   n_smarts_matches=0, n_distinct_sites=0,
                   sites_passing_beta_rule=0, max_tier="", pains_flag="",
                   MW="", cLogP="", HBD="", HBA="", TPSA="", ArRings="",
                   isomeric_smiles="")
        if mol is None:
            row.update(status="parse_failure",
                       reason="RDKit could not read or sanitise the record")
            counts["parse_failure"] += 1
            rows.append(row)
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        row["title"] = name
        m = re.search(r"(SANC\d+)", name)
        row["source_id"] = m.group(1) if m else name
        try:
            mol = clean(mol)
        except Exception as exc:
            row.update(status="parse_failure",
                       reason="standardisation failed: %s" % exc)
            counts["parse_failure"] += 1
            rows.append(row)
            continue

        matches, atoms, sites, tiers, n_smarts = [], [], set(), set(), 0
        for name_q, (q, tier, pairs) in queries.items():
            hits = mol.GetSubstructMatches(q, uniquify=True)
            if not hits:
                continue
            matches.append("%s:%d" % (name_q, len(hits)))
            tiers.add(tier)
            n_smarts += len(hits)
            for h in hits:
                atoms.append("%s[%s]" % (name_q, ",".join(str(x) for x in h)))
                for i, j in pairs:
                    sites.add(frozenset((h[i], h[j])))
        row["motif_matches"] = ";".join(matches)
        row["matched_atom_indices"] = ";".join(atoms)
        row["n_smarts_matches"] = n_smarts
        row["n_distinct_sites"] = len(sites)
        row["max_tier"] = max(tiers) if tiers else ""
        pains = _PAINS.GetFirstMatch(mol)
        row["pains_flag"] = pains.GetDescription() if pains else ""
        row["MW"] = round(Descriptors.MolWt(mol), 1)
        row["cLogP"] = round(Crippen.MolLogP(mol), 2)
        row["HBD"] = rdMolDescriptors.CalcNumHBD(mol)
        row["HBA"] = rdMolDescriptors.CalcNumHBA(mol)
        row["TPSA"] = round(rdMolDescriptors.CalcTPSA(mol), 1)
        row["ArRings"] = rdMolDescriptors.CalcNumAromaticRings(mol)
        row["isomeric_smiles"] = Chem.MolToSmiles(mol)

        if not matches:
            row.update(status="no_queried_motif",
                       reason="no queried electrophilic motif present")
            counts["no_motif"] += 1
            rows.append(row)
            continue

        try:
            key = Chem.MolToInchiKey(mol)
        except Exception:
            key = row["isomeric_smiles"]
        row["inchikey"] = key
        if key in seen:
            row.update(status="duplicate_motif_hit", duplicate_of=seen[key],
                       reason="identical InChIKey to an earlier motif-bearing record")
            counts["duplicate"] += 1
            rows.append(row)
            continue
        seen[key] = row["source_id"]

        passing = set()
        for name_q, q in beta_q.items():
            for h in mol.GetSubstructMatches(q, uniquify=True):
                i, j = MOTIFS[name_q][2][0]
                passing.add(frozenset((h[i], h[j])))
        for name_q, q in ring_q.items():
            for h in mol.GetSubstructMatches(q, uniquify=True):
                for i, j in MOTIFS[name_q][2]:
                    passing.add(frozenset((h[i], h[j])))
        row["sites_passing_beta_rule"] = len(passing & sites)

        if row["sites_passing_beta_rule"]:
            row.update(status="retained",
                       reason="at least one site passes the beta-hydrogen or retained-ring rule")
            counts["retained"] += 1
        else:
            row.update(status="failed_selection_rule",
                       reason="every matched site has a fully substituted beta-carbon "
                              "and no retained ring motif")
            counts["failed_selection"] += 1
        rows.append(row)

        for s in sorted(sites, key=lambda f: sorted(f)):
            i, j = sorted(s)
            sites_rows.append(dict(
                record=index, source_id=row["source_id"], site_atoms="%d,%d" % (i, j),
                site_elements="%s,%s" % (mol.GetAtomWithIdx(i).GetSymbol(),
                                         mol.GetAtomWithIdx(j).GetSymbol()),
                beta_rule_pass=s in passing,
                motifs=";".join(n for n, (q, _t, pairs) in queries.items()
                                for h in mol.GetSubstructMatches(q, uniquify=True)
                                for pi, pj in pairs if frozenset((h[pi], h[pj])) == s)))

    with (a.out / "Library_provenance_full.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with (a.out / "Library_electrophilic_sites.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sites_rows[0]))
        w.writeheader()
        w.writerows(sites_rows)

    motif_bearing = counts["duplicate"] + counts["failed_selection"] + counts["retained"]
    summary = dict(
        source_sdf=str(a.sdf),
        source_sha256=hashlib.sha256(a.sdf.read_bytes()).hexdigest(),
        records=counts["records"],
        parse_failures=counts["parse_failure"],
        parsed=counts["records"] - counts["parse_failure"],
        without_queried_motif=counts["no_motif"],
        motif_bearing_records=motif_bearing,
        duplicate_motif_hits=counts["duplicate"],
        unique_motif_bearing=motif_bearing - counts["duplicate"],
        failed_selection_rule=counts["failed_selection"],
        retained_parents=counts["retained"],
        distinct_sites_in_retained=sum(r["n_distinct_sites"] for r in rows if r["status"] == "retained"),
        smarts_matches_in_retained=sum(r["n_smarts_matches"] for r in rows if r["status"] == "retained"),
        beta_rule="topological substitution filter on the two-dimensional structure",
    )
    (a.out / "Library_provenance_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
