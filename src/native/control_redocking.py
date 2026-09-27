"""Covalent redocking controls, scored against the deposited product geometry.

Two controls are evaluated with the same procedure:

  3BWK_K11017   the vinyl sulfone K11017 redocked into its own falcipain-3
                complex at Cys51 - the target-specific control
  2OZ2_K11777   the vinyl sulfone K11777 redocked into cruzain at Cys25 - a
                cross-family control of the same reaction type

The docked ligand is compared with the deposited ligand after superposing the
two receptors on their chain-A Calpha atoms. The ligand is never fitted to the
reference: the reported RMSD is the residual in the receptor frame. Atom
correspondence comes from a maximum common substructure with exact bond-order
matching and ring constraints, and every symmetry-equivalent mapping is
enumerated so that the reported value is the best available correspondence
rather than an arbitrary one.

E64 is deliberately absent. Its epoxide-opening chemistry is not a Michael
addition, so reproducing it with a Michael-addition reaction definition would
not test the protocol used in this work.

  run.exe python3 control_redocking.py --control-3bwk <jobdir> \\
      --jobs <schrodinger-jobs> --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import hashlib
import json

import numpy as np
from schrodinger import structure


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ca_map(st):
    return {(a.chain.strip() or "A", a.resnum): np.array(a.xyz)
            for a in st.atom if a.pdbname.strip() == "CA"}


def kabsch(mobile, reference):
    mc, rc = mobile.mean(0), reference.mean(0)
    u, _, vt = np.linalg.svd((mobile - mc).T @ (reference - rc))
    rot = u @ vt
    if np.linalg.det(rot) < 0:
        u[:, -1] *= -1
        rot = u @ vt
    residual = float(np.sqrt(np.mean(np.sum(((mobile - mc) @ rot + rc - reference) ** 2, axis=1))))
    return rot, mc, rc, residual


PROTEIN = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR "
              "TRP TYR VAL HID HIE HIP CYX ASH GLH LYN HOH SO4 ACE NMA".split())


def ligand_atoms(st, resname=None, chain=None):
    """Heavy atoms of one covalent-ligand residue.

    The CovDock attachment property marks a single atom, so residue identity is
    used instead. When several copies are present - a multi-chain deposited
    structure, for example - only the requested chain is taken.
    """
    groups = {}
    for a in st.atom:
        if a.atomic_number <= 1:
            continue
        name = a.pdbres.strip()
        if name in PROTEIN:
            continue
        if resname and name != resname:
            continue
        key = (a.chain.strip(), a.resnum, name)
        groups.setdefault(key, []).append(a.index)
    if chain is not None:
        groups = {k: v for k, v in groups.items() if k[0] == chain.strip()}
    if not groups:
        return []
    # the covalent ligand is the largest non-protein residue present
    return sorted(max(groups.values(), key=len))


def best_rmsd(ref_st, ref_idx, dock_st, dock_idx, rot, mc, rc, attach_local=None):
    """Symmetry-aware correspondence via RDKit MCS on the two extracted ligands."""
    from rdkit import Chem
    from rdkit.Chem import rdFMCS

    def to_rdkit(st, idx, name):
        sub = st.extract(idx)
        path = Path(name)
        sub.write(str(path))
        mol = next(iter(Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)))
        path.unlink(missing_ok=True)
        mol = Chem.RemoveHs(mol, sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)
        return mol

    ref = to_rdkit(ref_st, ref_idx, "_ctrl_ref.sdf")
    dock = to_rdkit(dock_st, dock_idx, "_ctrl_dock.sdf")
    mcs = rdFMCS.FindMCS([ref, dock], timeout=60, ringMatchesRingOnly=True,
                         completeRingsOnly=True,
                         atomCompare=rdFMCS.AtomCompare.CompareElements,
                         bondCompare=rdFMCS.BondCompare.CompareOrderExact)
    query = Chem.MolFromSmarts(mcs.smartsString)
    rm = ref.GetSubstructMatches(query, uniquify=False, maxMatches=2000)
    dm = dock.GetSubstructMatches(query, uniquify=False, maxMatches=2000)
    rp = np.array(ref.GetConformer().GetPositions())
    dp = (np.array(dock.GetConformer().GetPositions()) - mc) @ rot + rc
    best = None
    for a in rm:
        for b in dm:
            value = float(np.sqrt(np.mean(np.sum((rp[list(a)] - dp[list(b)]) ** 2, axis=1))))
            if best is None or value < best[0]:
                best = (value, list(zip(a, b)))
    out = dict(reference_heavy_atoms=ref.GetNumAtoms(),
               docked_heavy_atoms=dock.GetNumAtoms(),
               common_atoms=mcs.numAtoms, mcs_cancelled=bool(mcs.canceled),
               ligand_receptor_frame_rmsd_A=None if best is None else best[0],
               mapping_zero_based=[] if best is None else best[1])
    if best is not None and attach_local is not None:
        # resolve the control at the reaction centre: the protocol has to place
        # the reacting end, whereas distal flexible groups are free to differ
        import collections
        adj = collections.defaultdict(set)
        for b in dock.GetBonds():
            adj[b.GetBeginAtomIdx()].add(b.GetEndAtomIdx())
            adj[b.GetEndAtomIdx()].add(b.GetBeginAtomIdx())
        depth, frontier, seen = {attach_local: 0}, [attach_local], {attach_local}
        while frontier:
            nxt = []
            for x in frontier:
                for y in adj[x]:
                    if y not in seen:
                        seen.add(y)
                        depth[y] = depth[x] + 1
                        nxt.append(y)
            frontier = nxt
        for cutoff in (4, 6):
            pairs = [(i, j) for i, j in best[1] if depth.get(j, 99) <= cutoff]
            if pairs:
                sq = np.mean([np.sum((rp[i] - dp[j]) ** 2) for i, j in pairs])
                out["reaction_centre_rmsd_A_within_%d_bonds" % cutoff] = float(np.sqrt(sq))
                out["reaction_centre_atoms_within_%d_bonds" % cutoff] = len(pairs)
        per_atom = sorted(float(np.linalg.norm(rp[i] - dp[j])) for i, j in best[1])
        out["median_atom_deviation_A"] = per_atom[len(per_atom) // 2]
        out["max_atom_deviation_A"] = per_atom[-1]
    return out


def best_record(path):
    """The record carrying a covalent ligand and the best available score."""
    best, best_score = None, None
    for st in structure.StructureReader(str(path)):
        if not ligand_atoms(st):
            continue
        score = st.property.get("r_i_docking_score",
                                st.property.get("r_i_sample_docking_score"))
        if best is None or (score is not None and (best_score is None or score < best_score)):
            best, best_score = st, score
    return best, best_score


def evaluate(reference_complex, docked_complex, label, reaction,
             ref_resname=None, ref_chain="A"):
    ref_st = list(structure.StructureReader(str(reference_complex)))[0]
    dock_st, dock_score = best_record(docked_complex)
    if dock_st is None:
        return dict(control=label, reaction=reaction, evaluated=False,
                    reason="no record in the docked file carries a covalent ligand")
    r, d = ca_map(ref_st), ca_map(dock_st)
    keys = sorted(set(r) & set(d))
    if len(keys) < 50:
        raise SystemExit("%s: only %d shared Calpha atoms" % (label, len(keys)))
    R = np.array([r[k] for k in keys])
    D = np.array([d[k] for k in keys])
    rot, mc, rc, ca_rmsd = kabsch(D, R)

    ref_idx = ligand_atoms(ref_st, resname=ref_resname, chain=ref_chain)
    dock_idx = ligand_atoms(dock_st)
    if not ref_idx or not dock_idx:
        return dict(control=label, reaction=reaction, evaluated=False,
                    reason=("the covalent ligand could not be identified by the mapping rule in "
                            "one of the two structures; reference heavy atoms=%d, docked "
                            "heavy atoms=%d" % (len(ref_idx), len(dock_idx))))

    out = dict(control=label, reaction=reaction,
               reference_complex=str(reference_complex),
               reference_sha256=sha(reference_complex),
               docked_complex=str(docked_complex),
               docked_sha256=sha(docked_complex),
               protein_CA_count=len(keys), protein_CA_rmsd_A=ca_rmsd,
               docked_score=dock_score, reference_chain=ref_chain,
               method=("chain-A Calpha superposition of the two receptors; element and exact "
                       "bond-order MCS with ring constraints; all symmetry mappings "
                       "enumerated; no ligand-only fit"))
    attach = [a.index for a in dock_st.atom
              if a.property.get("i_cdock_lig_attach") is not None]
    attach_local = dock_idx.index(attach[0]) if attach and attach[0] in dock_idx else None
    out["attachment_atom_index"] = attach[0] if attach else None
    out.update(best_rmsd(ref_st, ref_idx, dock_st, dock_idx, rot, mc, rc, attach_local))
    return out


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--control-3bwk", type=Path, required=True,
                   help="directory holding the 3BWK covalent redocking output")
    p.add_argument("--jobs", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    results = {}
    docked = sorted(a.control_3bwk.glob("*_LO-out.maegz")) or \
        sorted(a.control_3bwk.glob("*control-out.maegz")) or \
        sorted(a.control_3bwk.glob("*-out.maegz"))
    ref = a.jobs / "covalent_docking" / "3BWK_split_receptor1.mae"
    ref_complex = a.jobs / "covalent_docking" / "3BWK_prep.maegz"
    if docked and (ref_complex.is_file() or ref.is_file()):
        results["3BWK_K11017"] = evaluate(
            ref_complex if ref_complex.is_file() else ref, docked[0],
            "3BWK_K11017", "Michael addition at Cys51 of falcipain-3, CovDock enrichment mode", "C1P")

    k777 = a.jobs / "covdock_K777_control-out.maegz"
    k777_ref = a.jobs / "2OZ2_prep.maegz"
    if k777.is_file() and k777_ref.is_file():
        results["2OZ2_K11777"] = evaluate(
            k777_ref, k777, "2OZ2_K11777", "Michael addition at Cys25 of cruzain", None)

    results["E64"] = dict(
        control="E64", excluded=True,
        reason=("E64 forms its adduct by epoxide opening, not by Michael addition, so it "
                "cannot test a Michael-addition reaction definition and is not used as a "
                "covalent-docking control."))
    (a.out / "Control_redocking.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    for k, v in results.items():
        if v.get("evaluated") is False:
            print(k, "not evaluated:", v["reason"][:80])
        elif v.get("excluded"):
            print(k, "excluded:", v["reason"][:60], "...")
        else:
            print(k, "ligand RMSD %.2f A over %d atoms; CA RMSD %.3f A over %d atoms"
                  % (v["ligand_receptor_frame_rmsd_A"], v["common_atoms"],
                     v["protein_CA_rmsd_A"], v["protein_CA_count"]))


if __name__ == "__main__":
    main()
