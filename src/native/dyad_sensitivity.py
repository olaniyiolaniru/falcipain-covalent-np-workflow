"""Measure what the catalytic-dyad assignment contributes to a matched contrast.

Two FP-2 receptors are compared. They are the same prepared chain A, with
identical heavy-atom coordinates, and they differ by the position of a single
proton on the catalytic dyad: the standardized receptor carries a neutral Cys42
thiol with an S-gamma hydrogen and a neutral HIE His174, while the alternative
receptor carries the thiolate-imidazolium ion pair, with that hydrogen on the
His174 N-delta-1 instead. The grid box, force field, precision, prepared ligand
states and selection rule are identical, so re-docking the same twelve matched
states into the two receptors isolates the dyad term.

The script reports, for each receptor: the selected 4-OH and 4-Cl scores under
both selection rules, the Cl-minus-OH contrast, the selected prepared state, and
the residue-level receptor difference including partial charges at the dyad.

  run.exe python3 dyad_sensitivity.py \
      --standard ../jobs/matched_standard/FP2_matched_verified_pv.maegz \
      --alternative ../jobs/matched_alternative/FP2_matched_verified_pv.maegz \
      --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import json

import numpy as np
from schrodinger import structure

OH, CL = "Parent_SANC00867", "AnalogueA"
CYS, HIS = 42, 174


def dyad_profile(st):
    out = {}
    for res in st.residue:
        if res.resnum not in (CYS, HIS):
            continue
        atoms = [(a.pdbname.strip(), a.element, round(float(a.partial_charge), 4))
                 for a in res.atom if a.element in ("S", "N") or a.pdbname.strip() in
                 ("HG", "HD1", "HE2")]
        out["res%d" % res.resnum] = dict(
            resname=res.pdbres.strip(),
            atoms=sum(1 for _ in res.atom),
            formal_charge=sum(a.formal_charge for a in res.atom),
            key_atoms=sorted(atoms))
    return out


def read_run(path):
    records = list(structure.StructureReader(str(path)))
    receptor, poses = records[0], records[1:]
    rows = []
    for n, st in enumerate(poses, 1):
        rows.append(dict(record=n, title=st.title,
                         variant=st.property.get("s_lp_Variant", st.title),
                         DockingScore=float(st.property["r_i_docking_score"]),
                         GlideScore=float(st.property["r_i_glide_gscore"])))
    return receptor, rows


def contrast(rows, field):
    def best(title):
        block = [r for r in rows if r["title"] == title]
        return min(block, key=lambda r: r[field])
    a, b = best(OH), best(CL)
    return dict(OH_value=a[field], Cl_value=b[field],
                Cl_minus_OH=b[field] - a[field],
                OH_state=a["variant"], Cl_state=b["variant"],
                OH_record=a["record"], Cl_record=b["record"])


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--standard", type=Path, required=True)
    p.add_argument("--alternative", type=Path, required=True)
    p.add_argument("--reproduce", type=Path,
                   help="independent re-run of the alternative receptor, used as a "
                        "reproduction check on the reported shift")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    out = {}
    receptors = {}
    for label, path in (("standardized", a.standard), ("alternative", a.alternative)):
        receptor, rows = read_run(path)
        receptors[label] = receptor
        ds = contrast(rows, "DockingScore")
        gs = contrast(rows, "GlideScore")
        out[label] = dict(
            source="/".join(path.parts[-2:]),
            poses_returned=len(rows),
            receptor_atoms=receptor.atom_total,
            receptor_heavy_atoms=sum(1 for a in receptor.atom
                                     if a.atomic_number > 1),
            receptor_residues=sum(1 for _ in receptor.residue),
            receptor_formal_charge=receptor.formal_charge,
            dyad=dyad_profile(receptor),
            OH_DockingScore=ds["OH_value"], Cl_DockingScore=ds["Cl_value"],
            Cl_minus_OH_DockingScore=ds["Cl_minus_OH"],
            OH_GlideScore=gs["OH_value"], Cl_GlideScore=gs["Cl_value"],
            Cl_minus_OH_GlideScore=gs["Cl_minus_OH"],
            OH_selected_state=ds["OH_state"], Cl_selected_state=ds["Cl_state"])

    # heavy-atom agreement between the two receptor models
    # Heavy atoms are paired in file order, after checking that the two
    # receptors carry the same elements in the same sequence. Keying on chain,
    # residue number and atom name instead would merge the atoms of residues
    # that share a number, and silently compare fewer atoms than are present.
    def heavy(st):
        return [(at.atomic_number, (at.x, at.y, at.z))
                for at in st.atom if at.atomic_number > 1]
    ha, hb = heavy(receptors["standardized"]), heavy(receptors["alternative"])
    if [e for e, _ in ha] != [e for e, _ in hb]:
        raise SystemExit("the two receptors do not share a heavy-atom sequence")
    disp = [float(np.linalg.norm(np.array(x) - np.array(y)))
            for (_, x), (_, y) in zip(ha, hb)]
    shared = ha

    out["comparison"] = dict(
        heavy_atoms_compared=len(shared),
        max_heavy_atom_displacement_A=round(max(disp), 6) if disp else None,
        heavy_atom_rmsd_A=round(float(np.sqrt(np.mean(np.square(disp)))), 6) if disp else None,
        shift_DockingScore=round(out["alternative"]["Cl_minus_OH_DockingScore"]
                                 - out["standardized"]["Cl_minus_OH_DockingScore"], 4),
        shift_GlideScore=round(out["alternative"]["Cl_minus_OH_GlideScore"]
                               - out["standardized"]["Cl_minus_OH_GlideScore"], 4),
        OH_score_shift=round(out["alternative"]["OH_DockingScore"]
                             - out["standardized"]["OH_DockingScore"], 4),
        variable_changed=("position of one catalytic-dyad proton: Cys%d S-gamma hydrogen "
                          "in the standardized receptor, His%d N-delta-1 hydrogen in the "
                          "alternative receptor" % (CYS, HIS)),
        held_fixed=("prepared chain A, grid centre and box, force field, docking precision, "
                    "prepared ligand states, selection rule"))

    if a.reproduce and a.reproduce.is_file():
        _, rows = read_run(a.reproduce)
        ds = contrast(rows, "DockingScore")
        ref = out["alternative"]["Cl_minus_OH_DockingScore"]
        out["reproduction_check"] = dict(
            source=a.reproduce.name,
            description=("the alternative-dyad receptor re-docked as an independent job, "
                         "from the same prepared ligand states and grid archive"),
            Cl_minus_OH_DockingScore=ds["Cl_minus_OH"],
            difference_from_reported=round(ds["Cl_minus_OH"] - ref, 9),
            reproduces_exactly=bool(abs(ds["Cl_minus_OH"] - ref) < 1e-6))

    (a.out / "Dyad_sensitivity.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
