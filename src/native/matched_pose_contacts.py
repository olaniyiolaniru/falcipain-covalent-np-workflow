"""Structural account of the matched OH/Cl comparison at FP-2 and FP-3.

The four matched structures share one para-substituted cinnamoyl ester unit.
That unit is located by a single SMARTS query, so the same twelve heavy atoms
are identified in every member and the OH and Cl members are atom-mapped to
each other in the order the query returns.

For every returned matched pose the script records, in the receptor frame of
the pose viewer the pose came from:

  * protein residues within 4.0 A of any ligand heavy atom;
  * polar heavy-atom contacts (N/O/S pairs within 3.5 A);
  * the local environment of the para position itself - the phenol oxygen of
    the OH members, the chlorine of the Cl members - with the three nearest
    protein heavy atoms and, for chlorine, the C-Cl...X angle needed to judge
    whether the geometry is a halogen bond;
  * for the state selected by each scoring rule, the atom-mapped displacement
    of the shared cinnamoyl core between the OH and Cl members.

Pose viewer record 1 is the receptor; ligand records follow. Structure reading
uses the licensed Schrodinger Python; no calculation is started.

  run.exe python3 matched_pose_contacts.py --matched <matched_check dir> --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import csv
import json
import math

from schrodinger import structure
from schrodinger.structutils import analyze, measure

CONTACT = 4.0
HBOND = 3.5
POLAR = {"N", "O", "S"}

# Para-substituted cinnamoyl ester shared by all four matched structures.
# Match order: probe, 6 aromatic carbons, 2 vinyl carbons, acyl C, carbonyl O, ester O.
CORE_SMARTS = "[$([OX2H1]),$([Cl])]c1ccc(cc1)C=CC(=O)O"

SOURCES = [("FP2", "FP2_matched_verified_pv.maegz"),
           ("FP3", "FP3_matched_verified_pv.maegz")]

MEMBERS = {
    "Parent_SANC00867": ("glycosylated", "OH"),
    "AnalogueA": ("glycosylated", "Cl"),
    "AnalogueB": ("aglycone_ester", "OH"),
    "OPT1": ("aglycone_ester", "Cl"),
}

PAIRS = [("glycosylated", "Parent_SANC00867", "AnalogueA"),
         ("aglycone_ester", "AnalogueB", "OPT1")]


def residue_label(atom):
    return "{}:{}{}".format(atom.chain.strip() or "A", atom.pdbres.strip(), atom.resnum)


def core_match(st):
    """Single cinnamoyl-ester match; returns 1-based atom indices in query order."""
    hits = analyze.evaluate_smarts(st, CORE_SMARTS, unique_sets=True)
    if not hits:
        return None
    if len(hits) > 1:
        # keep the match whose probe atom carries the substituent of interest
        hits = sorted(hits, key=lambda h: (st.atom[h[0]].element != "Cl", h[0]))
    return hits[0]


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--matched", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    contacts, probes = [], []
    by_member = {}

    for target, fname in SOURCES:
        sts = list(structure.StructureReader(str(a.matched / fname)))
        receptor = sts[0]
        prot = [at for at in receptor.atom if at.atomic_number > 1]

        for rec, lig in enumerate(sts[1:], start=2):
            title = lig.title
            arch, group = MEMBERS.get(title, ("", ""))
            variant = lig.property.get("s_lp_Variant", title)
            dock = lig.property.get("r_i_docking_score")
            glide = lig.property.get("r_i_glide_gscore")
            core = core_match(lig)
            heavy = [at for at in lig.atom if at.atomic_number > 1]

            near, hb = {}, []
            for at in heavy:
                for pt in prot:
                    d = measure.measure_distance(at, pt)
                    if d <= CONTACT:
                        lab = residue_label(pt)
                        if lab not in near or d < near[lab]:
                            near[lab] = d
                        if d <= HBOND and at.element in POLAR and pt.element in POLAR:
                            hb.append((lab, at.element + str(at.index), pt.pdbname.strip(), round(d, 3)))
            contacts.append(dict(
                target=target, pose_record=rec, title=title, variant=variant,
                architecture=arch, para_group=group,
                DockingScore=dock, GlideScore=glide,
                formal_charge=lig.formal_charge,
                core_atoms_found=0 if core is None else len(core),
                contact_residues_4A=";".join("%s(%.2f)" % (k, v)
                                             for k, v in sorted(near.items(), key=lambda x: x[1])),
                n_contact_residues_4A=len(near),
                polar_contacts_3p5A=";".join("%s|%s-%s|%s" % t for t in sorted(hb, key=lambda x: x[3])),
                n_polar_contacts_3p5A=len(hb)))

            if core is None:
                continue
            at = lig.atom[core[0]]
            kind = "Cl" if at.element == "Cl" else "OH"
            bonded = [b.atom2 if b.atom1.index == at.index else b.atom1 for b in at.bond
                      if (b.atom2 if b.atom1.index == at.index else b.atom1).atomic_number > 1]
            row = dict(target=target, pose_record=rec, title=title, variant=variant,
                       architecture=arch, para_group=group, probe=kind,
                       probe_atom=at.element + str(at.index),
                       DockingScore=dock, GlideScore=glide)
            for n, (d, pt) in enumerate(sorted(((measure.measure_distance(at, pt), pt)
                                                for pt in prot), key=lambda x: x[0])[:3], 1):
                row["near%d_residue" % n] = residue_label(pt)
                row["near%d_atom" % n] = pt.pdbname.strip()
                row["near%d_distance_A" % n] = round(d, 3)
                row["near%d_C_X_angle_deg" % n] = (
                    round(measure.measure_bond_angle(bonded[0], at, pt), 2) if bonded else None)
                row["near%d_is_polar" % n] = pt.element in POLAR
            probes.append(row)

            key = (target, title)
            by_member.setdefault(key, []).append((rec, lig, core, dock, glide))

    core_rows = []
    for target, _ in SOURCES:
        for arch, oh_name, cl_name in PAIRS:
            for rule, field in (("minimum_DockingScore", 3), ("minimum_GlideScore", 4)):
                oh = by_member.get((target, oh_name))
                cl = by_member.get((target, cl_name))
                if not oh or not cl:
                    continue
                o = min(oh, key=lambda x: x[field])
                c = min(cl, key=lambda x: x[field])
                n = min(len(o[2]), len(c[2]))
                sq = 0.0
                for x, y in zip(o[2][:n], c[2][:n]):
                    ax, ay = o[1].atom[x], c[1].atom[y]
                    sq += (ax.x - ay.x) ** 2 + (ax.y - ay.y) ** 2 + (ax.z - ay.z) ** 2
                core_rows.append(dict(
                    target=target, architecture=arch, selection_rule=rule,
                    OH_member=oh_name, OH_pose_record=o[0], OH_DockingScore=o[3], OH_GlideScore=o[4],
                    Cl_member=cl_name, Cl_pose_record=c[0], Cl_DockingScore=c[3], Cl_GlideScore=c[4],
                    mapped_core_heavy_atoms=n,
                    core_displacement_rmsd_A=round(math.sqrt(sq / n), 3),
                    para_probe_separation_A=round(
                        measure.measure_distance(o[1].atom[o[2][0]], c[1].atom[c[2][0]]), 3),
                    Cl_minus_OH_DockingScore=round(c[3] - o[3], 6),
                    Cl_minus_OH_GlideScore=round(c[4] - o[4], 6)))

    def write(name, rows):
        with (a.out / name).open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    write("Matched_pose_contacts.csv", contacts)
    write("Matched_para_environment.csv", probes)
    write("Matched_core_displacement.csv", core_rows)

    summary = dict(poses=len(contacts), para_probes=len(probes),
                   core_comparisons=len(core_rows),
                   core_smarts=CORE_SMARTS,
                   contact_cutoff_A=CONTACT, polar_cutoff_A=HBOND,
                   note=("Contacts are geometric distance criteria in the pose-viewer "
                         "receptor frame; they are not interaction energies."))
    (a.out / "Matched_pose_contacts.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    for r in core_rows:
        print(r)


if __name__ == "__main__":
    main()
