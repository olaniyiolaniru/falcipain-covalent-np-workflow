"""Screen every chlorinated matched pose for halogen-bond geometry.

A halogen bond requires the acceptor to lie close to the extension of the C-X
bond: conventionally C-X...A >= 155 deg with X...A at or below the sum of the
van der Waals radii (Cl + O = 3.27 A, Cl + N = 3.30 A, Cl + S = 3.55 A).
Every protein N/O/S within 4.0 A of each ligand chlorine is reported with its
angle so that the criterion can be checked rather than assumed.
"""
from pathlib import Path
from argparse import ArgumentParser
import csv, json
from schrodinger import structure
from schrodinger.structutils import measure

VDW_SUM = {"O": 3.27, "N": 3.30, "S": 3.55}
ANGLE_MIN = 155.0
CUTOFF = 4.0
SOURCES = [("FP2", "FP2_matched_verified_pv.maegz"), ("FP3", "FP3_matched_verified_pv.maegz")]


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--matched", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    rows = []
    for target, fname in SOURCES:
        sts = list(structure.StructureReader(str(a.matched / fname)))
        prot = [x for x in sts[0].atom if x.element in VDW_SUM]
        for rec, lig in enumerate(sts[1:], start=2):
            for cl in (x for x in lig.atom if x.element == "Cl"):
                c = [b.atom2 if b.atom1.index == cl.index else b.atom1 for b in cl.bond
                     if (b.atom2 if b.atom1.index == cl.index else b.atom1).atomic_number > 1]
                if not c:
                    continue
                for pt in prot:
                    d = measure.measure_distance(cl, pt)
                    if d > CUTOFF:
                        continue
                    ang = measure.measure_bond_angle(c[0], cl, pt)
                    rows.append(dict(
                        target=target, pose_record=rec, title=lig.title,
                        DockingScore=lig.property.get("r_i_docking_score"),
                        acceptor="{}:{}{}".format(pt.chain.strip() or "A", pt.pdbres.strip(), pt.resnum),
                        acceptor_atom=pt.pdbname.strip(), acceptor_element=pt.element,
                        Cl_acceptor_distance_A=round(d, 3), C_Cl_acceptor_angle_deg=round(ang, 2),
                        vdW_sum_A=VDW_SUM[pt.element],
                        meets_distance=bool(d <= VDW_SUM[pt.element]),
                        meets_angle=bool(ang >= ANGLE_MIN),
                        halogen_bond=bool(d <= VDW_SUM[pt.element] and ang >= ANGLE_MIN)))
    with (a.out / "Matched_halogen_bond_check.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    summary = dict(
        contacts_screened=len(rows),
        chlorinated_poses=len({(r["target"], r["pose_record"]) for r in rows}),
        angle_criterion_deg=ANGLE_MIN, distance_criterion="element-specific vdW sum",
        contacts_meeting_distance=sum(r["meets_distance"] for r in rows),
        contacts_meeting_angle=sum(r["meets_angle"] for r in rows),
        halogen_bonds_found=sum(r["halogen_bond"] for r in rows),
        max_angle_deg=max(r["C_Cl_acceptor_angle_deg"] for r in rows),
        min_distance_A=min(r["Cl_acceptor_distance_A"] for r in rows))
    (a.out / "Matched_halogen_bond_check.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
