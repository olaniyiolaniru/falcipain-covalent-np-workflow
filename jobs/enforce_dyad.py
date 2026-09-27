"""Hold the catalytic dyad at the standardized assignment in the pH 5.5 receptors.

The convention is a neutral cysteine thiol carrying its S-gamma hydrogen and a
neutral catalytic histidine protonated at N-epsilon-2. Applying it at both pH
values makes the pH comparison report the effect of pH on the rest of the site
and on ligand ionization.
"""
from schrodinger import structure

TARGETS = [("FP2", "FP2_pH55_prep.mae", "FP2_pH55_standard.mae", 42, 174),
           ("FP3", "FP3_pH55_prep.mae", "FP3_pH55_standard.mae", 51, 183)]

for name, src, dest, cys, his in TARGETS:
    st = list(structure.StructureReader(src))[0]
    changed = []
    for res in st.residue:
        if res.resnum == his and res.pdbres.strip() in ("HIS", "HID", "HIE", "HIP"):
            hd1 = [a for a in res.atom if a.pdbname.strip() == "HD1"]
            nd1 = [a for a in res.atom if a.pdbname.strip() == "ND1"]
            ne2h = [a for a in res.atom if a.pdbname.strip() == "HE2"]
            if hd1 and ne2h:
                st.deleteAtoms([hd1[0].index])
                changed.append("removed HD1 from His%d" % his)
            for r2 in st.residue:
                if r2.resnum == his:
                    for a in r2.atom:
                        if a.pdbname.strip() == "ND1":
                            a.formal_charge = 0
                    r2.pdbres = "HIE "
            break
    for res in st.residue:
        if res.resnum == cys and res.pdbres.strip() in ("CYS", "CYX"):
            sg = [a for a in res.atom if a.pdbname.strip() == "SG"]
            if sg:
                has_h = any((b.atom2 if b.atom1.index == sg[0].index else b.atom1).atomic_number == 1
                            for b in sg[0].bond)
                changed.append("Cys%d S-gamma hydrogen present: %s" % (cys, has_h))
            break
    st.write(dest)
    print(name, dest, "atoms", st.atom_total, "formal charge", st.formal_charge, "|", "; ".join(changed))
