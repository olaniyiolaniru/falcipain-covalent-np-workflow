"""Headless render of the SANC00867 covalent-pose panels (Figure 3).

Run from this directory:  pymol -cq render_panels.py
Writes Figure3_FP2.png and Figure3_FP3.png from the deposited PDB coordinates.
"""
import os
from pymol import cmd

HERE = os.path.dirname(os.path.abspath(__file__))

def panel(pdb, cys, extra_resi, labels, out, label_offsets=None):
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.set("label_size", 22)
    cmd.set("label_color", "black")
    cmd.set("label_outline_color", "white")
    cmd.load(pdb, "cplx")
    cmd.remove("solvent")
    cmd.remove("hydro")
    cmd.hide("everything")
    cmd.show("cartoon", "polymer")
    cmd.color("grey80", "polymer")
    cmd.set("cartoon_transparency", 0.55)
    cmd.select("lig", "resn UNK")
    cmd.show("sticks", "lig")
    cmd.color("yellow", "lig and elem C")
    keys = "+".join(str(r) for r in ([cys] + extra_resi))
    cmd.select("key", f"(resi {keys}) and polymer")
    cmd.show("sticks", "key")
    cmd.color("cyan", "key and elem C")
    try:
        cmd.distance("cbond", "lig and name C3", f"resi {cys} and name SG")
        cmd.hide("labels", "cbond")
        cmd.set("dash_color", "black"); cmd.set("dash_width", 4)
    except Exception:
        pass
    for resi, txt in labels.items():
        cmd.label(f"resi {resi} and name CA", f'"{txt}"')
    if label_offsets:
        for resi, off in label_offsets.items():
            cmd.set("label_position", off, f"resi {resi} and name CA")
    cmd.util.cnc("lig"); cmd.util.cnc("key")
    cmd.set("stick_radius", 0.16)
    cmd.orient("lig")
    cmd.zoom("lig", 5)
    cmd.turn("y", 8)
    cmd.ray(2200, 1650)
    cmd.png(out, dpi=300)
    print("wrote", out)

# FP-2: Cys42 covalent; Gln36 oxyanion hole; His174 catalytic
panel(os.path.join(HERE, "SANC00867_FP2_covalent.pdb"), 42, [36, 174],
      {42: "Cys42", 36: "Gln36", 174: "His174"}, os.path.join(HERE, "Figure3_FP2.png"))
# FP-3: Cys51 covalent; Gln45 oxyanion hole; Asp44 subsite
panel(os.path.join(HERE, "SANC00867_FP3_covalent.pdb"), 51, [45, 44],
      {51: "Cys51", 45: "Gln45", 44: "Asp44"}, os.path.join(HERE, "Figure3_FP3.png"),
      label_offsets={44: (-3.0, 2.0, 0.0), 45: (1.5, -3.0, 0.0), 51: (2.5, 1.0, 0.0)})
print("DONE")
