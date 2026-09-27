"""Record the receptor treatment used by the covalent docking stage.

Covalent docking assigns the reacting cysteine itself, so the covalent panel is
prepared and scored under its own receptor convention rather than under the
recognition-screen convention. This script states that convention explicitly for
each covalent target, so the covalent stage can be read as the separate protocol
it is, and reports the SHA-256 of every receptor file the panel used.

  run.exe python3 covalent_receptor_record.py --jobs <schrodinger-jobs>/covalent_docking \
      --records ../data/Covalent_stage_records.csv --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import csv
import hashlib
import json

import pandas as pd
from schrodinger import structure

# target, receptor file, catalytic Cys, catalytic His, PDB
TARGETS = [
    ("FP2", "FP2_rec.mae", 42, 174, "3BPF"),
    ("FP3", "FP3_rec.maegz", 51, 183, "3BPM"),
    ("CatK", "CatK_rec.maegz", 25, 162, "1ATK"),
    ("CatL", "CatL_rec.maegz", 25, 163, "5MQY"),
]


def sha(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def dyad(path, cys_num, his_num):
    st = next(structure.StructureReader(str(path)))
    cys_state, his_state = "not found", "not found"
    for res in st.residue:
        name = res.pdbres.strip()
        if name in ("HIS", "HID", "HIE", "HIP") and res.resnum == his_num:
            his_state = name
        if name in ("CYS", "CYX", "CYT") and res.resnum == cys_num:
            sg = [a for a in res.atom if a.pdbname.strip() == "SG"]
            if sg:
                hs = [b for b in sg[0].bond
                      if (b.atom2 if b.atom1.index == sg[0].index else b.atom1).atomic_number == 1]
                cys_state = "neutral thiol (S-gamma H present)" if hs else "thiolate (no S-gamma H)"
    return dict(atoms=st.atom_total,
                heavy_atoms=sum(1 for a in st.atom if a.atomic_number > 1),
                residues=sum(1 for _ in st.residue),
                formal_charge=st.formal_charge,
                catalytic_cys="Cys" + str(cys_num) + ": " + cys_state,
                catalytic_his="His" + str(his_num) + ": " + his_state)


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--jobs", type=Path, required=True)
    p.add_argument("--records", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    stage = pd.read_csv(a.records)
    rows = []
    for name, filename, cys_num, his_num, pdb in TARGETS:
        path = a.jobs / filename
        block = stage[stage.target.eq(name)]
        row = dict(target=name, pdb=pdb, receptor_file=filename,
                   compounds_docked=int(block.compound.nunique()),
                   reaction_type=";".join(sorted(set(block.reaction))) if len(block) else "",
                   reactive_residue=";".join(sorted(set(str(x) for x in block.reactive_residue)))
                   if len(block) else "",
                   covdock_mode="enrichment",
                   protocol=("CovDock enrichment mode: the reacting cysteine is treated by the "
                             "covalent docking protocol itself, which forms the bond to the "
                             "specified residue and minimises the product complex, so this stage "
                             "carries its own receptor convention and its scores are not on the "
                             "same scale as the noncovalent recognition screen"))
        if path.is_file():
            row["receptor_sha256"] = sha(path)
            row.update(dyad(path, cys_num, his_num))
        rows.append(row)

    fields = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with (a.out / "Covalent_receptor_record.csv").open("w", newline="",
                                                       encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    (a.out / "Covalent_receptor_record.json").write_text(json.dumps(rows, indent=2),
                                                         encoding="utf-8")
    for r in rows:
        print("%-6s %-16s %-6s %s | %s" % (
            r["target"], r["receptor_file"], r.get("atoms", "-"),
            r.get("catalytic_cys", "-"), r.get("catalytic_his", "-")))


if __name__ == "__main__":
    main()
