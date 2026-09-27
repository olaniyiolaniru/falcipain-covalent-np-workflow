"""Residue-level difference between two prepared receptors.

Reports every residue whose assigned protonation state, formal charge or heavy
atom count differs, so that a controlled comparison can state exactly which
assignments change and which are held fixed.

The two inputs are the pose-viewer files written by the docking runs being
compared, so the receptors profiled here are the receptors that were scored
rather than a separately prepared copy of them.

  run.exe python3 receptor_diff.py --a <run_a_pv.maegz> --b <run_b_pv.maegz> \n      --label-a A --label-b B --out ../exports/FP2_dyad_receptor_difference.json
"""
from argparse import ArgumentParser
from pathlib import Path
import json

from schrodinger import structure


def profile(path):
    """Residue profile of the first structure in the file.

    A prepared receptor can carry the same residue number twice, so the key
    counts occurrences as well. Keying on chain and number alone merges those
    residues and undercounts both the residues and their atoms.
    """
    st = list(structure.StructureReader(str(path)))[0]
    out, seen = {}, {}
    for r in st.residue:
        base = "%s:%d" % (r.chain.strip() or "_", r.resnum)
        seen[base] = seen.get(base, 0) + 1
        key = base if seen[base] == 1 else "%s#%d" % (base, seen[base])
        heavy = sum(1 for a in r.atom if a.atomic_number > 1)
        charge = sum(a.formal_charge for a in r.atom)
        out[key] = dict(resname=r.pdbres.strip(), heavy=heavy,
                        atoms=sum(1 for _ in r.atom), formal_charge=charge)
    return st, out


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--a", type=Path, required=True)
    p.add_argument("--b", type=Path, required=True)
    p.add_argument("--label-a", default="A")
    p.add_argument("--label-b", default="B")
    p.add_argument("--out", type=Path)
    a = p.parse_args()

    sa, pa = profile(a.a)
    sb, pb = profile(a.b)
    diffs = []
    def order(k):
        chain, rest = k.split(":", 1)
        num, _, occ = rest.partition("#")
        return (chain, int(num), int(occ or 1))

    for key in sorted(set(pa) | set(pb), key=order):
        x, y = pa.get(key), pb.get(key)
        if x == y:
            continue
        diffs.append(dict(residue=key, a=x, b=y))
    summary = dict(
        label_a=a.label_a, label_b=a.label_b,
        # the two runs write files of the same name, so the run directory is
        # kept as well and the two inputs stay distinguishable in the record
        file_a="/".join(a.a.parts[-2:]), file_b="/".join(a.b.parts[-2:]),
        atoms_a=sa.atom_total, atoms_b=sb.atom_total,
        heavy_atoms_a=sum(1 for t in sa.atom if t.atomic_number > 1),
        heavy_atoms_b=sum(1 for t in sb.atom if t.atomic_number > 1),
        formal_charge_a=sa.formal_charge, formal_charge_b=sb.formal_charge,
        residues_a=len(pa), residues_b=len(pb),
        differing_residues=len(diffs), differences=diffs)
    print(json.dumps(summary, indent=2))
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
