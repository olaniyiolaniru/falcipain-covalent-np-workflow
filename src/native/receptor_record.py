"""Export the receptor preparation and grid record for every docked stage.

Every docking stage in the study is listed with the receptor it used, the grid
box it used, and the assigned protonation state of the catalytic dyad, so that
the recognition screen, the matched-analogue comparison, the catalytic-dyad
sensitivity experiment, the retrospective benchmark and the pH panel can each
be traced to one file.

The catalytic dyad is standardized across the two falcipain receptors used for
the reported screen and the reported matched comparison: the catalytic cysteine
carries a neutral thiol with an explicit S-gamma hydrogen, and the catalytic
histidine is the neutral HIE tautomer. The alternative assignment, in which the
cysteine is a thiolate and the histidine is the doubly protonated HIP form, is
retained as a sensitivity experiment on the identical grid box.

  run.exe python3 receptor_record.py --jobs <schrodinger-jobs> --work .. --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import csv
import hashlib
import json
import re
import zipfile

from schrodinger import structure

PROTEIN = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER "
              "THR TRP TYR VAL HID HIE HIP CYX CYT ASH GLH LYN ACE NMA NME".split())
WATER = {"HOH", "SPC", "T3P", "WAT"}

PREP = ("Protein Preparation Wizard: bond orders and hydrogens assigned, side chains "
        "completed, PROPKA hydrogen-bond assignment at pH 7.0, restrained minimisation "
        "to 0.30 A heavy-atom RMSD, OPLS_2005")


def sha(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def keyword(text, key):
    m = re.search(r"^\s*" + key + r"\s+(.+)$", text, re.M | re.I)
    return m.group(1).strip().strip('"') if m else ""


def grid_receptor(grid_zip, scratch):
    """Extract the receptor structure a Glide grid archive was built from."""
    with zipfile.ZipFile(str(grid_zip)) as z:
        names = [n for n in z.namelist() if n.endswith("_recep.mae")]
        if not names:
            return None
        scratch.mkdir(parents=True, exist_ok=True)
        return Path(z.extract(names[0], str(scratch)))


def describe(path, cys_num, his_num):
    st = next(structure.StructureReader(str(path)))
    het, waters = {}, 0
    for res in st.residue:
        name = res.pdbres.strip()
        if name in WATER:
            waters += 1
        elif name not in PROTEIN:
            het[name] = het.get(name, 0) + 1
    cys_state, his_state = "", ""
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
    dyad = ""
    if cys_state and his_state:
        if cys_state.startswith("neutral") and his_state == "HIE":
            dyad = "standardized neutral dyad"
        elif his_state == "HIP":
            dyad = "ion-pair dyad"
        else:
            dyad = "other"
    return dict(chains=";".join(sorted(set((a.chain.strip() or "?") for a in st.atom))),
                residues=sum(1 for _ in st.residue),
                atoms=st.atom_total,
                heavy_atoms=sum(1 for a in st.atom if a.atomic_number > 1),
                formal_charge=st.formal_charge,
                retained_waters=waters,
                non_standard_residues=";".join(k + " x" + str(v)
                                               for k, v in sorted(het.items())) or "none",
                catalytic_cys="Cys" + str(cys_num) + ": " + (cys_state or "not found"),
                catalytic_his="His" + str(his_num) + ": " + (his_state or "not found"),
                catalytic_dyad_convention=dyad)


def stages_for(jobs, work):
    return [
        ("primary recognition screen", "FP-2", "3BPF", 42, 174,
         work / "fp2_standard_screen/FP2_standard_xp.in",
         jobs / "FP3_grid_CYS42/FP3_grid_CYS42.in",
         work / "fp2_standard_screen/FP2_standard_grid.zip"),
        ("primary recognition screen", "FP-3", "3BPM", 51, 183,
         jobs / "XP_FP3_apo_185.in", jobs / "FP3_apo_grid.in", jobs / "FP3_apo_grid.zip"),
        ("matched-analogue comparison", "FP-2", "3BPF", 42, 174,
         None, jobs / "FP3_grid_CYS42/FP3_grid_CYS42.in",
         work / "fp2_standard_screen/FP2_standard_grid.zip"),
        ("matched-analogue comparison", "FP-3", "3BPM", 51, 183,
         None, jobs / "FP3_apo_grid.in", jobs / "FP3_apo_grid.zip"),
        ("catalytic-dyad sensitivity experiment", "FP-2", "3BPF", 42, 174,
         work / "matched_dyad_alt/FP2_alternative_xp.in",
         jobs / "FP3_grid_CYS42/FP3_grid_CYS42.in",
         jobs / "FP3_grid_CYS42/FP3_grid_CYS42.zip"),
        ("retrospective enrichment benchmark", "FP-2", "3BPF", 42, 174,
         work / "bench_fp2/sp.in", jobs / "FP3_grid_CYS42/FP3_grid_CYS42.in",
         work / "bench_fp2/FP2_standard_grid.zip"),
        ("pH 5.5 sensitivity panel", "FP-2", "3BPF", 42, 174,
         work / "ph_panel/FP2_pH55_xp.in", work / "ph_panel/FP2_pH55_grid.in",
         work / "ph_panel/FP2_pH55_grid.zip"),
        ("pH 5.5 sensitivity panel", "FP-3", "3BPM", 51, 183,
         work / "ph_panel/FP3_pH55_xp.in", work / "ph_panel/FP3_pH55_grid.in",
         work / "ph_panel/FP3_pH55_grid.zip"),
    ]


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--jobs", type=Path, required=True)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    scratch = a.out / "_grid_receptors"

    rows = []
    for stage, name, pdb, cys_num, his_num, job_in, grid_in, grid_zip in stages_for(
            a.jobs, a.work / "jobs"):
        row = dict(stage=stage, target=name, pdb=pdb,
                   preparation=PREP.replace("pH 7.0", "pH 5.5") if "5.5" in stage else PREP)
        if job_in is not None and Path(job_in).is_file():
            text = Path(job_in).read_text(errors="replace")
            row["docking_job_input"] = Path(job_in).name
            row["docking_job_sha256"] = sha(job_in)
            row["docking_precision"] = keyword(text, "PRECISION")
        if Path(grid_in).is_file():
            text = Path(grid_in).read_text(errors="replace")
            row["grid_centre_A"] = keyword(text, "GRID_CENTER")
            row["inner_box_A"] = keyword(text, "INNERBOX")
            row["outer_box_A"] = keyword(text, "OUTERBOX")
            row["grid_force_field"] = keyword(text, "FORCEFIELD")
        if Path(grid_zip).is_file():
            row["grid_archive"] = Path(grid_zip).name
            row["grid_archive_sha256"] = sha(grid_zip)
            tag = (name + "_" + stage).replace(" ", "_")
            rec = grid_receptor(grid_zip, scratch / tag)
            if rec is not None:
                row["receptor_from_grid"] = rec.name
                row["receptor_sha256"] = sha(rec)
                row.update(describe(rec, cys_num, his_num))
        rows.append(row)

    fields = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with (a.out / "Receptor_preparation_record.csv").open("w", newline="",
                                                          encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    (a.out / "Receptor_preparation_record.json").write_text(json.dumps(rows, indent=2),
                                                            encoding="utf-8")
    for r in rows:
        print("%-38s %-5s %-30s %-28s %s" % (
            r["stage"], r["target"], r.get("grid_centre_A", "-"),
            r.get("catalytic_dyad_convention", "-"),
            r.get("catalytic_cys", "") + " / " + r.get("catalytic_his", "")))


if __name__ == "__main__":
    main()
