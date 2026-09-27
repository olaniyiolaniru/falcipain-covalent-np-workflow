"""Export the complete electronic-structure evidence for the cinnamate series.

Produces, for every neutral/anion record used in the manuscript:

  Quantum_job_manifest.csv   one row per Jaguar job: input and output paths and
                             SHA-256, charge, multiplicity, functional, basis,
                             dispersion, solvation model and its dielectric and
                             probe radius, whether the geometry was optimised,
                             normal completion, final geometry gradient, ESP fit
                             RMS error, final <S**2>, orbital energies, total
                             energies and elapsed time.
  Quantum_atom_charges.csv   every atom of every job: label, element, Cartesian
                             coordinates, ESP and Mulliken charges, Mulliken
                             spin population and fragment assignment.
  Quantum_geometry_checks.csv alkene configuration and the neutral/anion
                             coordinate identity check for each pair.

Reading uses the licensed Schrodinger Python for the Jaguar .01.mae restart
files; no calculation is started.

  run.exe python3 quantum_evidence.py --jobs <schrodinger-jobs> \
      --extra ../jobs/fukui_clean_ext --series <Cinnamate_series.csv> \
      --check <Charge_scheme_check.csv> --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import csv
import hashlib
import json
import math
import re

from schrodinger import structure

HARTREE_EV = 27.211386245988

CORE = {"C1", "O2", "C3", "O4", "C5", "C6"}
PARA = {"OMe": {"O11", "C12"}, "OH": {"O11"}, "H": set(),
        "Cl": {"Cl11"}, "CN": {"C11", "N12"}, "NO2": {"N11", "O12", "O13"}}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def last(pattern, text, cast=float):
    hits = re.findall(pattern, text, re.I)
    return cast(hits[-1]) if hits else None


def region_map(atoms, sub):
    """Assign each atom to ester_alkene, aryl or para_substituent."""
    para = PARA[sub]
    region = {}
    heavy = [a for a in atoms if a["element"] != "H"]
    for a in heavy:
        lab = a["label"]
        region[lab] = ("ester_alkene" if lab in CORE
                       else "para_substituent" if lab in para else "aryl")
    for a in atoms:
        if a["element"] != "H":
            continue
        best, bd = None, 1e9
        for h in heavy:
            d = math.dist((a["x"], a["y"], a["z"]), (h["x"], h["y"], h["z"]))
            if d < bd:
                best, bd = h["label"], d
        assert bd < 1.35, (a["label"], bd)
        region[a["label"]] = region[best]
    return region


def read_job(mae, jobs_root, label_root):
    mae = Path(mae)
    base = mae.name.replace(".01.mae", "")
    inp, out = mae.with_name(base + ".in"), mae.with_name(base + ".out")
    it = inp.read_text(errors="replace") if inp.is_file() else ""
    ot = out.read_text(errors="replace") if out.is_file() else ""
    st = list(structure.StructureReader(str(mae)))[-1]
    rel = str(mae.relative_to(label_root)).replace("\\", "/")

    def key(name, default=""):
        m = re.search(r"\b" + name + r"\s*=\s*([^\s&]+)", it, re.I)
        return m.group(1) if m else default

    homo, lumo = st.property.get("r_j_HOMO"), st.property.get("r_j_LUMO")
    row = dict(
        job=rel, mae_sha256=sha(mae),
        input_file=inp.name, input_sha256=sha(inp) if inp.is_file() else "",
        output_file=out.name, output_sha256=sha(out) if out.is_file() else "",
        atoms=st.atom_total,
        molecular_charge=key("molchg"), multiplicity=key("multip"),
        functional=key("dftname"), basis=key("basis"),
        dispersion=("D3 (Grimme) via dftname" if "d3" in key("dftname").lower() else "none"),
        geometry_optimised=(key("igeopt") == "1"),
        solvation_keyword=key("isolv"), solvent=key("solvent"),
        solvation_model=("PBF Poisson-Boltzmann continuum"
                         if "Solvation energy will be computed using PBF" in ot else ""),
        continuum_dielectric=last(r"Continuum dielectric constant:\s*([\d.]+)", ot),
        solvent_probe_radius_A=last(r"Solvent probe molecule radius:\s*([\d.]+)", ot),
        normal_completion=bool(re.search(r"Job \S+ completed on", ot)),
        geometry_optimisation_complete=("Geometry optimization complete" in ot),
        solvation_converged=("stopping: solvation energy converged" in ot),
        final_gradient_max_au=last(r"gradient maximum:\s*([\d.Ee+-]+)", ot),
        final_gradient_rms_au=last(r"gradient rms:\s*([\d.Ee+-]+)", ot),
        esp_fit_rms_error_au=last(r"RMS Error\s*([\d.Ee+-]+) hartrees", ot),
        final_S2=last(r"<S\*\*2>\s*\.\.\.\s*([\d.]+)", ot),
        gas_phase_energy_au=st.property.get("r_j_Gas_Phase_Energy"),
        solution_phase_energy_au=st.property.get("r_j_Solution_Phase_Energy"),
        final_energy_au=st.property.get("r_j_Final_Energy"),
        solvation_energy_kcal=st.property.get("r_j_Solvation_Energy_(kcal/mol)"),
        HOMO_au=homo, LUMO_au=lumo,
        HOMO_eV=None if homo is None else homo * HARTREE_EV,
        LUMO_eV=None if lumo is None else lumo * HARTREE_EV,
        elapsed_seconds=last(r"Total elapsed time:\s*([\d.]+)", ot),
    )
    atoms = [dict(job=rel, index=a.index,
                  label=a.property.get("s_m_atom_name", "").strip(),
                  element=a.element, x=a.x, y=a.y, z=a.z,
                  ESP_charge=a.property.get("r_j_ESP_Charges"),
                  Mulliken_charge=a.property.get("r_j_Mulliken_Charges"),
                  Mulliken_spin=a.property.get("r_j_Mulliken_Spin_Populations"))
             for a in st.atom]
    return row, atoms, st


def alkene_check(st):
    """C5=C6 configuration of the cinnamoyl unit from the stored geometry."""
    by = {a.property.get("s_m_atom_name", "").strip(): a for a in st.atom}
    need = ("C3", "C5", "C6", "C7")
    if not all(k in by for k in need):
        return None, None
    import numpy as np
    p = [np.array([by[k].x, by[k].y, by[k].z]) for k in need]
    b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    ang = math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))
    return ang, ("E" if abs(ang) > 90 else "Z")


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--jobs", type=Path, required=True)
    p.add_argument("--extra", type=Path, default=None)
    p.add_argument("--series", type=Path, required=True)
    p.add_argument("--check", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    with a.series.open() as fh:
        series = list(csv.DictReader(fh))
    with a.check.open() as fh:
        check = list(csv.DictReader(fh))

    pairs = []
    for r in series:
        pairs.append(("primary", r["substituent"], a.jobs / r["neutral_source"],
                      a.jobs / r["anion_source"]))
    for r in check:
        pairs.append(("fixed_geometry_check", r["substituent"], a.jobs / r["neutral_source"],
                      a.jobs / r["anion_source"]))
    if a.extra:
        for sub in ("OH", "H"):
            n = a.extra / ("fc_%s_N.01.mae" % sub)
            an = a.extra / ("fc_%s_A.01.mae" % sub)
            if n.is_file() and an.is_file():
                pairs.append(("fixed_geometry_check", sub, n, an))

    manifest, charges, geoms, frag_rows, pair_rows = [], [], [], [], []
    seen = set()
    for dataset, sub, npath, apath in pairs:
        recs = {}
        for role, path in (("neutral", npath), ("anion", apath)):
            root = a.jobs if str(path).startswith(str(a.jobs)) else path.parent
            row, atoms, st = read_job(path, a.jobs, root)
            row.update(dataset=dataset, substituent=sub, role=role)
            if row["job"] not in seen:
                seen.add(row["job"])
                manifest.append(row)
                charges.extend(atoms)
            recs[role] = (row, atoms, st)
            ang, cfg = alkene_check(st)
            geoms.append(dict(dataset=dataset, substituent=sub, role=role, job=row["job"],
                              C3_C5_C6_C7_dihedral_deg=None if ang is None else round(ang, 2),
                              alkene_configuration=cfg,
                              geometry_optimised=row["geometry_optimised"]))

        n_at = {x["label"]: x for x in recs["neutral"][1]}
        a_at = {x["label"]: x for x in recs["anion"][1]}
        assert set(n_at) == set(a_at)
        max_shift = max(math.dist((n_at[k]["x"], n_at[k]["y"], n_at[k]["z"]),
                                  (a_at[k]["x"], a_at[k]["y"], a_at[k]["z"])) for k in n_at)
        region = region_map(recs["neutral"][1], sub)
        for scheme, col in (("ESP", "ESP_charge"), ("Mulliken", "Mulliken_charge")):
            if any(n_at[k][col] is None or a_at[k][col] is None for k in n_at):
                continue
            f = {k: n_at[k][col] - a_at[k][col] for k in n_at}
            pair_rows.append(dict(
                dataset=dataset, substituent=sub, scheme=scheme,
                neutral_job=recs["neutral"][0]["job"], anion_job=recs["anion"][0]["job"],
                qN_sum=sum(n_at[k][col] for k in n_at),
                qA_sum=sum(a_at[k][col] for k in a_at),
                fplus_sum=sum(f.values()), fplus_beta_C6=f.get("C6"),
                fplus_alpha_C5=f.get("C5"),
                max_atom=max(f, key=f.get), max_atom_fplus=max(f.values()),
                min_atom=min(f, key=f.get), min_atom_fplus=min(f.values()),
                negative_atom_count=sum(1 for v in f.values() if v < 0),
                max_neutral_anion_coordinate_shift_A=max_shift))
            for reg in ("ester_alkene", "aryl", "para_substituent"):
                labs = [k for k in f if region[k] == reg]
                frag_rows.append(dict(dataset=dataset, substituent=sub, scheme=scheme,
                                      region=reg, fplus_sum=sum(f[k] for k in labs),
                                      atom_labels=";".join(sorted(labs))))

    def write(name, rows):
        with (a.out / name).open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    write("Quantum_job_manifest.csv", manifest)
    write("Quantum_atom_charges.csv", charges)
    write("Quantum_geometry_checks.csv", geoms)
    write("Quantum_pair_response.csv", pair_rows)
    write("Quantum_fragment_response.csv", frag_rows)

    s2 = [r["final_S2"] for r in manifest if r["final_S2"] is not None]
    summary = dict(
        jobs=len(manifest),
        pairs=len(pairs),
        normal_completion=sum(bool(r["normal_completion"]) for r in manifest),
        solvation_converged=sum(bool(r["solvation_converged"]) for r in manifest),
        anion_S2_range=[min(s2), max(s2)] if s2 else None,
        dielectrics=sorted({r["continuum_dielectric"] for r in manifest}),
        probe_radii_A=sorted({r["solvent_probe_radius_A"] for r in manifest}),
        alkene_configurations=sorted({g["alkene_configuration"] for g in geoms
                                      if g["alkene_configuration"]}),
        max_neutral_anion_coordinate_shift_A=max(r["max_neutral_anion_coordinate_shift_A"]
                                                 for r in pair_rows),
        charge_normalisation_max_error=max(
            max(abs(r["qN_sum"]), abs(r["qA_sum"] + 1), abs(r["fplus_sum"] - 1))
            for r in pair_rows),
        schemes_by_dataset={d: sorted({r["scheme"] for r in pair_rows if r["dataset"] == d})
                            for d in {r["dataset"] for r in pair_rows}},
        substituents_by_dataset={d: sorted({r["substituent"] for r in pair_rows if r["dataset"] == d})
                                 for d in {r["dataset"] for r in pair_rows}},
    )
    (a.out / "Quantum_evidence_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
