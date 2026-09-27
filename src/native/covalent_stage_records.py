"""Build the separated covalent-stage record table and a geometry audit.

Every covalent result is kept as its own named field with its own status, so
that enrichment sampling, Prime REAL_MIN optimisation, graph preservation,
minimisation convergence and post-reaction XP scoring cannot be confused:

  enrichment_*       CovDock enrichment (pre-reaction sampling) stage
  prime_*            Prime REAL_MIN covalent-complex optimisation stage
  postreaction_*     CovDock post-reaction XP scoring of the optimised complex
  apparent_affinity  vendor score = 0.5*(enrichment DockingScore + post-reaction DockingScore)

Structure reading uses the licensed Schrodinger Python; no calculation is run.

  run.exe python3 covalent_stage_records.py --base .. --polish ... --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import csv
import json
import math
import re

from schrodinger import structure
from schrodinger.structutils import measure

FINAL_G = re.compile(r"^\s*final\s+(-?[\d.]+)\s+([\d.Ee+-]+)", re.M)
CONVERGED = re.compile(r"Minimization (not yet converged|converged), dE \(\s*([\d.Ee+-]+)")
ITER = re.compile(r"Completed Iteration\s+(\d+)\s+of\s+(\d+)")
SOLV = re.compile(r"energy param solv (\S+)")
OPLS = re.compile(r"Using OPLS version (\S+)")
TARGET_G = re.compile(r"Desired final RMS gradient:\s*([\d.Ee+-]+)")


def parse_prime_log(path):
    txt = path.read_text(errors="replace")
    grads = FINAL_G.findall(txt)
    conv = CONVERGED.findall(txt)
    iters = ITER.findall(txt)
    solv = SOLV.search(txt)
    opls = OPLS.search(txt)
    tgt = TARGET_G.search(txt)
    return dict(
        prime_solvation_model=solv.group(1) if solv else "",
        prime_force_field=opls.group(1) if opls else "",
        prime_target_rms_gradient=float(tgt.group(1)) if tgt else None,
        prime_final_rms_gradient=float(grads[-1][1]) if grads else None,
        prime_macro_iterations_run=int(iters[-1][0]) if iters else None,
        prime_macro_iterations_requested=int(iters[-1][1]) if iters else None,
        prime_energy_converged=(conv[-1][0] == "converged") if conv else None,
        prime_final_dE_kcal=float(conv[-1][1]) if conv else None,
    )


def geometry_audit(path, sg_index, cb_index):
    st = list(structure.StructureReader(str(path)))[0]
    sg, cb = st.atom[sg_index], st.atom[cb_index]
    bonded_cb = [b.atom2 if b.atom1.index == cb_index else b.atom1 for b in cb.bond]
    heavy_nb = [x for x in bonded_cb if x.atomic_number > 1]
    ca = None
    for x in heavy_nb:
        if x.index != sg_index and x.element == "C":
            ca = x
            break
    out = dict(
        S_C_bond_A=measure.measure_distance(sg, cb),
        Cbeta_bonded_heavy_atoms=len(heavy_nb),
        Cbeta_total_bonded_atoms=len(bonded_cb),
        Cbeta_is_sp3_by_valence=(len(bonded_cb) == 4),
        Ca_Cbeta_S_angle_deg=None,
        Cbeta_S_CB_angle_deg=None,
        Ca_Cbeta_S_CB_dihedral_deg=None,
        Cbeta_improper_deg=None,
    )
    if ca is not None:
        out["Ca_Cbeta_S_angle_deg"] = measure.measure_bond_angle(ca, cb, sg)
        cg = [b.atom2 if b.atom1.index == sg_index else b.atom1 for b in sg.bond]
        cg = [x for x in cg if x.index != cb_index and x.atomic_number > 1]
        if cg:
            out["Cbeta_S_CB_angle_deg"] = measure.measure_bond_angle(cb, sg, cg[0])
            out["Ca_Cbeta_S_CB_dihedral_deg"] = measure.measure_dihedral_angle(ca, cb, sg, cg[0])
    if len(bonded_cb) >= 3:
        a, b, c = bonded_cb[0], bonded_cb[1], bonded_cb[2]
        out["Cbeta_improper_deg"] = measure.measure_dihedral_angle(a, cb, b, c)

    lig_idx = set(at.index for at in st.atom
                  if at.property.get("i_cdock_lig_attach") is not None)
    if not lig_idx:
        lig_idx = set(at.index for at in st.atom if at.pdbres.strip() in ("UNK", "LIG"))
    worst, count = math.inf, 0
    if lig_idx:
        bonded_pairs = set()
        for bd in st.bond:
            bonded_pairs.add((bd.atom1.index, bd.atom2.index))
            bonded_pairs.add((bd.atom2.index, bd.atom1.index))
        prot = [at for at in st.atom if at.index not in lig_idx and at.atomic_number > 1]
        for at in (st.atom[i] for i in sorted(lig_idx)):
            if at.atomic_number <= 1:
                continue
            for pt in prot:
                if (at.index, pt.index) in bonded_pairs:
                    continue
                d = measure.measure_distance(at, pt)
                if d < worst:
                    worst = d
                if d < 2.2:
                    count += 1
    out["min_nonbonded_heavy_contact_A"] = None if worst is math.inf else worst
    out["nonbonded_heavy_contacts_below_2p2A"] = count
    out["ligand_heavy_atoms"] = sum(1 for i in lig_idx if st.atom[i].atomic_number > 1)
    return out


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--base", type=Path, required=True,
                   help="directory holding the covalent docking stage output")
    p.add_argument("--polish", type=Path, required=True,
                   help="directory holding the covalent stage audit records")
    p.add_argument("--recovered-logs", type=Path, default=None,
                   help="directory of re-run Prime logs for complexes whose archived log "
                        "lost its minimisation section")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    with (a.polish / "native_audit/uniform_covdock_verified.csv").open() as fh:
        enrich = {r["target"] + "_" + r["compound"]: r for r in csv.DictReader(fh)}
    scored = {r["target"] + "_" + r["compound"]: r
              for r in json.loads((a.base / "scored_panel/all_rescore_results.json").read_text())}

    rows, audits = [], []
    for key in sorted(enrich):
        e = enrich[key]
        ref_dir = a.base / "refined_panel" / key
        refined = json.loads((ref_dir / "result.json").read_text())
        log_path = ref_dir / "Prime_refinement.log"
        log_source = "archived job log"
        if "STARTING TO WORK ON MINIMIZATION" not in log_path.read_text(errors="replace"):
            if a.recovered_logs is None:
                raise SystemExit("archived log for %s has no minimisation section and no "
                                 "--recovered-logs directory was given" % key)
            recovered = a.recovered_logs / key / "Prime_refinement.log"
            check = json.loads((a.recovered_logs / key / "recheck.json").read_text())
            if not check.get("reproduces_archive"):
                raise SystemExit("recovered log for %s does not reproduce the archive" % key)
            log_path = recovered
            log_source = ("re-run from the archived input; reproduces the archived Prime "
                          "energy and reaction-centre geometry exactly")
        log = parse_prime_log(log_path)
        log["prime_log_source"] = log_source
        s = scored.get(key, {})
        ig = refined.get("input_geometry", [])
        og = refined.get("output_geometry", [])
        sc_in = next((b["length_A"] for b in ig if "S" in b["element1"] + b["element2"]), None)
        sc_out = next((b["length_A"] for b in og if "S" in b["element1"] + b["element2"]), None)
        pre = float(e["score"])
        post = s.get("postreaction_docking_score")
        row = dict(
            target=e["target"], compound=e["compound"],
            reaction="Michael Addition", reactive_residue=e["SG_residue"],
            SG_atom_index=int(e["SG_index"]), Cbeta_atom_index=int(e["Cbeta_index"]),
            enrichment_source=e["source"], enrichment_pose_record=int(e["record"]),
            enrichment_DockingScore=pre, enrichment_GlideScore=float(e["sample_gscore"]),
            enrichment_minimisation_converged=int(e["minimization_converged"]),
            enrichment_rms_derivative=float(e["rms_derivative"]),
            enrichment_S_Cbeta_A=float(e["SG_Cbeta_A"]),
            prime_status="completed" if refined.get("success") else "failed",
            prime_graph_preserved=bool(refined.get("success")),
            prime_selected_atoms=refined.get("selected_atoms"),
            prime_energy_kcal=refined.get("Prime_energy"),
            prime_S_Cbeta_input_A=sc_in, prime_S_Cbeta_output_A=sc_out,
            prime_elapsed_s=refined.get("elapsed_seconds"),
        )
        row.update(log)
        row.update(dict(
            postreaction_status=s.get("status", "not run"),
            postreaction_DockingScore=post,
            postreaction_GlideScore=s.get("postreaction_glidescore"),
            apparent_affinity_score=s.get("apparent_affinity_score"),
            apparent_affinity_formula="0.5*(enrichment_DockingScore + postreaction_DockingScore)",
            apparent_affinity_recomputed=(None if post is None else 0.5 * (pre + post)),
            postreaction_grid=s.get("matching_grid"),
            postreaction_grid_sha256=s.get("grid_sha256"),
            postreaction_input_sha256=s.get("refined_input_sha256"),
            postreaction_precision="XP", postreaction_docking_method="optandscore",
            postreaction_capping_force_field=s.get("hydrogen_capping_forcefield"),
            postreaction_elapsed_s=s.get("elapsed_seconds"),
            covdock_commandline=e["commandline"],
        ))
        rows.append(row)

        audit = dict(target=e["target"], compound=e["compound"])
        audit.update(geometry_audit(ref_dir / "refined.maegz", int(e["SG_index"]), int(e["Cbeta_index"])))
        audits.append(audit)

    with (a.out / "Covalent_stage_records.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with (a.out / "Covalent_geometry_audit.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(audits[0]))
        w.writeheader()
        w.writerows(audits)

    ok = [r for r in rows if r["postreaction_status"] == "completed"]
    mism = [r["target"] + "_" + r["compound"] for r in ok
            if abs(r["apparent_affinity_score"] - r["apparent_affinity_recomputed"]) > 1e-9]
    ang = [x["Ca_Cbeta_S_angle_deg"] for x in audits if x["Ca_Cbeta_S_angle_deg"] is not None]
    contacts = [x["min_nonbonded_heavy_contact_A"] for x in audits
                if x["min_nonbonded_heavy_contact_A"] is not None]
    summary = dict(
        records=len(rows),
        prime_completed=sum(r["prime_status"] == "completed" for r in rows),
        prime_graph_preserved=sum(bool(r["prime_graph_preserved"]) for r in rows),
        prime_energy_converged=sum(bool(r["prime_energy_converged"]) for r in rows),
        prime_final_rms_gradient_max=max(r["prime_final_rms_gradient"] for r in rows
                                         if r["prime_final_rms_gradient"] is not None),
        prime_target_rms_gradient=rows[0]["prime_target_rms_gradient"],
        prime_solvation_models=sorted({r["prime_solvation_model"] for r in rows}),
        prime_force_fields=sorted({r["prime_force_field"] for r in rows}),
        postreaction_completed=len(ok),
        apparent_affinity_formula_mismatches=mism,
        S_C_input_range_A=[min(r["prime_S_Cbeta_input_A"] for r in rows),
                           max(r["prime_S_Cbeta_input_A"] for r in rows)],
        S_C_output_range_A=[min(r["prime_S_Cbeta_output_A"] for r in rows),
                            max(r["prime_S_Cbeta_output_A"] for r in rows)],
        postreaction_score_range=[min(r["postreaction_DockingScore"] for r in ok),
                                  max(r["postreaction_DockingScore"] for r in ok)],
        apparent_affinity_range=[min(r["apparent_affinity_score"] for r in ok),
                                 max(r["apparent_affinity_score"] for r in ok)],
        min_nonbonded_contact_A=min(contacts) if contacts else None,
        complexes_with_contact_below_2p2A=sum(x["nonbonded_heavy_contacts_below_2p2A"] > 0 for x in audits),
        Cbeta_sp3_by_valence=sum(bool(x["Cbeta_is_sp3_by_valence"]) for x in audits),
        Ca_Cbeta_S_angle_range_deg=[min(ang), max(ang)] if ang else None,
    )
    (a.out / "Covalent_stage_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
