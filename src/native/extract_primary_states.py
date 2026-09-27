"""Export every returned primary ligand state from the standardized-receptor screen.

Both falcipain receptors are prepared under one catalytic-dyad convention: a
neutral cysteine thiol with its S-gamma hydrogen and a neutral catalytic
histidine protonated at N-epsilon-2. The FP-2 pose viewer is the standardized
rerun of all prepared states; the FP-3 pose viewer is the primary screen, whose
receptor already carried that convention.

Requires the licensed Schrodinger installation only for structure reading; no
calculation is started. Produces the complete state-level record requested for
the primary falcipain screen: one row per returned pose record, plus a separate
inventory of prepared states that returned no pose.

  run.exe python3 extract_primary_states.py --jobs <schrodinger-jobs> --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import csv, hashlib, json

from schrodinger import structure
from schrodinger.structutils import analyze

SOURCES = [
    ("FP2", "FP2_standard_xp_pv.maegz", "primary_FP2_standardised"),
    ("FP3", "XP_FP3_apo_185_pv.maegz", "primary_FP3_standardised"),
]
PREPARED = "lp185.maegz"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def smiles(st):
    try:
        return analyze.generate_smiles(st, unique=True, stereo=analyze.STEREO_FROM_GEOMETRY)
    except Exception:
        try:
            return analyze.generate_smiles(st)
        except Exception:
            return ""


def chiral_tags(st):
    """Readable stereo descriptor list, e.g. 'C12:R;C15:S'."""
    out = []
    try:
        for idx, code in analyze.get_chiral_atoms(st).items():
            out.append(f"{st.atom[idx].element}{idx}:{code}")
    except Exception:
        pass
    return ";".join(out)


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--jobs", type=Path, required=True)
    p.add_argument("--fp2-dir", type=Path, default=None,
                   help="directory holding the standardized FP-2 pose viewer")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    prepared = {}
    prep_path = a.jobs / PREPARED
    for n, st in enumerate(structure.StructureReader(str(prep_path)), 1):
        variant = st.property.get("s_lp_Variant", st.title)
        prepared[variant] = dict(
            prepared_state_id=variant,
            parent_id=st.title.split("_")[0],
            prepared_record=n,
            formal_charge=st.formal_charge,
            atoms=st.atom_total,
            heavy_atoms=sum(at.atomic_number > 1 for at in st.atom),
            stereo_descriptors=chiral_tags(st),
            isomeric_smiles=smiles(st),
            epik_state_penalty=st.property.get("r_epik_State_Penalty"),
            force_field=st.property.get("s_lp_Force_Field"),
        )

    rows, seen = [], {t: set() for t, _, _ in SOURCES}
    for target, fname, protocol in SOURCES:
        path = (a.fp2_dir / fname) if (a.fp2_dir and (a.fp2_dir / fname).is_file()) else (a.jobs / fname)
        digest = sha(path)
        for n, st in enumerate(structure.StructureReader(str(path)), 1):
            if n == 1:
                continue  # record 1 is the receptor in a Glide pose viewer
            variant = st.property.get("s_lp_Variant", "")
            meta = prepared.get(variant, {})
            parent = meta.get("parent_id") or st.title.split("_")[0]
            dock = st.property.get("r_i_docking_score")
            glide = st.property.get("r_i_glide_gscore")
            seen[target].add(variant)
            rows.append(dict(
                target=target, protocol=protocol, parent_id=parent,
                prepared_state_id=variant,
                prepared_record=meta.get("prepared_record"),
                formal_charge=meta.get("formal_charge"),
                stereo_descriptors=meta.get("stereo_descriptors", ""),
                isomeric_smiles=meta.get("isomeric_smiles", ""),
                heavy_atoms=meta.get("heavy_atoms"),
                pose_record_index=n,
                DockingScore=dock,
                GlideScore=glide,
                effective_state_penalty=(None if dock is None or glide is None else dock - glide),
                epik_state_penalty=meta.get("epik_state_penalty"),
                glide_emodel=st.property.get("r_i_glide_emodel"),
                glide_energy=st.property.get("r_i_glide_energy"),
                source_job=fname, source_sha256=digest,
            ))

    with (a.out / "Primary_state_records.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    missing = []
    for target, fname, _ in SOURCES:
        for variant, meta in prepared.items():
            if variant not in seen[target]:
                missing.append(dict(target=target, parent_id=meta["parent_id"],
                                    prepared_state_id=variant,
                                    prepared_record=meta["prepared_record"],
                                    formal_charge=meta["formal_charge"],
                                    epik_state_penalty=meta["epik_state_penalty"],
                                    source_job=fname,
                                    outcome="no pose returned",
                                    reason="state produced no pose that passed Glide XP scoring/rejection filters"))
    with (a.out / "Primary_state_failures.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(missing[0]))
        w.writeheader()
        w.writerows(missing)

    summary = {
        "prepared_states": len(prepared),
        "prepared_parents": len({m["parent_id"] for m in prepared.values()}),
        "returned_state_records": {t: sum(r["target"] == t for r in rows) for t, _, _ in SOURCES},
        "parents_with_pose": {t: len({r["parent_id"] for r in rows if r["target"] == t}) for t, _, _ in SOURCES},
        "states_without_pose": {t: sum(m["target"] == t for m in missing) for t, _, _ in SOURCES},
        "source_sha256": dict(
            [(f, sha((a.fp2_dir / f) if (a.fp2_dir and (a.fp2_dir / f).is_file()) else (a.jobs / f)))
             for _, f, _ in SOURCES] + [(PREPARED, sha(prep_path))]),
        "catalytic_dyad_convention": ("neutral cysteine thiol (S-gamma H present) and neutral "
                                      "catalytic histidine (HIE, N-epsilon-2 H) at both targets"),
    }
    (a.out / "Primary_state_extraction.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
