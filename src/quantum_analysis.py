"""Recompute every reported electronic descriptor from the supplied charge table.

Portable: pandas/numpy only. Consumes Quantum_job_manifest.csv and
Quantum_atom_charges.csv and reproduces, without any quantum calculation:

  * global electrophilicity omega = mu^2/(2 eta) with mu = (eHOMO + eLUMO)/2 and
    the full-gap hardness eta = eLUMO - eHOMO, in eV;
  * the condensed response to electron addition f+_k = q_k(N) - q_k(N+1) at the
    beta carbon C6 and at every other atom;
  * the derived local philicity omega * f+_beta;
  * the fragment partition of the unit response over the ester-alkene unit, the
    aryl ring and the para substituent, with hydrogens assigned to the fragment
    of their nearest bonded heavy atom;
  * the same-geometry ESP/Mulliken comparison between the primary optimisation
    jobs and the fixed-geometry single-point check jobs.

  python quantum_analysis.py --data ../exports --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import json
import numpy as np
import pandas as pd

HARTREE_EV = 27.211386245988
CORE = ["C1", "O2", "C3", "O4", "C5", "C6"]
PARA = {"OMe": ["O11", "C12"], "OH": ["O11"], "H": [],
        "Cl": ["Cl11"], "CN": ["C11", "N12"], "NO2": ["N11", "O12", "O13"]}
BETA = "C6"


def regions(atoms, sub):
    heavy = atoms[atoms.element.ne("H")]
    reg = {}
    for r in heavy.itertuples():
        reg[r.label] = ("ester_alkene" if r.label in CORE
                        else "para_substituent" if r.label in PARA[sub] else "aryl")
    hxyz = heavy[["x", "y", "z"]].to_numpy()
    hlab = heavy.label.tolist()
    for r in atoms[atoms.element.eq("H")].itertuples():
        d = np.linalg.norm(hxyz - np.array([r.x, r.y, r.z]), axis=1)
        j = int(np.argmin(d))
        if d[j] >= 1.35:
            raise ValueError("hydrogen %s is %.3f A from every heavy atom" % (r.label, d[j]))
        reg[r.label] = reg[hlab[j]]
    return reg


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--data", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(a.data / "Quantum_job_manifest.csv")
    charges = pd.read_csv(a.data / "Quantum_atom_charges.csv")
    meta = manifest.set_index("job")

    pairs = (manifest[["dataset", "substituent", "role", "job"]]
             .pivot_table(index=["dataset", "substituent"], columns="role",
                          values="job", aggfunc="first"))

    desc, frag, atom = [], [], []
    for (dataset, sub), row in pairs.iterrows():
        n = charges[charges.job.eq(row.neutral)].set_index("label", drop=False)
        an = charges[charges.job.eq(row.anion)].set_index("label", drop=False)
        if not n.index.equals(an.index):
            raise ValueError("atom labels differ between %s and %s" % (row.neutral, row.anion))
        reg = regions(n, sub)
        m = meta.loc[row.neutral]
        homo, lumo = m.HOMO_eV, m.LUMO_eV
        eta = lumo - homo
        mu = (homo + lumo) / 2
        omega = mu * mu / (2 * eta)
        for scheme, col in (("ESP", "ESP_charge"), ("Mulliken", "Mulliken_charge")):
            if n[col].isna().any() or an[col].isna().any():
                continue
            f = n[col] - an[col]
            desc.append(dict(dataset=dataset, substituent=sub, scheme=scheme,
                             HOMO_eV=homo, LUMO_eV=lumo,
                             hardness_fullgap_eV=eta, chemical_potential_eV=mu,
                             omega_eV=omega,
                             fplus_beta=float(f.loc[BETA]),
                             local_philicity_eV=float(omega * f.loc[BETA]),
                             fplus_alpha_C5=float(f.loc["C5"]),
                             qN_sum=float(n[col].sum()), qA_sum=float(an[col].sum()),
                             fplus_sum=float(f.sum()),
                             neutral_job=row.neutral, anion_job=row.anion))
            for reg_name in ("ester_alkene", "aryl", "para_substituent"):
                labs = [k for k in f.index if reg[k] == reg_name]
                frag.append(dict(dataset=dataset, substituent=sub, scheme=scheme,
                                 region=reg_name, fplus_sum=float(f.loc[labs].sum()),
                                 atom_labels=";".join(labs)))
            for lab, v in f.items():
                atom.append(dict(dataset=dataset, substituent=sub, scheme=scheme,
                                 atom=lab, element=n.loc[lab, "element"],
                                 region=reg[lab], fplus=float(v),
                                 qN=float(n.loc[lab, col]), qA=float(an.loc[lab, col])))

    d = pd.DataFrame(desc)
    d.to_csv(a.out / "Series_electronic_descriptors.csv", index=False)
    pd.DataFrame(frag).to_csv(a.out / "Series_fragment_response.csv", index=False)
    pd.DataFrame(atom).to_csv(a.out / "Series_atom_response.csv", index=False)

    esp = d[d.scheme.eq("ESP")]
    prim = esp[esp.dataset.eq("primary")].set_index("substituent")
    chk = esp[esp.dataset.eq("fixed_geometry_check")].set_index("substituent")
    shared = sorted(set(prim.index) & set(chk.index))
    cmp_rows = []
    for s in shared:
        mull = d[(d.dataset.eq("fixed_geometry_check")) & (d.substituent.eq(s))
                 & (d.scheme.eq("Mulliken"))]
        cmp_rows.append(dict(
            substituent=s,
            fplus_ESP_primary=prim.loc[s, "fplus_beta"],
            fplus_ESP_check=chk.loc[s, "fplus_beta"],
            fplus_Mulliken_check=float(mull.fplus_beta.iloc[0]) if len(mull) else None,
            same_scheme_difference=chk.loc[s, "fplus_beta"] - prim.loc[s, "fplus_beta"],
            same_scheme_relative_percent=100 * (chk.loc[s, "fplus_beta"] - prim.loc[s, "fplus_beta"])
            / prim.loc[s, "fplus_beta"]))
    cmp = pd.DataFrame(cmp_rows)
    cmp.to_csv(a.out / "Series_charge_scheme_comparison.csv", index=False)

    nitro = {r.scheme + "_" + r.region: r.fplus_sum
             for r in pd.DataFrame(frag).itertuples()
             if r.substituent == "NO2" and r.dataset == "fixed_geometry_check"}
    nitro_esp = {r.region: r.fplus_sum for r in pd.DataFrame(frag).itertuples()
                 if r.substituent == "NO2" and r.dataset == "primary" and r.scheme == "ESP"}
    summary = dict(
        descriptor_records=len(d),
        substituents_primary=sorted(prim.index),
        substituents_checked=sorted(chk.index),
        omega_range_eV=[float(prim.omega_eV.min()), float(prim.omega_eV.max())],
        omega_NO2_over_OH=float(prim.loc["NO2", "omega_eV"] / prim.loc["OH", "omega_eV"]),
        fplus_beta_NO2_over_OH=float(prim.loc["NO2", "fplus_beta"] / prim.loc["OH", "fplus_beta"]),
        local_philicity_eV={s: float(prim.loc[s, "local_philicity_eV"]) for s in prim.index},
        charge_normalisation_max_error=float(
            max(d.qN_sum.abs().max(), (d.qA_sum + 1).abs().max(), (d.fplus_sum - 1).abs().max())),
        nitro_fragment_response=dict(primary_ESP=nitro_esp, check=nitro),
        max_same_scheme_difference=float(cmp.same_scheme_difference.abs().max()),
        hardness_convention="full gap, eta = eLUMO - eHOMO",
    )
    (a.out / "Quantum_analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
