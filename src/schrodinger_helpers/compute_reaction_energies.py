#!/usr/bin/env python3
"""Compute thia-Michael (methanethiolate) conjugate-addition reaction energies from Jaguar outputs.

For each capped warhead model, the reaction energy is

    dE = E(adduct anion) - [E(neutral warhead) + E(CH3S- anion)]

using B3LYP-D3/6-31+G(d,p) PBF(water) solution-phase energies parsed from the Jaguar `.out` files.
Point --dft-dir at the folder holding the Jaguar outputs and supply label:warhead pairs; the adduct
output for each warhead is assumed to be named `adduct_<warhead>`.
"""
import re
import os
import csv
import sys
import argparse

H2KCAL = 627.5094740631


def solution_phase_energy(path):
    """Return the last DFT solution-phase energy (hartree) from a finished Jaguar .out, else None."""
    if not os.path.exists(path):
        return None
    text = open(path, errors="ignore").read()
    if "Total elapsed time" not in text:
        return None  # calculation not finished
    matches = re.findall(r"Solution phase energy:\s*DFT\S*\s*(-?\d+\.\d+)", text)
    return float(matches[-1]) if matches else None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dft-dir", default=os.environ.get("DFT_DIR", "dft"),
                    help="directory containing the Jaguar .out files (default: ./dft or $DFT_DIR)")
    ap.add_argument("--thiolate", default="MeS_anion_water",
                    help="basename of the methanethiolate anion output")
    ap.add_argument("--pairs", nargs="+",
                    default=["OMe:wh_OMe_water", "Cl:wh_Cl_water", "CN:wh_CN_water"],
                    help="label:warhead_basename pairs (adduct assumed named adduct_<warhead>)")
    ap.add_argument("--out", default="reaction_energies.csv")
    args = ap.parse_args()

    energy = lambda base: solution_phase_energy(os.path.join(args.dft_dir, base + ".out"))
    mes = energy(args.thiolate)
    if mes is None:
        sys.exit(f"methanethiolate energy not found in {args.dft_dir!r}")

    rows = []
    for pair in args.pairs:
        label, warhead = pair.split(":")
        e_wh, e_adduct = energy(warhead), energy("adduct_" + warhead)
        if None in (e_wh, e_adduct):
            print(f"[skip] {label}: warhead or adduct output missing/unfinished")
            continue
        dE = (e_adduct - (e_wh + mes)) * H2KCAL
        rows.append({"substituent": label, "dE_CH3S_kcal": round(dE, 1)})
        print(f"{label}: dE(CH3S-) = {dE:.1f} kcal/mol")

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["substituent", "dE_CH3S_kcal"])
        writer.writeheader()
        writer.writerows(rows)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
