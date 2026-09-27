"""Measure the catalytic-dyad effect across the whole screened library.

The matched-pair experiment measures the dyad term on one pair of compounds.
This script measures it on every parent that the screen docked, using the two
FP-2 screens that differ in exactly one variable: the same prepared ligand file,
the same grid centre and box, the same force field, precision and selection
rule, and two receptors whose heavy-atom coordinates are identical and whose
catalytic dyad differs by the position of a single proton.

Reported: the per-parent score change, its distribution, the rank agreement
between the two receptors, how many parents change their selected prepared
state, the overlap of the leading sets, and the share of parents whose score
moves by more than a stated reference effect.

  run.exe python3 dyad_library_effect.py \
      --standard ../jobs/fp2_standard_screen/FP2_standard_xp_pv.maegz \
      --alternative <schrodinger-jobs>/XP_FP2_185_pv.maegz \
      --reference-effect 0.766 --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import csv
import hashlib
import json

import numpy as np
from scipy.stats import spearmanr
from schrodinger import structure


def bootstrap_rho(x, y, replicates=2000, seed=1):
    """Percentile bootstrap interval for a Spearman coefficient.

    The coefficient is reported on a modest number of compounds, so the interval
    says how much of the ordering is pinned down by the data rather than leaving
    the reader to guess.
    """
    import numpy as np
    rng = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    if n < 4:
        return [None, None]
    vals = []
    for _ in range(replicates):
        idx = rng.integers(0, n, n)
        if len(set(x[idx])) < 2 or len(set(y[idx])) < 2:
            continue
        r = spearmanr(x[idx], y[idx])
        vals.append(float(getattr(r, "statistic", r[0])))
    if not vals:
        return [None, None]
    vals.sort()
    return [round(vals[int(0.025 * len(vals))], 4),
            round(vals[int(0.975 * len(vals))], 4)]


def sha(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def best_per_parent(path):
    """Minimum DockingScore per parent, with the prepared state that gave it."""
    out = {}
    for n, st in enumerate(structure.StructureReader(str(path)), 1):
        if n == 1 and st.atom_total > 1000:
            continue
        score = st.property.get("r_i_docking_score")
        if score is None:
            continue
        variant = st.property.get("s_lp_Variant", st.title)
        # the prepared-state titles carry the source file name of the parent
        parent = st.title.split("-")[0].replace("_minRM1.pdb", "")
        rec = dict(parent=parent, variant=variant, score=float(score),
                   glide=float(st.property.get("r_i_glide_gscore", float("nan"))),
                   record=n)
        if parent not in out or rec["score"] < out[parent]["score"]:
            out[parent] = rec
    return out


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--standard", type=Path, required=True)
    p.add_argument("--alternative", type=Path, required=True)
    p.add_argument("--reference-effect", type=float, default=0.766,
                   help="substituent effect the dyad term is compared against")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    std = best_per_parent(a.standard)
    alt = best_per_parent(a.alternative)
    shared = sorted(set(std) & set(alt))
    if not shared:
        raise SystemExit("the two screens share no parent")

    rows = []
    for parent in shared:
        s, t = std[parent], alt[parent]
        rows.append(dict(compound=parent,
                         standardized_DockingScore=round(s["score"], 4),
                         alternative_DockingScore=round(t["score"], 4),
                         change=round(t["score"] - s["score"], 4),
                         standardized_state=s["variant"],
                         alternative_state=t["variant"],
                         state_changes=bool(s["variant"] != t["variant"])))

    delta = np.array([r["change"] for r in rows])
    x = np.array([std[c]["score"] for c in shared])
    y = np.array([alt[c]["score"] for c in shared])
    rho = spearmanr(x, y)
    by_std = [c for c in sorted(shared, key=lambda c: std[c]["score"])]
    by_alt = [c for c in sorted(shared, key=lambda c: alt[c]["score"])]

    summary = dict(
        experiment=("the same prepared ligand file docked into two FP-2 receptors whose "
                    "heavy-atom coordinates are identical and whose catalytic dyad differs "
                    "by the position of one proton, at the same grid centre and box"),
        standardized_file=a.standard.name, standardized_sha256=sha(a.standard),
        alternative_file=a.alternative.name, alternative_sha256=sha(a.alternative),
        parents_compared=len(shared),
        parents_standardized_only=sorted(set(std) - set(alt)),
        parents_alternative_only=sorted(set(alt) - set(std)),
        median_absolute_change=round(float(np.median(np.abs(delta))), 4),
        mean_change=round(float(delta.mean()), 4),
        iqr_absolute_change=[round(float(np.percentile(np.abs(delta), 25)), 4),
                             round(float(np.percentile(np.abs(delta), 75)), 4)],
        max_absolute_change=round(float(np.abs(delta).max()), 4),
        largest_mover=rows[int(np.abs(delta).argmax())]["compound"],
        rank_agreement_spearman=round(float(getattr(rho, "statistic", rho[0])), 4),
        rank_agreement_ci=bootstrap_rho(x, y),
        parents_changing_selected_state=int(sum(r["state_changes"] for r in rows)),
        reference_effect=a.reference_effect,
        parents_moving_more_than_reference=int((np.abs(delta) > a.reference_effect).sum()),
        fraction_moving_more_than_reference=round(
            float((np.abs(delta) > a.reference_effect).mean()), 4),
        top10_overlap=len(set(by_std[:10]) & set(by_alt[:10])),
        top20_overlap=len(set(by_std[:20]) & set(by_alt[:20])),
        leading_compound_standardized=by_std[0],
        leading_compound_alternative=by_alt[0])

    with (a.out / "Dyad_library_effect.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: -abs(r["change"])))
    (a.out / "Dyad_library_effect.json").write_text(json.dumps(summary, indent=2),
                                                     encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
