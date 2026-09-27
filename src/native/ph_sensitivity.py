"""pH sensitivity of the reported panel.

Falcipains act in the acidic food vacuole, so the panel is re-docked with both
the receptor and the ligands prepared at pH 5.5 and compared with the pH 7.0
protocol used for the reported screen. The catalytic dyad is held at the
standardized assignment at both pH values, so the comparison reports the effect
of pH on the remainder of the site and on ligand ionization rather than on the
dyad model.

Reported per target: the rank agreement between the two pH conditions over the
panel, the per-compound score change, and whether the matched 4-OH/4-Cl contrast
keeps its sign.

  run.exe python3 ph_sensitivity.py --ph55 FP2:<pv> FP3:<pv> \
      --reference ../exports/Primary_state_selection.csv \
      --matched ../exports/Matched_effects_by_rule.csv --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import csv
import json

from schrodinger import structure
from scipy.stats import spearmanr

MATCHED = {"Parent_SANC00867": ("glycosylated", "OH"),
           "AnalogueA": ("glycosylated", "Cl"),
           "AnalogueB": ("aglycone_ester", "OH"),
           "OPT1": ("aglycone_ester", "Cl")}
PAIRS = [("glycosylated", "Parent_SANC00867", "AnalogueA"),
         ("aglycone_ester", "AnalogueB", "OPT1")]


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


def best_per_title(path):
    """Minimum DockingScore per parent title in one pose-viewer file."""
    out = {}
    for n, st in enumerate(structure.StructureReader(str(path))):
        if n == 0 and st.atom_total > 1000:
            continue
        title = st.title
        ds = st.property.get("r_i_docking_score")
        gs = st.property.get("r_i_glide_gscore")
        if ds is None:
            continue
        parent = title.split("-")[0] if title not in MATCHED else title
        rec = dict(parent=parent, title=title, DockingScore=float(ds),
                   GlideScore=float(gs) if gs is not None else None,
                   variant=st.property.get("s_lp_Variant", title))
        if parent not in out or rec["DockingScore"] < out[parent]["DockingScore"]:
            out[parent] = rec
    return out


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--ph55", nargs="+", required=True, help="TARGET:path pairs")
    p.add_argument("--reference", type=Path, required=True,
                   help="Primary_state_selection.csv from the pH 7.0 screen")
    p.add_argument("--matched", type=Path, required=True,
                   help="Matched_effects_by_rule.csv from the pH 7.0 matched run")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    reference = {}
    for r in csv.DictReader(a.reference.open(encoding="utf-8")):
        reference[(r["target"], r["parent_id"])] = float(r["DockingScore"])
    matched_ref = {}
    for r in csv.DictReader(a.matched.open(encoding="utf-8")):
        if r["rule"] == "minimum_DockingScore":
            matched_ref[(r["target"], r["architecture"])] = float(r["Cl_minus_OH"])

    rows, summary = [], dict(
        protocol=("receptor prepared with the Protein Preparation Wizard at PROPKA pH 5.5 "
                  "and ligands prepared with LigPrep and Epik at pH 5.5 +/- 2.0; the "
                  "catalytic dyad is held at the standardized neutral assignment at both "
                  "pH values"),
        reference="pH 7.0 protocol used for the reported screen",
        rank_agreement={}, matched_contrast={}, panel_members={})

    ph55 = {}
    for spec in a.ph55:
        target, path = spec.split(":", 1)
        path = Path(path)
        if not path.is_file():
            print("missing", path)
            continue
        # a target may be supplied as more than one pose-viewer file, so that a
        # panel member docked in a later job joins the same comparison
        block = ph55.setdefault(target, {})
        for parent, rec in best_per_title(path).items():
            if parent not in block or rec["DockingScore"] < block[parent]["DockingScore"]:
                block[parent] = rec

    for target, block in sorted(ph55.items()):
        pairs = []
        for parent, rec in sorted(block.items()):
            ref = reference.get((target, parent))
            row = dict(target=target, compound=parent,
                       pH55_DockingScore=round(rec["DockingScore"], 4),
                       pH55_GlideScore=(round(rec["GlideScore"], 4)
                                        if rec["GlideScore"] is not None else ""),
                       pH55_state=rec["variant"],
                       pH70_DockingScore=round(ref, 4) if ref is not None else "",
                       change=round(rec["DockingScore"] - ref, 4) if ref is not None else "",
                       role=("matched-set member" if parent in MATCHED
                             else "panel member"))
            rows.append(row)
            if ref is not None:
                pairs.append((ref, rec["DockingScore"]))
        summary["panel_members"][target] = len(block)
        compared = [(parent, reference[(target, parent)], rec["DockingScore"])
                    for parent, rec in sorted(block.items())
                    if (target, parent) in reference]
        if len(compared) >= 3:
            rho = spearmanr([q[1] for q in compared], [q[2] for q in compared])
            summary["rank_agreement"][target] = round(
                float(getattr(rho, "statistic", rho[0])), 4)
            summary.setdefault("rank_agreement_ci", {})[target] = bootstrap_rho(
                [q[1] for q in compared], [q[2] for q in compared])
            summary.setdefault("compounds_compared", {})[target] = len(compared)
            by70 = [q[0] for q in sorted(compared, key=lambda q: q[1])]
            by55 = [q[0] for q in sorted(compared, key=lambda q: q[2])]
            summary.setdefault("leading_compound", {})[target] = dict(
                pH70=by70[0], pH55=by55[0], unchanged=bool(by70[0] == by55[0]))
            summary.setdefault("top_k_overlap", {})[target] = dict(
                ("top%d" % k, len(set(by70[:k]) & set(by55[:k])))
                for k in (3, 5, 10) if k <= len(compared))
            changes = [abs(q[2] - q[1]) for q in compared]
            summary.setdefault("absolute_score_change", {})[target] = dict(
                median=round(float(sorted(changes)[len(changes) // 2]), 3),
                maximum=round(float(max(changes)), 3))

        for arch, oh, cl in PAIRS:
            if oh in block and cl in block:
                d = block[cl]["DockingScore"] - block[oh]["DockingScore"]
                ref = matched_ref.get((target, arch))
                summary["matched_contrast"]["%s_%s" % (target, arch)] = dict(
                    pH55=round(d, 4),
                    pH70=round(ref, 4) if ref is not None else None,
                    sign_preserved=(bool(d * ref > 0) if ref is not None else None))

    with (a.out / "pH_sensitivity.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["target", "compound", "role",
                                           "pH70_DockingScore", "pH55_DockingScore",
                                           "change", "pH55_GlideScore", "pH55_state"])
        w.writeheader()
        w.writerows(rows)

    preserved = [v["sign_preserved"] for v in summary["matched_contrast"].values()
                 if v["sign_preserved"] is not None]
    summary["matched_pair_note"] = (
        "The matched 4-OH/4-Cl contrast keeps its sign in %d of the %d architecture and "
        "target combinations compared." % (sum(1 for x in preserved if x), len(preserved))
        if preserved else "")
    (a.out / "pH_sensitivity_summary.json").write_text(json.dumps(summary, indent=2),
                                                        encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
