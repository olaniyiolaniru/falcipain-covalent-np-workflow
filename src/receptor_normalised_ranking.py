"""Parent selection and receptor-normalized dual-target ranking.

Raw Glide scores from two different receptors are protocol-relative outputs of
separate grids. This script therefore builds the dual-target analysis on
within-receptor normalization and reports the raw composites alongside it:

  raw_mean            arithmetic mean of the two raw DockingScores
  raw_worse           the larger (less favourable) of the two raw scores
  mean_percentile     mean of the within-target percentile ranks
  worse_percentile    the worse of the two within-target percentile ranks
  mean_z              mean of the within-target z-scores
  rank_sum            sum of the within-target ordinal ranks

Percentile is computed within each target over all parents that returned a pose
at that target, so a value of 100 is the most favourable score at that receptor.
A Spearman agreement matrix over all six composites quantifies how much the
dual-target conclusion depends on the normalization chosen.

Size diagnostics are reported in the same table: heavy-atom count, score per
heavy atom, and the size-scaled rank, so that large glycosides and smaller
comparators can be read against each other.

Tie rule: minimum score per parent; ties within 1e-9 resolved by the lower
prepared-record index.

  python receptor_normalised_ranking.py --states ../exports/Primary_state_records.csv \
      --failures ../exports/Primary_state_failures.csv \
      --descriptors ../exports/Library_descriptors.csv --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import json
import itertools

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

TIE_TOL = 1e-9
COMPOSITES = ["raw_mean", "raw_worse", "mean_percentile", "worse_percentile",
              "mean_z", "rank_sum"]
# lower is better for raw scores and rank_sum; higher is better for percentiles
BETTER_LOW = {"raw_mean": True, "raw_worse": True, "mean_percentile": False,
              "worse_percentile": False, "mean_z": True, "rank_sum": True}


def select(frame, field):
    out, ties = [], 0
    for parent, block in frame.groupby("parent_id", sort=True):
        best = block[field].min()
        tied = block[np.isclose(block[field], best, rtol=0.0, atol=TIE_TOL)]
        if len(tied) > 1:
            ties += 1
        out.append(tied.sort_values("prepared_record").iloc[0])
    return pd.DataFrame(out).set_index("parent_id"), ties


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--states", type=Path, required=True)
    p.add_argument("--failures", type=Path, required=True)
    p.add_argument("--descriptors", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    states = pd.read_csv(a.states)
    failures = pd.read_csv(a.failures)
    desc = pd.read_csv(a.descriptors).set_index("compound")

    summary, selected, per_target = {}, [], {}
    for target in ("FP2", "FP3"):
        block = states[states.target.eq(target)]
        pen, pen_ties = select(block, "DockingScore")
        unp, unp_ties = select(block, "GlideScore")
        unp = unp.loc[pen.index]
        changed = (pen.prepared_state_id != unp.prepared_state_id)
        # within-target normalization over every parent that posed at this target
        s = pen.DockingScore.astype(float)
        pct = 100.0 * (1.0 - s.rank(method="average", ascending=True) / len(s))
        z = (s - s.mean()) / s.std(ddof=0)
        rank = s.rank(method="average", ascending=True)
        per_target[target] = pd.DataFrame(
            {"score": s, "percentile": pct, "z": z, "rank": rank})
        for parent in pen.index:
            x, y = pen.loc[parent], unp.loc[parent]
            selected.append(dict(
                target=target, parent_id=parent,
                DockingScore=x.DockingScore, GlideScore_of_selected_state=x.GlideScore,
                minimum_GlideScore=y.GlideScore,
                within_target_percentile=float(pct[parent]),
                within_target_z=float(z[parent]),
                within_target_rank=float(rank[parent]),
                DockingScore_state=x.prepared_state_id, GlideScore_state=y.prepared_state_id,
                DockingScore_pose_record=int(x.pose_record_index),
                GlideScore_pose_record=int(y.pose_record_index),
                selected_formal_charge=int(x.formal_charge),
                effective_state_penalty=x.effective_state_penalty,
                state_changes_under_GlideScore_rule=bool(changed[parent]),
                states_with_pose=int(len(block[block.parent_id.eq(parent)])),
                states_prepared=int(len(block[block.parent_id.eq(parent)])
                                    + len(failures[failures.target.eq(target)
                                                   & failures.parent_id.eq(parent)])),
                source_job=x.source_job))
        summary[target] = dict(
            returned_state_records=int(len(block)),
            parents_with_pose=int(len(pen)),
            states_without_pose=int(len(failures[failures.target.eq(target)])),
            changed_state_if_GlideScore_minimised=int(changed.sum()),
            rank_rho_DockingScore_vs_minimum_GlideScore=float(
                spearmanr(pen.DockingScore, unp.GlideScore).statistic),
            parents_with_ties=pen_ties + unp_ties,
            score_mean=float(s.mean()), score_sd=float(s.std(ddof=0)))

    sel = pd.DataFrame(selected)
    sel.to_csv(a.out / "Primary_state_selection.csv", index=False)

    paired = sel.pivot(index="parent_id", columns="target", values="DockingScore").dropna()
    for name, col in (("percentile", "within_target_percentile"),
                      ("z", "within_target_z"), ("rank", "within_target_rank")):
        w = sel.pivot(index="parent_id", columns="target", values=col)
        paired["FP2_" + name] = w.FP2
        paired["FP3_" + name] = w.FP3
    paired = paired.join(desc.drop(columns=[c for c in ("SMILES",) if c in desc.columns]))

    paired["raw_mean"] = paired[["FP2", "FP3"]].mean(axis=1)
    paired["raw_worse"] = paired[["FP2", "FP3"]].max(axis=1)
    paired["mean_percentile"] = paired[["FP2_percentile", "FP3_percentile"]].mean(axis=1)
    paired["worse_percentile"] = paired[["FP2_percentile", "FP3_percentile"]].min(axis=1)
    paired["mean_z"] = paired[["FP2_z", "FP3_z"]].mean(axis=1)
    paired["rank_sum"] = paired[["FP2_rank", "FP3_rank"]].sum(axis=1)
    paired["percentile_difference_FP3_minus_FP2"] = paired.FP3_percentile - paired.FP2_percentile
    paired["score_per_heavy_atom"] = paired.raw_mean / paired.heavy_atoms
    paired["size_scaled_rank"] = paired.score_per_heavy_atom.rank(method="average")
    for c in COMPOSITES:
        paired[c + "_rank"] = paired[c].rank(method="average",
                                             ascending=BETTER_LOW[c])
    paired["n_dominating_parents"] = [
        int(((paired[["FP2", "FP3"]] <= r[["FP2", "FP3"]]).all(axis=1)
             & (paired[["FP2", "FP3"]] < r[["FP2", "FP3"]]).any(axis=1)).sum())
        for _, r in paired.iterrows()]
    paired = paired.sort_values("mean_percentile", ascending=False)
    paired.reset_index().rename(columns={"parent_id": "compound"}).to_csv(
        a.out / "Recognition_objectives_paired.csv", index=False)

    agreement = {}
    for x, y in itertools.combinations(COMPOSITES, 2):
        rx = paired[x] * (1 if BETTER_LOW[x] else -1)
        ry = paired[y] * (1 if BETTER_LOW[y] else -1)
        agreement["%s_vs_%s" % (x, y)] = round(float(spearmanr(rx, ry).statistic), 4)

    top10 = {c: list(paired.nsmallest(10, c + "_rank").index) for c in COMPOSITES}
    top10["size_scaled"] = list(paired.nsmallest(10, "size_scaled_rank").index)
    overlap = {}
    for x, y in itertools.combinations(sorted(top10), 2):
        overlap["%s_vs_%s" % (x, y)] = len(set(top10[x]) & set(top10[y]))

    summary["paired_parents"] = int(len(paired))
    summary["composite_agreement_spearman"] = agreement
    summary["top10_by_composite"] = top10
    summary["top10_overlap"] = overlap
    summary["score_size_rho"] = {t: float(spearmanr(paired[t], paired.heavy_atoms).statistic)
                                 for t in ("FP2", "FP3", "raw_mean")}
    summary["pareto_front"] = list(paired[paired.n_dominating_parents.eq(0)].index)
    summary["tie_rule"] = ("minimum score per parent; ties within 1e-9 resolved by the "
                           "lower prepared-record index")
    summary["normalisation"] = ("percentile, z-score and ordinal rank computed within each "
                                "target over the parents that returned a pose at that target")
    (a.out / "Recognition_ranking_summary.json").write_text(json.dumps(summary, indent=2),
                                                            encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
