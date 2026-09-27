"""Compare retrospective enrichment between the two catalytic-dyad receptors.

The benchmark was run through the same funnel on both FP-2 receptors, which
differ only in the position of one catalytic-dyad proton. The two runs score
different numbers of compounds, because the receptor changes which ligands
return a pose, so a direct comparison of their headline metrics would confound
the dyad with the composition of the scored set.

This script removes that confound by restricting both runs to the compounds
scored by both, and recomputing every metric on that common set. The result is
a like-for-like statement of what the dyad assignment does to the retrospective
performance of the protocol.

  python benchmark_dyad_comparison.py --standard <scores.csv> --alternative <scores.csv> \
      --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import csv
import json

from enrichment_metrics import metrics, bootstrap_ci


def read(path):
    rows = {}
    for r in csv.DictReader(Path(path).open(encoding="utf-8")):
        rows[r["id"]] = dict(id=r["id"], label=int(r["label"]),
                             cls=r.get("class", ""), score=float(r["score"]))
    return rows


def subsets(rows):
    return dict(all_actives=list(rows),
                covalent_only=[r for r in rows
                               if r["label"] == 0 or r["cls"] == "covalent"])


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--standard", type=Path, required=True)
    p.add_argument("--alternative", type=Path, required=True)
    p.add_argument("--replicates", type=int, default=2000)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    std, alt = read(a.standard), read(a.alternative)
    common = sorted(set(std) & set(alt))
    if not common:
        raise SystemExit("the two runs share no scored compound")

    out = dict(
        experiment=("the same benchmark funnel run on two FP-2 receptors whose catalytic "
                    "dyad differs by the position of one proton, restricted to the "
                    "compounds that both receptors scored"),
        standard_file=a.standard.name, alternative_file=a.alternative.name,
        scored_by_standard=len(std), scored_by_alternative=len(alt),
        scored_by_both=len(common))

    for label, rows in (("standardized", std), ("alternative", alt)):
        block = [rows[c] for c in common]
        out[label] = {}
        for name, subset in subsets(block).items():
            m = metrics(subset)
            ci = bootstrap_ci(subset, a.replicates)
            m.update(dict(roc_ci=ci["roc_auc"], ef1_ci=ci["ef1"], ef5_ci=ci["ef5"]))
            out[label][name] = m

    # The difference between the two receptors is the quantity of interest, and
    # it has to be estimated on the paired data. Comparing each receptor
    # separately against random and noting that one clears the threshold while
    # the other does not is not evidence that the two differ, so the interval
    # below resamples compounds jointly and is reported alongside the values.
    import numpy as np
    rng = np.random.default_rng(1)
    out["difference"] = {}
    for name in ("all_actives", "covalent_only"):
        s, t = out["standardized"][name], out["alternative"][name]
        sset = subsets([std[c] for c in common])[name]
        aset = subsets([alt[c] for c in common])[name]
        n = len(sset)
        draws = []
        for _ in range(a.replicates):
            idx = rng.integers(0, n, n)
            draws.append(metrics([sset[i] for i in idx])["roc_auc"]
                         - metrics([aset[i] for i in idx])["roc_auc"])
        draws.sort()
        out["difference"][name] = dict(
            roc_auc=round(s["roc_auc"] - t["roc_auc"], 4),
            roc_auc_ci=[round(draws[int(0.025 * len(draws))], 4),
                        round(draws[int(0.975 * len(draws))], 4)],
            roc_auc_interval_excludes_zero=bool(
                draws[int(0.025 * len(draws))] > 0 or draws[int(0.975 * len(draws))] < 0),
            pr_auc=round(s["pr_auc"] - t["pr_auc"], 4),
            ef1=round(s["ef1"] - t["ef1"], 4),
            ef5=round(s["ef5"] - t["ef5"], 4),
            bedroc20=round(s["bedroc20"] - t["bedroc20"], 4),
            pairing=("compounds resampled jointly across the two receptors, so the "
                     "interval is on the difference itself"))

    rows_out = []
    for label in ("standardized", "alternative"):
        for name in ("all_actives", "covalent_only"):
            m = out[label][name]
            rows_out.append(dict(
                receptor=label, subset=name, n=m["n"], actives=m["actives"],
                roc_auc=round(m["roc_auc"], 4),
                roc_ci_low=round(m["roc_ci"][0], 4), roc_ci_high=round(m["roc_ci"][1], 4),
                pr_auc=round(m["pr_auc"], 4), ef1=round(m["ef1"], 4),
                ef5=round(m["ef5"], 4), bedroc20=round(m["bedroc20"], 4)))
    with (a.out / "Benchmark_dyad_comparison.csv").open("w", newline="",
                                                        encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_out[0]))
        w.writeheader()
        w.writerows(rows_out)
    (a.out / "Benchmark_dyad_comparison.json").write_text(json.dumps(out, indent=2),
                                                           encoding="utf-8")
    print(json.dumps(out["difference"], indent=2))
    for r in rows_out:
        print("%-13s %-14s n=%d actives=%d ROC=%.3f [%.3f, %.3f] EF5=%.2f BEDROC=%.3f"
              % (r["receptor"], r["subset"], r["n"], r["actives"], r["roc_auc"],
                 r["roc_ci_low"], r["roc_ci_high"], r["ef5"], r["bedroc20"]))


if __name__ == "__main__":
    main()
