"""Separate the effect of docking precision from the effect of decoy-pool size.

The retrospective benchmark was run twice on the same receptor and grid: once
under the recognition protocol itself, Glide XP over every prepared state on a
sampled decoy set, and once through an HTVS-then-SP funnel over the whole decoy
pool. Comparing those two runs directly would change precision and decoy count
at the same time, so neither could be assigned a size.

This script changes one variable at a time. The precision term is measured on
the compounds both runs scored, so the decoy set is held fixed. The decoy-pool
term is measured on the funnel alone, so the precision is held fixed. It also
reports the rank correlation between the two protocols, which is the quantity
that shows how far the two orderings differ even where the areas agree.

  python benchmark_precision_comparison.py --xp <scores.csv> --funnel <scores.csv> \
      --labels <subset labels.csv> --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import csv
import json
import math

from enrichment_metrics import metrics, bootstrap_ci


def read(path):
    return {r["id"]: dict(id=r["id"], label=int(r["label"]), cls=r.get("class", ""),
                          score=float(r["score"]))
            for r in csv.DictReader(Path(path).open(encoding="utf-8"))}


def spearman(u, v):
    def rank(z):
        order = sorted(range(len(z)), key=lambda i: z[i])
        out = [0.0] * len(z)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and z[order[j + 1]] == z[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                out[order[k]] = avg
            i = j + 1
        return out
    a, b = rank(u), rank(v)
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))
    return num / den if den else float("nan")


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--xp", type=Path, required=True,
                   help="protocol-matched scores, Glide XP over every prepared state")
    p.add_argument("--funnel", type=Path, required=True,
                   help="HTVS-then-SP scores over the full decoy pool")
    p.add_argument("--labels", type=Path, required=True,
                   help="labels of the sampled subset the XP run used")
    p.add_argument("--replicates", type=int, default=2000)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    xp, funnel = read(a.xp), read(a.funnel)
    sampled = {r["id"] for r in csv.DictReader(a.labels.open(encoding="utf-8"))}
    funnel_sampled = dict((k, v) for k, v in funnel.items() if k in sampled)
    common = sorted(set(xp) & set(funnel))

    def block(rows):
        m = metrics(rows)
        ci = bootstrap_ci(rows, a.replicates)
        return dict(n=m["n"], actives=m["actives"], roc_auc=round(m["roc_auc"], 4),
                    roc_ci=[round(ci["roc_auc"][0], 4), round(ci["roc_auc"][1], 4)])

    out = dict(
        design=("one variable at a time: the precision term holds the decoy set fixed by "
                "comparing the two protocols on the compounds both scored, and the "
                "decoy-pool term holds the precision fixed by comparing the funnel over "
                "the sampled and the full decoy pool"),
        xp_sampled_decoys=block(list(xp.values())),
        funnel_sampled_decoys=block(list(funnel_sampled.values())),
        funnel_full_pool=block(list(funnel.values())),
        xp_common=block([xp[c] for c in common]),
        funnel_common=block([funnel[c] for c in common]),
        compounds_scored_by_both=len(common))

    out["precision_effect"] = dict(
        definition="Glide XP minus HTVS/SP funnel, over the compounds both scored",
        delta_roc_auc=round(out["xp_common"]["roc_auc"] - out["funnel_common"]["roc_auc"], 4),
        spearman_rho=round(spearman([xp[c]["score"] for c in common],
                                    [funnel[c]["score"] for c in common]), 3))
    out["decoy_pool_effect"] = dict(
        definition="full decoy pool minus sampled decoys, both under the HTVS/SP funnel",
        delta_roc_auc=round(out["funnel_full_pool"]["roc_auc"]
                            - out["funnel_sampled_decoys"]["roc_auc"], 4))

    (a.out / "Benchmark_precision_comparison.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({k: out[k] for k in
                      ("precision_effect", "decoy_pool_effect",
                       "compounds_scored_by_both")}, indent=2))


if __name__ == "__main__":
    main()
