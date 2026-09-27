"""Retrospective enrichment metrics with bootstrap confidence intervals.

Consumes a scored actives-and-decoys table and reports ROC-AUC, PR-AUC, EF1%,
EF5% and BEDROC(alpha = 20) for the full active set and for the covalent-active
subset, each with a bootstrap 95% confidence interval obtained by resampling
actives and decoys separately. The ROC and precision-recall curves are exported
so the reported areas can be read against the curves that produced them.

Scores are Glide docking scores, so a more negative value ranks earlier.

  python enrichment_metrics.py --scores scores.csv --labels labels.csv \
      --receptor "FP-2 standardized dyad" --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import csv
import json
import math
import random

BOOTSTRAP = 2000
SEED = 1
ALPHA = 20.0


def rank_best_first(rows):
    return sorted(rows, key=lambda r: r["score"])


def roc_points(ordered):
    P = sum(r["label"] for r in ordered)
    N = len(ordered) - P
    tp = fp = 0
    pts = [(0.0, 0.0)]
    for r in ordered:
        if r["label"] == 1:
            tp += 1
        else:
            fp += 1
        pts.append((fp / N if N else 0.0, tp / P if P else 0.0))
    return pts


def pr_points(ordered):
    P = sum(r["label"] for r in ordered)
    tp = fp = 0
    pts = []
    for r in ordered:
        if r["label"] == 1:
            tp += 1
        else:
            fp += 1
        pts.append((tp / P if P else 0.0, tp / (tp + fp)))
    return pts


def roc_auc(ordered):
    P = sum(r["label"] for r in ordered)
    N = len(ordered) - P
    if P == 0 or N == 0:
        return float("nan")
    tp = 0
    area = 0.0
    for r in ordered:
        if r["label"] == 1:
            tp += 1
        else:
            area += tp
    return area / (P * N)


def pr_auc(ordered):
    P = sum(r["label"] for r in ordered)
    if P == 0:
        return float("nan")
    tp = fp = 0
    prev = 0.0
    area = 0.0
    for r in ordered:
        if r["label"] == 1:
            tp += 1
        else:
            fp += 1
        recall = tp / P
        area += (tp / (tp + fp)) * (recall - prev)
        prev = recall
    return area


def ef(ordered, frac):
    P = sum(r["label"] for r in ordered)
    n = max(1, int(round(frac * len(ordered))))
    hits = sum(r["label"] for r in ordered[:n])
    return (hits / n) / (P / len(ordered)) if P else float("nan")


def bedroc(ordered, alpha=ALPHA):
    """Truchon and Bayly (2007), best-scoring first, one-based ranks."""
    N = len(ordered)
    n = sum(r["label"] for r in ordered)
    if n == 0 or n == N:
        return float("nan")
    ra = n / N
    s = sum(math.exp(-alpha * (i + 1) / N)
            for i, r in enumerate(ordered) if r["label"] == 1)
    rie = (s / n) / ((1.0 / N) * (1 - math.exp(-alpha)) / (math.exp(alpha / N) - 1))
    return (rie * (ra * math.sinh(alpha / 2))
            / (math.cosh(alpha / 2) - math.cosh(alpha / 2 - alpha * ra))
            + 1.0 / (1 - math.exp(alpha * (1 - ra))))


def ranking_rule(rows):
    """Describe the ranking from the stages that actually supplied the scores.

    A run scored entirely at one precision is described by that precision; a
    funnel that mixes stages is described as the mixture it is. Writing a fixed
    sentence here would let the record drift away from the calculation.
    """
    stages = sorted(set(r.get("stage", "") for r in rows) - {""})
    if stages == ["XP"]:
        return ("each compound is ranked by its minimum Glide XP DockingScore "
                "across its prepared states, as the recognition screen ranks it")
    if stages == ["SP"]:
        return ("each compound is ranked by its minimum Glide SP DockingScore "
                "across its prepared states")
    if set(stages) == {"HTVS", "SP"}:
        return ("each compound is ranked by its SP score where the funnel produced "
                "one and by its HTVS score otherwise, which is the ranking a "
                "prospective run of this funnel produces")
    return "compounds ranked by score; stages present: " + ", ".join(stages or ["unlabelled"])


def metrics(rows):
    o = rank_best_first(rows)
    P = sum(r["label"] for r in o)
    return dict(n=len(o), actives=P, prevalence=P / len(o) if o else float("nan"),
                roc_auc=roc_auc(o), pr_auc=pr_auc(o),
                ef1=ef(o, 0.01), ef5=ef(o, 0.05), bedroc20=bedroc(o))


def bootstrap_ci(rows, replicates=BOOTSTRAP, seed=SEED):
    random.seed(seed)
    acts = [r for r in rows if r["label"] == 1]
    decs = [r for r in rows if r["label"] == 0]
    keys = ["roc_auc", "pr_auc", "ef1", "ef5", "bedroc20"]
    samples = dict((k, []) for k in keys)
    for _ in range(replicates):
        bs = ([random.choice(acts) for _ in acts]
              + [random.choice(decs) for _ in decs])
        m = metrics(bs)
        for k in keys:
            samples[k].append(m[k])
    ci = {}
    for k in keys:
        s = sorted(v for v in samples[k] if v == v)
        ci[k] = [s[int(0.025 * len(s))], s[int(0.975 * len(s))]] if s else [None, None]
    return ci


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--scores", type=Path, required=True,
                   help="csv with id and score columns; best score is most negative")
    p.add_argument("--labels", type=Path, required=True,
                   help="csv with id, label and class columns")
    p.add_argument("--receptor", default="", help="receptor the scores came from")
    p.add_argument("--replicates", type=int, default=BOOTSTRAP)
    p.add_argument("--tag", default="", help="suffix for the output file names")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    labels = {}
    for r in csv.DictReader(a.labels.open(encoding="utf-8")):
        labels[r["id"]] = (int(r["label"]), r.get("class", ""))

    rows, seen = [], set()
    for r in csv.DictReader(a.scores.open(encoding="utf-8")):
        cid = r["id"]
        if cid not in labels or cid in seen:
            continue
        label, cls = labels[cid]
        seen.add(cid)
        rows.append(dict(id=cid, label=label, cls=cls, score=float(r["score"]),
                         stage=r.get("stage", "")))
    if not rows:
        raise SystemExit("no scored rows matched the label table")

    # A compound that returns no pose is not scored, and how those compounds are
    # handled changes the answer. Dropping them assumes they would have ranked at
    # random; ranking them last assumes a docking failure is informative, which is
    # what happens in practice because such a compound is never bought. Both
    # treatments are reported so the reader can see the size of the choice.
    unscored = [cid for cid in labels if cid not in seen]
    worst = max(r["score"] for r in rows) + 1.0
    ranked_last = list(rows) + [dict(id=cid, label=labels[cid][0], cls=labels[cid][1],
                                     score=worst, stage="unscored")
                               for cid in unscored]

    subsets = dict(all_actives=rows,
                   covalent_only=[r for r in rows
                                  if r["label"] == 0 or r["cls"] == "covalent"],
                   all_actives_unscored_ranked_last=ranked_last,
                   covalent_only_unscored_ranked_last=[
                       r for r in ranked_last
                       if r["label"] == 0 or r["cls"] == "covalent"])
    summary = dict(
        receptor=a.receptor,
        labelled_compounds=len(labels),
        scored_compounds=len(rows),
        unscored_compounds=len(labels) - len(rows),
        unscored_actives=sum(1 for cid in labels
                             if cid not in seen and labels[cid][0] == 1),
        unscored_decoys=sum(1 for cid in labels
                            if cid not in seen and labels[cid][0] == 0),
        scored_by_stage=dict((st, sum(1 for r in rows if r["stage"] == st))
                             for st in sorted(set(r["stage"] for r in rows)) if st),
        ranking_rule=ranking_rule(rows),
        unscored_treatment=("metrics are reported both over the scored compounds alone "
                            "and over the whole labelled set with the unscored compounds "
                            "ranked last"),
        bootstrap_replicates=a.replicates,
        bootstrap_scheme="actives and decoys resampled separately with replacement",
        bedroc_alpha=ALPHA,
        score_field="Glide docking score, more negative ranks earlier",
        scores_file=a.scores.name, labels_file=a.labels.name)

    curves = []
    for name, subset in subsets.items():
        m = metrics(subset)
        ci = bootstrap_ci(subset, a.replicates)
        m.update(dict(roc_ci=ci["roc_auc"], pr_ci=ci["pr_auc"], ef1_ci=ci["ef1"],
                      ef5_ci=ci["ef5"], bedroc20_ci=ci["bedroc20"]))
        summary[name] = m
        o = rank_best_first(subset)
        roc, pr = roc_points(o), pr_points(o)
        step = max(1, len(o) // 800)
        for i in range(0, len(roc), step):
            fpr, tpr = roc[i]
            j = min(max(i - 1, 0), len(pr) - 1)
            recall, precision = pr[j]
            curves.append(dict(subset=name, index=i, fpr=round(fpr, 6),
                               tpr=round(tpr, 6), recall=round(recall, 6),
                               precision=round(precision, 6)))
        fpr, tpr = roc[-1]
        recall, precision = pr[-1]
        curves.append(dict(subset=name, index=len(roc) - 1, fpr=round(fpr, 6),
                           tpr=round(tpr, 6), recall=round(recall, 6),
                           precision=round(precision, 6)))

    tag = ("_" + a.tag) if a.tag else ""
    with (a.out / ("Benchmark_curves%s.csv" % tag)).open("w", newline="",
                                                         encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["subset", "index", "fpr", "tpr",
                                           "recall", "precision"])
        w.writeheader()
        w.writerows(curves)
    with (a.out / ("Benchmark_scores%s.csv" % tag)).open("w", newline="",
                                                         encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "label", "cls", "score", "stage"])
        w.writeheader()
        w.writerows(rows)
    (a.out / ("Benchmark_summary%s.json" % tag)).write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
