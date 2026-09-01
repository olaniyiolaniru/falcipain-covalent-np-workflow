#!/usr/bin/env python3
"""Compute retrospective enrichment metrics with bootstrap 95% CIs from a scored
actives+decoys set. Input: scores.csv with columns [id, label(1=active,0=decoy), score]
where 'score' is a Glide docking score (more negative = better). Emits EF1%, EF5%,
ROC-AUC, PR-AUC, BEDROC(alpha=20) and bootstrap CIs. Pure-stdlib (no numpy needed)."""
import csv, math, random, sys

def load(fn):
    rows = []
    for r in csv.DictReader(open(fn)):
        rows.append((r["id"], int(r["label"]), float(r["score"])))
    return rows

def rank_desc(rows):
    # better score = more negative -> ascending sort puts best first
    return sorted(rows, key=lambda x: x[2])

def roc_auc(ordered):
    P = sum(l for _, l, _ in ordered); Ntot = len(ordered); Nn = Ntot - P
    if P == 0 or Nn == 0: return float("nan")
    tp = fp = 0; auc = 0.0; prev_fp = 0
    for _, l, _ in ordered:
        if l == 1: tp += 1
        else:
            auc += tp  # area increment per negative encountered
            fp += 1
    return auc / (P * Nn)

def pr_auc(ordered):
    P = sum(l for _, l, _ in ordered)
    if P == 0: return float("nan")
    tp = 0; fp = 0; prev_recall = 0.0; area = 0.0
    for _, l, _ in ordered:
        if l == 1: tp += 1
        else: fp += 1
        recall = tp / P; prec = tp / (tp + fp)
        area += prec * (recall - prev_recall); prev_recall = recall
    return area

def ef(ordered, frac):
    P = sum(l for _, l, _ in ordered); Ntot = len(ordered)
    n = max(1, int(round(frac * Ntot)))
    hits = sum(l for _, l, _ in ordered[:n])
    return (hits / n) / (P / Ntot) if P else float("nan")

def bedroc_clean(ordered, alpha=20.0):
    """Truchon & Bayly (2007) BEDROC. ordered = best-scoring first; r_i = 1-based rank."""
    N = len(ordered); n = sum(l for _, l, _ in ordered)
    if n == 0 or n == N: return float("nan")
    ra = n / N
    sum_exp = sum(math.exp(-alpha * (i + 1) / N)
                  for i, (_, l, _) in enumerate(ordered) if l == 1)
    # Robust Initial Enhancement
    rie = (sum_exp / n) / ((1.0 / N) * (1 - math.exp(-alpha)) / (math.exp(alpha / N) - 1))
    bed = rie * (ra * math.sinh(alpha / 2)) / \
        (math.cosh(alpha / 2) - math.cosh(alpha / 2 - alpha * ra)) \
        + 1.0 / (1 - math.exp(alpha * (1 - ra)))
    return bed

def metrics(rows):
    o = rank_desc(rows)
    return {"N": len(o), "n_active": sum(l for _, l, _ in o),
            "EF1%": ef(o, 0.01), "EF5%": ef(o, 0.05),
            "ROC_AUC": roc_auc(o), "PR_AUC": pr_auc(o), "BEDROC20": bedroc_clean(o)}

def bootstrap(rows, B=2000, seed=1):
    random.seed(seed)
    acts = [r for r in rows if r[1] == 1]; decs = [r for r in rows if r[1] == 0]
    keys = ["EF1%", "EF5%", "ROC_AUC", "PR_AUC", "BEDROC20"]
    samples = {k: [] for k in keys}
    for _ in range(B):
        bs = [random.choice(acts) for _ in acts] + [random.choice(decs) for _ in decs]
        m = metrics(bs)
        for k in keys: samples[k].append(m[k])
    ci = {}
    for k in keys:
        s = sorted(v for v in samples[k] if v == v)
        if not s: ci[k] = (float("nan"), float("nan")); continue
        ci[k] = (s[int(0.025 * len(s))], s[int(0.975 * len(s))])
    return ci

if __name__ == "__main__":
    fn = sys.argv[1] if len(sys.argv) > 1 else "scores.csv"
    rows = load(fn)
    m = metrics(rows); ci = bootstrap(rows)
    print(f"Set: {fn}  N={m['N']}  actives={m['n_active']}")
    for k in ["EF1%", "EF5%", "ROC_AUC", "PR_AUC", "BEDROC20"]:
        lo, hi = ci[k]
        print(f"  {k:9s} = {m[k]:.3f}   95% CI [{lo:.3f}, {hi:.3f}]")
