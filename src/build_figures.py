"""Generate every manuscript figure and the table-of-contents graphic from the
supplied tables. Panel letters and panel content here are the definitive
reference for the captions.

  Figure 1  receptor-normalized recognition landscape
  Figure 2  cinnamate electronic response
  Figure 3  matched analogue structures and scores
  Figure 4  structural basis and receptor dependence of the matched contrast
  Figure S1 retrospective enrichment benchmark
  TOC graphic

  python build_figures.py --data ../exports --structures ../data/structures --out ../figures
"""
from pathlib import Path
from argparse import ArgumentParser
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image

TEAL, BLUE, ORANGE, VIOLET = "#287b66", "#245e80", "#b96921", "#765590"
GREY, INK = "#a6b8c4", "#22323c"
SUBS = ["OMe", "OH", "H", "Cl", "CN", "NO2"]
SUBLAB = ["OMe", "OH", "H", "Cl", "CN", "NO$_2$"]
COMPOSITES = ["raw_mean", "raw_worse", "mean_percentile", "worse_percentile",
              "mean_z", "rank_sum"]
COMPLAB = ["raw mean", "raw worse", "mean pct", "worse pct", "mean z", "rank sum"]

plt.rcParams.update({
    "font.family": ["Arial", "DejaVu Sans"], "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#5b6a72", "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK,
    "pdf.fonttype": 42, "savefig.facecolor": "white",
})


def save(fig, out, stem):
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(str(out / (stem + "." + ext)), dpi=450, bbox_inches="tight")
    plt.close(fig)


def save_exact(fig, out, stem):
    """Save at the figure's declared size, opaque and with no alpha channel.

    The table-of-contents graphic has to arrive at the size the journal asks
    for, so the bounding box is not trimmed to the drawn content, and the PNG
    is flattened onto white because a transparent background can print black.
    """
    out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "eps"):
        fig.savefig(str(out / (stem + "." + ext)), dpi=450, facecolor="white")
    plt.close(fig)
    png = out / (stem + ".png")
    im = Image.open(png)
    if im.mode in ("RGBA", "LA", "P"):
        im = im.convert("RGBA")
        flat = Image.new("RGB", im.size, (255, 255, 255))
        flat.paste(im, mask=im.split()[-1])
        flat.save(png, dpi=(450, 450))
    # TIF at 300 dpi is the format the journal asks for
    im = Image.open(png).convert("RGB")
    im.save(out / (stem + ".tif"), format="TIFF", compression="tiff_lzw",
            dpi=(450, 450))


def panel(ax, letter, title):
    ax.set_title(letter + "  " + title, loc="left", fontweight="bold", fontsize=9.5, pad=6)


# --------------------------------------------------------------------------- 1
def figure1(data, out):
    p = pd.read_csv(data / "Recognition_objectives_paired.csv").set_index("compound")
    summary = json.loads((data / "Recognition_ranking_summary.json").read_text(encoding="utf-8"))
    keys = ["SANC00370", "SANC01049", "SANC00867", "SANC00853"]
    names = ["Acteoside", "Aloeresin C", "Tribuloside", "Rosmarinic acid"]
    colours = [TEAL, BLUE, ORANGE, VIOLET]

    fig = plt.figure(figsize=(7.2, 6.6))
    gs = fig.add_gridspec(2, 2, hspace=.46, wspace=.46)

    ax = fig.add_subplot(gs[0, 0])
    ax.scatter(p.FP2, p.FP3, s=11, color=GREY, alpha=.6, linewidths=0, zorder=2)
    lim = (-11.0, 5.0)
    ax.plot(lim, lim, color="#b1b8bc", lw=.8, ls="--", zorder=1)
    # label offsets are chosen so no two labels overlap at this figure size
    offs = {"SANC00370": (9, -10), "SANC01049": (0, 32),
            "SANC00867": (10, 8), "SANC00853": (9, -10)}
    for k, nm, c in zip(keys, names, colours):
        r = p.loc[k]
        o = offs[k]
        ax.scatter(r.FP2, r.FP3, s=30, color=c, edgecolor="white", linewidth=.5, zorder=4)
        ax.annotate(nm, (r.FP2, r.FP3), xytext=o, textcoords="offset points",
                    fontsize=7.2, color=c, ha="left",
                    arrowprops=dict(arrowstyle="-", color=c, lw=.5, shrinkA=0, shrinkB=2))
    ax.set_xlim(*lim)
    ax.set_ylim(*lim)
    ax.set_xlabel("FP-2 XP DockingScore")
    ax.set_ylabel("FP-3 XP DockingScore")
    panel(ax, "A", "Paired recognition, 174 parents")

    # B: percentile-percentile with the top-decile-at-both region shaded
    ax = fig.add_subplot(gs[0, 1])
    ax.add_patch(Rectangle((90, 90), 12, 12, facecolor="#e6f0ea", edgecolor="none", zorder=1))
    ax.scatter(p.FP2_percentile, p.FP3_percentile, s=11, color=GREY, alpha=.6,
               linewidths=0, zorder=2)
    front = summary.get("pareto_front", [])
    poffs = {"SANC00370": (-10, -16), "SANC01049": (-10, -30)}
    for k in front:
        r = p.loc[k]
        c = TEAL if k == "SANC00370" else BLUE
        ax.scatter(r.FP2_percentile, r.FP3_percentile, s=34, color=c,
                   edgecolor="white", linewidth=.5, zorder=4)
        ax.annotate(names[keys.index(k)] if k in keys else k,
                    (r.FP2_percentile, r.FP3_percentile), xytext=poffs.get(k, (-10, -16)),
                    textcoords="offset points", fontsize=7.2, color=c, ha="right",
                    arrowprops=dict(arrowstyle="-", color=c, lw=.5, shrinkA=0, shrinkB=2))
    ax.text(50, 94.5, "top decile at both", ha="center", fontsize=6.8, color=TEAL, zorder=3)
    ax.annotate("", xy=(89, 94.5), xytext=(70, 94.5),
                arrowprops=dict(arrowstyle="->", color=TEAL, lw=.7))
    ax.set_xlim(0, 102)
    ax.set_ylim(0, 102)
    ax.set_xlabel("FP-2 within-target percentile")
    ax.set_ylabel("FP-3 within-target percentile")
    panel(ax, "B", "Receptor-normalized profile")

    # C: composite agreement
    ax = fig.add_subplot(gs[1, 0])
    agree = summary["composite_agreement_spearman"]
    n = len(COMPOSITES)
    m = np.ones((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            k1 = COMPOSITES[i] + "_vs_" + COMPOSITES[j]
            k2 = COMPOSITES[j] + "_vs_" + COMPOSITES[i]
            m[i, j] = agree.get(k1, agree.get(k2, np.nan))
    im = ax.imshow(m, cmap="BuGn", vmin=.90, vmax=1.0)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, "%.2f" % m[i, j], ha="center", va="center", fontsize=6.4,
                    color="white" if m[i, j] > .97 else INK)
    ax.set_xticks(range(n))
    ax.set_xticklabels(COMPLAB, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(n))
    ax.set_yticklabels(COMPLAB, fontsize=7)
    ax.set_xticks(np.arange(-.5, n, 1), minor=True)
    ax.set_yticks(np.arange(-.5, n, 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.2)
    ax.tick_params(which="minor", length=0)
    cb = fig.colorbar(im, ax=ax, fraction=.040, pad=.02, shrink=.86)
    cb.set_label(r"Spearman $\rho$", fontsize=7.0, labelpad=2)
    cb.ax.tick_params(labelsize=6.2)
    panel(ax, "C", "Composites agree")

    # D: size axis
    ax = fig.add_subplot(gs[1, 1])
    ax.scatter(p.heavy_atoms, p.raw_mean, s=11, color=GREY, alpha=.6, linewidths=0, zorder=2)
    for k, n_, c in zip(keys, names, colours):
        r = p.loc[k]
        ax.scatter(r.heavy_atoms, r.raw_mean, s=30, color=c, edgecolor="white",
                   linewidth=.5, zorder=4)
    rho = summary["score_size_rho"]["raw_mean"]
    ax.text(.03, .06, r"$\rho$ = %.3f" % rho, transform=ax.transAxes, fontsize=7.4,
            color="#5b6a72")
    ax.set_xlabel("Heavy atoms")
    ax.set_ylabel("Mean raw DockingScore")
    panel(ax, "D", "Size is an independent axis")

    save(fig, out, "Figure_1_Recognition")


# --------------------------------------------------------------------------- 2
def figure2(data, structures, out):
    d = pd.read_csv(data / "Series_electronic_descriptors.csv")
    prim = d[(d.dataset == "primary") & (d.scheme == "ESP")].set_index("substituent")
    frag = pd.read_csv(data / "Series_fragment_response.csv")
    chk = frag[frag.dataset == "fixed_geometry_check"]
    checked = [s for s in SUBS if s in set(chk.substituent)]

    fig = plt.figure(figsize=(7.2, 7.2))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.05], hspace=.62, wspace=.34)

    for idx, (letter, col, ylab, title, fmt) in enumerate([
            ("A", "omega_eV", r"$\omega$ (eV)", "Global electrophilicity", "%.2f"),
            ("B", "fplus_beta", r"$f^{+}_{\beta}$ (ESP)",
             r"Response at the $\beta$-carbon (C6)", "%.3f")]):
        ax = fig.add_subplot(gs[0, idx])
        vals = [prim.loc[s, col] for s in SUBS]
        cols = [ORANGE if s == "NO2" else BLUE for s in SUBS]
        ax.bar(range(6), vals, color=cols, width=.66)
        for i, v in enumerate(vals):
            ax.text(i, v + max(vals) * .03, fmt % v, ha="center", fontsize=7.4)
        ax.set_xticks(range(6))
        ax.set_xticklabels(SUBLAB)
        ax.set_ylim(0, max(vals) * 1.22)
        ax.set_ylabel(ylab)
        ax.set_xlabel("para substituent")
        panel(ax, letter, title)

    chkd = d[d.dataset == "fixed_geometry_check"]
    e = chkd[chkd.scheme == "ESP"].set_index("substituent").fplus_beta
    m = chkd[chkd.scheme == "Mulliken"].set_index("substituent").fplus_beta
    ax = fig.add_subplot(gs[1, 0])
    offs = {"OMe": (7, 2), "OH": (-8, 1), "H": (0, 8), "Cl": (7, -9),
            "CN": (8, -2), "NO2": (8, -2)}
    for s, lab_s in zip(SUBS, SUBLAB):
        c = ORANGE if s in ("NO2", "CN") else BLUE
        ax.scatter(e[s], m[s], s=38, color=c, edgecolor="white", linewidth=.5, zorder=3)
        dx, dy = offs[s]
        ax.annotate(lab_s, (e[s], m[s]), xytext=(dx, dy), textcoords="offset points",
                    fontsize=7.4, color=c,
                    ha="center" if dx == 0 else ("left" if dx > 0 else "right"))
    ax.plot([e["OH"], e["Cl"]], [m["OH"], m["Cl"]], color="#b1b8bc", lw=1.0, ls="--", zorder=1)
    ax.set_xlim(0, .26)
    ax.set_ylim(0, .16)
    ax.set_xlabel(r"$f^{+}_{\beta}$ (ESP)")
    ax.set_ylabel(r"$f^{+}_{\beta}$ (Mulliken)")
    ax.text(.03, .955, r"Spearman $\rho$ = %.2f"
            % spearmanr([e[s] for s in SUBS], [m[s] for s in SUBS]).statistic,
            transform=ax.transAxes, fontsize=7.2, color="#5b6a72", va="top")
    panel(ax, "C", "OH/Cl order depends on the partition")

    ax = fig.add_subplot(gs[1, 1])
    order = ["ester_alkene", "aryl", "para_substituent"]
    lab = ["Ester-alkene", "Aryl", "para substituent"]
    cmap = [TEAL, "#8d9aa2", ORANGE]
    width = .34
    for j, s in enumerate(checked):
        for k, scheme in enumerate(("ESP", "Mulliken")):
            base, x = 0.0, j + (k - .5) * width
            for reg, c in zip(order, cmap):
                v = float(chk[(chk.substituent == s) & (chk.scheme == scheme)
                              & (chk.region == reg)].fplus_sum.iloc[0])
                ax.bar(x, v, bottom=base, width=width * .9, color=c,
                       hatch="///" if scheme == "Mulliken" else None,
                       edgecolor="white", linewidth=.4)
                base += v
            ax.text(x, 1.04, "E" if scheme == "ESP" else "M", ha="center",
                    fontsize=6.6, color="#5b6a72")
    ax.set_xticks(range(len(checked)))
    ax.set_xticklabels([SUBLAB[SUBS.index(s)] for s in checked])
    ax.set_ylim(0, 1.16)
    ax.set_ylabel(r"$\sum_A f^{+}_A$")
    ax.set_xlabel("para substituent")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in cmap]
    ax.legend(handles, lab, frameon=False, fontsize=6.6, loc="upper center",
              bbox_to_anchor=(.5, -.30), ncol=3, handlelength=1.1, columnspacing=1.0)
    panel(ax, "D", "Fragment partition of the unit response")

    ax = fig.add_subplot(gs[2, 0])
    img = structures / "NO2_fragment_diagram.png"
    if img.is_file():
        ax.imshow(Image.open(str(img)))
    ax.axis("off")
    panel(ax, "E", "Fragment definition, methyl 4-nitrocinnamate")

    esp = {r.region: r.fplus_sum
           for r in chk[(chk.substituent == "NO2") & (chk.scheme == "ESP")].itertuples()}
    mul = {r.region: r.fplus_sum
           for r in chk[(chk.substituent == "NO2") & (chk.scheme == "Mulliken")].itertuples()}
    ax = fig.add_subplot(gs[2, 1])
    ax.axis("off")
    ax.text(0, .92,
            "Methyl 4-nitrocinnamate, fixed-geometry check jobs\n\n"
            "Nitro fragment       %.3f (ESP)   %.3f (Mulliken)\n"
            "Aryl ring                %.3f (ESP)   %.3f (Mulliken)\n"
            "Ester-alkene         %.3f (ESP)   %.3f (Mulliken)\n\n"
            "Fragment sums include attached hydrogens and total\n"
            "one electron. They are signed charge-response\n"
            "partitions, not orbital populations or probabilities."
            % (esp["para_substituent"], mul["para_substituent"],
               esp["aryl"], mul["aryl"],
               esp["ester_alkene"], mul["ester_alkene"]),
            transform=ax.transAxes, fontsize=7.8, va="top", color=INK,
            family=["DejaVu Sans Mono", "monospace"])
    save(fig, out, "Figure_2_Electronic_Response")


# --------------------------------------------------------------------------- 3
def figure3(data, structures, out):
    eff = pd.read_csv(data / "Matched_effects_by_rule.csv")
    eff = eff[eff.rule.isin(["minimum_DockingScore", "minimum_GlideScore"])]
    records = json.loads((structures / "analogue_structures.json").read_text(encoding="utf-8"))

    fig = plt.figure(figsize=(7.2, 7.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 1.15], hspace=.18, wspace=.14)
    letters = "ABCD"
    for i, r in enumerate(records):
        ax = fig.add_subplot(gs[i // 2, i % 2])
        img = structures / (r["compound"] + ".png")
        if img.is_file():
            ax.imshow(Image.open(str(img)))
        ax.axis("off")
        ax.set_title(letters[i] + "  " + r["label"], fontsize=8.4, fontweight="bold",
                     color=BLUE if i < 2 else ORANGE, loc="left", pad=2)

    for j, t in enumerate(["FP2", "FP3"]):
        ax = fig.add_subplot(gs[2, j])
        for arch, label, c in (("glycosylated", "Glycosylated", BLUE),
                               ("aglycone_ester", "Aglycone ester", ORANGE)):
            row = eff[(eff.target == t) & (eff.architecture == arch)
                      & (eff.rule == "minimum_DockingScore")].iloc[0]
            ax.plot([0, 1], [row.OH_value, row.Cl_value], "o-", color=c, lw=1.8, ms=5,
                    label=label)
            dy = 9 if arch == "aglycone_ester" else -15
            for x, y in ((0, row.OH_value), (1, row.Cl_value)):
                ax.annotate("%.2f" % y, (x, y), xytext=(0, dy), textcoords="offset points",
                            ha="center", fontsize=7.4, color=c)
            row2 = eff[(eff.target == t) & (eff.architecture == arch)
                       & (eff.rule == "minimum_GlideScore")].iloc[0]
            ax.plot([0, 1], [row2.OH_value, row2.Cl_value], "o--", color=c, lw=1.0,
                    ms=3.4, alpha=.55)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["4-OH", "4-Cl"])
        ax.set_xlim(-.25, 1.25)
        ax.set_ylim(-9.1, -2.4)
        ax.set_ylabel("XP score (lower is more favourable)")
        panel(ax, "EF"[j], t.replace("FP", "FP-"))
        ax.grid(axis="y", alpha=.12)
    handles = [plt.Line2D([], [], color=BLUE, marker="o", lw=1.8, ms=5),
               plt.Line2D([], [], color=ORANGE, marker="o", lw=1.8, ms=5),
               plt.Line2D([], [], color="#414b53", marker="o", lw=1.6, ms=4.4),
               plt.Line2D([], [], color="#414b53", marker="o", ls="--", lw=1.0, ms=3.4,
                          alpha=.55)]
    fig.legend(handles,
               ["Glycosylated", "Aglycone ester",
                "minimum DockingScore", "minimum GlideScore"],
               frameon=False, fontsize=7.4, ncol=4, loc="lower center",
               bbox_to_anchor=(.5, -.015), handlelength=1.8, columnspacing=1.6)
    save(fig, out, "Figure_3_Analogue_Comparison")


# --------------------------------------------------------------------------- 4
def figure4(data, out):
    core = pd.read_csv(data / "Matched_core_displacement.csv")
    core = core[core.selection_rule == "minimum_DockingScore"]
    hal = pd.read_csv(data / "Matched_halogen_bond_check.csv")
    con = pd.read_csv(data / "Matched_pose_contacts.csv")
    fam = pd.read_csv(data / "Matched_pose_family_summary.csv")
    eff = pd.read_csv(data / "Matched_effects_by_rule.csv")
    dyad_path = data / "Dyad_sensitivity.json"
    dyad = json.loads(dyad_path.read_text(encoding="utf-8")) if dyad_path.is_file() else None

    ncol = 4 if dyad else 3
    widths = [1.16, 1.12, 1.0, 1.05][:ncol]
    fig, axs = plt.subplots(1, ncol, figsize=(2.10 * ncol + 1.2, 3.05),
                            gridspec_kw={"width_ratios": widths})

    # A: selected-pose displacement against the ensemble distribution
    ax = axs[0]
    labels, vals, cols, lo, med, hi = [], [], [], [], [], []
    for r in core.itertuples():
        key = fam[(fam.target == r.target) & (fam.architecture == r.architecture)]
        labels.append(r.target.replace("FP", "FP-") + "\n"
                      + ("glyco." if r.architecture == "glycosylated" else "aglyc."))
        vals.append(r.core_displacement_rmsd_A)
        cols.append(BLUE if r.architecture == "glycosylated" else ORANGE)
        if len(key):
            lo.append(float(key.min_pairwise_core_rmsd_A.iloc[0]))
            med.append(float(key.median_pairwise_core_rmsd_A.iloc[0]))
            hi.append(float(key.max_pairwise_core_rmsd_A.iloc[0]))
        else:
            lo.append(np.nan); med.append(np.nan); hi.append(np.nan)
    x = np.arange(len(vals))
    for i in range(len(vals)):
        if not np.isnan(lo[i]):
            ax.vlines(x[i], lo[i], hi[i], color="#c3ced5", lw=5, zorder=1,
                      capstyle="round")
            ax.scatter(x[i], med[i], marker="_", s=110, color="#5b6a72", zorder=2)
    ax.scatter(x, vals, s=44, color=cols, edgecolor="white", linewidth=.6, zorder=3)
    for i, v in enumerate(vals):
        ax.text(i + .18, v, "%.1f" % v, fontsize=7.2, va="center")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_xlim(-.6, len(vals) - .2)
    ax.set_ylabel("Mapped-core RMSD (Å)")
    ax.set_ylim(0, max([h for h in hi if not np.isnan(h)] + vals) * 1.12)
    panel(ax, "A", "Selected pose vs ensemble")

    ax = axs[1]
    ax.add_patch(Rectangle((2.2, 155), 1.15, 27, facecolor="#dff0e6", edgecolor="none", zorder=1))
    ax.text(2.78, 168, "halogen-bond\ncriterion", ha="center", va="center",
            fontsize=6.6, color=TEAL, zorder=2)
    marks = {"O": "o", "N": "s", "S": "^"}
    for r in hal.itertuples():
        ax.scatter(r.Cl_acceptor_distance_A, r.C_Cl_acceptor_angle_deg,
                   marker=marks.get(r.acceptor_element, "o"), s=34,
                   color=ORANGE, edgecolor="white", linewidth=.5, zorder=3)
    ax.set_xlim(2.2, 4.4)
    ax.set_ylim(60, 190)
    ax.set_xlabel("Cl···acceptor distance (Å)")
    ax.set_ylabel("C–Cl···acceptor angle (°)")
    panel(ax, "B", "No halogen-bond geometry")

    ax = axs[2]
    sel = {("FP2", "Parent_SANC00867"): BLUE, ("FP2", "AnalogueA"): ORANGE,
           ("FP3", "Parent_SANC00867"): BLUE, ("FP3", "AnalogueA"): ORANGE}
    xs, heights, cols2, labels2 = [], [], [], []
    pos = 0
    for target in ("FP2", "FP3"):
        for title in ("Parent_SANC00867", "AnalogueA"):
            block = con[(con.target == target) & (con.title == title)]
            row = block.loc[block.DockingScore.idxmin()]
            xs.append(pos)
            heights.append(row.n_polar_contacts_3p5A)
            cols2.append(sel[(target, title)])
            labels2.append(target.replace("FP", "FP-") + "\n"
                           + ("4-OH" if "Parent" in title else "4-Cl"))
            pos += 1
        pos += .4
    ax.bar(xs, heights, color=cols2, width=.66)
    for xx, v in zip(xs, heights):
        ax.text(xx, v + .15, str(int(v)), ha="center", fontsize=7.4)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels2, fontsize=7)
    ax.set_ylabel("Polar contacts ≤ 3.5 Å")
    ax.set_ylim(0, max(heights) * 1.25)
    panel(ax, "C", "Glycosylated pair")

    if dyad:
        ax = axs[3]
        fp3 = float(eff[(eff.target == "FP3") & (eff.architecture == "glycosylated")
                        & (eff.rule == "minimum_DockingScore")].Cl_minus_OH.iloc[0])
        fp3g = float(eff[(eff.target == "FP3") & (eff.architecture == "glycosylated")
                         & (eff.rule == "minimum_GlideScore")].Cl_minus_OH.iloc[0])
        std = float(dyad["standardized"]["Cl_minus_OH_DockingScore"])
        stdg = float(dyad["standardized"]["Cl_minus_OH_GlideScore"])
        alt = float(dyad["alternative"]["Cl_minus_OH_DockingScore"])
        altg = float(dyad["alternative"]["Cl_minus_OH_GlideScore"])
        names = ["FP-3\nneutral\ndyad", "FP-2\nneutral\ndyad", "FP-2\nion-pair\ndyad"]
        vals3 = [fp3, std, alt]
        gvals = [fp3g, stdg, altg]
        cols3 = [TEAL, BLUE, ORANGE]
        ax.axhline(0, color="#b1b8bc", lw=.8, zorder=1)
        ax.bar(range(3), vals3, color=cols3, width=.58, zorder=2)
        for i, (v, g) in enumerate(zip(vals3, gvals)):
            ax.hlines(g, i - .29, i + .29, color="#5b6a72", lw=1.1, ls="--", zorder=3)
            ax.annotate("%+.2f" % v, (i, v), xytext=(0, 5 if v >= 0 else -11),
                        textcoords="offset points", ha="center", fontsize=7.0)
        ax.annotate("", xy=(2.42, alt), xytext=(2.42, std),
                    arrowprops=dict(arrowstyle="<->", color=INK, lw=.9))
        ax.text(2.52, (alt + std) / 2, "%.2f\nunits" % abs(alt - std), fontsize=6.8,
                color=INK, va="center", ha="left")
        ax.set_xticks(range(3))
        ax.set_xticklabels(names, fontsize=6.6)
        ax.set_xlim(-.62, 3.15)
        ax.set_ylim(min(vals3) * 1.35, max(vals3) * 1.30)
        ax.set_ylabel("Cl − OH contrast (score units)", fontsize=8)
        panel(ax, "D", "Dyad assignment")

    fig.tight_layout(w_pad=1.7)
    save(fig, out, "Figure_4_Structural_Basis")


# ------------------------------------------------------------------ benchmark
def figure_benchmark(data, out):
    """ROC, precision-recall and early recognition at both reported receptors.

    The falcipain-3 covalent subset holds three actives and is left off the
    figure, as it is left out of every statement in the article.
    """
    path = data / "Benchmark_curves.csv"
    summary_path = data / "Benchmark_summary.json"
    if not path.is_file() or not summary_path.is_file():
        print("benchmark curves not available; skipping Figure S1")
        return
    c = pd.read_csv(path)
    s = json.loads(summary_path.read_text(encoding="utf-8"))

    # the protocol-matched runs lead; the funnel is shown as the lower-precision
    # comparison over the full decoy pool
    def load(tag):
        cp = data / ("Benchmark_curves_%s.csv" % tag)
        sp = data / ("Benchmark_summary_%s.json" % tag)
        if cp.is_file() and sp.is_file():
            return pd.read_csv(cp), json.loads(sp.read_text(encoding="utf-8"))
        return None, None

    cx2, sx2 = load("FP2_XP")
    cx3, sx3 = load("FP3_XP")
    series = []
    if sx2 is not None:
        series.append(("FP-2 XP, all actives", cx2, sx2, "all_actives", BLUE))
        series.append(("FP-2 XP, covalent", cx2, sx2, "covalent_only", ORANGE))
    series.append(("FP-2 HTVS/SP, all actives", c, s, "all_actives", TEAL))
    if sx3 is not None:
        series.append(("FP-3 XP, all actives", cx3, sx3, "all_actives", VIOLET))

    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.7))

    ax = axs[0]
    for label, curves, summary, key, col in series:
        b = curves[curves.subset.eq(key)]
        if not len(b):
            continue
        ax.plot(b.fpr, b.tpr, color=col, lw=1.6,
                label="%s (%.3f)" % (label.replace(", all actives", ", all"),
                                     summary[key]["roc_auc"]))
    ax.plot([0, 1], [0, 1], color="#b1b8bc", lw=.8, ls="--")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.legend(frameon=False, fontsize=6.4, loc="lower right")
    panel(ax, "A", "ROC")

    ax = axs[1]
    for label, curves, summary, key, col in series:
        b = curves[curves.subset.eq(key)]
        if not len(b):
            continue
        ax.plot(b.recall, b.precision, color=col, lw=1.6,
                label="%s (%.3f)" % (label.replace(" actives", ""),
                                     summary[key]["pr_auc"]))
        ax.axhline(summary[key]["prevalence"], color=col, lw=.8, ls=":")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.legend(frameon=False, fontsize=6.6, loc="upper right")
    panel(ax, "B", "Precision–recall")

    # EF 1% is omitted for falcipain-3: the top 1% of that list is two compounds
    ax = axs[2]
    metrics = ["ef1", "ef5", "bedroc20"]
    labels = ["EF 1%", "EF 5%", "BEDROC$_{20}$"]
    xs = np.arange(len(metrics))
    width = .84 / len(series)
    for k, (label, _curves, summary, key, col) in enumerate(series):
        vals = [summary[key][m] for m in metrics]
        off = (k - (len(series) - 1) / 2.0) * width
        ax.bar(xs + off, vals, width=width * .92, color=col, label=label)
    ax.axhline(1.0, color="#b1b8bc", lw=.8, ls="--")
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=7.6)
    ax.set_ylabel("Value")
    ax.set_ylim(0, max(v for _, _, sm, ky, _ in series
                       for v in (sm[ky]["ef1"], sm[ky]["ef5"])) * 1.42)
    ax.legend(frameon=False, fontsize=6.6, loc="upper right")
    panel(ax, "C", "Early recognition")

    fig.tight_layout(w_pad=1.6)
    save(fig, out, "Figure_S1_Benchmark")


# ---------------------------------------------------------------------- TOC
def toc(data, structures, out):
    """ACS table-of-contents graphic, 3.25 x 1.75 in, legible at final size."""
    p = pd.read_csv(data / "Recognition_objectives_paired.csv").set_index("compound")
    d = pd.read_csv(data / "Series_electronic_descriptors.csv")
    prim = d[(d.dataset == "primary") & (d.scheme == "ESP")].set_index("substituent")
    eff = pd.read_csv(data / "Matched_effects_by_rule.csv")

    fig = plt.figure(figsize=(3.25, 1.75))
    # the canvas is not trimmed to content, so the axes are placed to fill it
    gs = fig.add_gridspec(1, 3, wspace=.90, left=.115, right=.975, top=.90, bottom=.16)
    tick = dict(labelsize=5.2, length=1.8, pad=1.5)

    # 1: prioritize
    ax = fig.add_subplot(gs[0, 0])
    ax.add_patch(Rectangle((90, 90), 14, 14, facecolor="#e6f0ea", edgecolor="none", zorder=1))
    ax.scatter(p.FP2_percentile, p.FP3_percentile, s=3.2, color=GREY, alpha=.65,
               linewidths=0, zorder=2)
    for k, c in (("SANC00370", TEAL), ("SANC01049", BLUE)):
        r = p.loc[k]
        ax.scatter(r.FP2_percentile, r.FP3_percentile, s=15, color=c,
                   edgecolor="white", linewidth=.35, zorder=3)
    ax.set_xticks([0, 50, 100])
    ax.set_yticks([0, 50, 100])
    ax.set_xlim(-4, 104)
    ax.set_ylim(-4, 104)
    ax.tick_params(**tick)
    ax.set_xlabel("FP-2 percentile", fontsize=5.4, labelpad=1)
    ax.set_ylabel("FP-3 percentile", fontsize=5.4, labelpad=1)
    ax.set_title("Prioritize", fontsize=6.8, fontweight="bold", color=INK, pad=3)

    # 2: tune - omega up, beta-carbon response down
    ax = fig.add_subplot(gs[0, 1])
    omega = np.array([prim.loc[s, "omega_eV"] for s in SUBS])
    fbeta = np.array([prim.loc[s, "fplus_beta"] for s in SUBS])
    x = np.arange(6)
    ax.bar(x, omega / omega.max(), width=.72,
           color=["#c8d4da"] * 5 + [ORANGE], zorder=2)
    ax.plot(x, fbeta / fbeta.max(), "o-", color=TEAL, lw=1.2, ms=2.4, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(SUBLAB, fontsize=4.8, rotation=45, ha="right")
    ax.set_yticks([0, .5, 1])
    ax.set_yticklabels(["0", "", "max"], fontsize=5.2)
    ax.tick_params(**tick)
    ax.set_ylim(0, 1.30)
    ax.text(.02, 1.22, r"$\omega$", color="#6d7a82", fontsize=6.0, ha="left")
    ax.text(5.0, 1.22, r"$f^{+}_{\beta}$", color=TEAL, fontsize=6.0, ha="right")
    ax.set_ylabel("scaled response", fontsize=5.4, labelpad=1)
    ax.set_title("Tune", fontsize=6.8, fontweight="bold", color=INK, pad=3)

    # 3: test
    ax = fig.add_subplot(gs[0, 2])
    row = eff[(eff.target == "FP3") & (eff.architecture == "glycosylated")
              & (eff.rule == "minimum_DockingScore")].iloc[0]
    ax.plot([0, 1], [row.OH_value, row.Cl_value], "o-", color=BLUE, lw=1.5, ms=3.4)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["4-OH", "4-Cl"], fontsize=5.2)
    ax.tick_params(**tick)
    ax.set_yticks([])
    ax.set_xlim(-.35, 1.35)
    ax.set_ylim(row.Cl_value - .35, row.OH_value + .35)
    ax.set_ylabel("FP-3 score", fontsize=5.4, labelpad=3)
    # the contrast is a pH 7.0 result and the graphic says so
    ax.text(0.98, 0.04, "pH 7.0", transform=ax.transAxes, fontsize=5.0,
            color="#6d7a82", ha="right", va="bottom")
    ax.set_title("Test", fontsize=6.8, fontweight="bold", color=INK, pad=3)

    save_exact(fig, out, "TOC_Graphic")


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--data", type=Path, required=True)
    p.add_argument("--structures", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--only", nargs="*")
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    stages = [("1", lambda: figure1(a.data, a.out)),
              ("2", lambda: figure2(a.data, a.structures, a.out)),
              ("3", lambda: figure3(a.data, a.structures, a.out)),
              ("4", lambda: figure4(a.data, a.out)),
              ("benchmark", lambda: figure_benchmark(a.data, a.out)),
              ("toc", lambda: toc(a.data, a.structures, a.out))]
    for name, fn in stages:
        if a.only and name not in a.only:
            continue
        fn()
        print("figure", name, "written")


if __name__ == "__main__":
    main()
