"""Binding-mode family analysis of the matched OH/Cl pose ensemble.

The matched set was docked under three protocols: the reported single-pose XP
run, a wider XP run retaining up to ten poses per prepared state, and an SP run
with four-fold enhanced sampling. This script uses every pose those runs
returned to ask whether the 4-OH and 4-Cl members occupy the same binding-mode
family, rather than resting the structural interpretation on one selected pose.

Method. The shared para-substituted cinnamoyl ester is located in each ligand by
one SMARTS query, which returns the same twelve heavy atoms in the same order in
every member, so poses are compared on an atom-mapped common substructure with
no ligand superposition. Poses of one compound at one receptor are clustered by
average-linkage on the pairwise mapped-core RMSD at a stated cutoff. For each
cluster the script reports population, best score, and the mapped-core
displacement to every cluster of the partner compound, so that the reported
OH-to-Cl core displacement can be read against the full distribution rather than
as a single number.

  run.exe python3 pose_ensemble_clustering.py --runs <dirs...> --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import csv
import itertools
import json
import math

import numpy as np
from schrodinger import structure
from schrodinger.structutils import analyze

CORE_SMARTS = "[$([OX2H1]),$([Cl])]c1ccc(cc1)C=CC(=O)O"
MEMBERS = {"Parent_SANC00867": ("glycosylated", "OH"),
           "AnalogueA": ("glycosylated", "Cl"),
           "AnalogueB": ("aglycone_ester", "OH"),
           "OPT1": ("aglycone_ester", "Cl")}
PAIRS = [("glycosylated", "Parent_SANC00867", "AnalogueA"),
         ("aglycone_ester", "AnalogueB", "OPT1")]


def core_coords(st):
    hits = analyze.evaluate_smarts(st, CORE_SMARTS, unique_sets=True)
    if not hits:
        return None
    if len(hits) > 1:
        hits = sorted(hits, key=lambda h: (st.atom[h[0]].element != "Cl", h[0]))
    idx = hits[0]
    return np.array([[st.atom[i].x, st.atom[i].y, st.atom[i].z] for i in idx])


def rmsd(a, b):
    n = min(len(a), len(b))
    return float(np.sqrt(np.mean(np.sum((a[:n] - b[:n]) ** 2, axis=1))))


def average_linkage(coords, cutoff):
    """Agglomerative average-linkage clustering on the mapped-core RMSD."""
    clusters = [[i] for i in range(len(coords))]
    if len(coords) < 2:
        return clusters
    d = np.zeros((len(coords), len(coords)))
    for i, j in itertools.combinations(range(len(coords)), 2):
        d[i, j] = d[j, i] = rmsd(coords[i], coords[j])
    while len(clusters) > 1:
        best, pair = None, None
        for x, y in itertools.combinations(range(len(clusters)), 2):
            m = float(np.mean([d[i, j] for i in clusters[x] for j in clusters[y]]))
            if best is None or m < best:
                best, pair = m, (x, y)
        if best is None or best > cutoff:
            break
        x, y = pair
        clusters[x] = clusters[x] + clusters[y]
        clusters.pop(y)
    return sorted(clusters, key=len, reverse=True)


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--runs", nargs="+", required=True,
                   help="pose-viewer files as TARGET:RUN:PATH")
    p.add_argument("--cutoff", type=float, default=2.0,
                   help="mapped-core RMSD cutoff in angstrom (default 2.0)")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    poses = {}
    pose_rows = []
    for spec in a.runs:
        target, run, path = spec.split(":", 2)
        path = Path(path)
        if not path.is_file():
            print("missing", path)
            continue
        for n, st in enumerate(structure.StructureReader(str(path)), 1):
            if n == 1 and st.atom_total > 1000:
                continue  # receptor entry of a pose viewer
            title = st.title
            if title not in MEMBERS:
                continue
            c = core_coords(st)
            if c is None:
                continue
            key = (target, title)
            poses.setdefault(key, []).append(
                dict(run=run, record=n, coords=c,
                     variant=st.property.get("s_lp_Variant", title),
                     dock=st.property.get("r_i_docking_score"),
                     glide=st.property.get("r_i_glide_gscore")))

    cluster_rows, family_rows = [], []
    clusters_by_key = {}
    for key, plist in sorted(poses.items()):
        target, title = key
        coords = [x["coords"] for x in plist]
        groups = average_linkage(coords, a.cutoff)
        clusters_by_key[key] = (plist, groups)
        for ci, members in enumerate(groups, 1):
            scores = [plist[i]["dock"] for i in members if plist[i]["dock"] is not None]
            centroid = np.mean([coords[i] for i in members], axis=0)
            spread = float(np.mean([rmsd(coords[i], centroid) for i in members]))
            cluster_rows.append(dict(
                target=target, compound=title,
                architecture=MEMBERS[title][0], para_group=MEMBERS[title][1],
                cluster=ci, poses=len(members),
                fraction_of_ensemble=round(len(members) / len(plist), 3),
                best_DockingScore=min(scores) if scores else None,
                mean_DockingScore=float(np.mean(scores)) if scores else None,
                within_cluster_core_rmsd_A=round(spread, 3),
                runs=";".join(sorted({plist[i]["run"] for i in members})),
                contains_reported_pose=any(plist[i]["run"] == "reported" for i in members)))

    for target in sorted({t for t, _ in poses}):
        for arch, oh, cl in PAIRS:
            ka, kb = (target, oh), (target, cl)
            if ka not in clusters_by_key or kb not in clusters_by_key:
                continue
            pa, ga = clusters_by_key[ka]
            pb, gb = clusters_by_key[kb]
            best = None
            for i, ma in enumerate(ga, 1):
                ca = np.mean([pa[k]["coords"] for k in ma], axis=0)
                for j, mb in enumerate(gb, 1):
                    cb = np.mean([pb[k]["coords"] for k in mb], axis=0)
                    d = rmsd(ca, cb)
                    family_rows.append(dict(
                        target=target, architecture=arch,
                        OH_cluster=i, OH_poses=len(ma),
                        Cl_cluster=j, Cl_poses=len(mb),
                        centroid_core_rmsd_A=round(d, 3),
                        OH_best_DockingScore=min(
                            (pa[k]["dock"] for k in ma if pa[k]["dock"] is not None), default=None),
                        Cl_best_DockingScore=min(
                            (pb[k]["dock"] for k in mb if pb[k]["dock"] is not None), default=None)))
                    if best is None or d < best:
                        best = d
            # displacement distribution over every OH/Cl pose pair
            all_d = [rmsd(x["coords"], y["coords"]) for x in pa for y in pb]
            family_rows[-1]["_"] = None
            fam = dict(target=target, architecture=arch,
                       OH_poses=len(pa), Cl_poses=len(pb),
                       OH_clusters=len(ga), Cl_clusters=len(gb),
                       min_pairwise_core_rmsd_A=round(min(all_d), 3),
                       median_pairwise_core_rmsd_A=round(float(np.median(all_d)), 3),
                       max_pairwise_core_rmsd_A=round(max(all_d), 3),
                       closest_cluster_centroid_rmsd_A=round(best, 3),
                       shared_family_at_cutoff=bool(min(all_d) <= a.cutoff))
            cluster_rows.append({})  # spacer removed below
            cluster_rows.pop()
            family_rows.append(dict(summary=json.dumps(fam)))

    summaries = [json.loads(r["summary"]) for r in family_rows if "summary" in r]
    family_rows = [r for r in family_rows if "summary" not in r]
    for r in family_rows:
        r.pop("_", None)

    def write(name, rows):
        if not rows:
            return
        with (a.out / name).open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    write("Matched_pose_clusters.csv", cluster_rows)
    write("Matched_pose_family_comparison.csv", family_rows)
    write("Matched_pose_family_summary.csv", summaries)

    out = dict(cutoff_A=a.cutoff, core_smarts=CORE_SMARTS,
               poses_analysed=sum(len(v) for v in poses.values()),
               compounds=len(poses), summaries=summaries,
               method=("average-linkage clustering on atom-mapped cinnamoyl-core RMSD; "
                       "no ligand superposition"))
    (a.out / "Matched_pose_clustering.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
