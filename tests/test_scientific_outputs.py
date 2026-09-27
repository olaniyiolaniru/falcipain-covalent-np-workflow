"""Checks on the released scientific outputs.

These assert the values the article reports, against the released tables. A
failure here means the release and the article disagree.
"""
from pathlib import Path
import json

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"


def load_json(name):
    path = DATA / name
    if not path.is_file():
        pytest.skip("not in this build: " + name)
    return json.loads(path.read_text(encoding="utf-8"))


def load_csv(name):
    path = DATA / name
    if not path.is_file():
        pytest.skip("not in this build: " + name)
    return pd.read_csv(path)


def test_library_triage_counts_reconcile():
    s = load_json("Library_provenance_summary.json")
    assert s["records"] == 1012
    assert s["parsed"] + s["parse_failures"] == s["records"]
    assert s["motif_bearing_records"] - s["duplicate_motif_hits"] == s["unique_motif_bearing"]
    assert s["unique_motif_bearing"] - s["failed_selection_rule"] == s["retained_parents"]
    assert s["retained_parents"] == 185


def test_paired_parents_and_pareto_front():
    s = load_json("Recognition_ranking_summary.json")
    assert s["paired_parents"] == 174
    assert s["pareto_front"] == ["SANC00370", "SANC01049"]


def test_every_composite_pair_agrees():
    s = load_json("Recognition_ranking_summary.json")
    assert min(s["composite_agreement_spearman"].values()) >= 0.90


def test_size_scaling_reorders_the_list_completely():
    s = load_json("Recognition_ranking_summary.json")
    for key, value in s["top10_overlap"].items():
        if "size_scaled" in key:
            assert value == 0


def test_both_receptors_carry_the_standardized_dyad():
    rec = load_csv("Receptor_preparation_record.csv")
    reported = rec[rec.stage.isin(["primary recognition screen",
                                   "matched-analogue comparison",
                                   "retrospective enrichment benchmark"])]
    assert len(reported) >= 4
    assert set(reported.catalytic_dyad_convention) == {"standardized neutral dyad"}
    alt = rec[rec.stage.eq("catalytic-dyad sensitivity experiment")]
    assert set(alt.catalytic_dyad_convention) == {"ion-pair dyad"}


def test_dyad_experiment_changes_only_the_dyad():
    d = load_json("Dyad_sensitivity.json")
    c = d["comparison"]
    assert c["heavy_atom_rmsd_A"] == 0.0
    assert c["max_heavy_atom_displacement_A"] == 0.0
    assert c["heavy_atoms_compared"] == 1877
    assert round(c["shift_DockingScore"], 2) == 3.89


def test_the_two_dyad_receptors_carry_the_same_net_charge():
    """The experiment moves one proton, so the net charge cannot change.

    Both released records describe the same two receptors and are built from
    the pose-viewer files of the two docking runs, so they are required to
    agree with each other as well as with the experiment.
    """
    d = load_json("Dyad_sensitivity.json")
    std, alt = d["standardized"], d["alternative"]
    assert std["receptor_formal_charge"] == alt["receptor_formal_charge"] == -11
    assert std["receptor_atoms"] == alt["receptor_atoms"] == 3666
    assert std["receptor_heavy_atoms"] == alt["receptor_heavy_atoms"] == 1877
    assert std["receptor_residues"] == alt["receptor_residues"] == 240
    # the proton moves off the cysteine and onto the histidine
    assert std["dyad"]["res42"]["formal_charge"] == 0
    assert std["dyad"]["res174"]["formal_charge"] == 0
    assert alt["dyad"]["res42"]["formal_charge"] == -1
    assert alt["dyad"]["res174"]["formal_charge"] == 1

    r = load_json("FP2_dyad_receptor_difference.json")
    assert r["formal_charge_a"] == r["formal_charge_b"] == -11
    assert r["heavy_atoms_a"] == r["heavy_atoms_b"] == 1877
    assert r["residues_a"] == r["residues_b"] == 240
    assert r["differing_residues"] == 2
    assert r["file_a"] != r["file_b"]


def test_dyad_experiment_reproduces():
    d = load_json("Dyad_sensitivity.json")
    if "reproduction_check" not in d:
        pytest.skip("reproduction run not in this build")
    assert d["reproduction_check"]["reproduces_exactly"]


def test_nitro_relocates_the_response_in_both_partitions():
    q = load_json("Quantum_analysis_summary.json")
    esp = q["nitro_fragment_response"]["check"]["ESP_para_substituent"]
    mul = q["nitro_fragment_response"]["check"]["Mulliken_para_substituent"]
    assert esp > 0.60 and mul > 0.60
    assert abs(esp - mul) < 0.02
    assert q["omega_NO2_over_OH"] > 1.6
    assert q["fplus_beta_NO2_over_OH"] < 0.2


def test_charge_normalisation():
    q = load_json("Quantum_analysis_summary.json")
    assert q["charge_normalisation_max_error"] < 1e-9


def test_fp3_glycosylated_contrast_keeps_its_sign_under_every_protocol():
    """The magnitude is invariant within XP; across SP only the sign transfers."""
    eff = load_csv("Matched_effects_by_rule.csv")
    block = eff[(eff.target == "FP3") & (eff.architecture == "glycosylated")]
    assert len(block) >= 2
    assert (block.Cl_minus_OH < 0).all()
    xp = block[~block.rule.str.contains("SPes4")]
    assert xp.Cl_minus_OH.max() - xp.Cl_minus_OH.min() < 1e-6
    assert round(xp.Cl_minus_OH.iloc[0], 2) == -0.77


def test_no_halogen_bond_geometry_anywhere():
    hal = load_csv("Matched_halogen_bond_check.csv")
    assert int(hal.halogen_bond.sum()) == 0
    assert hal.C_Cl_acceptor_angle_deg.max() < 155.0


def test_pose_ensembles_redistribute_and_stay_connected():
    """The reported 8.7 A is a median over cross-ensemble pose pairs.

    The closest family centroids lie inside the clustering cutoff, so the two
    ensembles share a family and the article does not describe them as disjoint.
    """
    fam = load_csv("Matched_pose_family_summary.csv")
    fp3 = fam[(fam.target == "FP3") & (fam.architecture == "glycosylated")].iloc[0]
    assert round(fp3.median_pairwise_core_rmsd_A, 1) == 8.7
    assert fp3.min_pairwise_core_rmsd_A < fp3.median_pairwise_core_rmsd_A
    assert fp3.OH_clusters >= 1 and fp3.Cl_clusters >= 1
    c = load_json("Matched_pose_clustering.json")
    s3 = [x for x in c["summaries"]
          if x["target"] == "FP3" and x["architecture"] == "glycosylated"][0]
    assert s3["closest_cluster_centroid_rmsd_A"] < c["cutoff_A"]
    assert s3["shared_family_at_cutoff"] is True


def test_benchmark_is_reported_on_the_reported_receptor():
    b = load_json("Benchmark_summary.json")
    assert "standardized" in b["receptor"]
    for name in ("all_actives", "covalent_only"):
        m = b[name]
        assert m["n"] > 1000 and m["actives"] > 50
        assert 0.0 <= m["roc_auc"] <= 1.0
        assert m["roc_ci"][0] <= m["roc_auc"] <= m["roc_ci"][1]
    # the ranking carries a resolvable signal: the interval excludes random
    assert b["all_actives"]["roc_ci"][0] > 0.5
    # and the covalent actives enrich in the leading 5% above chance
    assert b["covalent_only"]["ef5_ci"][0] > 1.0


def test_protocol_matched_benchmark_is_xp_only():
    """The headline benchmark must be the protocol that produces the ranking.

    Every score in it comes from Glide XP over a prepared state, as in the
    recognition screen. A funnel score appearing here would mean the article was
    quoting a different calculation from the one it reports.
    """
    b = load_json("Benchmark_summary_FP2_XP.json")
    assert set(b["scored_by_stage"]) == {"XP"}
    assert "Glide XP" in b["receptor"]
    assert b["all_actives"]["actives"] == 263
    # the interval excludes random over the full active set and over the covalent subset
    assert b["all_actives"]["roc_ci"][0] > 0.5
    assert b["covalent_only"]["roc_ci"][0] > 0.5
    assert b["covalent_only"]["roc_auc"] > b["all_actives"]["roc_auc"]


def test_xp_summaries_describe_an_xp_ranking():
    """The recorded ranking rule must match the stages that supplied the scores.

    A protocol-matched run is scored entirely at XP, so its record must not
    describe a funnel. This catches metadata drifting away from the calculation.
    """
    for tag in ("FP2_XP", "FP3_XP"):
        b = load_json("Benchmark_summary_%s.json" % tag)
        assert set(b["scored_by_stage"]) == {"XP"}
        rule = b["ranking_rule"]
        assert "XP" in rule
        assert "HTVS" not in rule and "SP score" not in rule and "funnel" not in rule
    # and the funnel record must still describe the funnel
    f = load_json("Benchmark_summary.json")
    assert set(f["scored_by_stage"]) == {"HTVS", "SP"}
    assert "funnel" in f["ranking_rule"]


def test_precision_and_decoy_pool_effects_are_separated():
    """Each comparison changes one variable, so each term has a size.

    The precision term is measured on the compounds both protocols scored and
    the decoy-pool term on one protocol across two decoy sets, so neither is a
    mixture of the two.
    """
    j = load_json("Benchmark_precision_comparison.json")
    assert j["xp_common"]["n"] == j["funnel_common"]["n"] == j["compounds_scored_by_both"]
    assert j["funnel_sampled_decoys"]["actives"] == j["funnel_full_pool"]["actives"]
    assert abs(j["precision_effect"]["delta_roc_auc"]) < 0.02
    assert abs(j["decoy_pool_effect"]["delta_roc_auc"]) < 0.02
    # the areas agree while the orderings do not
    assert j["precision_effect"]["spearman_rho"] < 0.6


def test_falcipain3_benchmark_keeps_its_stated_interval():
    """The falcipain-3 set is quoted with the interval its size supports.

    Twenty-two actives give an interval about three times the width of the
    falcipain-2 one, and the article quotes it as such. This test fails if the
    released data ever stop matching that description.
    """
    b = load_json("Benchmark_summary_FP3_XP.json")
    assert set(b["scored_by_stage"]) == {"XP"}
    m = b["all_actives"]
    assert m["actives"] == 22
    assert m["roc_ci"][0] < 0.5 < m["roc_ci"][1]


def test_dyad_comparison_uses_one_common_scored_set():
    j = load_json("Benchmark_dyad_comparison.json")
    assert j["scored_by_both"] <= min(j["scored_by_standard"],
                                      j["scored_by_alternative"])
    for name in ("all_actives", "covalent_only"):
        assert (j["standardized"][name]["n"] == j["alternative"][name]["n"])
        assert (j["standardized"][name]["actives"]
                == j["alternative"][name]["actives"])


def test_covalent_stage_statuses_are_separate():
    cov = load_csv("Covalent_stage_records.csv")
    assert len(cov) == 64
    assert set(cov.prime_status) == {"completed"}
    assert int((cov.postreaction_status == "completed").sum()) == 62


def test_covalent_receptors_use_the_thiolate_by_design():
    rec = load_csv("Covalent_receptor_record.csv")
    assert all("thiolate" in s for s in rec.catalytic_cys)


def test_panel_is_defined_by_the_stated_rule():
    s = load_json("Panel_property_summary.json")
    paired = load_csv("Recognition_objectives_paired.csv")
    expected = set(paired[paired.worse_percentile >= 90].compound)
    assert set(s["panel_members"]) == expected
    assert not s["members_without_qikprop"]


def test_property_panel_regenerates_all_ten_members():
    """Run the properties stage and check the regenerated table, not the shipped one.

    The tenth panel member is covered by a separate QikProp run, so this test
    fails if that input stops being released or stops being passed to the stage.
    """
    import subprocess
    import sys
    import tempfile

    expected = set(load_json("Panel_property_summary.json")["panel_members"])
    assert len(expected) == 10
    with tempfile.TemporaryDirectory() as tmp:
        r = subprocess.run(
            [sys.executable, str(ROOT / "src" / "panel_properties.py"),
             "--paired", str(DATA / "Recognition_objectives_paired.csv"),
             "--workbook", str(DATA / "Panel_property_predictions.xlsx"),
             "--names", str(DATA / "Library_compound_names.csv"),
             "--extra", str(DATA / "Panel_property_extra_qikprop.csv"),
             "--out", tmp],
            capture_output=True, text=True)
        assert r.returncode == 0, r.stderr[-2000:]
        got = json.loads((Path(tmp) / "Panel_property_summary.json")
                         .read_text(encoding="utf-8"))
    assert set(got["panel_members"]) == expected
    assert not got["members_without_qikprop"]


def test_manifest_lists_every_reported_stage():
    man = pd.read_csv(ROOT / "MANIFEST.csv")
    assert set(man.manuscript_section) >= {"2.1", "2.2", "2.3", "2.4", "2.5"}
    absent = man[man.output_sha256 == "not present in this build"]
    assert absent.empty, ("the release is missing the output for: "
                          + ", ".join(absent.output_record))


def test_requirements_cover_every_imported_package():
    """Anything a released script imports must be installable from requirements.

    A package that is imported but unlisted makes the documented install
    incomplete, and the failure only appears when someone else runs the code.
    """
    import ast
    STDLIB = set("""argparse ast collections csv datetime functools gzip hashlib io
        itertools json math os pathlib random re shutil subprocess sys tempfile time
        traceback typing warnings zipfile statistics textwrap glob copy"""
        .split())
    local = {f.stem for f in ROOT.rglob("*.py")}
    used = set()
    for f in ROOT.rglob("*.py"):
        try:
            tree = ast.parse(f.read_text(encoding="utf-8"))
        except SyntaxError:
            continue
        for n in ast.walk(tree):
            if isinstance(n, ast.Import):
                used.update(a.name.split(".")[0] for a in n.names)
            elif isinstance(n, ast.ImportFrom) and n.module and n.level == 0:
                used.add(n.module.split(".")[0])
    # the Schrodinger modules are licensed and installed with the suite
    external = {m for m in used if m not in STDLIB and m not in local
                and m != "schrodinger"}
    alias = {"PIL": "pillow", "docx": "python-docx"}
    listed = {l.split(">")[0].split("=")[0].strip().lower()
              for l in (ROOT / "requirements.txt").read_text(encoding="utf-8").splitlines()
              if l.strip() and not l.startswith("#")}
    missing = sorted(m for m in external
                     if alias.get(m, m).lower() not in listed)
    assert not missing, "imported but not in requirements.txt: %s" % missing


def test_checksums_cover_the_tree_and_are_readable_by_sha256sum():
    """Every file is listed, and the list is in the format the tool expects.

    A carriage return at the end of a line becomes part of the file name, so
    "sha256sum -c" then fails to open every file it is asked to check.
    """
    raw = (ROOT / "SHA256SUMS.txt").read_bytes()
    assert chr(13).encode() not in raw, "SHA256SUMS.txt must use LF line endings"
    sums = raw.decode("utf-8").splitlines()
    listed = {line.split("  ", 1)[1] for line in sums if line.strip()}
    assert all(not n.strip() != n for n in listed)
    present = {f.relative_to(ROOT).as_posix() for f in ROOT.rglob("*")
               if f.is_file() and f.name != "SHA256SUMS.txt"
               and "__pycache__" not in f.parts and ".pytest_cache" not in f.parts
               and "outputs" not in f.parts and ".git" not in f.parts}
    assert present <= listed
