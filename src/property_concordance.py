"""Reconcile the QikProp and ADMETlab property predictions.

The two platforms are reported together, so this script states their coverage
explicitly, puts the shared permeability endpoint on one scale, and reports
where the predictions agree and where they differ. ADMETlab returns apparent
Caco-2 permeability as log10(cm/s); QikProp returns QPPCaco in nm/s. The
conversion used here is

    QPPCaco_equivalent (nm/s) = 10**(caco2_logcm_s) * 1e7

because 1 cm/s = 1e7 nm/s. Predicted liability probabilities are reported as
model-predicted flags on the platform's own scale.

  python property_concordance.py --workbook <xlsx> --out ../exports
"""
from pathlib import Path
from argparse import ArgumentParser
import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

LIABILITY = ["hERG_blocker_prob", "DILI_prob", "Ames_prob", "Genotoxicity_prob",
             "Carcinogenicity_prob", "H_HT_prob", "CYP3A4_inh",
             "Aggregator_prob", "Reactive_prob"]


def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--workbook", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    admet = pd.read_excel(a.workbook, sheet_name="ADMETlab_predictions")
    qp = pd.read_excel(a.workbook, sheet_name="QikProp_shortlist")
    qp["compound"] = qp.molecule.astype(str).str.replace("_minRM1.pdb", "", regex=False)

    admet["caco2_nm_s_equivalent"] = (10 ** admet.caco2_logcm_s) * 1e7
    rows = []
    for r in admet.itertuples():
        block = qp[qp.compound.eq(r.SANC_id)]
        if block.empty:
            continue
        rows.append(dict(
            compound=r.SANC_id,
            qikprop_records=int(len(block)),
            QPPCaco_min_nm_s=float(block.QPPCaco.min()),
            QPPCaco_max_nm_s=float(block.QPPCaco.max()),
            ADMETlab_caco2_log_cm_s=float(r.caco2_logcm_s),
            ADMETlab_caco2_nm_s_equivalent=round(float(r.caco2_nm_s_equivalent), 2),
            same_order_of_magnitude=bool(
                abs(np.log10(max(r.caco2_nm_s_equivalent, 1e-9))
                    - np.log10(max(block.QPPCaco.mean(), 1e-9))) < 1.0),
            QikProp_oral_absorption_min=float(block.PercentHumanOralAbsorption.min()),
            QikProp_oral_absorption_max=float(block.PercentHumanOralAbsorption.max()),
            QikProp_RuleOfFive=int(block.RuleOfFive.max())))
    frame = pd.DataFrame(rows)
    frame.to_csv(a.out / "Property_platform_concordance.csv", index=False)

    flags = admet[["SANC_id"] + [c for c in LIABILITY if c in admet.columns]].copy()
    flags.to_csv(a.out / "ADMETlab_predicted_flags.csv", index=False)

    rho = None
    if len(frame) >= 3:
        rho = float(spearmanr(frame.QPPCaco_min_nm_s, frame.ADMETlab_caco2_nm_s_equivalent).statistic)

    summary = dict(
        admetlab_compounds=sorted(admet.SANC_id.tolist()),
        admetlab_coverage=("ADMETlab predictions were obtained for %d of the %d compounds "
                           "carried through the property stage" % (len(admet),
                                                                   qp.compound.nunique())),
        qikprop_coverage=("QikProp predictions were obtained for all %d compounds carried "
                          "through the property stage" % qp.compound.nunique()),
        shared_endpoint="apparent Caco-2 permeability",
        caco2_rank_agreement_spearman=rho,
        compounds_same_order_of_magnitude=int(frame.same_order_of_magnitude.sum()),
        compounds_compared=int(len(frame)),
        liability_endpoints=LIABILITY,
        reporting_rule=("ADMETlab outputs are reported as model-predicted liability flags on "
                        "the platform scale; QikProp outputs are reported as ranges over every "
                        "prepared-state record"))
    (a.out / "Property_platform_concordance.json").write_text(json.dumps(summary, indent=2),
                                                              encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print(frame.to_string(index=False))


if __name__ == "__main__":
    main()
