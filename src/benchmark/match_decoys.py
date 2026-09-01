#!/usr/bin/env python3
"""DUD-E-style property-matched decoy selection.

For each active, pick N decoys from the background pool that match six physicochemical
properties within progressively widened DUD-E windows, are not themselves falcipain
actives, and are assigned to only one active. Topological dissimilarity to the active
is approximated by requiring a different Bemis-Murcko scaffold hash where available
(computed separately with Schrodinger); here we enforce property matching + uniqueness.

DUD-E matching properties: MW, cLogP, #HBD, #HBA, #rotatable bonds, net charge.
"""
import csv, json, os

BM = r"data/validation"
os.chdir(BM)
N_PER_ACTIVE = 50
POTENCY_CUT_nM = 10000        # actives defined as <=10 uM

pool = json.load(open("decoy_pool.json"))
aprops = json.load(open("actives_props.json"))

def charge_num(species):
    return {"ACID": -1, "BASE": 1, "NEUTRAL": 0, "ZWITTERION": 0}.get(species, 0)

def num(x):
    try: return float(x)
    except: return None

# Build pool feature table (drop entries with missing core props)
pool_feats = {}
for cid, p in pool.items():
    mw, lp = num(p.get("mw")), num(p.get("alogp"))
    hbd, hba, rtb = num(p.get("hbd")), num(p.get("hba")), num(p.get("rtb"))
    if None in (mw, lp, hbd, hba, rtb): continue
    pool_feats[cid] = (mw, lp, hbd, hba, rtb, charge_num(p.get("charge")), p.get("smi"))

# DUD-E widening windows (mw, logp, hbd, hba, rtb)
WINDOWS = [(20,0.5,0,0,1),(35,1.0,1,1,2),(50,1.5,1,2,3),(75,2.0,2,3,4),(100,3.0,3,4,5)]

for grp in ("FP2", "FP3"):
    actives = [r for r in csv.DictReader(open(f"actives_{grp}_classified.csv"))
               if num(r["gmean_nM"]) is not None and num(r["gmean_nM"]) <= POTENCY_CUT_nM]
    used = set()
    assignments = []
    unmatched = 0
    for a in actives:
        ap = aprops.get(a["molecule_chembl_id"])
        if not ap:
            unmatched += 1; continue
        amw, alp = num(ap.get("mw")), num(ap.get("alogp"))
        ahbd, ahba, artb = num(ap.get("hbd")), num(ap.get("hba")), num(ap.get("rtb"))
        ach = charge_num(ap.get("charge"))
        if None in (amw, alp, ahbd, ahba, artb):
            unmatched += 1; continue
        picked = []
        for (dmw, dlp, dhbd, dhba, drtb) in WINDOWS:
            for cid, (mw, lp, hbd, hba, rtb, ch, smi) in pool_feats.items():
                if cid in used or cid in aprops: continue
                if abs(mw-amw) <= dmw and abs(lp-alp) <= dlp and abs(hbd-ahbd) <= dhbd \
                   and abs(hba-ahba) <= dhba and abs(rtb-artb) <= drtb and ch == ach:
                    picked.append((cid, smi)); used.add(cid)
                    if len(picked) >= N_PER_ACTIVE: break
            if len(picked) >= N_PER_ACTIVE: break
        for cid, smi in picked:
            assignments.append((cid, smi, a["molecule_chembl_id"]))
    with open(f"decoys_{grp}.csv", "w", newline="") as o:
        w = csv.writer(o); w.writerow(["decoy_chembl_id","smiles","matched_active"])
        w.writerows(assignments)
    print(f"{grp}: {len(actives)} actives, {len(assignments)} decoys "
          f"({len(assignments)/max(len(actives),1):.1f}/active), {unmatched} actives lacked props")
