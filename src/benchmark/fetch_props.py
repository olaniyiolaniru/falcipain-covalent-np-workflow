#!/usr/bin/env python3
"""Resumable fetch of ChEMBL molecule_properties for actives + a random drug-like
background decoy pool. Saves incrementally so it survives interruptions."""
import csv, json, os, time, urllib.request, random

BM = r"data/validation"
os.chdir(BM)

def fetch(url, tries=5):
    for _ in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return json.load(r)
        except Exception:
            time.sleep(3)
    return None

# ---- 1. actives properties (resumable) ----
props = {}
if os.path.exists("actives_props.json"):
    props = json.load(open("actives_props.json"))
ids = []
for grp in ("FP2", "FP3"):
    for r in csv.DictReader(open(f"actives_{grp}_classified.csv")):
        ids.append(r["molecule_chembl_id"])
ids = sorted(set(ids))
todo = [i for i in ids if i not in props]
print(f"actives: {len(ids)} total, {len(todo)} to fetch", flush=True)
for i in range(0, len(todo), 30):
    chunk = todo[i:i+30]
    d = fetch("https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_chembl_id__in="
              + ",".join(chunk) + "&limit=30")
    if not d:
        print("chunk fail", i, flush=True); continue
    for m in d.get("molecules", []):
        mp = m.get("molecule_properties") or {}
        props[m["molecule_chembl_id"]] = {
            "smi": (m.get("molecule_structures") or {}).get("canonical_smiles"),
            "mw": mp.get("full_mwt"), "alogp": mp.get("alogp"), "hbd": mp.get("hbd"),
            "hba": mp.get("hba"), "rtb": mp.get("rtb"),
            "charge": mp.get("molecular_species"), "arom": mp.get("aromatic_rings")}
    json.dump(props, open("actives_props.json", "w"))
    print(f"actives props {len(props)}/{len(ids)}", flush=True)
print("ACTIVES_PROPS_DONE", flush=True)

# ---- 2. background decoy pool: random drug-like ChEMBL molecules with properties ----
pool = {}
if os.path.exists("decoy_pool.json"):
    pool = json.load(open("decoy_pool.json"))
TARGET_POOL = 20000            # universe to draw property-matched decoys from
# ChEMBL molecule endpoint filtered to drug-like space; page via offset over a shuffled range.
offsets = list(range(0, 2_000_000, 1000))
random.seed(42); random.shuffle(offsets)
oi = 0
while len(pool) < TARGET_POOL and oi < len(offsets):
    off = offsets[oi]; oi += 1
    url = ("https://www.ebi.ac.uk/chembl/api/data/molecule.json?"
           "molecule_properties__full_mwt__gte=250&molecule_properties__full_mwt__lte=600&"
           "molecule_properties__alogp__gte=-2&molecule_properties__alogp__lte=6&"
           f"limit=1000&offset={off}")
    d = fetch(url)
    if not d:
        continue
    for m in d.get("molecules", []):
        ms = m.get("molecule_structures") or {}
        mp = m.get("molecule_properties") or {}
        smi = ms.get("canonical_smiles")
        if not smi:
            continue
        pool[m["molecule_chembl_id"]] = {
            "smi": smi, "mw": mp.get("full_mwt"), "alogp": mp.get("alogp"),
            "hbd": mp.get("hbd"), "hba": mp.get("hba"), "rtb": mp.get("rtb"),
            "charge": mp.get("molecular_species"), "arom": mp.get("aromatic_rings")}
    json.dump(pool, open("decoy_pool.json", "w"))
    print(f"decoy pool {len(pool)}/{TARGET_POOL}", flush=True)
print("DECOY_POOL_DONE", len(pool), flush=True)
