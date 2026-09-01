import csv, re
from rdkit import Chem

acc = set()
with open("sancdb_hits_accessible.csv") as fh:
    for r in csv.DictReader(fh):
        acc.add(r["id"])
accpar = set(re.search(r'(SANC\d+)', x).group(1) for x in acc)

supp = Chem.SDMolSupplier("sancdb_hits.sdf", removeHs=False)
w = Chem.SDWriter("accessible_185.sdf")
kept = set(); bad = 0
for m in supp:
    if m is None:
        bad += 1; continue
    name = m.GetProp("_Name")
    pid = re.search(r'(SANC\d+)', name)
    if pid and (name in acc or pid.group(1) in accpar):
        w.write(m); kept.add(pid.group(1))
w.close()
print("accessible parents wanted:", len(accpar))
print("records written (unique parents):", len(kept))
print("parents NOT found in 267 SDF:", sorted(accpar - kept))
print("unparseable SDF records skipped:", bad)
