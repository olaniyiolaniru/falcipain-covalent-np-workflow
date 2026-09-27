"""Take the prepared states of the benchmark subset from the existing preparation.

The states are not re-prepared: they are the same LigPrep and Epik output the
full benchmark used, so the subset differs from it only in which compounds it
contains.
"""
import csv
from schrodinger import structure

keep = {r["id"] for r in csv.DictReader(open("bench_FP2_xp1k_labels.csv"))}
n = kept = 0
seen = set()
with structure.StructureWriter("prep.maegz") as out:
    for st in structure.StructureReader("../bench_fp2/prep.maegz"):
        n += 1
        if st.title in keep:
            out.append(st); kept += 1; seen.add(st.title)
print("states read %d, states kept %d, compounds represented %d of %d"
      % (n, kept, len(seen), len(keep)))
missing = sorted(keep - seen)
print("compounds with no prepared state:", len(missing))
