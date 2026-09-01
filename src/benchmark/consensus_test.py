#!/usr/bin/env python3
"""Test whether reactivity-weighting improves retrospective enrichment over docking alone.
Classifies all benchmark ligands by the paper's warhead SMARTS tiers, then compares
docking-alone vs reactivity-alone vs consensus (z(-dock)+z(tier)) enrichment with bootstrap CIs.
Run with Schrodinger python (RDKit)."""
import csv, math, random, statistics as st, io, sys, os
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding='utf-8',errors='replace')

BM=r"data/validation"
SCORES=r"data/validation/scores_FP2_full.csv"

WARHEADS={  # (smarts, tier) - from michael_filter.py
 "ab_unsat_carbonyl":("[CX3]=[CX3][CX3]=[OX1]",1),
 "exo_methylene_carbonyl":("[CH2]=[CX3][CX3]=[OX1]",1),
 "vinyl_nitrile":("[CX3]=[CX3][CX2]#[NX1]",1),
 "vinyl_sulfone":("[CX3]=[CX3][SX4](=[OX1])(=[OX1])",2),
 "maleimide":("O=C1C=CC(=O)N1",2),
 "para_quinone":("O=C1C=CC(=O)C=C1",2),
 "ortho_quinone":("O=C1C=CC=CC1=O",2),
}
Q=[(Chem.MolFromSmarts(s),t) for s,t in WARHEADS.values()]
_lfc=rdMolStandardize.LargestFragmentChooser(); _norm=rdMolStandardize.Normalizer()
def tier(smiles):
    if not smiles: return 0
    m=Chem.MolFromSmiles(smiles)
    if m is None: return 0
    try: m=_norm.normalize(_lfc.choose(m))
    except Exception: pass
    t=0
    for q,tr in Q:
        if q is not None and m.HasSubstructMatch(q): t=max(t,tr)
    return t

# --- SMILES for every benchmark id ---
smi={}
for r in csv.DictReader(open(BM+"/actives_FP2_classified.csv")):
    smi[r["molecule_chembl_id"]]=r["canonical_smiles"]
for r in csv.DictReader(open(BM+"/decoys_FP2.csv")):
    smi[r["decoy_chembl_id"]]=r["smiles"]

# --- scored set: id,label,score,class ---
rows=[]
for r in csv.DictReader(open(SCORES)):
    rid=r["id"]; lab=int(r["label"]); score=float(r["score"]); cls=r["class"]
    rows.append([rid,lab,score,cls,tier(smi.get(rid,""))])
print(f"scored ligands: {len(rows)}  actives={sum(1 for x in rows if x[1]==1)}  "
      f"covalent-actives={sum(1 for x in rows if x[1]==1 and x[3]=='covalent')}")
# warhead prevalence
for grp,pred in [("actives",lambda x:x[1]==1),("decoys",lambda x:x[1]==0)]:
    sub=[x for x in rows if pred(x)]
    t12=sum(1 for x in sub if x[4]>=1)
    print(f"  {grp}: {t12}/{len(sub)} bear a warhead ({100*t12/len(sub):.0f}%)")

# --- scoring vectors (higher = better) ---
def zscore(vals):
    m=st.mean(vals); s=st.pstdev(vals) or 1.0
    return [(v-m)/s for v in vals]
dock=[-x[2] for x in rows]         # higher better (more negative dock score)
react=[float(x[4]) for x in rows]  # tier
zc=zscore(dock); zr=zscore(react)
consensus=[a+b for a,b in zip(zc,zr)]

def auc(scores, labels):  # Mann-Whitney, tie-aware (0.5 for ties)
    pos=[s for s,l in zip(scores,labels) if l==1]
    neg=[s for s,l in zip(scores,labels) if l==0]
    if not pos or not neg: return float('nan')
    # rank-based
    allv=sorted(pos+neg)
    import bisect
    # average-rank U
    U=0.0
    negs=sorted(neg)
    for p in pos:
        lo=bisect.bisect_left(negs,p); hi=bisect.bisect_right(negs,p)
        U+= lo + 0.5*(hi-lo)   # negs strictly below + half of ties
    return U/(len(pos)*len(neg))

def ef(scores,labels,frac):
    order=sorted(range(len(scores)),key=lambda i:-scores[i])
    n=max(1,int(round(frac*len(order)))); P=sum(labels)
    hits=sum(labels[i] for i in order[:n])
    return (hits/n)/(P/len(labels)) if P else float('nan')

def boot_auc(scores,labels,B=1000,seed=1):
    random.seed(seed)
    ia=[i for i,l in enumerate(labels) if l==1]; idd=[i for i,l in enumerate(labels) if l==0]
    out=[]
    for _ in range(B):
        sa=[random.choice(ia) for _ in ia]; sd=[random.choice(idd) for _ in idd]
        idx=sa+sd
        out.append(auc([scores[i] for i in idx],[labels[i] for i in idx]))
    out=sorted(out); return out[int(.025*B)],out[int(.975*B)]

def evaluate(name, active_pred):
    # build subset: chosen actives + all decoys
    keep=[i for i,x in enumerate(rows) if (x[1]==0) or active_pred(x)]
    labs=[rows[i][1] for i in keep]
    print(f"\n=== {name}  (actives={sum(labs)}, decoys={len(labs)-sum(labs)}) ===")
    for tag,vec in [("docking-alone",dock),("reactivity-alone",react),("consensus (z_dock+z_react)",consensus)]:
        s=[vec[i] for i in keep]
        a=auc(s,labs); lo,hi=boot_auc(s,labs)
        e1=ef(s,labs,0.01)
        print(f"  {tag:28} ROC-AUC={a:.3f} [{lo:.3f}, {hi:.3f}]   EF1%={e1:.2f}")

evaluate("ALL actives vs decoys", lambda x:x[1]==1)
evaluate("COVALENT actives vs decoys", lambda x:x[1]==1 and x[3]=='covalent')
evaluate("NONCOVALENT actives vs decoys", lambda x:x[1]==1 and x[3]!='covalent')
