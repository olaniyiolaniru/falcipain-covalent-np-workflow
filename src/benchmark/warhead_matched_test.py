#!/usr/bin/env python3
"""Stringent warhead-matched benchmark: covalent actives vs warhead-BEARING decoys (both carry a
Michael warhead). Tests whether enrichment survives once the trivial warhead-presence signal is
removed. Uses existing docking scores only."""
import csv, math, random, statistics as st, io, sys, bisect
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
RDLogger.DisableLog('rdApp.*')
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding='utf-8',errors='replace')
BM=r"data/validation"; SCORES=r"data/validation/scores_FP2_full.csv"
WARHEADS=[("[CX3]=[CX3][CX3]=[OX1]",1),("[CH2]=[CX3][CX3]=[OX1]",1),("[CX3]=[CX3][CX2]#[NX1]",1),
 ("[CX3]=[CX3][SX4](=[OX1])(=[OX1])",2),("O=C1C=CC(=O)N1",2),("O=C1C=CC(=O)C=C1",2),("O=C1C=CC=CC1=O",2)]
Q=[(Chem.MolFromSmarts(s),t) for s,t in WARHEADS]
_lfc=rdMolStandardize.LargestFragmentChooser(); _norm=rdMolStandardize.Normalizer()
def tier(sm):
    if not sm: return 0
    m=Chem.MolFromSmiles(sm)
    if m is None: return 0
    try: m=_norm.normalize(_lfc.choose(m))
    except: pass
    t=0
    for q,tr in Q:
        if q is not None and m.HasSubstructMatch(q): t=max(t,tr)
    return t
smi={}
for r in csv.DictReader(open(BM+"/actives_FP2_classified.csv")): smi[r["molecule_chembl_id"]]=r["canonical_smiles"]
for r in csv.DictReader(open(BM+"/decoys_FP2.csv")): smi[r["decoy_chembl_id"]]=r["smiles"]
rows=[]
for r in csv.DictReader(open(SCORES)):
    rows.append([r["id"],int(r["label"]),float(r["score"]),r["class"],tier(smi.get(r["id"],""))])

cov_act=[x for x in rows if x[1]==1 and x[3]=="covalent" and x[4]>=1]
wh_dec =[x for x in rows if x[1]==0 and x[4]>=1]
print(f"covalent actives (warhead-bearing): {len(cov_act)}")
print(f"warhead-bearing decoys: {len(wh_dec)}")
def tdist(s):
    from collections import Counter; return dict(Counter(x[4] for x in s))
print("active tiers:",tdist(cov_act)," decoy tiers:",tdist(wh_dec))

# tier-stratified matching: cap decoys so tier ratio matches actives (stringent)
random.seed(7)
from collections import defaultdict
dec_by=defaultdict(list)
for x in wh_dec: dec_by[x[4]].append(x)
# target: 8 decoys per active within same tier (as available)
matched=[]
act_by=defaultdict(list)
for x in cov_act: act_by[x[4]].append(x)
for tr,acts in act_by.items():
    pool=dec_by.get(tr,[])[:]; random.shuffle(pool)
    need=min(len(pool), 8*len(acts))
    matched+=pool[:need]
print(f"\ntier-matched decoys selected: {len(matched)} (ratio {len(matched)/len(cov_act):.1f}:1)")

def build(active_list, decoy_list):
    s=[(-x[2]) for x in active_list+decoy_list]      # docking (higher=better)
    r=[float(x[4]) for x in active_list+decoy_list]   # reactivity tier
    lab=[1]*len(active_list)+[0]*len(decoy_list)
    def z(v):
        m=st.mean(v); sd=st.pstdev(v) or 1
        return [(a-m)/sd for a in v]
    cons=[a+b for a,b in zip(z(s),z(r))]
    return s,r,cons,lab
def auc(sc,lab):
    pos=[a for a,l in zip(sc,lab) if l==1]; neg=sorted(a for a,l in zip(sc,lab) if l==0)
    if not pos or not neg: return float('nan')
    U=0.0
    for p in pos:
        lo=bisect.bisect_left(neg,p); hi=bisect.bisect_right(neg,p); U+=lo+0.5*(hi-lo)
    return U/(len(pos)*len(neg))
def ef(sc,lab,frac):
    o=sorted(range(len(sc)),key=lambda i:-sc[i]); n=max(1,int(round(frac*len(o)))); P=sum(lab)
    return (sum(lab[i] for i in o[:n])/n)/(P/len(lab)) if P else float('nan')
def boot(sc,lab,B=1000,seed=3):
    random.seed(seed); ia=[i for i,l in enumerate(lab) if l==1]; idd=[i for i,l in enumerate(lab) if l==0]
    o=[]
    for _ in range(B):
        idx=[random.choice(ia) for _ in ia]+[random.choice(idd) for _ in idd]
        o.append(auc([sc[i] for i in idx],[lab[i] for i in idx]))
    o.sort(); return o[int(.025*B)],o[int(.975*B)]

for name,dec in [("ALL warhead-bearing decoys",wh_dec),("tier-matched decoys (8:1)",matched)]:
    s,r,cons,lab=build(cov_act,dec)
    print(f"\n=== Covalent actives vs {name} (actives={sum(lab)}, decoys={len(lab)-sum(lab)}) ===")
    for tag,vec in [("docking-alone",s),("reactivity-alone",r),("consensus",cons)]:
        a=auc(vec,lab); lo,hi=boot(vec,lab)
        print(f"  {tag:18} ROC-AUC={a:.3f} [{lo:.3f}, {hi:.3f}]   EF1%={ef(vec,lab,0.01):.2f}")
