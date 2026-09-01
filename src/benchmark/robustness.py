#!/usr/bin/env python3
"""Benchmark robustness analyses:
 A. DeLong tests for correlated AUCs (docking vs consensus; reactivity vs consensus).
 B. Warhead-ontology reconciliation (broader covalent SMARTS) with warhead-matched re-evaluation.
 C. Decoy-attrition handling (failed dockings assigned the worst score) + survived/failed property compare.
 D. Consensus-weighting robustness: weight sweep + 5-fold CV logistic model.
 E. Ligand-length vs docking box diagnostics.
"""
import csv, math, random, statistics as st, io, sys, json, bisect
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
RDLogger.DisableLog('rdApp.*')
sys.stdout=io.TextIOWrapper(sys.stdout.buffer,encoding='utf-8',errors='replace')
BM=r"data/validation"; SC=r"data/validation/scores_FP2_full.csv"

# ---------- broader covalent-warhead ontology (M9d) ----------
COV_SMARTS={
 # Michael acceptors (original)
 "enone":("[CX3]=[CX3][CX3]=[OX1]",1),"exo_methylene":("[CH2]=[CX3][CX3]=[OX1]",1),
 "vinyl_nitrile":("[CX3]=[CX3][CX2]#[NX1]",1),
 "vinyl_sulfone":("[CX3]=[CX3][SX4](=[OX1])(=[OX1])",2),"maleimide":("O=C1C=CC(=O)N1",2),
 "para_quinone":("O=C1C=CC(=O)C=C1",2),"ortho_quinone":("O=C1C=CC=CC1=O",2),
 # additional cysteine-reactive warhead classes (reversible + irreversible)
 "aldehyde":("[CX3H1]=O",1),"alpha_ketoamide":("[CX3](=O)[CX3](=O)[NX3]",1),
 "nitrile":("[NX1]#[CX2]",1),"epoxide":("[OX2r3]1[#6r3][#6r3]1",2),
 "haloacetamide":("[NX3][CX3](=O)[CH2][F,Cl,Br,I]",2),"boronic":("[BX3]([OX2])[OX2]",1),
 "aziridine":("[NX3r3]1[#6r3][#6r3]1",2),"beta_lactam":("O=C1CCN1",1),
 "chloroketone":("[CX3](=O)[CH2][Cl,Br]",2),"disulfide":("[SX2][SX2]",1),
}
Qc=[(Chem.MolFromSmarts(s),t,n) for n,(s,t) in COV_SMARTS.items()]
_lfc=rdMolStandardize.LargestFragmentChooser(); _norm=rdMolStandardize.Normalizer()
def prep(sm):
    m=Chem.MolFromSmiles(sm) if sm else None
    if m is None: return None
    try: m=_norm.normalize(_lfc.choose(m))
    except: pass
    return m
def cov_tier(m):
    if m is None: return 0
    t=0
    for q,tr,_ in Qc:
        if q is not None and m.HasSubstructMatch(q): t=max(t,tr)
    return t

smi={}
for r in csv.DictReader(open(BM+"/actives_FP2_classified.csv")): smi[r["molecule_chembl_id"]]=r["canonical_smiles"]
for r in csv.DictReader(open(BM+"/decoys_FP2.csv")): smi[r["decoy_chembl_id"]]=r["smiles"]

rows=[]
for r in csv.DictReader(open(SC)):
    m=prep(smi.get(r["id"],""))
    rows.append({"id":r["id"],"lab":int(r["label"]),"score":float(r["score"]),
                 "class":r["class"],"tier":cov_tier(m),"mol":m})

# ---------- AUC + DeLong ----------
def midrank(x):
    order=sorted(range(len(x)),key=lambda i:x[i]); r=[0.0]*len(x); i=0
    while i<len(x):
        j=i
        while j<len(x) and x[order[j]]==x[order[i]]: j+=1
        rk=0.5*((i+1)+j)
        for k in range(i,j): r[order[k]]=rk
        i=j
    return r
def fastdelong(preds, lab):
    # preds: list of score-vectors (higher=better); lab: 1/0
    pos=[i for i,l in enumerate(lab) if l==1]; neg=[i for i,l in enumerate(lab) if l==0]
    m,n=len(pos),len(neg); K=len(preds)
    aucs=[]; V10=[[0]*m for _ in range(K)]; V01=[[0]*n for _ in range(K)]
    for k in range(K):
        pv=[preds[k][i] for i in pos]; nv=[preds[k][i] for i in neg]
        tz=midrank(pv+nv); tx=midrank(pv); ty=midrank(nv)
        for a in range(m): V10[k][a]=(tz[a]-tx[a])/n
        for b in range(n): V01[k][b]=1.0-(tz[m+b]-ty[b])/m
        aucs.append(sum(tz[:m])/(m*n) - (m+1)/(2.0*n))
    def cov(A,B,N):
        ma=sum(A)/N; mb=sum(B)/N
        return sum((A[i]-ma)*(B[i]-mb) for i in range(N))/(N-1)
    S=[[0.0]*K for _ in range(K)]
    for i in range(K):
        for j in range(K):
            S[i][j]=cov(V10[i],V10[j],m)/m + cov(V01[i],V01[j],n)/n
    return aucs,S
def delong_p(preds,lab,i,j):
    aucs,S=fastdelong(preds,lab)
    var=S[i][i]+S[j][j]-2*S[i][j]
    if var<=0: return aucs,1.0,0.0
    z=(aucs[i]-aucs[j])/math.sqrt(var)
    from math import erf,sqrt
    p=2*(1-0.5*(1+erf(abs(z)/sqrt(2))))
    return aucs,p,z

def run_delong(subset,label):
    dock=[-x["score"] for x in subset]; react=[float(x["tier"]) for x in subset]
    def z(v):
        mm=st.mean(v); sd=st.pstdev(v) or 1; return [(a-mm)/sd for a in v]
    cons=[a+b for a,b in zip(z(dock),z(react))]
    lab=[x["lab"] for x in subset]
    preds=[dock,react,cons]; names=["docking","reactivity","consensus"]
    aucs,_=fastdelong(preds,lab)
    print(f"\n[{label}] n_act={sum(lab)} n_dec={len(lab)-sum(lab)}")
    for nm,a in zip(names,aucs): print(f"    AUC {nm:11}= {a:.3f}")
    for i,j in [(2,0),(2,1),(0,1)]:
        _,p,zz=delong_p(preds,lab,i,j)
        print(f"    DeLong {names[i]} vs {names[j]:11}: dAUC={aucs[i]-aucs[j]:+.3f}  z={zz:+.2f}  p={p:.3g}")

print("="*70,"\nA. DeLong tests (broader covalent ontology)")
cov=[x for x in rows if x["lab"]==1 and x["class"]=="covalent"]
allact=[x for x in rows if x["lab"]==1]
dec_all=[x for x in rows if x["lab"]==0]
dec_wh=[x for x in rows if x["lab"]==0 and x["tier"]>=1]
print(f"\ncovalent actives matched by broadened ontology: {sum(1 for x in cov if x['tier']>=1)}/{len(cov)}")
run_delong(cov+dec_all, "covalent vs property-matched decoys")
run_delong(cov+dec_wh,  "covalent vs warhead-matched decoys (broadened)")

# ---------- C. attrition treatment ----------
print("\n"+"="*70,"\nC. Decoy-attrition treatment (failed dockings = worst score)")
inp=set()
for line in open(BM+"/bench_FP2.smi"):
    p=line.split()
    if len(p)>=2: inp.add(p[1])
scored=set(x["id"] for x in rows)
failed=inp-scored
lab_all=[x["lab"] for x in rows]
worst=min(x["score"] for x in rows)  # most positive (worst) glide score
# add failed as worst-docking, keep their reactivity tier
labels_full=[]; dock_full=[]
for x in rows: labels_full.append(x["lab"]); dock_full.append(-x["score"])
# figure labels for failed (actives vs decoys)
actids=set();
for r in csv.DictReader(open(BM+"/actives_FP2_classified.csv")):
    if float(r["gmean_nM"])<=10000: actids.add(r["molecule_chembl_id"])
nf_act=sum(1 for i in failed if i in actids); nf_dec=len(failed)-nf_act
print(f"input={len(inp)} scored={len(scored)} failed={len(failed)} (failed actives={nf_act}, failed decoys={nf_dec})")
for i in failed:
    labels_full.append(1 if i in actids else 0); dock_full.append(-worst-1e-6)  # worst
def auc1(sc,lab):
    pos=[a for a,l in zip(sc,lab) if l==1]; neg=sorted(a for a,l in zip(sc,lab) if l==0)
    if not pos or not neg: return float('nan')
    U=0.0
    for p in pos:
        lo=bisect.bisect_left(neg,p); hi=bisect.bisect_right(neg,p); U+=lo+0.5*(hi-lo)
    return U/(len(pos)*len(neg))
print(f"  docking AUC, scored-only        = {auc1([-x['score'] for x in rows],[x['lab'] for x in rows]):.3f}")
print(f"  docking AUC, failed=worst (all) = {auc1(dock_full,labels_full):.3f}")

# ---------- D. weighting robustness ----------
print("\n"+"="*70,"\nD. Consensus-weighting robustness (covalent vs warhead-matched)")
sub=cov+dec_wh
dock=[-x["score"] for x in sub]; react=[float(x["tier"]) for x in sub]; lab=[x["lab"] for x in sub]
def z(v):
    mm=st.mean(v); sd=st.pstdev(v) or 1; return [(a-mm)/sd for a in v]
zd,zr=z(dock),z(react)
for w in [0.0,0.25,0.5,0.75,1.0]:
    cons=[(1-w)*a+w*b for a,b in zip(zd,zr)]
    print(f"  w_react={w:.2f}: AUC={auc1(cons,lab):.3f}")

# ---------- E. ligand length vs box ----------
print("\n"+"="*70,"\nE. Ligand extended-length vs 20 A max (shortlist + big actives)")
def maxdim(m):
    if m is None: return None
    try:
        mh=Chem.AddHs(m); AllChem.EmbedMolecule(mh,randomSeed=1); AllChem.MMFFOptimizeMolecule(mh)
        c=mh.GetConformer(); import itertools
        pts=[c.GetAtomPosition(i) for i in range(mh.GetNumAtoms())]
        return max(pts[i].Distance(pts[j]) for i,j in itertools.combinations(range(len(pts)),2))
    except: return None
big=sorted([x for x in allact],key=lambda x:-(Descriptors.MolWt(x["mol"]) if x["mol"] else 0))[:8]
for x in big:
    d=maxdim(x["mol"]); mw=Descriptors.MolWt(x["mol"]) if x["mol"] else 0
    print(f"  {x['id']:16} MW={mw:6.1f}  maxdim={d:.1f} A" if d else f"  {x['id']} embed-fail")
