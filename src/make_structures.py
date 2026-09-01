# RDKit-only: molecular structure figures (no matplotlib in this process)
import os, csv, re
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image

OUT = "figures"
import os as _os; _os.makedirs("figures", exist_ok=True)
os.makedirs(OUT, exist_ok=True)
SM={}
with open("data/sancdb_hits_accessible.csv") as fh:
    for r in csv.DictReader(fh):
        m=re.search(r'(SANC\d+)',r["id"])
        if m: SM[m.group(1)]=r["SMILES"]
SM["OPT1"]="O=c1c(OC(=O)/C=C/c2ccc(Cl)cc2)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12"

def draw(smi,legend,size=(380,310)):
    m=Chem.MolFromSmiles(smi)
    patt=Chem.MolFromSmarts("[C;!$(C=O)]=[C;!$(C=O)][CX3]=[OX1]")
    hit=[a for mt in m.GetSubstructMatches(patt) for a in mt]
    hb=[b.GetIdx() for b in m.GetBonds() if b.GetBeginAtomIdx() in hit and b.GetEndAtomIdx() in hit]
    d=rdMolDraw2D.MolDraw2DCairo(*size)
    o=d.drawOptions(); o.legendFontSize=20; o.bondLineWidth=2
    col={a:(0.90,0.62,0.0) for a in hit}
    rdMolDraw2D.PrepareAndDrawMolecule(d,m,legend=legend,highlightAtoms=hit,highlightBonds=hb,
        highlightAtomColors=col,highlightBondColors={b:(0.90,0.62,0.0) for b in hb})
    d.FinishDrawing()
    fn=os.path.join(OUT,"tmp_%s.png"%re.sub(r'\W','',legend)[:14]); open(fn,"wb").write(d.GetDrawingText())
    return fn

# Figure 2: grid of top leads
leads=[("SANC00867","top dual lead"),("SANC00964","dual, cinnamate"),
       ("SANC00370","strong dual"),("SANC01067","strong recognition"),
       ("SANC00711","recognition binder"),("SANC01026","tier-2 quinone")]
ims=[Image.open(draw(SM[s],f"{s}  ({t})")) for s,t in leads]
w,h=ims[0].size; cols,rows,pad=3,2,12
canvas=Image.new("RGB",(cols*w+(cols+1)*pad,rows*h+(rows+1)*pad),"white")
for i,im in enumerate(ims):
    r,c=divmod(i,cols); canvas.paste(im,(pad+c*(w+pad),pad+r*(h+pad)))
canvas.save(os.path.join(OUT,"Figure2_lead_structures.png"),dpi=(300,300))

# structures for Figure 4 (composited later with the chart)
Image.open(draw(SM["SANC00867"],"SANC00867 (parent)",size=(430,330))).save(os.path.join(OUT,"struct_parent.png"))
Image.open(draw(SM["OPT1"],"OPT1 (optimized)",size=(430,330))).save(os.path.join(OUT,"struct_opt1.png"))
for f in os.listdir(OUT):
    if f.startswith("tmp_"): os.remove(os.path.join(OUT,f))
print("structures done")
