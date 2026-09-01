import os, csv, re
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image

OUT = "figures"
import os as _os; _os.makedirs("figures", exist_ok=True)
os.makedirs(OUT, exist_ok=True)

# ---- global publication style (Okabe-Ito, CVD-safe) ----
BLUE, ORANGE, GREEN, VERM, PURPLE, SKY = "#0072B2","#E69F00","#009E73","#D55E00","#CC79A7","#56B4E9"
INK, MUTED, GRID = "#1a1a1a", "#555555", "#d9d9d9"
mpl.rcParams.update({
    "font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
    "font.size":9,"axes.titlesize":10,"axes.labelsize":9.5,"axes.edgecolor":MUTED,
    "axes.linewidth":0.8,"xtick.color":INK,"ytick.color":INK,"text.color":INK,
    "axes.labelcolor":INK,"axes.titlecolor":INK,"figure.dpi":300,"savefig.dpi":300,
    "axes.spines.top":False,"axes.spines.right":False,"xtick.direction":"out",
    "ytick.direction":"out","axes.grid":False,
})

def smiles_map():
    d={}
    with open("data/sancdb_hits_accessible.csv") as fh:
        for r in csv.DictReader(fh):
            m=re.search(r'(SANC\d+)',r["id"])
            if m: d[m.group(1)]=r["SMILES"]
    return d
SM=smiles_map()

# =========================================================
# FIGURE 1 - compound filtering funnel
# =========================================================
def fig1():
    stages=["SANCDB\n(raw)","Valid\n(parsed)","Warhead-bearing\n(unique)","\u03b2-accessible\n(working set)"]
    vals=[1012,995,267,185]
    fig,ax=plt.subplots(figsize=(4.6,3.0))
    y=np.arange(len(vals))[::-1]
    maxv=max(vals)
    for yi,v,s in zip(y,vals,stages):
        w=v/maxv
        ax.barh(yi,w,height=0.62,color=BLUE,edgecolor="white",linewidth=0,zorder=3,
                left=(1-w)/2)
        ax.text(0.5,yi,f"{v:,}",ha="center",va="center",color="white",
                fontweight="bold",fontsize=10,zorder=4)
        ax.text(-0.02,yi,s,ha="right",va="center",fontsize=8.2,color=INK)
    ax.set_xlim(-0.28,1.02); ax.set_ylim(-0.6,len(vals)-0.4)
    ax.axis("off")
    ax.set_title("Covalent-warhead filtering cascade",fontweight="bold",loc="left",x=-0.28)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT,"Figure1_funnel.png"),bbox_inches="tight")
    plt.close(fig)

# =========================================================
# FIGURE 2 - structures of top dual covalent leads (warhead highlighted)
# =========================================================
def draw_mol(smi,legend,size=(360,300)):
    m=Chem.MolFromSmiles(smi)
    patt=Chem.MolFromSmarts("[C;!$(C=O)]=[C;!$(C=O)][CX3]=[OX1]")
    hit=[a for match in m.GetSubstructMatches(patt) for a in match]
    hb=[b.GetIdx() for b in m.GetBonds() if b.GetBeginAtomIdx() in hit and b.GetEndAtomIdx() in hit]
    d=rdMolDraw2D.MolDraw2DCairo(*size)
    o=d.drawOptions(); o.legendFontSize=18; o.bondLineWidth=2
    hlcol={a:(0.90,0.62,0.0) for a in hit}  # orange warhead
    rdMolDraw2D.PrepareAndDrawMolecule(d,m,legend=legend,highlightAtoms=hit,
        highlightBonds=hb,highlightAtomColors=hlcol,
        highlightBondColors={b:(0.90,0.62,0.0) for b in hb})
    d.FinishDrawing()
    fn=os.path.join(OUT,"_m_%s.png"%re.sub(r'\W','',legend)[:12]); open(fn,"wb").write(d.GetDrawingText())
    return fn

def fig2():
    leads=[("SANC00867","top dual lead"),("SANC00964","dual, cinnamate"),
           ("SANC00370","strong dual"),("SANC01067","strong recognition"),
           ("SANC00711","recognition binder"),("SANC01026","tier-2 quinone")]
    imgs=[draw_mol(SM[s],f"{s}\n({t})") for s,t in leads]
    ims=[Image.open(f) for f in imgs]
    w,h=ims[0].size; cols=3; rows=2; pad=10
    canvas=Image.new("RGB",(cols*w+(cols+1)*pad,rows*h+(rows+1)*pad),"white")
    for i,im in enumerate(ims):
        r,c=divmod(i,cols); canvas.paste(im,(pad+c*(w+pad),pad+r*(h+pad)))
    canvas.save(os.path.join(OUT,"Figure2_lead_structures.png"),dpi=(300,300))

# =========================================================
# FIGURE 3 - DFT warhead reactivity tuning map (small multiples)
# =========================================================
def fig3():
    rows=list(csv.DictReader(open("data/leadopt_warhead_tuning.csv")))
    lab=[r["para_substituent"].replace("(parent)","") for r in rows]
    sig=[float(r["hammett_sigma_p"]) for r in rows]
    om=[float(r["omega_eV"]) for r in rows]
    fk=[float(r["fukui_plus_betaC"]) for r in rows]
    dE=[float(r["dE_CH3S_kcal"]) if r["dE_CH3S_kcal"] not in ("","-") else None for r in rows]
    order=np.argsort(sig); sig=np.array(sig)[order]
    def rr(a): return [a[i] for i in order]
    lab=rr(lab); om=rr(om); fk=rr(fk); dE=rr(dE)
    parent_i=[i for i,l in enumerate(lab) if "4-OH" in l][0]
    cl_i=[i for i,l in enumerate(lab) if l.strip()=="4-Cl"][0]

    fig,axes=plt.subplots(1,3,figsize=(7.2,2.8))
    panels=[("Global electrophilicity","\u03c9 (eV)",om,BLUE),
            ("Local Fukui at \u03b2-C","f\u207a",fk,GREEN),
            ("Thiolate reaction energy","\u0394E (kcal mol\u207b\u00b9)",dE,VERM)]
    for ax,(ttl,yl,yv,col) in zip(axes,panels):
        xs=[s for s,v in zip(sig,yv) if v is not None]
        ys=[v for v in yv if v is not None]
        ax.plot(xs,ys,"-",color=col,lw=1.6,zorder=2)
        ax.scatter(xs,ys,s=34,color=col,zorder=3,edgecolor="white",linewidth=0.8)
        # emphasize parent (4-OH) and recommended (4-Cl)
        for idx,mk,mc in [(parent_i,"D","#000000"),(cl_i,"*","#000000")]:
            if yv[idx] is not None:
                ax.scatter([sig[idx]],[yv[idx]],s=90 if mk=="*" else 46,marker=mk,
                           facecolor="none",edgecolor=mc,linewidth=1.3,zorder=4)
        ax.set_title(ttl,fontsize=9,fontweight="bold")
        ax.set_ylabel(yl); ax.set_xlabel("Hammett \u03c3$_p$")
        ax.tick_params(length=3)
    axes[0].annotate("4-NO$_2$",(sig[-1],om[-1]),textcoords="offset points",
                     xytext=(-4,-12),fontsize=7.5,color=MUTED,ha="right")
    axes[1].annotate("4-NO$_2$\ncollapses",(sig[-1],fk[-1]),textcoords="offset points",
                     xytext=(-2,8),fontsize=7.5,color=VERM,ha="right",fontweight="bold")
    fig.suptitle("DFT-guided warhead electronic tuning  (\u25c6 parent 4-OH   \u2605 recommended 4-Cl)",
                 fontsize=9.5,fontweight="bold",y=1.04,x=0.5)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT,"Figure3_DFT_tuning.png"),bbox_inches="tight")
    plt.close(fig)

# =========================================================
# FIGURE 4 - lead optimization before/after (structures + property panels)
# =========================================================
def fig4():
    parent_smi=SM["SANC00867"]
    opt_smi="O=c1c(OC(=O)/C=C/c2ccc(Cl)cc2)c(-c2ccc(O)cc2)oc2cc(O)cc(O)c12"
    p_img=draw_mol(parent_smi,"SANC00867 (parent)",size=(420,320))
    o_img=draw_mol(opt_smi,"OPT1 (optimized)",size=(420,320))
    pim,oim=Image.open(p_img),Image.open(o_img)

    fig=plt.figure(figsize=(7.2,4.6))
    gs=fig.add_gridspec(2,4,height_ratios=[1.15,1.0],hspace=0.55,wspace=0.55)
    axp=fig.add_subplot(gs[0,0:2]); axp.imshow(pim); axp.axis("off")
    axo=fig.add_subplot(gs[0,2:4]); axo.imshow(oim); axo.axis("off")
    axp.annotate("",xy=(1.06,0.5),xytext=(0.98,0.5),xycoords="axes fraction",
                 arrowprops=dict(arrowstyle="-|>",color=INK,lw=1.6))

    metrics=[("FP-2 XP\n(kcal mol\u207b\u00b9)",-8.24,-4.35,"lower=better",False),
             ("Caco-2\n(nm s\u207b\u00b9)",2,80,"higher=better",True),
             ("Oral abs.\n(%)",0,83,"higher=better",True),
             ("RuleOf5\nviolations",3,0,"lower=better",False)]
    for j,(name,pv,ov,note,up) in enumerate(metrics):
        ax=fig.add_subplot(gs[1,j])
        ax.bar([0,1],[abs(pv),abs(ov)],color=[BLUE,ORANGE],width=0.62,zorder=3,edgecolor="white")
        for x,v in zip([0,1],[pv,ov]):
            ax.text(x,abs(v),f"{v:g}",ha="center",va="bottom",fontsize=8,fontweight="bold")
        ax.set_xticks([0,1]); ax.set_xticklabels(["parent","OPT1"],fontsize=7.5)
        ax.set_title(name,fontsize=8,fontweight="bold")
        ax.set_ylim(0,abs(max([pv,ov],key=abs))*1.28+0.5)
        ax.set_yticks([]); ax.tick_params(length=0)
        for s in ["left"]: ax.spines[s].set_visible(False)
    fig.suptitle("DFT-guided lead optimization: SANC00867 \u2192 OPT1  (ADMET recovered; FP-2 binding cost)",
                 fontsize=9.5,fontweight="bold",y=0.99)
    fig.savefig(os.path.join(OUT,"Figure4_lead_optimization.png"),bbox_inches="tight")
    plt.close(fig)

fig1(); print("Figure 1 done")
fig2(); print("Figure 2 done")
fig3(); print("Figure 3 done")
fig4(); print("Figure 4 done")
# cleanup temp mol pngs
for f in os.listdir(OUT):
    if f.startswith("_m_"): os.remove(os.path.join(OUT,f))
print("all figures ->", OUT)
