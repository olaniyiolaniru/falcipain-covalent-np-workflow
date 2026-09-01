# matplotlib-only: chart figures + Fig4 composite (PIL). No rdkit in this process.
import os, csv
import numpy as np
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

OUT="figures"
import os as _os; _os.makedirs("figures", exist_ok=True)
BLUE,ORANGE,GREEN,VERM = "#0072B2","#E69F00","#009E73","#D55E00"
INK,MUTED = "#1a1a1a","#555555"
mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","DejaVu Sans"],
 "font.size":9,"axes.edgecolor":MUTED,"axes.linewidth":0.8,"text.color":INK,"axes.labelcolor":INK,
 "xtick.color":INK,"ytick.color":INK,"axes.spines.top":False,"axes.spines.right":False,
 "figure.dpi":300,"savefig.dpi":300})

# ---- Figure 1: funnel ----
fig,ax=plt.subplots(figsize=(4.6,3.0))
stages=["SANCDB (raw)","Valid (parsed)","Warhead-bearing (unique)","\u03b2-accessible (working set)"]
vals=[1012,995,267,185]; y=np.arange(len(vals))[::-1]; mx=max(vals)
for yi,v,s in zip(y,vals,stages):
    w=v/mx; ax.barh(yi,w,height=0.62,color=BLUE,left=(1-w)/2,zorder=3)
    ax.text(0.5,yi,f"{v:,}",ha="center",va="center",color="white",fontweight="bold",fontsize=10,zorder=4)
    ax.text(-0.03,yi,s,ha="right",va="center",fontsize=8.2)
ax.set_xlim(-0.5,1.02); ax.set_ylim(-0.6,len(vals)-0.4); ax.axis("off")
ax.set_title("Covalent-warhead filtering cascade",fontweight="bold",loc="left",x=-0.5,pad=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT,"Figure1_funnel.png"),bbox_inches="tight"); plt.close(fig)

# ---- Figure 3: DFT tuning (small multiples) ----
rows=list(csv.DictReader(open("data/leadopt_warhead_tuning.csv")))
sig=np.array([float(r["hammett_sigma_p"]) for r in rows]); od=np.argsort(sig); sig=sig[od]
lab=[rows[i]["para_substituent"].replace("(parent)","").strip() for i in od]
om=[float(rows[i]["omega_eV"]) for i in od]
fk=[float(rows[i]["fukui_plus_betaC"]) for i in od]
dE=[float(rows[i]["dE_CH3S_kcal"]) if rows[i]["dE_CH3S_kcal"] not in("","-") else None for i in od]
pi=[i for i,l in enumerate(lab) if l=="4-OH"][0]; ci=[i for i,l in enumerate(lab) if l=="4-Cl"][0]
fig,ax=plt.subplots(1,3,figsize=(7.2,2.7))
for a,(ttl,yl,yv,col) in zip(ax,[("Global electrophilicity","$\\omega$ (eV)",om,BLUE),
        ("Local Fukui at $\\beta$-carbon","$f^{+}$",fk,GREEN),
        ("Thiolate reaction energy","$\\Delta E$ (kcal mol$^{-1}$)",dE,VERM)]):
    xs=[s for s,v in zip(sig,yv) if v is not None]; ys=[v for v in yv if v is not None]
    a.plot(xs,ys,"-",color=col,lw=1.6,zorder=2); a.scatter(xs,ys,s=32,color=col,zorder=3,edgecolor="white",lw=0.8)
    for idx,mk,sz in [(pi,"D",44),(ci,"*",95)]:
        if yv[idx] is not None: a.scatter([sig[idx]],[yv[idx]],s=sz,marker=mk,facecolor="none",edgecolor="black",lw=1.3,zorder=4)
    a.set_title(ttl,fontsize=9,fontweight="bold"); a.set_ylabel(yl); a.set_xlabel("Hammett $\\sigma_{p}$"); a.tick_params(length=3)
ax[1].annotate("4-NO$_2$ collapses",(sig[-1],fk[-1]),xytext=(-3,10),textcoords="offset points",
               fontsize=7.5,color=VERM,ha="right",fontweight="bold")
from matplotlib.lines import Line2D
leg=[Line2D([0],[0],marker="D",color="none",markerfacecolor="none",markeredgecolor="black",markersize=7,label="parent (4-OH)"),
     Line2D([0],[0],marker="*",color="none",markerfacecolor="none",markeredgecolor="black",markersize=11,label="recommended (4-Cl)")]
fig.legend(handles=leg,ncol=2,loc="upper center",bbox_to_anchor=(0.5,1.02),frameon=False,fontsize=8.5,handletextpad=0.3,columnspacing=1.5)
fig.suptitle("DFT-guided warhead electronic tuning",fontsize=10,fontweight="bold",y=1.11)
fig.tight_layout(); fig.savefig(os.path.join(OUT,"Figure3_DFT_tuning.png"),bbox_inches="tight"); plt.close(fig)

# ---- Figure 4: property-change panel, then composite with structures ----
fig,ax=plt.subplots(1,4,figsize=(7.2,1.9))
metrics=[("FP-2 XP\n(kcal mol$^{-1}$)",-8.24,-4.35),("Caco-2\n(nm s$^{-1}$)",2,80),
         ("Oral abs. (%)",0,83),("RuleOf5\nviolations",3,0)]
for a,(name,pv,ov) in zip(ax,metrics):
    a.bar([0,1],[abs(pv),abs(ov)],color=[BLUE,ORANGE],width=0.62,zorder=3)
    for x,v in zip([0,1],[pv,ov]): a.text(x,abs(v),f"{v:g}",ha="center",va="bottom",fontsize=8,fontweight="bold")
    a.set_xticks([0,1]); a.set_xticklabels(["parent","OPT1"],fontsize=7.5)
    a.set_title(name,fontsize=8,fontweight="bold"); a.set_ylim(0,abs(max([pv,ov],key=abs))*1.3+0.5)
    a.set_yticks([]); a.tick_params(length=0); a.spines["left"].set_visible(False)
fig.tight_layout(); fig.savefig(os.path.join(OUT,"tmp_fig4bars.png"),bbox_inches="tight",dpi=300); plt.close(fig)

# composite: [parent | opt1] over [bars]
p=Image.open(os.path.join(OUT,"struct_parent.png")); o=Image.open(os.path.join(OUT,"struct_opt1.png"))
b=Image.open(os.path.join(OUT,"tmp_fig4bars.png"))
rowW=p.width+o.width+30; scale=rowW/b.width; b=b.resize((rowW,int(b.height*scale)))
canvas=Image.new("RGB",(rowW,p.height+b.height+20),"white")
canvas.paste(p,(0,0)); canvas.paste(o,(p.width+30,0)); canvas.paste(b,(0,p.height+20))
canvas.save(os.path.join(OUT,"Figure4_lead_optimization.png"),dpi=(300,300))
os.remove(os.path.join(OUT,"tmp_fig4bars.png"))
for f in ["struct_parent.png","struct_opt1.png"]:
    fp=os.path.join(OUT,f);
    if os.path.exists(fp): os.remove(fp)
print("charts + Fig4 composite done")
