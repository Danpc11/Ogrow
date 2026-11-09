"""Figure 4 — a self-averaging criterion separates specifiable from irreducible."""
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle

from organgrow import figstyle as F
from organgrow.config import ROOT

F.apply()
RES = ROOT / "results"
OUTDIR = ROOT / "figures" / "out"
OUTDIR.mkdir(parents=True, exist_ok=True)
K = F.C
MM = F.MM



A0 = json.load(open(RES / 'scaling_nodisorder.json'))
A1 = json.load(open(RES / 'scaling_disorder.json'))
N0 = np.array([r['N'] for r in A0])

OBS = [("cv_N", "Unit number, $N$", K["pancreas"]),
       ("cv_mean_unit", "Mean unit mass", K["grey"]),
       ("cv_mass", "Total mass", K["liver"])]
NSA = [("vessel_identity", "Vessel identity", K["sim"]),
       ("density_discord", "Local density pattern", "#9b59b6")]

fig = plt.figure(figsize=(F.WIDTH_FULL, 104 * F.MM))
gs = fig.add_gridspec(2, 3, height_ratios=[0.92, 1.0], hspace=0.50, wspace=0.50,
                      left=0.065, right=0.985, top=0.945, bottom=0.095)




def slope(x, y):
    g = np.polyfit(np.log(x), np.log(y), 1)
    return g[0], np.exp(g[1])


# ------------------------------------------------------------ a  schematic
axa = fig.add_subplot(gs[0, 0])
axa.set_xlim(0, 100); axa.set_ylim(0, 100); axa.axis("off")
rng = np.random.default_rng(3)
for i in range(3):
    x0 = 2 + i * 10
    axa.add_patch(FancyBboxPatch((x0, 66), 8.6, 26, boxstyle="round,pad=0.3,rounding_size=1",
                                 fc="white", ec=K["grey"], lw=0.6))
    for _ in range(16):
        axa.add_patch(Circle((x0 + rng.uniform(1.3, 7.3), 67.5 + rng.uniform(1.3, 23)),
                             0.8, fc=K["grey"], ec="none", alpha=0.75))
axa.text(17, 96, "replicates", ha="center", fontsize=6.2, color=K["grey"])
axa.text(17, 60, "each a realisation\nof quenched disorder",
         ha="center", va="top", fontsize=6.0)
axa.annotate("", xy=(44, 84), xytext=(32, 82),
             arrowprops=dict(arrowstyle="-|>", lw=0.7, color=K["good"]))
axa.annotate("", xy=(44, 36), xytext=(30, 52),
             arrowprops=dict(arrowstyle="-|>", lw=0.7, color=K["sim"]))
axa.add_patch(FancyBboxPatch((45, 72), 53, 24, boxstyle="round,pad=0.4,rounding_size=1.4",
                             fc=K["box"], ec=K["good"], lw=0.9))
axa.text(71.5, 90, "Self-averaging", ha="center", fontsize=6.3,
         fontweight="bold", color=K["good"])
axa.text(71.5, 79, "sums over the volume\n$\\mathrm{CV}\\propto N^{-a}$, $a>0$",
         ha="center", va="center", fontsize=6.0)
axa.add_patch(FancyBboxPatch((45, 22), 53, 24, boxstyle="round,pad=0.4,rounding_size=1.4",
                             fc=K["box"], ec=K["sim"], lw=0.9))
axa.text(71.5, 40, "Not self-averaging", ha="center", fontsize=6.3,
         fontweight="bold", color=K["sim"])
axa.text(71.5, 29, "properties of the realisation\n$a\\approx0$",
         ha="center", va="center", fontsize=6.0)

# ------------------------------------------------- b  CV vs N, no param floor
axb = fig.add_subplot(gs[0, 1])
for key, name, col in OBS:
    y = np.array([r[key] for r in A0])
    axb.loglog(N0, y, "o", color=col, ms=3.2, mew=0)
    a, c = slope(N0, y)
    xx = np.array([N0[0] * 0.85, N0[-1] * 1.18])
    axb.loglog(xx, c * xx ** a, "-", color=col, lw=0.9,
               label="%s  $a$=%.2f" % (name, -a))
for key, name, col in NSA:
    y = np.array([r[key] for r in A0])
    axb.loglog(N0, y, "s", color=col, ms=3.0, mew=0, alpha=0.85)
    a, _ = slope(N0, y)
    axb.loglog(N0, y, "--", color=col, lw=0.8, label="%s  $a$=%.2f" % (name, -a))
xr = np.array([70, 900])
axb.loglog(xr, 0.33 * xr ** -0.5, ":", color=K["ink"], lw=0.8)
axb.text(700, 0.010, "$N^{-1/2}$", fontsize=6.2, color=K["ink"])
axb.set_xlabel("Units, $N$"); axb.set_ylabel("CV across replicates")
axb.legend(loc="upper left", fontsize=F.SZ["small"] - 0.5, ncol=1,
           labelspacing=0.18, handlelength=1.0, borderpad=0.2)
axb.set_ylim(5e-4, 300)

# ---------------------------------------------------- c  exponents with CI
axc = fig.add_subplot(gs[0, 2])
items = OBS + NSA
names, exps, errs, cols = [], [], [], []
for key, name, col in items:
    y = np.array([r[key] for r in A0])
    lx, ly = np.log(N0), np.log(y)
    g, cov = np.polyfit(lx, ly, 1, cov=True)
    names.append(name); exps.append(-g[0])
    errs.append(np.sqrt(cov[0, 0])); cols.append(col)
yy = np.arange(len(names))[::-1]
axc.errorbar(exps, yy, xerr=errs, fmt="o", ms=3.4, mew=0, lw=0,
             elinewidth=0.8, ecolor=K["grey"])
for x, y, c in zip(exps, yy, cols):
    axc.scatter([x], [y], s=16, c=c, lw=0, zorder=3)
axc.axvline(0.5, color=K["ink"], ls=":", lw=0.8)
axc.axvline(0.0, color=K["grey"], lw=0.6)
axc.text(0.54, len(names) - 0.35, "$a=1/2$", fontsize=6.0, color=K["ink"])
axc.set_yticks(yy)
axc.set_yticklabels([n.replace(" ", "\n", 1) if len(n) > 11 else n
                     for n in names], fontsize=F.SZ["small"], linespacing=0.95)
axc.set_xlabel("Fitted exponent, $a$ in $\\mathrm{CV}\\propto N^{-a}$")
axc.set_xlim(-0.45, 1.75)
axc.spines["left"].set_visible(False); axc.tick_params(axis="y", length=0)

# --------------------------------------------- d  the parameter floor masks it
axd = fig.add_subplot(gs[1, 0])
for key, name, col in OBS:
    y0 = np.array([r[key] for r in A0])
    y1 = np.array([r[key] for r in A1])
    axd.loglog(N0, y0, "o-", color=col, ms=3.0, mew=0, lw=0.9)
    axd.loglog(N0, y1, "s--", color=col, ms=3.0, mew=0, lw=0.9, alpha=0.75)
axd.axhline(0.04, color=K["ink"], lw=0.7, ls=(0, (3, 2)))
axd.text(0.97, 0.95, "$\\sigma_{param}=0.04$", transform=axd.transAxes,
         fontsize=F.SZ["small"], color=K["ink"], ha="right", va="top")
axd.plot([], [], "o-", color=K["grey"], ms=3.0, mew=0, lw=0.9, label="no parameter disorder")
axd.plot([], [], "s--", color=K["grey"], ms=3.0, mew=0, lw=0.9, label="with parameter disorder")
axd.set_xlabel("Units, $N$"); axd.set_ylabel("CV across replicates")
axd.legend(frameon=False, handlelength=1.5, borderpad=0.1, labelspacing=0.22,
           handletextpad=0.45, loc="lower left", fontsize=5.8)

# ------------------------------------- e  correlation between the invariants
axe = fig.add_subplot(gs[1, 1])
REF = np.load(RES/RES / 'ref_cohort.npy')
INV = ["Unit no.", "Mean mass", "Mass spread", "Mass/vol",
       "Spacing", "Spacing sd", "Supply dist.", "Segment len."]
C = np.corrcoef(REF.T)
im = axe.imshow(C, cmap="RdBu_r", vmin=-1, vmax=1)
axe.set_xticks(np.arange(len(INV))); axe.set_yticks(np.arange(len(INV)))
axe.set_xticklabels(INV, rotation=55, ha="right", fontsize=5.2)
axe.set_yticklabels(INV, fontsize=5.2)
r = np.abs(C[np.triu_indices(len(INV), 1)]).mean()
deff = len(INV) / (1 + (len(INV) - 1) * r)
axe.set_title("mean $|r|$ = %.2f, effective axes = %.1f" % (r, deff),
              fontsize=5.9, pad=3)
cb = fig.colorbar(im, ax=axe, fraction=0.040, pad=0.04)
cb.set_label("Pearson $r$", fontsize=5.6)
cb.ax.tick_params(labelsize=5.3); cb.outline.set_linewidth(0.5)

# ---------------------------------------------- f  what can be specified
axf = fig.add_subplot(gs[1, 2])
axf.set_xlim(0, 100); axf.set_ylim(0, 100); axf.axis("off")
spec = ["Unit number", "Total mass", "Mean unit mass", "Exchange surface"]
nspec = ["Vessel identity", "Local density pattern"]
axf.text(2, 99, "Specifiable", fontsize=6.3, fontweight="bold", color=K["good"], va="top")
y = 92
for name in spec:
    axf.add_patch(FancyBboxPatch((2, y - 8), 96, 8.6,
                                 boxstyle="round,pad=0.2,rounding_size=1",
                                 fc=K["box"], ec=K["good"], lw=0.7))
    axf.text(6, y - 3.8, name, fontsize=6.1, va="center")
    axf.text(94, y - 3.8, "set a target", fontsize=5.8, va="center",
             ha="right", color=K["good"])
    y -= 10.2
axf.text(2, y - 1, "Not specifiable", fontsize=6.3, fontweight="bold",
         color=K["sim"], va="top")
y -= 10
for name in nspec:
    axf.add_patch(FancyBboxPatch((2, y - 8), 96, 8.6,
                                 boxstyle="round,pad=0.2,rounding_size=1",
                                 fc=K["box"], ec=K["sim"], lw=0.7))
    axf.text(6, y - 3.8, name, fontsize=6.1, va="center")
    axf.text(94, y - 3.8, "report only", fontsize=5.8, va="center",
             ha="right", color=K["sim"])
    y -= 10.2

F.panel_grid(fig, [[axa, axb, axc], [axd, axe, axf]],
             ["a", "b", "c", "d", "e", "f"])
fig.savefig(str(OUTDIR / "Fig4.pdf"), dpi=600, bbox_inches="tight")
fig.savefig(str(OUTDIR / "Fig4.png"), dpi=400, bbox_inches="tight")
print("ok")
