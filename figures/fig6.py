"""Figure 6 — size belongs to the individual, shape to the surroundings."""
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch

from organgrow import figstyle as F
from organgrow.config import ROOT

F.apply()
RES = ROOT / "results"
OUTDIR = ROOT / "figures" / "out"
OUTDIR.mkdir(parents=True, exist_ok=True)
K = F.C
MM = F.MM



P = json.load(open(RES / "partition.json"))
D = json.load(open(RES / "trusted_split.json"))
PRETTY = {"volume_cm3": "Volume", "surface_cm2": "Surface",
          "sphericity": "Sphericity", "elong": "Elongation",
          "flat": "Flatness", "surf_vol": "Surface / volume"}
keys = list(P["per_invariant"])

fig = plt.figure(figsize=(F.WIDTH_FULL, 100 * F.MM))
gs = fig.add_gridspec(2, 3, height_ratios=[0.92, 1.0], hspace=0.44, wspace=0.52,
                      left=0.075, right=0.975, top=0.945, bottom=0.105)




# ------------------------------------------------------------- a  the design
axa = fig.add_subplot(gs[0, 0])
axa.set_xlim(0, 100); axa.set_ylim(0, 108); axa.axis("off")
axa.text(50, 104, "96 kidneys, 48 donors", ha="center", fontsize=6.4,
         fontweight="bold", va="top")
for i, x in enumerate((16, 58)):
    axa.add_patch(FancyBboxPatch((x, 56), 26, 24,
                                 boxstyle="round,pad=0.5,rounding_size=2",
                                 fc=K["box"], ec=K["grey"], lw=0.7))
    axa.text(x + 13, 76, "donor %d" % (i + 1), ha="center", fontsize=6.0,
             color=K["grey"])
    for j, s in enumerate("LR"):
        axa.add_patch(Circle((x + 8 + j * 11, 65), 4.4, fc=K["kidney"], ec="none",
                             alpha=0.50 + 0.30 * j))
        axa.text(x + 8 + j * 11, 65, s, ha="center", va="center", fontsize=5.8,
                 color="white", fontweight="bold")
axa.add_patch(FancyArrowPatch((24, 57), (35, 57), arrowstyle="<|-|>",
                              mutation_scale=5, lw=0.8, color=K["liver"]))
axa.text(29, 51, "within", ha="center", fontsize=6.0, color=K["liver"])
axa.add_patch(FancyArrowPatch((29, 82), (71, 82),
                              connectionstyle="arc3,rad=-0.28",
                              arrowstyle="<|-|>", mutation_scale=5, lw=0.8,
                              color=K["kidney"]))
axa.text(50, 88, "between", ha="center", fontsize=6.0, color=K["kidney"])
axa.text(50, 38, "what the individual fixes appears\nbetween donors; what the "
                 "organ's\nsurroundings fix appears within one",
         ha="center", va="center", fontsize=6.1)
axa.text(50, 16, "ratio %.2f    $P$ = %.4f" % (P["ratio"], P["p_permutation"]),
         ha="center", fontsize=6.4, color=K["ink"], fontweight="bold")
axa.text(50, 7, "permutation null 1.00 $\\pm$ %.2f" % P["null_sd"],
         ha="center", fontsize=5.9, color=K["grey"])

# --------------------------------------------- b  where each invariant lives
axb = fig.add_subplot(gs[0, 1:])
frac = np.array([P["per_invariant"][k]["between_frac"] for k in keys]) * 100
order = np.argsort(-frac)
for i, o in enumerate(order):
    axb.barh(i, frac[o], height=0.62, color=K["kidney"], lw=0)
    axb.barh(i, 100 - frac[o], left=frac[o], height=0.62, color=K["liver"],
             alpha=0.78, lw=0)
    axb.text(frac[o] / 2, i, "%.0f%%" % frac[o], va="center", ha="center",
             fontsize=5.9, color="white", fontweight="bold")
axb.set_yticks(range(len(keys)))
axb.set_yticklabels([PRETTY[keys[o]] for o in order], fontsize=6.0)
axb.invert_yaxis()
axb.axvline(50, color=K["ink"], lw=0.7, ls=(0, (2, 2)))
axb.set_xlabel("Share of variance (%)"); axb.set_xlim(0, 100)
axb.spines["left"].set_visible(False); axb.tick_params(axis="y", length=0)
axb.text(0.16, 1.05, "between donors", transform=axb.transAxes, fontsize=6.2,
         color=K["kidney"], ha="center", fontweight="bold")
axb.text(0.84, 1.05, "within donor", transform=axb.transAxes, fontsize=6.2,
         color=K["liver"], ha="center", fontweight="bold")

# ------------------------------------- c  the within-donor part is systematic
axc = fig.add_subplot(gs[1, 0])
rel = np.array([P["per_invariant"][k]["lr_rel_diff"] for k in keys])
nL = np.array([P["per_invariant"][k]["n_left_greater"] for k in keys])
npair = P["per_invariant"][keys[0]]["n_pairs"]
o2 = np.argsort(rel)
for i, o in enumerate(o2):
    col = K["liver"] if rel[o] > 0 else K["kidney"]
    axc.barh(i, rel[o], height=0.62, color=col, lw=0)
    p = P["per_invariant"][keys[o]]["p_wilcoxon"]
    ptxt = "$P$ < 0.0001" if p < 1e-4 else "$P$ = %.4f" % p
    sig = p < P["alpha_bonferroni"]
    axc.text(rel[o] + 0.6, i, "%d/%d, %s%s" % (nL[o], npair, ptxt,
             "  *" if sig else ""), va="center", ha="left",
             fontsize=F.SZ["small"] - 0.5,
             color=K["ink"] if sig else K["grey"])
axc.set_yticks(range(len(keys)))
axc.set_yticklabels([PRETTY[keys[o]] for o in o2], fontsize=6.0)
axc.axvline(0, color=K["ink"], lw=0.7)
axc.set_xlabel("Left $-$ right (% of mean)"); axc.set_xlim(-30, 24)
axc.spines["left"].set_visible(False); axc.tick_params(axis="y", length=0)
axc.text(0.5, -0.30, "paired Wilcoxon; * survives Bonferroni "
         "($\\alpha$ = %.4f); 3 of 6" % P["alpha_bonferroni"],
         transform=axc.transAxes, ha="center", fontsize=F.SZ["small"],
         color=K["grey"], style="italic")

# ------------------------------------------------ d  the paired measurements
axd = fig.add_subplot(gs[1, 1])
pid = {}
for k, v in D.items():
    pid.setdefault(v["pid"], {})[v["side"]] = v
pr = [(d["L"], d["R"]) for d in pid.values() if "L" in d and "R" in d]
L = np.array([a["flat"] for a, _ in pr])
R = np.array([b["flat"] for _, b in pr])
axd.scatter(R, L, s=11, c=K["liver"], lw=0, alpha=0.65)
lim = [min(L.min(), R.min()) * 0.92, max(L.max(), R.max()) * 1.06]
axd.plot(lim, lim, color=K["ink"], lw=0.8, ls=(0, (3, 2)))
sl, ic = np.polyfit(R, L, 1)
xx = np.linspace(lim[0], lim[1], 50)
axd.plot(xx, sl * xx + ic, color=K["liver"], lw=1.1)
r = np.corrcoef(R, L)[0, 1]
axd.text(0.04, 0.88, "slope %.2f, $r$ = %.2f" % (sl, r),
         transform=axd.transAxes, fontsize=F.SZ["small"], va="top",
         color=K["liver"])
axd.set_xlim(lim); axd.set_ylim(lim)
axd.set_xlabel("Right kidney, flatness")
axd.set_ylabel("Left kidney, flatness")
axd.text(0.04, 0.96, "%d of %d above the line" % (int((L > R).sum()), len(L)),
         transform=axd.transAxes, fontsize=6.0, va="top", color=K["liver"])
axd.set_aspect("equal", adjustable="box")

# ------------------------------------------------------ e  what this implies
axe = fig.add_subplot(gs[1, 2])
axe.set_xlim(0, 100); axe.set_ylim(0, 100); axe.axis("off")
axe.text(2, 99, "Three sources", fontsize=6.4, fontweight="bold", va="top")
rows = [("Systemic", "size", "specify it", K["kidney"]),
        ("Positional", "mean shape", "set by where it sits", K["liver"]),
        ("Contingent", "the remainder", "report, do not target", K["sim"])]
y = 84
for name, what, act, col in rows:
    axe.add_patch(FancyBboxPatch((2, y - 24), 96, 24,
                                 boxstyle="round,pad=0.3,rounding_size=1.2",
                                 fc=K["box"], ec=col, lw=0.8))
    axe.text(7, y - 6.5, name, fontsize=F.SZ["annot"], va="center", color=col,
             fontweight="bold")
    axe.text(7, y - 14, "governs " + what, fontsize=F.SZ["small"],
             va="center")
    axe.text(7, y - 21, act, fontsize=F.SZ["small"], va="center",
             color=col, style="italic")
    y -= 28
axe.text(50, -6, "a construct can reach its mass target and\n"
                 "still not match native proportions",
         ha="center", fontsize=5.9, color=K["grey"], style="italic",
         clip_on=False)

F.panel_grid(fig, [[axa, axb], [axc, axd, axe]],
             ["a", "b", "c", "d", "e"])
fig.savefig(str(OUTDIR / "Fig6.pdf"), dpi=600, bbox_inches="tight")
fig.savefig(str(OUTDIR / "Fig6.png"), dpi=400, bbox_inches="tight")
print("wrote figures/out/Fig6.pdf")
