"""Figure 5 — structural identifiability determines what must be measured."""
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from organgrow import figstyle as F
from organgrow.config import ROOT

F.apply()
RES = ROOT / "results"
OUTDIR = ROOT / "figures" / "out"
OUTDIR.mkdir(parents=True, exist_ok=True)
K = F.C
MM = F.MM



D = json.load(open(RES/RES / 'ident.json'))
PAR = D["params"]
PRETTY = {"nu": "$\\nu$", "c_set": "$c_{set}$", "lam": "$\\lambda$",
          "beta": "$\\beta$", "C0": "$C_0$"}
SETS = list(D["sets"].keys())
COLS = [K["pancreas"], K["liver"], K["kidney"]]

fig = plt.figure(figsize=(F.WIDTH_FULL, 108 * F.MM))
gs = fig.add_gridspec(2, 3, height_ratios=[0.88, 1.0], hspace=0.52, wspace=0.52,
                      left=0.070, right=0.985, top=0.945, bottom=0.095)




# ------------------------------------------------------- a  observable sets
axa = fig.add_subplot(gs[0, 0])
axa.set_xlim(0, 100); axa.set_ylim(0, 100); axa.axis("off")
axa.text(2, 99, "What is measured", fontsize=6.4, fontweight="bold", va="top")
rows = [("Final $N$ only", "1 value", K["pancreas"]),
        ("$N$ trajectory", "6 values", K["liver"]),
        ("$N$ + per-unit supply", "12 values", K["kidney"])]

y = 86
for name, note, col in rows:
    axa.add_patch(FancyBboxPatch((2, y - 11), 96, 11.5,
                                 boxstyle="round,pad=0.25,rounding_size=1",
                                 fc=K["box"], ec=col, lw=0.8))
    axa.text(6, y - 5.6, name, fontsize=6.2, va="center")
    axa.text(94, y - 5.6, note, fontsize=5.8, va="center", ha="right", color=K["grey"])
    y -= 15.5
axa.text(50, 8, "five parameters to recover:\n" + "  ".join(PRETTY[p] for p in PAR),
         ha="center", va="center", fontsize=6.1)

# ------------------------------------------- b  sensitivity map (where info is)
axb = fig.add_subplot(gs[0, 1])
SZ = np.load(RES/RES / 'sens.npz')
J = SZ["J"]
keys = [str(k) for k in SZ["keys"]]
days = [int(d) for d in SZ["days"]]
PNAMES = {"N": "$N$", "cu": "$c_u$", "nv": "vessels"}
v = np.nanmax(np.abs(J))
im = axb.imshow(J, cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
axb.set_yticks(np.arange(len(PAR)))
axb.set_yticklabels([PRETTY[p] for p in PAR], fontsize=6.4)
ticks, labels = [], []
for i, k in enumerate(keys):
    ticks.append(i * len(days) + len(days) / 2 - 0.5); labels.append(PNAMES[k])
    if i:
        axb.axvline(i * len(days) - 0.5, color="white", lw=1.2)
axb.set_xticks(ticks); axb.set_xticklabels(labels, fontsize=6.2)
axb.set_xlabel("Observable (day 10 to 70 within each block)", labelpad=2,
               fontsize=6.0)
cb = fig.colorbar(im, ax=axb, fraction=0.030, pad=0.02, aspect=18)
cb.ax.set_title("sens.", fontsize=F.SZ["small"], pad=3)
cb.ax.tick_params(labelsize=5.3); cb.outline.set_linewidth(0.5)

# ------------------------------------------------------- c  condition number
axc = fig.add_subplot(gs[0, 2])
cond = [D["sets"][s]["cond"] for s in SETS]
rank = [D["sets"][s]["rank"] for s in SETS]
yy = np.arange(len(SETS))[::-1]
for y, c, r, col in zip(yy, cond, rank, COLS):
    if r < len(PAR):
        axc.barh(y, 1e4, color=K["pale"], height=0.6)
        axc.text(1.3e2, y, "rank %d/%d" % (r, len(PAR)), va="center",
                 fontsize=6.0, color=K["bad"])
    else:
        axc.barh(y, c, color=col, height=0.6)
        axc.text(c * 1.25, y, "%.0f" % c, va="center", fontsize=6.0)
axc.set_xscale("log")
axc.set_yticks(yy)
axc.set_yticklabels([s.replace(" + ", "\n+ ") for s in SETS], fontsize=5.8)
axc.text(0.98, 0.04, "cond. varies 5 to 9 across seeds",
         transform=axc.transAxes, ha="right", fontsize=F.SZ["small"] - 0.5,
         color=K["grey"], style="italic")
axc.set_xlabel("Condition number")
axc.set_xlim(1, 2e5)
axc.spines["left"].set_visible(False); axc.tick_params(axis="y", length=0)

# ---------------------------------------------------- d  replicates needed
axd = fig.add_subplot(gs[1, 0])
w = 0.20
for i, (lbl, col) in enumerate(zip(SETS[1:], COLS[1:])):
    n = np.array(D["sets"][lbl]["nrep"])
    axd.bar(np.arange(len(PAR)) + (i - 1) * w, np.maximum(n, 1),
            width=w, color=col, lw=0, label=lbl)
axd.set_yscale("log")
axd.set_xticks(np.arange(len(PAR)))
axd.set_xticklabels([PRETTY[p] for p in PAR])
axd.set_ylabel("Replicates for $\\pm$20%")
axd.axhline(10, color=K["grey"], lw=0.6, ls=(0, (2, 2)))
axd.text(0.99, 0.03, "feasible below", transform=axd.transAxes,
         fontsize=F.SZ["small"], color=K["grey"], ha="right")
axd.legend(frameon=False, handlelength=1.2, borderpad=0.1, labelspacing=0.22,
           handletextpad=0.4, loc="upper left", fontsize=5.5)

# ------------------------------------------- e  the flat direction in detail
axe = fig.add_subplot(gs[1, 1])
se = np.array(D["sets"]["N trajectory"]["se"])
se2 = np.array(D["sets"]["N + per-unit supply"]["se"])
yy = np.arange(len(PAR))[::-1]
axe.barh(yy + 0.18, se, height=0.34, color=K["liver"], lw=0, label="N trajectory")
axe.barh(yy - 0.18, se2, height=0.34, color=K["kidney"], lw=0, label="N + per-unit supply")
axe.set_xscale("log")
axe.set_yticks(yy); axe.set_yticklabels([PRETTY[p] for p in PAR])
axe.set_xlabel("Standard error of $\\log$ parameter")
axe.axvline(0.20, color=K["ink"], lw=0.7, ls=(0, (3, 2)))
axe.text(0.23, len(PAR) - 0.6, "$\\pm$20%", fontsize=6.0, color=K["ink"])
axe.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2,
           fontsize=F.SZ["small"])
axe.spines["left"].set_visible(False); axe.tick_params(axis="y", length=0)

# ------------------------------------------------- f  what this prescribes
axf = fig.add_subplot(gs[1, 2])
axf.set_xlim(0, 100); axf.set_ylim(0, 100); axf.axis("off")
axf.text(2, 99, "Prescription", fontsize=6.4, fontweight="bold", va="top")
msgs = [("Count units over time", "rank 1 to 5 of 5", K["good"]),
        ("Add per-unit perfusion", "condition 85 to 10.1", K["good"]),
        
        ("$\\beta$ is the costliest", "3 replicates vs 1 for $\\nu$", K["warn"])]
y = 85
for name, note, col in msgs:
    axf.add_patch(FancyBboxPatch((2, y - 12), 96, 12.5,
                                 boxstyle="round,pad=0.25,rounding_size=1",
                                 fc=K["box"], ec=col, lw=0.8))
    axf.text(6, y - 3.6, name, fontsize=6.2, va="center")
    axf.text(6, y - 9.0, note, fontsize=5.8, va="center", color=col)
    y -= 16.5

F.panel_grid(fig, [[axa, axb, axc], [axd, axe, axf]],
             ["a", "b", "c", "d", "e", "f"])
fig.savefig(str(OUTDIR / "Fig5.pdf"), dpi=600, bbox_inches="tight")
fig.savefig(str(OUTDIR / "Fig5.png"), dpi=400, bbox_inches="tight")
print("ok")
