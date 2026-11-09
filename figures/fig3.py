"""Figure 3 — window closure determines whether unit-number variance contracts."""
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

from organgrow import figstyle as F
from organgrow.config import ROOT

F.apply()
RES = ROOT / "results"
OUTDIR = ROOT / "figures" / "out"
OUTDIR.mkdir(parents=True, exist_ok=True)
K = F.C
MM = F.MM


D0, D1 = 0, 69                     # day 1 and day 70 (0-based)

ens = np.load(RES/RES / 'f3_main.npz')
swL = np.load(RES / 'f3_sweep_liver.npz'); swK = np.load(RES / 'f3_sweep_kidney.npz')
sig = swL["sig"]


def cv(a, axis=0):
    return a.std(axis=axis, ddof=1) / a.mean(axis=axis)


fig = plt.figure(figsize=(F.WIDTH_FULL, 104 * F.MM))
gs = fig.add_gridspec(2, 4, height_ratios=[0.90, 1.0], hspace=0.46, wspace=0.52,
                      left=0.065, right=0.985, top=0.945, bottom=0.085)




# ------------------------------------------------------------ a  design
axa = fig.add_subplot(gs[0, 0])
axa.set_xlim(0, 100); axa.set_ylim(0, 100); axa.axis("off")
axa.add_patch(FancyBboxPatch((4, 72), 92, 22, boxstyle="round,pad=0.6,rounding_size=2",
                             fc=K["box"], ec=K["grey"], lw=0.6))
axa.text(50, 88, "Identical input", ha="center", fontsize=6.4, fontweight="bold")
axa.text(50, 79, "same seeds · $\\sigma_{init}=0.6$ · $\\sigma_{param}=0$",
         ha="center", fontsize=6.0, color=K["grey"])
for i, (nm, mech, col) in enumerate([("Saturation", r"$\mathbb{C}/N \to c_{set}$", K["liver"]),
                                     ("Exhaustion", "$P_i\\to 0$", K["kidney"])]):
    x = 8 + i * 48
    axa.add_patch(FancyArrowPatch((50, 71), (x + 18, 56),
                                  connectionstyle="arc3,rad=%.2f" % (0.2 - 0.4 * i),
                                  arrowstyle="-|>", mutation_scale=6, lw=0.7, color=K["grey"]))
    axa.add_patch(FancyBboxPatch((x, 34), 36, 21, boxstyle="round,pad=0.6,rounding_size=2",
                                 fc="none", ec=col, lw=1.0))
    axa.text(x + 18, 49, nm, ha="center", fontsize=6.4, fontweight="bold", color=col)
    axa.text(x + 18, 40, mech, ha="center", fontsize=6.4)
axa.text(50, 62.5, "one parameter differs", ha="center", fontsize=6.2,
         color=K["ink"], style="italic")
axa.add_patch(FancyBboxPatch((4, 6), 92, 20, boxstyle="round,pad=0.6,rounding_size=2",
                             fc=K["box"], ec=K["grey"], lw=0.6))
axa.text(50, 20, "Readout", ha="center", fontsize=6.4, fontweight="bold")
axa.text(50, 11.5, "absolute s.d. of $N$ across replicates, day 1 vs day 70",
         ha="center", fontsize=6.0, color=K["grey"])

# -------------------------------------------------- b  per-replicate values
axb = fig.add_subplot(gs[0, 1])
rng = np.random.default_rng(0)
pos = {"liver": (0.8, 2.0), "kidney": (3.4, 4.6)}
for kind, col in [("liver", K["liver"]), ("kidney", K["kidney"])]:
    M = ens[kind]
    for j, d in enumerate([D0, D1]):
        v = M[:, d] - np.mean(M[:, d])
        x = pos[kind][j] + rng.normal(0, 0.055, len(v))
        axb.scatter(x, v, s=5, c=col, lw=0, alpha=0.70)
        axb.plot([pos[kind][j] - 0.22, pos[kind][j] + 0.22], [0, 0],
                 color=K["ink"], lw=1.1)
        axb.text(pos[kind][j], 21.5, "%.1f" % M[:, d].std(ddof=1),
                 ha="center", va="bottom", fontsize=5.6, color=col)
axb.text(0.05, 24.5, "s.d. of $N$", ha="left", va="bottom", fontsize=5.4, color=K["grey"])
axb.set_xlim(0.0, 5.3); axb.set_ylim(-20, 29)
axb.set_xticks([0.8, 2.0, 3.4, 4.6])
axb.set_xticklabels(["day 1", "day 70", "day 1", "day 70"],
                    fontsize=F.SZ["small"])
axb.set_ylabel("$N$ $-$ mean")
axb.text(1.5, -18, "Saturation", ha="center", fontsize=6.2, color=K["liver"])
axb.text(3.9, -18, "Exhaustion", ha="center", fontsize=6.2, color=K["kidney"])

# ------------------------------------------------------------- c  CV(t)
axc = fig.add_subplot(gs[0, 2])
t = np.arange(1, 71)
for kind, name, col in [("liver", "Saturation", K["liver"]),
                        ("kidney", "Exhaustion (kidney)", K["kidney"]),
                        ("lung", "Exhaustion (lung)", K["lung"]),
                        ("pancreas", "Limiting pool", K["pancreas"])]:
    A = ens[kind]
    axc.plot(t, A.std(axis=0, ddof=1) / A[:, 0].std(ddof=1), color=col, lw=1.1,
             label=name)
axc.axhline(1, color=K["grey"], lw=0.6, ls=(0, (2, 2)))
axc.text(2, 1.05, "", fontsize=5.8, color=K["grey"], ha="left")
axc.set_xlabel("Day")
axc.set_yscale("log")
axc.set_ylabel("s.d. of $N$, relative to day 1")
axc.legend(loc="lower right", fontsize=F.SZ["small"] - 0.5)
axc.set_ylim(0.05, 20)

# ------------------------------------------- d  contraction vs growth ratio
axd = fig.add_subplot(gs[1, 0])
names, meas, lin, quad, cols = [], [], [], [], []
for kind, name, col in [("liver", "Saturation", K["liver"]),
                        ("kidney", "Exhaustion\n(kidney)", K["kidney"]),
                        ("lung", "Exhaustion\n(lung)", K["lung"]),
                        ("pancreas", "Limiting\npool", K["pancreas"])]:
    A = ens[kind]
    r = A[:, D1].mean() / A[:, D0].mean()
    names.append(name); cols.append(col)
    meas.append(cv(A[:, D0]) / cv(A[:, D1])); lin.append(r); quad.append(r ** 2)
x = np.arange(len(names)); w = 0.26
axd.bar(x - w, lin, width=w, color=K["pale"], lw=0, label="dilution, $N_f/N_i$")
axd.bar(x, meas, width=w, color=cols, lw=0, label="measured")
axd.bar(x + w, quad, width=w, color="none", lw=0.9, edgecolor=K["ink"],
        label="bound, $(N_f/N_i)^2$")
axd.set_yscale("log"); axd.set_xticks(x)
short = {"Saturation": "Sat.", "Exhaustion (kidney)": "Exh.\n(kidney)",
         "Exhaustion (lung)": "Exh.\n(lung)", "Limiting pool": "Pool"}
axd.set_xticklabels(["Sat.", "Exh.\n(kidney)", "Exh.\n(lung)", "Pool"][:len(names)],
                    fontsize=F.SZ["small"] - 1)
axd.set_ylim(top=axd.get_ylim()[1] * 8)
axd.set_ylabel("Contraction of CV")
axd.legend(frameon=False, handlelength=1.1, borderpad=0.1, labelspacing=0.22,
           handletextpad=0.4, loc="upper left", fontsize=5.3)

# --------------------------------------------- e  sweep of parameter disorder
axe = fig.add_subplot(gs[1, 1])
for sw, col, name in [(swL, K["liver"], "Saturation"), (swK, K["kidney"], "Exhaustion")]:
    A = sw["N"]
    c = np.array([A[i][:, D1].std(ddof=1) / A[i][:, D0].std(ddof=1)
                  for i in range(len(sig))])
    axe.plot(sig, c, "o-", color=col, ms=3.2, mew=0, lw=1.1, label=name)
axe.axhline(1, color=K["grey"], lw=0.6)
axe.text(0.98, 0.06, "", transform=axe.transAxes,
         fontsize=F.SZ["small"], color=K["grey"], ha="right", va="bottom")
axe.set_xlabel("Parameter disorder, $\\sigma_{param}$")
axe.set_ylabel("Final / initial s.d. of $N$")
axe.legend(loc="upper left")

# ------------------------------------------------------------- f  the floor
axf = fig.add_subplot(gs[1, 2])
for sw, col, name in [(swL, K["liver"], "Saturation"), (swK, K["kidney"], "Exhaustion")]:
    A = sw["N"]
    cf = np.array([cv(A[i][:, D1]) for i in range(len(sig))])  # terminal CV
    axf.plot(sig, cf, "o-", color=col, ms=3.2, mew=0, lw=1.1, label=name)
# the dashed line is the parameter disorder itself, for reference: the terminal
# CV equals it only when the window has closed on the parameter that carries it
axf.plot(sig, sig, ls=(0, (3, 2)), color=K["ink"], lw=0.7)
axf.text(0.055, 0.082, "$\\sigma_{param}$ (reference)", fontsize=6.0,
         color=K["ink"], rotation=36)
axf.set_xlabel("Parameter disorder, $\\sigma_{param}$")
axf.set_ylabel("Terminal CV of $N$")
axf.legend(frameon=False, handlelength=1.4, borderpad=0.15, labelspacing=0.25,
           handletextpad=0.45, loc="upper left")

# ------------------------------------- d  the tree injects its own disorder
vd_path = RES / "vascular_disorder.json"
if not vd_path.exists():
    raise SystemExit("Fig3 needs results/vascular_disorder.json; "
                     "run scripts/11_vascular_disorder.py first")
VD = json.load(open(vd_path))
axg = fig.add_subplot(gs[0, 3])
sw = VD["eta_sweep"]
eta = np.array([r["rate"] for r in sw])
cvg = np.array([r["cv"] for r in sw])
lo = np.array([r.get("cv_lo", r["cv"]) for r in sw])
hi = np.array([r.get("cv_hi", r["cv"]) for r in sw])
axg.fill_between(eta, lo, hi, color=K["kidney"], alpha=0.18, lw=0)
axg.plot(eta, cvg, "o-", color=K["kidney"], ms=3.0, mew=0, lw=1.1)
axg.set_xlabel("Calibre relaxation rate, $\\eta$")
axg.set_ylabel("CV of supply across replicates")
axg.text(0.03, 0.95, "fixed $N$, no seeded disorder:\nall spread is the tree's\n"
         "band: 90% bootstrap over replicates",
         transform=axg.transAxes, fontsize=F.SZ["small"], va="top",
         color=K["grey"], style="italic")
axg.set_ylim(0, max(hi) * 1.25)

# ------------------------------------------ h  it settles, it does not diverge
axh = fig.add_subplot(gs[1, 3])
ot = VD["over_time"]
dd = np.array([r["day"] for r in ot])
axh.plot(dd, [r["cv_C"] for r in ot], "o-", color=K["kidney"], ms=2.6, mew=0,
         lw=1.1, label="supply, $\\mathbb{C}$")
axh.plot(dd, [r["cv_N"] for r in ot], "s-", color=K["liver"], ms=2.6, mew=0,
         lw=1.1, label="unit number, $N$")
axh.set_xlabel("Day"); axh.set_ylabel("CV across replicates")
axh.legend(loc="upper left")
g = VD["gamma"]["liver"]["mean"]
axh.text(0.97, 0.60, "$\\mathbb{C}\\sim N^{%.2f}$, so no runaway" % g,
         transform=axh.transAxes, ha="right", fontsize=F.SZ["small"],
         color=K["grey"], style="italic")

axc.text(0.03, 0.96, "only saturation reads $\\mathbb{C}$,\nso only it is exposed",
         transform=axc.transAxes, ha="left", va="top",
         fontsize=F.SZ["small"], color=K["liver"], style="italic")

F.panel_grid(fig, [[axa, axb, axc, axg], [axd, axe, axf, axh]],
             ["a", "b", "c", "d", "e", "f", "g", "h"])
fig.savefig(str(OUTDIR / "Fig3.pdf"), dpi=600, bbox_inches="tight")
fig.savefig(str(OUTDIR / "Fig3.png"), dpi=400, bbox_inches="tight")
print("ok")
