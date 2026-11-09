"""Figure 1 — a generative engine in which organ identity is one parameter."""
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle

from organgrow import figstyle as F
from organgrow.config import ROOT

F.apply()
RES = ROOT / "results"
OUTDIR = ROOT / "figures" / "out"
OUTDIR.mkdir(parents=True, exist_ok=True)
K, SZ = F.C, F.SZ


def sdf2(kind, x, y):
    """Cross-section of the tissue domain; negative inside."""
    if kind == "liver":
        a = np.hypot((x - 30) / 1.45, (y + 5) / 0.88) - 105
        b = np.hypot((x + 92) / 0.95, (y + 12) / 0.70) - 62
        d = np.minimum(a, b)
        notch = np.hypot((x + 22) / 0.62, (y + 96) / 1.05) - 54
        return np.maximum(np.maximum(d, -notch), y - 62)
    if kind == "kidney":
        e = np.hypot(x / 0.78, y / 1.30) - 100
        cut = np.hypot((x + 92) / 0.80, y / 0.92) - 72
        return np.maximum(e, -cut)
    if kind == "lung":
        e = np.hypot(x / 0.92, (y + 10) / 1.45) - 115
        hil = np.hypot((x + 86) / 0.70, (y + 20) / 0.85) - 58
        return np.maximum(np.maximum(e, -hil), -(y + 150))
    t = np.clip(x / 165.0, -1, 1)
    cy = 34 * np.sin(t * 1.35) - 6
    rad = 50 * (1 - 0.52 * (t * 0.5 + 0.5))
    return np.maximum(np.abs(y - cy) - rad, np.abs(x) - 165)


ORG = [("liver", "Liver", K["liver"]), ("kidney", "Kidney", K["kidney"]),
       ("lung", "Lung", K["lung"]), ("pancreas", "Pancreas", K["pancreas"])]

fig = plt.figure(figsize=(F.WIDTH_FULL, 106 * F.MM))
gs = fig.add_gridspec(3, 4, height_ratios=[0.60, 0.86, 1.0],
                      hspace=0.34, wspace=0.40,
                      left=0.055, right=0.985, top=0.925, bottom=0.095)

# ------------------------------------------------------------ a  the loop
axa = fig.add_subplot(gs[0, :])
axa.set_xlim(-1, 101); axa.set_ylim(-9.5, 15); axa.axis("off")
steps = [("Update\nsupply", "$\\mathbb{C},\\ c_u$"),
         ("Update\nprogenitor pool", "$P$"),
         ("Branch\nprimary tree", "topology"),
         ("Nucleate units", "$N$"),
         ("Extend\nsecond tree", "vessels"),
         ("Grow units", "$m_j$")]
w, h, y0 = 13.6, 9.4, 3.2
for i, (lab, out) in enumerate(steps):
    x = 1.0 + i * 16.5
    hot = (i == 1)
    axa.add_patch(FancyBboxPatch((x, y0), w, h,
                                 boxstyle="round,pad=0.28,rounding_size=0.9",
                                 fc="#fdf6e3" if hot else K["box"],
                                 ec=K["pancreas"] if hot else K["grey"],
                                 lw=1.0 if hot else 0.6, zorder=2))
    axa.text(x + w / 2, y0 + h * 0.64, lab, ha="center", va="center",
             fontsize=SZ["annot"], zorder=3)
    axa.text(x + w / 2, y0 + h * 0.20, out, ha="center", va="center",
             fontsize=SZ["small"], color=K["grey"], zorder=3)
    if i < 5:
        axa.add_patch(FancyArrowPatch((x + w + 0.3, y0 + h / 2),
                                      (x + 16.5 - 0.3, y0 + h / 2),
                                      arrowstyle="-|>", mutation_scale=6,
                                      lw=0.7, color=K["ink"], zorder=1))
xr, xl = 1.0 + 5 * 16.5 + w / 2, 1.0 + w / 2
axa.plot([xr, xr, xl, xl], [y0 - 0.3, -3.2, -3.2, y0 - 0.9], color=K["grey"],
         lw=0.7, ls=(0, (3, 2)), zorder=1, clip_on=False)
axa.add_patch(FancyArrowPatch((xl, -3.2 + 0.6), (xl, y0 - 0.9),
                              arrowstyle="-|>", mutation_scale=6, lw=0.7,
                              color=K["grey"], zorder=1))
axa.text(50, -6.2, "accumulated demand feeds back into supply", ha="center",
         va="top", fontsize=SZ["small"], color=K["grey"])
axa.text(1.0 + 16.5 + w / 2, y0 + h + 1.0, "organ identity enters here",
         ha="center", va="bottom", fontsize=SZ["annot"],
         color=K["pancreas"], fontweight="bold")

# ---------------------------------------------------------- b  the domains
axb = []
for j, (kind, name, col) in enumerate(ORG):
    ax = fig.add_subplot(gs[1, j]); axb.append(ax)
    X, Y = np.meshgrid(np.linspace(-215, 215, 400), np.linspace(-150, 150, 300))
    D = sdf2(kind, X, Y)
    ax.contourf(X, Y, np.clip(D, -90, 90), levels=18, cmap="Greys_r", alpha=0.26)
    ax.contour(X, Y, D, levels=[0], colors=[col], linewidths=1.3)
    ax.set_xlim(-215, 215); ax.set_ylim(-150, 150)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_color(K["pale"]); sp.set_linewidth(0.5)
    ax.text(0.5, -0.05, name, transform=ax.transAxes, ha="center", va="top",
            fontsize=SZ["annot"], color=col)
axb[0].text(0.03, 0.05, "$\\phi(\\mathbf{x})<0$", transform=axb[0].transAxes,
            fontsize=SZ["small"])

# -------------------------------------------------------- c  the local rules
axc = fig.add_subplot(gs[2, 0:2])
axc.set_xlim(0, 100); axc.set_ylim(0, 100); axc.axis("off")
axc.text(0, 100, "Nucleation", fontsize=SZ["annot"], fontweight="bold", va="top")
axc.text(0, 90, "$\\dot N=\\nu\\,(\\mathbb{C}/N-c_{set})$", fontsize=SZ["label"],
         va="top")
ax1 = axc.inset_axes([0.00, 0.44, 0.27, 0.26])
cu = np.linspace(0, 1, 120)
ax1.plot(cu, np.maximum(cu - 0.42, 0), color=K["ink"], lw=1.0)
ax1.axvline(0.42, color=K["grey"], lw=0.6, ls=(0, (2, 2)))
ax1.text(0.46, 0.38, "$c_{set}$", fontsize=SZ["small"], color=K["grey"])
ax1.set_xticks([]); ax1.set_yticks([])
ax1.set_xlabel("$\\mathbb{C}/N$", fontsize=SZ["small"], labelpad=1)
ax1.set_ylabel("$\\dot N$", fontsize=SZ["small"], labelpad=1)
for sp in ax1.spines.values():
    sp.set_linewidth(0.5)

axc.text(52, 100, "Spatial exclusion", fontsize=SZ["annot"], fontweight="bold",
         va="top")
cx, cy = 72, 72
rng = np.random.default_rng(4)
pts = []
while len(pts) < 9:
    p = rng.uniform(-13, 13, 2)
    if all((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 > 46 for q in pts):
        pts.append(p)
for p in pts:
    axc.add_patch(Circle((cx + p[0], cy + p[1] * 0.7), 2.5, fc=K["grey"],
                         ec="none", alpha=0.8))
axc.add_patch(Circle((cx + 16, cy + 10), 2.5, fc="none", ec=K["grey"], lw=0.7,
                     ls=(0, (2, 1.6))))
axc.text(cx + 16, cy + 15, "rejected", fontsize=SZ["small"], color=K["grey"],
         ha="center")
_p, _q = np.array(pts[0]), np.array(pts[1])
_a = np.array([cx + _p[0], cy + _p[1] * 0.7])
_b = np.array([cx + _q[0], cy + _q[1] * 0.7])
_u = (_b - _a) / np.linalg.norm(_b - _a)
axc.annotate("", xy=tuple(_b - _u * 2.7), xytext=tuple(_a + _u * 2.7),
             arrowprops=dict(arrowstyle="<|-|>", lw=0.8, color=K["ink"],
                             mutation_scale=5, shrinkA=0, shrinkB=0))
axc.text(cx, cy - 22, "$\\lambda$, minimum unit spacing", fontsize=SZ["small"],
         ha="center")

axc.text(0, 30, "Supply attenuation", fontsize=SZ["annot"], fontweight="bold",
         va="top")
axc.text(0, 20, "$c_{u,j}=c_u\\,e^{-d_j/\\ell_a}$",
         fontsize=SZ["label"], va="top")
ax2 = axc.inset_axes([0.50, 0.01, 0.46, 0.24])
dd = np.linspace(0, 200, 200)
ax2.plot(dd, np.exp(-dd / 58.0), color=K["kidney"], lw=1.0)
ax2.axvline(58.0, color=K["grey"], lw=0.6, ls=(0, (2, 2)))
ax2.text(62, 0.80, "$\\ell_a$", fontsize=SZ["small"], color=K["grey"])
ax2.set_xticks([]); ax2.set_yticks([])
ax2.set_xlabel("distance to nearest vessel", fontsize=SZ["small"], labelpad=1)
ax2.set_ylabel("$c_{u,j}/c_u$", fontsize=SZ["small"], labelpad=1)
for sp in ax2.spines.values():
    sp.set_linewidth(0.5)

# ------------------------------------------------------ d  identity as a knob
axd = fig.add_subplot(gs[2, 2])
axd.set_xlim(0, 100); axd.set_ylim(-8, 100); axd.axis("off")
axd.add_patch(FancyBboxPatch((4, 84), 92, 14,
                             boxstyle="round,pad=0.3,rounding_size=1",
                             fc="#fdf6e3", ec=K["pancreas"], lw=0.9))
axd.text(50, 91, "Window closure, $\\sigma(t)$", ha="center", va="center",
         fontsize=SZ["annot"], fontweight="bold")
axd.plot([4, 96], [76, 76], color=K["pale"], lw=0.7)
rows = [("Liver", "capacity saturation", K["liver"], "contracts"),
        ("Kidney", "progenitor exhaustion", K["kidney"], "no"),
        ("Lung", "progenitor exhaustion", K["lung"], "no"),
        ("Pancreas", "limiting pool", K["pancreas"], "no")]
y = 62
for nm, mech, col, contr in rows:
    axd.add_patch(plt.Rectangle((4, y - 2), 2.4, 11, fc=col, ec="none"))
    axd.text(10, y + 6, nm, fontsize=SZ["annot"], va="center", color=col,
             fontweight="bold")
    axd.text(10, y + 0.5, mech, fontsize=SZ["small"], va="center",
             color=K["grey"])
    axd.text(97, y + 6, contr, fontsize=SZ["small"], va="center", ha="right",
             color=K["good"] if contr == "contracts" else K["grey"],
             fontweight="bold" if contr == "contracts" else "normal")
    y -= 16
axd.text(50, -4, "one parameter; everything else is shared", ha="center",
         fontsize=SZ["small"], color=K["grey"], style="italic")

# ------------------------------------------------------------ e  the cost
axe = fig.add_subplot(gs[2, 3])
cost = json.load(open(RES / "cost.json"))
E = np.array([c["E"] for c in cost], float)
tt = np.array([c["tree"] for c in cost]) * 1e3
ts = np.array([c["sparse"] for c in cost]) * 1e3
m = ~np.isnan(ts)
axe.loglog(E, tt, "o-", color=K["good"], ms=3.0, mew=0, lw=1.0,
           label="Tree traversal")
axe.loglog(E[m], ts[m], "s--", color=K["grey"], ms=3.0, mew=0, lw=1.0,
           label="Sparse solve")
ref = E / E[0] * tt[0]
axe.loglog(E, ref, ":", color=K["ink"], lw=0.7)
axe.text(E[-2] * 1.2, ref[-2] * 0.40, "$\\mathcal{O}(E)$", fontsize=SZ["small"])
axe.set_xlabel("Supply-tree edges, $E$")
axe.set_ylabel("Cost of one solve (ms)")
axe.legend(loc="upper left")
axe.grid(True, which="major", lw=0.3, color=K["pale"], alpha=0.8)
F.tidy(axe)

F.panel_grid(fig, [[axa], axb, [axc, axd, axe]],
             ["a", "b", "", "", "", "c", "d", "e"])
axe.text(0.02, 0.97, "benchmark; not part of the loop",
         transform=axe.transAxes, fontsize=SZ["small"] - 0.5, va="top",
         color=K["grey"], style="italic")
fig.savefig(str(OUTDIR / "Fig1.pdf"))
fig.savefig(str(OUTDIR / "Fig1.png"), dpi=400)
print("wrote Fig1")
