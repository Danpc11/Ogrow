"""Figure 1 (framework and regimes) and the mechanism figure of the paper."""
import glob
import json

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from scipy import stats

from organgrow import figstyle as F
from organgrow.config import ROOT

F.apply()
RES, OUT = ROOT / "results", ROOT / "figures" / "out"
OUT.mkdir(parents=True, exist_ok=True)
K = F.C
REG = [("liver", "Liver", K["liver"]), ("kidney", "Kidney", K["kidney"]),
       ("pancreas", "Pancreas", K["pancreas"])]
if not (RES / "paper" / "snapshot_liver.npz").exists():
    raise SystemExit("framework figures skipped: run scripts/17_paper_runs.py snapshot")
SD = ROOT / "source_data"

# ------------------------------------------------------------------ Figure 1
d = np.load(RES / "paper" / "snapshot_liver.npz")
fig = plt.figure(figsize=(F.WIDTH_FULL, 64 * F.MM))
gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.0, 0.95], wspace=0.32,
                      left=0.02, right=0.985, top=0.9, bottom=0.17)

ax = fig.add_subplot(gs[0, 0])
scale = 7.4                                   # unit spacing: tree units to tissue units
st, en, r = d["starts"] * scale, d["ends"] * scale, d["radius"]
w = 0.25 + 2.2 * (r - r.min()) / max(np.ptp(r), 1e-9)
for a, b, lw in zip(st, en, w):
    ax.plot([a[0], b[0]], [a[1], b[1]], color=K["ink"], lw=lw, alpha=0.55, solid_capstyle="round")
U = d["U"]
ax.scatter(U[:, 0], U[:, 1], s=2.5, color=K["liver"], alpha=0.35, lw=0)
tip_xy = d["ends"][d["tips"]] * scale
stranded = (d["q_tip"] <= d["q_set"]) & (d["P_tip"] > 0)
ax.scatter(tip_xy[stranded, 0], tip_xy[stranded, 1], s=7, color=K["kidney"], lw=0)
ax.set_aspect("equal"); ax.axis("off")
ax.text(0.0, 1.02, "units", color=K["liver"], transform=ax.transAxes,
        fontsize=F.SZ["small"], va="bottom")
ax.text(0.22, 1.02, "tips holding progenitors below the threshold", color=K["kidney"],
        transform=ax.transAxes, fontsize=F.SZ["small"], va="bottom")

ax = fig.add_subplot(gs[0, 1]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
def box(x, y, txt, col):
    ax.add_patch(FancyBboxPatch((x - 0.48, y - 0.09), 0.96, 0.18, boxstyle="round,pad=0.01",
                                fc="white", ec=col, lw=0.8))
    ax.text(x, y, txt, ha="center", va="center", fontsize=F.SZ["small"])
box(0.5, 0.86, "tip $i$\nsupply per unit $q_i$, progenitors $P_i$", K["ink"])
box(0.5, 0.53, r"nucleation $\lambda_i = \nu P_i\,(q_i - q_{\mathrm{set}})_+$", K["ink"])
box(0.5, 0.2, "stops when $q_i \\leq q_{\\mathrm{set}}$ (supply bound)\nor $P_i = 0$ (progenitor bound)", K["ink"])
for y0, y1 in [(0.77, 0.62), (0.44, 0.29)]:
    ax.annotate("", xy=(0.5, y1), xytext=(0.5, y0),
                arrowprops=dict(arrowstyle="-|>", lw=0.8, color=K["ink"]))

ax = fig.add_subplot(gs[0, 2])
rows = []
for reg, name, col in REG:
    fs = sorted(glob.glob(str(RES / "transmission" / f"{reg}_baseline_k*.json")))
    x = [json.load(open(f)) for f in fs]
    pl = np.array([v["pool_final"] / v["P0"] for v in x]); qs = np.array([v["q_final"] / v["q_set"] for v in x])
    ax.scatter(pl, qs, s=10, color=col, lw=0, alpha=0.85)
    lab = {"liver": (0.62, 0.86, "right"), "kidney": (-0.02, 1.135, "left"), "pancreas": (0.24, 1.33, "left")}[reg]
    ax.text(lab[0], lab[1], name, color=col, ha=lab[2], fontsize=F.SZ["small"])
    rows += [[reg, i, a, b] for i, (a, b) in enumerate(zip(pl, qs))]
ax.axhline(1, color=K["grey"], lw=0.6, ls=(0, (2, 2)))
ax.set_xlabel("Progenitor pool left / initial"); ax.set_ylabel("Supply per unit / threshold")
ax.set_xlim(-0.05, 1.0)
F.panel_grid(fig, [[fig.axes[0], fig.axes[1], ax]], ["a", "b", "c"])
fig.savefig(str(OUT / "Fig1_closed.pdf")); fig.savefig(str(OUT / "Fig1_closed.png"), dpi=300)
plt.close(fig)
import csv
with open(SD / "Fig1closed_regime_map.csv", "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(["regime", "replicate", "pool_left_fraction", "supply_over_threshold"]); wr.writerows(rows)

# ------------------------------------------------------- mechanism figure
E = json.load(open(RES / "elasticity.json"))
QS = {"ref": 0.6, "q_lo": 0.5, "q_hi": 0.7, "nu_lo": 0.6, "nu_hi": 0.6}
def predicted(reg, st):
    pr = []
    for f in glob.glob(str(RES / "elasticity" / f"{reg}_{st}_P0*_k*.npz")):
        z = np.load(f); q, P, N = z["q_tip"], z["P_tip"], z["N"][-1]; qs = QS[st]
        pr.append(-qs * P[(q <= qs) & (q > qs - 0.05)].sum() / 0.05 / N)
    return float(np.mean(pr))

fig, ax = plt.subplots(1, 3, figsize=(F.WIDTH_FULL, 60 * F.MM))
fig.subplots_adjust(left=0.07, right=0.985, bottom=0.2, top=0.88, wspace=0.42)
bins = np.linspace(0, 2.5, 41)
for reg, name, col in REG[::2]:
    q, P = [], []
    for f in glob.glob(str(RES / "elasticity" / f"{reg}_ref_P0*_k*.npz")):
        z = np.load(f); q += list(z["q_tip"] / 0.6); P += list(z["P_tip"])
    q, P = np.array(q), np.array(P)
    ax[0].hist(q, bins=bins, weights=P / max(P.sum(), 1e-9), histtype="step", color=col, lw=1.1)
    ax[0].text(0.97, 0.92 - 0.1 * (reg == "pancreas"), name, color=col, ha="right",
               transform=ax[0].transAxes, fontsize=F.SZ["small"])
ax[0].axvline(1, color=K["ink"], lw=0.7)
ax[0].axvspan(1 - 0.05 / 0.6, 1, color=K["grey"], alpha=0.25, lw=0)
ax[0].set_xlabel("Tip supply / threshold"); ax[0].set_ylabel("Share of remaining progenitors")

mrows = []
for reg, name, col in REG:
    pr, me, se = [], [], []
    for st in QS:
        k = f"{reg}_{st}"
        if k not in E:
            continue
        pr.append(predicted(reg, st)); me.append(E[k]["q_set"]["mean"]); se.append(E[k]["q_set"]["se"])
        mrows.append([reg, st, pr[-1], me[-1], se[-1]])
    pr, me, se = map(np.array, (pr, me, se))
    ax[1].errorbar(pr, me, yerr=se, fmt="o", ms=3.5, color=col, lw=0.8, capsize=0)
    if reg != "kidney":
        g = np.polyfit(pr, me, 1); xx = np.linspace(pr.min(), pr.max(), 10)
        ax[1].plot(xx, np.polyval(g, xx), color=col, lw=0.8)
        rr, pp = stats.pearsonr(pr, me)
        ax[1].text(0.03, 0.97 - 0.09 * (reg == "pancreas"), f"{name}: r = {rr:.2f}, P = {pp:.3f}",
                   color=col, transform=ax[1].transAxes, va="top", fontsize=F.SZ["small"])
ax[1].text(0.03, 0.79, "Kidney: no stranded pool", color=K["kidney"], transform=ax[1].transAxes,
           va="top", fontsize=F.SZ["small"])
lim = [-3.6, 0.3]
ax[1].plot(lim, lim, color=K["grey"], lw=0.6, ls=(0, (2, 2)))
ax[1].set_xlim(lim); ax[1].set_ylim(lim)
ax[1].set_xlabel(r"Predicted $\varepsilon_{q_{\mathrm{set}}}$"); ax[1].set_ylabel(r"Measured $\varepsilon_{q_{\mathrm{set}}}$")

for i, (reg, name, col) in enumerate(REG):
    v = [E[f"{reg}_{st}"]["P0"] for st in QS if f"{reg}_{st}" in E]
    xs = i + np.linspace(-0.22, 0.22, len(v))
    ax[2].errorbar(xs, [a["mean"] for a in v], yerr=[a["se"] for a in v], fmt="o", ms=3.5,
                   color=col, lw=0.8, capsize=0)
for y in (0, 1):
    ax[2].axhline(y, color=K["grey"], lw=0.6, ls=(0, (2, 2)))
ax[2].set_xticks(range(3)); ax[2].set_xticklabels([r[1] for r in REG])
ax[2].set_ylabel(r"$\varepsilon_{P_0}$, five settings"); ax[2].set_ylim(-0.1, 1.15)
F.panel_grid(fig, [list(ax)], ["a", "b", "c"])
fig.savefig(str(OUT / "FigMechanism_closed.pdf")); fig.savefig(str(OUT / "FigMechanism_closed.png"), dpi=300)
plt.close(fig)
with open(SD / "FigMechanism_elasticity.csv", "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(["regime", "setting", "predicted_eps_setpoint", "measured_eps_setpoint", "se"]); wr.writerows(mrows)
print("wrote Fig1_closed and FigMechanism_closed")
