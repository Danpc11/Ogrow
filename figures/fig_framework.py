"""Figure 1 (framework and regimes) and the mechanism figure of the paper."""
import glob
import json

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

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

# ------------------------------------------------- Figure 5: the pool at closure
import glob as _glob
E = json.load(open(RES / "elasticity.json"))
Rb = json.load(open(RES / "robustness.json"))
QS = ["ref", "nu_lo", "nu_hi", "q_lo", "q_hi"]
fig, ax = plt.subplots(1, 3, figsize=(F.WIDTH_FULL, 62 * F.MM),
                       gridspec_kw=dict(width_ratios=[0.8, 0.9, 0.9]))
fig.subplots_adjust(left=0.07, right=0.93, bottom=0.27, top=0.88, wspace=0.95)
# the same four colours Figure 1e uses for the same four shares
parts = [("spent", "spent", K["ink"]), ("supply_held", "held by supply", K["supply"]),
         ("space_held", "held by room", K["grey"]), ("active", "active", K["pale"])]
crow = []
for i, (reg, name, col) in enumerate(REG):
    runs = [json.load(open(f)) for f in sorted(_glob.glob(str(RES / "robustness" / f"reference_{reg}_s*_0.json")))]
    bottom = 0.0
    for key, lab, c in parts:
        v = float(np.mean([r[key] for r in runs]))
        ax[0].bar(i, v, bottom=bottom, color=c, width=0.62, lw=0)
        bottom += v
    crow += [[reg, k, r["spent"], r["supply_held"], r["space_held"], r["active"]] for k, r in enumerate(runs)]
ax[0].set_xticks(range(3)); ax[0].set_xticklabels([r[1] for r in REG])
ax[0].set_ylabel("Share of the initial progenitor pool"); ax[0].set_ylim(0, 1.05)
for j, (key, lab, c) in enumerate(parts):
    ax[0].text(0.5 + 0.0 * j, -0.22 - 0.075 * j, lab, color=c, transform=ax[0].transAxes,
               ha="center", fontsize=F.SZ["small"])

for i, (reg, name, col) in enumerate(REG):
    v = [E[f"{reg}_{st}"]["P0"] for st in QS if f"{reg}_{st}" in E]
    xs = i + np.linspace(-0.22, 0.22, len(v))
    ax[1].errorbar(xs, [a["mean"] for a in v], yerr=[a["se"] for a in v], fmt="o", ms=3.5,
                   color=col, lw=0.8, capsize=0)
for y in (0, 1):
    ax[1].axhline(y, color=K["grey"], lw=0.6, ls=(0, (2, 2)))
ax[1].set_xticks(range(3)); ax[1].set_xticklabels([r[1] for r in REG])
ax[1].set_ylabel(r"$\varepsilon_{P_0}$, five settings"); ax[1].set_ylim(-0.15, 1.2)

order = ["reference", "root_surface", "root_random", "no_deficit", "no_load",
         "no_crowding", "no_mesenchyme", "territory", "binomial", "binomial_half_day",
         "binomial_quarter_day", "competing_fates", "no_self_renewal"]
labels = ["reference", "root at surface", "root at random", "no deficit sprouting",
          "no load sprouting", "no crowding sprouting", "no mesenchymal signal",
          "territory assignment", "binomial nucleation", "binomial, half-day step",
          "binomial, quarter-day step", "competing fates",
          "liver pool consumable"]
rrow = []
for j, v in enumerate(order):
    L, Pn = Rb[v]["liver"], Rb[v]["pancreas"]
    c = Rb[v].get("contrast")
    d = c["mean"] if c else Pn["eps_pool"] - L["eps_pool"]
    se = c["se"] if c else float(np.hypot(Pn["eps_se"], L["eps_se"]))
    ax[2].errorbar(d, j, xerr=se, fmt="o", ms=3.5, color=K["ink"], lw=0.8, capsize=0)
    ax[2].text(1.03, j, f"z = {d / se:.1f}", va="center", fontsize=F.SZ["small"] - 0.5,
               transform=ax[2].get_yaxis_transform())
    row = [v, Pn["eps_pool"], Pn["eps_se"], L["eps_pool"], L["eps_se"],
           d, se, d / se, c["n"] if c else ""]
    for a, b in (("kidney", "liver"), ("pancreas", "kidney")):
        e = Rb[v].get(f"{a}_minus_{b}")
        row += [e["mean"], e["se"], e["z"]] if e else ["", "", ""]
    rrow.append(row)
ax[2].axvline(0, color=K["grey"], lw=0.6, ls=(0, (2, 2)))
ax[2].set_yticks(range(len(order))); ax[2].set_yticklabels(labels, fontsize=F.SZ["small"])
ax[2].invert_yaxis(); ax[2].set_xlim(-0.5, 2.1)
ax[2].set_xlabel(r"$\varepsilon_{P_0}$, pancreas minus liver")
F.panel_grid(fig, [list(ax)], ["a", "b", "c"])
fig.savefig(str(OUT / "Fig5_pool.pdf")); fig.savefig(str(OUT / "Fig5_pool.png"), dpi=300)
plt.close(fig)
with open(SD / "Fig5_closure_shares.csv", "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(["regime", "seed", "spent", "held_by_supply", "held_by_room", "active"]); wr.writerows(crow)
with open(SD / "Fig5_robustness_contrast.csv", "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(["variant", "eps_pancreas", "se_pancreas", "eps_liver", "se_liver",
                                     "paired_difference", "paired_se", "z", "n_seeds",
                                     "kidney_minus_liver", "se", "z",
                                     "pancreas_minus_kidney", "se", "z"]); wr.writerows(rrow)
print("wrote Fig1_closed and Fig5_pool")
