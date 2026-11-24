"""Figure: injury after closure, and the two things recovery needs."""
import csv
import glob
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from organgrow import figstyle as F
from organgrow.config import ROOT

F.apply()
RES, OUT, SD = ROOT / "results", ROOT / "figures" / "out", ROOT / "source_data"
OUT.mkdir(parents=True, exist_ok=True); SD.mkdir(parents=True, exist_ok=True)
K = dict(ink=F.C["ink"], grey=F.C["grey"], liver=F.C["liver"],
         kidney=F.C["kidney"], warm=F.C["cells"])
CLOSE, FRAC = 140, 0.3
CONDS = ("as grown", "reactivated", "reactivated, fixed inflow")
SHORT = ("as grown", "reactivated", "reactivated,\nfixed inflow")


def load(reg, cond, frac=FRAC, kind="parenchyma"):
    t = f"{reg}_{kind}_{cond.replace(' ', '').replace(',', '')}_f{frac}"
    return [json.load(open(f)) for f in sorted(glob.glob(str(RES / "injury" / f"{t}_s*.json")))]


fig = plt.figure(figsize=(F.WIDTH_FULL, 104 * F.MM))
gs = fig.add_gridspec(2, 2, left=0.075, right=0.985, bottom=0.08, top=0.93,
                      wspace=0.32, hspace=0.60)
ax = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]),
      fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]

# A: the liver's mass through the injury, under the three conditions
rows_a = []
for cond, col in zip(CONDS, (K["grey"], K["warm"], K["liver"])):
    Y = np.array([r["series"]["mass"] for r in load("liver", cond)], float)
    Y = Y / Y[:, CLOSE - 2][:, None]
    y, lo, hi = Y.mean(0), Y.min(0), Y.max(0)
    t = np.arange(1, Y.shape[1] + 1)
    ax[0].fill_between(t, lo, hi, color=col, alpha=0.18, lw=0)
    ax[0].plot(t, y, "-", lw=1.1, color=col)
    rows_a += [["liver", cond, int(a), float(b), float(c), float(d)]
               for a, b, c, d in zip(t, y, lo, hi)]
ax[0].axvline(CLOSE, color=K["ink"], lw=0.6, ls=(0, (2, 2)))
ax[0].axhline(1, color=K["grey"], lw=0.5, ls=(0, (1, 2)))
ax[0].text(CLOSE + 3, 0.63, "30 % lost", fontsize=5.5, color=K["ink"])
ax[0].set_xlabel("Day"); ax[0].set_ylabel("Mass, relative to before injury")
ax[0].set_xlim(100, 220); ax[0].set_ylim(0.6, 2.1)
ax[0].set_title("Liver", fontsize=7, color=K["liver"])
ax[0].legend(handles=[Patch(facecolor=c, label=l) for l, c in
                      zip(SHORT, (K["grey"], K["warm"], K["liver"]))],
             frameon=False, fontsize=5.5, loc="upper left", handlelength=1.1,
             handleheight=0.8, borderpad=0.1)

# B: where each readout ends, by condition
I = json.load(open(RES / "injury.json"))
rows_b = []
for i, key in enumerate(("units", "mass", "mass_per_unit", "perfused")):
    for j, cond in enumerate(CONDS):
        d = I[f"liver_{cond}"][key]
        x = i + (j - 1) * 0.26
        ax[1].bar(x, d["mean"], width=0.24,
                  color=[K["grey"], K["warm"], K["liver"]][j], lw=0)
        ax[1].plot(np.full(len(d["seeds"]), x), d["seeds"], ".", ms=1.4,
                   color=K["ink"], alpha=0.5, mew=0)
        rows_b += [["liver", cond, key, float(v)] for v in d["seeds"]]
ax[1].axhline(1, color=K["ink"], lw=0.6, ls=(0, (2, 2)))
ax[1].set_xticks(range(4))
ax[1].set_xticklabels(["units", "mass", "mass per\nunit", "perfused\nmass"], fontsize=5.5)
ax[1].set_ylabel("Day 220, relative to before injury"); ax[1].set_ylim(0, 2.25)
ax[1].set_title("Liver, by condition", fontsize=7, color=K["liver"])

# C: the kidney cannot make nephrons, so it enlarges the ones it kept
rows_c = []
for i, frac in enumerate((0.3, 0.5)):
    d = I[f"kidney_compensation_{frac}"]
    for j, (lab, key) in enumerate((("nephrons kept", "units"),
                                    ("mass per nephron", "mass_per_unit"),
                                    ("full compensation", "full"))):
        x = i + (j - 1) * 0.26
        seeds = d[key]["seeds"] if key != "full" else None
        v = d[key]["mean"] if key != "full" else d["full"]
        ax[2].bar(x, v, width=0.24,
                  color=[K["kidney"], K["warm"], K["grey"]][j], lw=0,
                  hatch="" if j < 2 else "///", edgecolor=K["grey"])
        if seeds:
            ax[2].plot(np.full(len(seeds), x), seeds, ".", ms=1.4,
                       color=K["ink"], alpha=0.5, mew=0)
            rows_c += [["kidney", frac, lab, float(w)] for w in seeds]
        else:
            rows_c.append(["kidney", frac, lab, float(v)])
ax[2].axhline(1, color=K["ink"], lw=0.6, ls=(0, (2, 2)))
ax[2].set_xticks(range(2)); ax[2].set_xticklabels(["30 % lost", "50 % lost"], fontsize=6)
ax[2].set_ylim(0, 2.25); ax[2].set_ylabel("Day 220, relative to before injury")
ax[2].set_title("Kidney: no new nephrons", fontsize=7, color=K["kidney"])
ax[2].legend(handles=[Patch(facecolor=c, label=l, hatch=h, edgecolor=K["grey"]) for l, c, h in
                      (("nephrons kept", K["kidney"], ""),
                       ("mass per nephron", K["warm"], ""),
                       ("full compensation", K["grey"], "///"))],
             frameon=False, fontsize=5.5, loc="upper left", handlelength=1.1,
             handleheight=0.8, borderpad=0.1)

# D: where a construct lands when the delivery it is given is varied, which
#    is the part of the fixed-inflow result that is not arithmetic
rows_d = []
_law = RES / "inflow_law.json"
if _law.exists():
    L = json.load(open(_law))["liver"]
    xs = np.array([r["scale"] for r in L["rows"]], float)
    for r in L["rows"]:
        ax[3].plot(np.full(len(r["seeds"]), r["scale"]), r["seeds"], ".",
                   ms=1.8, color=K["ink"], alpha=0.45, mew=0)
        rows_d += [["liver", r["scale"], float(v)] for v in r["seeds"]]
    ys = np.array([r["mass_rel"] for r in L["rows"]], float)
    es = np.array([r["mass_rel_sd"] for r in L["rows"]], float)
    ax[3].errorbar(xs, ys, yerr=es, fmt="o-", ms=3.2, lw=1.1, color=K["liver"],
                   elinewidth=0.7, capsize=1.5, mew=0)
    ax[3].plot([xs.min(), xs.max()], [xs.min(), xs.max()], ls=(0, (2, 2)),
               lw=0.8, color=K["grey"])
    ax[3].text(xs.max(), xs.max() * 1.03, "budget $Q/\\phi$", fontsize=5.5,
               color=K["grey"], ha="right", va="bottom")
    ax[3].text(0.04, 0.93, f"mass $\\propto Q^{{{L['slope']:.2f}}}$",
               transform=ax[3].transAxes, fontsize=6, color=K["liver"])
    ax[3].set_xlabel("Delivered flow, relative to before injury")
    ax[3].set_ylabel("Day 220 mass, relative\nto before injury")
    ax[3].set_title("Liver: the mass follows the supply", fontsize=7,
                    color=K["liver"])
else:
    ax[3].axis("off")
    ax[3].text(0.5, 0.5, "run scripts/32_inflow_law.py", ha="center",
               va="center", fontsize=6, color=K["grey"])

F.panel_grid(fig, [ax[:2], ax[2:]], ["a", "b", "c", "d"])
F.audit(fig, "Fig8_injury")
fig.savefig(str(OUT / "Fig8_injury.pdf")); fig.savefig(str(OUT / "Fig8_injury.png"), dpi=300)
plt.close(fig)
for name, hdr, rows in (("Fig8_trajectories", ["regime", "condition", "day", "mean_relative_mass",
                                             "min_over_seeds", "max_over_seeds"], rows_a),
                        ("Fig8_liver_recovery", ["regime", "condition", "readout", "relative_per_seed"], rows_b),
                        ("Fig8_kidney_compensation",
                         ["regime", "fraction_lost", "readout", "relative_per_seed"], rows_c),
                        ("Fig8_inflow_law",
                         ["regime", "delivered_flow_relative", "mass_relative_per_seed"], rows_d)):
    with open(SD / f"{name}.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(hdr); w.writerows(rows)
print("wrote Fig8_injury")
