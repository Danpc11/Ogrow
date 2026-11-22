"""Figure: injury after closure, and the two things recovery needs."""
import csv
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
K = dict(ink="#2b2b2b", grey="#9a9a9a", liver="#c0392b", kidney="#2e6da4",
         warm="#b8860b")
CLOSE, SEEDS, FRAC = 140, (300, 301, 302), 0.3
CONDS = ("as grown", "reactivated", "reactivated, fixed inflow")
SHORT = ("as grown", "reactivated", "reactivated,\nfixed inflow")


def load(reg, cond, frac=FRAC, kind="parenchyma"):
    t = f"{reg}_{kind}_{cond.replace(' ', '').replace(',', '')}_f{frac}"
    return [json.load(open(RES / "injury" / f"{t}_s{s}.json")) for s in SEEDS]


fig, ax = plt.subplots(1, 3, figsize=(F.WIDTH_FULL, 54 * F.MM),
                       gridspec_kw=dict(width_ratios=[1.15, 1.0, 1.0]))
fig.subplots_adjust(left=0.07, right=0.985, bottom=0.3, top=0.86, wspace=0.42)

# A: the liver's mass through the injury, under the three conditions
rows_a = []
for cond, col in zip(CONDS, (K["grey"], K["warm"], K["liver"])):
    y = np.mean([r["series"]["mass"] for r in load("liver", cond)], axis=0)
    y = y / y[CLOSE - 2]
    t = np.arange(1, len(y) + 1)
    ax[0].plot(t, y, "-", lw=1.1, color=col)
    rows_a += [["liver", cond, int(a), float(b)] for a, b in zip(t, y)]
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
        v = I[f"liver_{cond}"][key]
        ax[1].bar(i + (j - 1) * 0.26, v, width=0.24,
                  color=[K["grey"], K["warm"], K["liver"]][j], lw=0)
        rows_b.append(["liver", cond, key, v])
ax[1].axhline(1, color=K["ink"], lw=0.6, ls=(0, (2, 2)))
ax[1].set_xticks(range(4))
ax[1].set_xticklabels(["units", "mass", "mass per\nunit", "perfused\nmass"], fontsize=5.5)
ax[1].set_ylabel("Day 220, relative to before injury"); ax[1].set_ylim(0, 2.25)
ax[1].set_title("Liver, by condition", fontsize=7, color=K["liver"])

# C: the kidney cannot make nephrons, so it enlarges the ones it kept
rows_c = []
for i, frac in enumerate((0.3, 0.5)):
    d = I[f"kidney_compensation_{frac}"]
    for j, (lab, v) in enumerate((("nephrons kept", d["units"]),
                                  ("mass per nephron", d["mass_per_unit"]),
                                  ("full compensation", d["full"]))):
        ax[2].bar(i + (j - 1) * 0.26, v, width=0.24,
                  color=[K["kidney"], K["warm"], K["grey"]][j], lw=0,
                  hatch="" if j < 2 else "///", edgecolor=K["grey"])
        rows_c.append(["kidney", frac, lab, v])
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

F.panel_grid(fig, [list(ax)], ["a", "b", "c"])
fig.savefig(str(OUT / "Fig8_injury.pdf")); fig.savefig(str(OUT / "Fig8_injury.png"), dpi=300)
plt.close(fig)
for name, hdr, rows in (("Fig8_trajectories", ["regime", "condition", "day", "relative_mass"], rows_a),
                        ("Fig8_liver_recovery", ["regime", "condition", "readout", "relative"], rows_b),
                        ("Fig8_kidney_compensation", ["regime", "fraction_lost", "readout", "relative"], rows_c)):
    with open(SD / f"{name}.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(hdr); w.writerows(rows)
print("wrote Fig8_injury")
