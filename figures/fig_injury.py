"""Figure: injury after closure, and whether function follows structure."""
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
CLOSE, SEEDS = 140, (300, 301, 302)
LAB = dict(N="units", mass="mass", perfused="perfused mass", Q="root flow")


def load(reg, kind):
    return [json.load(open(RES / "injury" / f"{reg}_{kind}_s{s}.json")) for s in SEEDS]


fig, ax = plt.subplots(1, 3, figsize=(F.WIDTH_FULL, 52 * F.MM),
                       gridspec_kw=dict(width_ratios=[1.1, 1.0, 1.0]))
fig.subplots_adjust(left=0.07, right=0.985, bottom=0.24, top=0.86, wspace=0.42)

# A: the liver's two trajectories through a parenchymal injury
rows_a = []
for key, style, lab in (("N", "-", "units"), ("perfused", "--", "perfused mass")):
    for kind, col in (("control", K["grey"]), ("parenchyma", K["liver"])):
        y = np.mean([r["series"][key] for r in load("liver", kind)], axis=0)
        y = y / y[CLOSE - 2]
        t = np.arange(1, len(y) + 1)
        ax[0].plot(t, y, style, lw=1.1, color=col)
        rows_a += [[kind, key, int(tt), float(vv)] for tt, vv in zip(t, y)]
ax[0].axvline(CLOSE, color=K["ink"], lw=0.6, ls=(0, (2, 2)))
ax[0].text(CLOSE + 3, 0.49, "injury", fontsize=5.5, color=K["ink"])
ax[0].set_xlabel("Day"); ax[0].set_ylabel("Relative to the day before injury")
ax[0].set_xlim(100, 220); ax[0].set_ylim(0.45, 1.15)
ax[0].set_title("Liver, 30 % of units lost", fontsize=7, color=K["liver"])
for lab, st, c in (("units", "-", K["liver"]), ("perfused mass", "--", K["liver"]),
                   ("uninjured", "-", K["grey"])):
    ax[0].plot([], [], st, color=c, lw=1.1, label=lab)
ax[0].legend(frameon=False, fontsize=5.5, loc="lower left", handlelength=1.6)

# B: what is recovered, liver, by quantity and injury
I = json.load(open(RES / "injury.json"))
kinds = ("parenchyma", "vascular", "combined")
rows_b = []
for i, key in enumerate(("N", "mass", "perfused", "Q")):
    for j, kind in enumerate(kinds):
        v = I[f"liver_{kind}"][key]
        x = i + (j - 1) * 0.26
        val = v["of_control"]
        ax[1].bar(x, val, width=0.24, color=[K["liver"], K["warm"], K["ink"]][j], lw=0)
        rows_b.append([kind, key, v["recovered"], v["se"], val])
ax[1].axhline(1, color=K["grey"], lw=0.6, ls=(0, (2, 2)))
ax[1].set_xticks(range(4)); ax[1].set_xticklabels([LAB[k] for k in ("N", "mass", "perfused", "Q")],
                                                  fontsize=5.5, rotation=20, ha="right")
ax[1].set_ylabel("Day 220, relative to uninjured"); ax[1].set_ylim(0, 1.48)
ax[1].set_title("Liver", fontsize=7, color=K["liver"])
ax[1].legend(handles=[Patch(facecolor=c, label=l) for l, c in
                      zip(kinds, [K["liver"], K["warm"], K["ink"]])],
             frameon=False, fontsize=5.5, loc="upper left", handlelength=1.1,
             handleheight=0.8, borderpad=0.1)

# C: liver against kidney, same injury
rows_c = []
for j, reg in enumerate(("liver", "kidney")):
    for i, key in enumerate(("N", "mass", "perfused", "Q")):
        v = I[f"{reg}_parenchyma"][key]["of_control"]
        ax[2].bar(i + (j - 0.5) * 0.34, v, width=0.3, color=K[reg], lw=0)
        rows_c.append([reg, key, v])
ax[2].axhline(1, color=K["grey"], lw=0.6, ls=(0, (2, 2)))
ax[2].axhline(0.7, color=K["grey"], lw=0.6, ls=(0, (1, 2)))
ax[2].text(-0.42, 0.645, "no recovery", fontsize=5, color=K["grey"], ha="left")
ax[2].set_xticks(range(4)); ax[2].set_xticklabels([LAB[k] for k in ("N", "mass", "perfused", "Q")],
                                                  fontsize=5.5, rotation=20, ha="right")
ax[2].set_ylim(0, 1.48); ax[2].set_ylabel("Day 220, relative to uninjured")
ax[2].set_title("Same injury, two organs", fontsize=7)
ax[2].legend(handles=[Patch(facecolor=K[r], label=r) for r in ("liver", "kidney")],
             frameon=False, fontsize=5.5, loc="upper left", handlelength=1.1,
             handleheight=0.8, borderpad=0.1)

F.panel_grid(fig, [list(ax)], ["a", "b", "c"])
fig.savefig(str(OUT / "Fig8_injury.pdf")); fig.savefig(str(OUT / "Fig8_injury.png"), dpi=300)
plt.close(fig)
for name, hdr, rows in (("Fig8_trajectories", ["condition", "quantity", "day", "relative"], rows_a),
                        ("Fig8_liver_recovery", ["injury", "quantity", "recovered", "se", "of_control"], rows_b),
                        ("Fig8_two_organs", ["regime", "quantity", "of_control"], rows_c)):
    with open(SD / f"{name}.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(hdr); w.writerows(rows)
print("wrote Fig8_injury")
