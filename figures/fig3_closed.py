"""Closed model: which individual spread reaches the final unit count."""
import json

import numpy as np
import matplotlib.pyplot as plt

from organgrow import figstyle as F
from organgrow.config import ROOT

F.apply()
RES = ROOT / "results"
OUTDIR = ROOT / "figures" / "out"
OUTDIR.mkdir(parents=True, exist_ok=True)
K = F.C

f = RES / "transmission.json"
if not f.exists():
    raise SystemExit("Fig3_closed skipped: run make closed first")
T = json.load(open(f))

REG = [("liver", "Liver\nsupply-bound", K["liver"]),
       ("kidney", "Kidney\npool-bound, loss", K["kidney"]),
       ("pancreas", "Pancreas\npool-bound, no loss", K["pancreas"])]
CON = [("baseline", "none"), ("pool", "pool"), ("setpoint", "setpoint")]

fig = plt.figure(figsize=(F.WIDTH_FULL, 68 * F.MM))
gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 1.0], wspace=0.32,
                      left=0.07, right=0.985, top=0.86, bottom=0.2)

# A: final unit count per replicate, relative to the regime's baseline mean
axa = fig.add_subplot(gs[0, 0])
x, ticks, labs = 0, [], []
rng = np.random.default_rng(0)
for reg, name, col in REG:
    base = np.mean(T[f"{reg}_baseline"]["N"])
    for cond, cname in CON:
        n = np.asarray(T[f"{reg}_{cond}"]["N"], float) / base
        jit = rng.uniform(-0.17, 0.17, len(n))
        axa.scatter(x + jit, n, s=9, color=col, alpha=0.55 if cond == "baseline" else 0.9,
                    lw=0)
        axa.plot([x - 0.25, x + 0.25], [n.mean()] * 2, color=K["ink"], lw=0.8)
        ticks.append(x); labs.append(cname)
        x += 1
    x += 0.6
axa.set_xticks(ticks)
axa.set_xticklabels(labs, rotation=90)
axa.set_ylabel("Final unit count\n(relative to no added spread)")
axa.axhline(1, color=K["grey"], lw=0.5, ls=(0, (2, 2)))
for j, (reg, name, col) in enumerate(REG):
    axa.text(j * 3.6 + 1, 1.04, name, transform=axa.get_xaxis_transform(),
             ha="center", va="bottom", fontsize=F.SZ["small"], color=col)
axa.set_xlabel("Individual spread added", labelpad=4)

# B: transmission of each spread, per regime
axb = fig.add_subplot(gs[0, 1])
w = 0.36
for j, (reg, name, col) in enumerate(REG):
    for i, (cond, cname) in enumerate(CON[1:]):
        v = T[f"{reg}_{cond}"]["transmission"]
        xx = j + (i - 0.5) * w
        axb.bar(xx, v, width=w * 0.92, color=col, alpha=0.45 if cond == "pool" else 1.0,
                lw=0)
        axb.text(xx, v + 0.03, f"{v:.2f}", ha="center", va="bottom",
                 fontsize=F.SZ["small"] - 0.5)
axb.axhline(1, color=K["grey"], lw=0.5, ls=(0, (2, 2)))
axb.set_xticks(range(3))
axb.set_xticklabels([r[1].split("\n")[0] for r in REG])
axb.set_ylabel("Transmission to unit count")
axb.set_ylim(0, 1.75)
from matplotlib.patches import Patch
axb.legend(handles=[Patch(color=K["grey"], alpha=0.45, lw=0, label="pool spread"),
                    Patch(color=K["grey"], lw=0, label="setpoint spread")],
           loc="upper right", fontsize=F.SZ["small"], frameon=False)

F.panel_grid(fig, [[axa, axb]], ["a", "b"])
fig.savefig(str(OUTDIR / "Fig3_closed.pdf"))
fig.savefig(str(OUTDIR / "Fig3_closed.png"), dpi=300)
print("wrote figures/out/Fig3_closed.pdf")
