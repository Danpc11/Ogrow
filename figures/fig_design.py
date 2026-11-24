"""Figure 5 of the manuscript: what domain size costs.

The scaling sweep. Its point is not the exponent but the spread:
the per-seed counts are drawn, not only their means, because past an
intermediate domain the same rules and parameters stop giving the same number
of units, and a mean would hide that.
"""
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
OUT.mkdir(parents=True, exist_ok=True)
SD.mkdir(parents=True, exist_ok=True)
K = F.C
REG = ("liver", "kidney", "pancreas")
COND = (("none", "unchanged"), ("supply", "supply relieved"), ("room", "room relieved"))
SHADE = (0.45, 0.75, 1.0)        # one regime colour, three tints by condition


def tint(hexcol, f):
    """Lighten a colour towards white by 1 - f, so one hue carries a regime."""
    c = np.array([int(hexcol[i:i + 2], 16) for i in (1, 3, 5)], float) / 255
    return tuple(1 - f * (1 - c))


# ========================================================= Figure 5: scaling
Sc = json.load(open(RES / "scaling.json")) if (RES / "scaling.json").exists() else None
if Sc:
    MODES = [("geometry", "domain alone", K["sim"]),
             ("allometric", "with the allometries", K["liver"])]
    fig, ax = plt.subplots(1, 3, figsize=(F.WIDTH_FULL, 56 * F.MM))
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.30, top=0.88, wspace=0.40)
    rows6 = []

    # A: every seed, not only the mean, against domain volume
    hand = []
    for mode, lab, col in MODES:
        rows = Sc[mode]["rows"]
        V = np.array([r["S"] ** 3 for r in rows], float)
        N = np.array([r["N"] for r in rows], float)
        ax[0].plot(V, N, "o", ms=2.2, mew=0, color=col, alpha=0.5)
        b, c = np.polyfit(np.log(V), np.log(N), 1)
        vv = np.array([V.min(), V.max()])
        ax[0].plot(vv, np.exp(c) * vv ** b, "-", lw=1.1, color=col)
        hand.append(plt.Line2D([0], [0], color=col, lw=1.1,
                               label=f"{lab}, $V^{{{b:.2f}}}$"))
        rows6 += [[mode, r["S"], r["seed"], r["N"]] for r in rows]
    Vr = np.array([(0.4) ** 3, 1.0])
    n0 = np.mean([r["N"] for r in Sc["geometry"]["rows"] if r["S"] == 0.4])
    ax[0].plot(Vr, n0 * (Vr / Vr[0]) ** 1.0, ls=(0, (2, 2)), lw=0.8, color=K["grey"])
    hand.append(plt.Line2D([0], [0], color=K["grey"], lw=0.8, ls=(0, (2, 2)),
                           label="$V^{1}$, conserved density"))
    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    ax[0].set_xlabel("Domain volume (relative)"); ax[0].set_ylabel("Units at arrest")
    ax[0].set_xlim(0.045, 1.9)
    ax[0].set_title("One point per seed", fontsize=F.SZ["label"])
    ax[0].legend(handles=hand, loc="upper center", bbox_to_anchor=(0.5, -0.34),
                 ncol=1)

    S = sorted({r["S"] for r in Sc["geometry"]["rows"]})

    def per_scale(mode, fn):
        return [fn([r for r in Sc[mode]["rows"] if r["S"] == s]) for s in S]

    # B: which bound holds the pool, against size
    for key, lab, col, mk in (("supply_held", "held by perfusion", K["kidney"], "o"),
                              ("space_held", "held by room", K["pancreas"], "s")):
        y = per_scale("geometry", lambda g, k=key: float(np.mean([r[k] for r in g])))
        ax[1].plot(S, y, mk + "-", ms=3, lw=1.0, color=col, mew=0)
    ax[1].set_xlabel("Domain scale"); ax[1].set_ylabel("Share of the pool")
    ax[1].set_ylim(0, 1.0)
    ax[1].set_title("Which bound binds", fontsize=F.SZ["label"])
    ax[1].legend(handles=[Patch(facecolor=K["kidney"], label="held by perfusion"),
                          Patch(facecolor=K["pancreas"], label="held by room")],
                 loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1)

    # C: the spread of the unit count, which is the result
    for mode, lab, col in MODES:
        cv = per_scale(mode, lambda g: float(np.std([r["N"] for r in g], ddof=1)
                                             / max(np.mean([r["N"] for r in g]), 1e-9)))
        ax[2].plot(S, cv, "o-", ms=3, lw=1.0, color=col, mew=0)
    ax[2].set_xlabel("Domain scale")
    ax[2].set_ylabel("CV of the unit count")
    ax[2].set_ylim(0, None)
    ax[2].set_title("Reproducibility at size", fontsize=F.SZ["label"])
    ax[2].legend(handles=[Patch(facecolor=c, label=l) for _, l, c in MODES],
                 loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1)

    for a in ax:
        F.tidy(a)
    F.panel_grid(fig, [list(ax)], ["a", "b", "c"])
    F.audit(fig, "Fig6_scaling")
    fig.savefig(str(OUT / "Fig6_scaling.pdf"))
    fig.savefig(str(OUT / "Fig6_scaling.png"), dpi=300)
    plt.close(fig)
    with open(SD / "Fig6_scaling.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["sweep", "domain_scale", "seed", "units"])
        wr.writerows(rows6)
    print("wrote Fig6_scaling")
