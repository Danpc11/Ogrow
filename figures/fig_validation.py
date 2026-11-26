"""Figure 6 of the manuscript: three numbers of the model beside measurements.

Filled markers are the model, open markers are what was measured, and nothing
is fitted. Panel A is the comparison the model passes and panel B the one it
passes only in part, and they are drawn on the same axes so that the second
cannot be read more kindly than the first.
"""
import csv
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from organgrow import figstyle as F
from organgrow.config import ROOT

F.apply()
RES, OUT, SD = ROOT / "results", ROOT / "figures" / "out", ROOT / "source_data"
OUT.mkdir(parents=True, exist_ok=True)
SD.mkdir(parents=True, exist_ok=True)
K = F.C

V = json.load(open(RES / "validation.json")) if (RES / "validation.json").exists() else None
if V is None:
    print("no results/validation.json; run scripts/36_validation.py first")
else:
    fig, ax = plt.subplots(1, 3, figsize=(F.WIDTH_FULL, 60 * F.MM))
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.34, top=0.88, wspace=0.42)
    rows = []

    def model(a, reg, key, col, label):
        """The model's curve: what is left against the share of the pool seeded."""
        pts = sorted((float(k), v[key]) for k, v in V["ablation"][reg].items())
        x = np.array([p[0] for p in pts] + [1.0])
        y = np.array([p[1]["value"] for p in pts] + [1.0])
        lo = np.array([p[1]["lo"] for p in pts] + [1.0])
        hi = np.array([p[1]["hi"] for p in pts] + [1.0])
        a.plot(x, y, "-", lw=1.0, color=col, zorder=2)
        a.errorbar(x, y, yerr=[y - lo, hi - y], fmt="o", ms=3.2, mew=0, color=col,
                   elinewidth=0.8, capsize=0, zorder=3, label=label)
        rows.extend([reg, "model", key, xi, yi, li, hi_] for xi, yi, li, hi_
                    in zip(x, y, lo, hi))

    def measured(a, reg, col, marker, study, label, dx=0.0):
        """One study's measurements, open, so they cannot be taken for the model."""
        pts = [r for r in V["published"][reg] if r["study"] == study]
        x = np.array([r["kept"] + dx for r in pts])
        y = np.array([r["left"] for r in pts])
        err = None
        if all("lo" in r for r in pts):
            err = [y - np.array([r["lo"] for r in pts]),
                   np.array([r["hi"] for r in pts]) - y]
        a.errorbar(x, y, yerr=err, fmt=marker, ms=4.6, mfc="white", mec=col, mew=1.0,
                   ecolor=col, elinewidth=0.8, capsize=1.6, zorder=4, label=label)
        rows.extend([reg, study, r["at"], r["kept"], r["left"], r.get("lo", ""),
                     r.get("hi", "")] for r in pts)

    def frame(a, title):
        a.plot([0, 1], [0, 1], ls=(0, (2, 2)), lw=0.8, color=K["grey"], zorder=1)
        a.set_xlim(0.15, 1.05); a.set_ylim(0, 1.12)
        a.set_xlabel("Share of the progenitor pool left")
        a.set_title(title, fontsize=F.SZ["label"])

    # A: pancreas and liver, the contrast of one experiment
    model(ax[0], "pancreas", "mass", K["pancreas"], "pancreas-like, model")
    model(ax[0], "liver", "mass", K["liver"], "liver-like, model")
    measured(ax[0], "pancreas", K["pancreas"], "s", "Stanger 2007", "mouse pancreas")
    measured(ax[0], "liver", K["liver"], "s", "Stanger 2007", "mouse liver")
    frame(ax[0], "Pancreas and liver")
    ax[0].set_ylabel("Share of the organ that results")
    ax[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2,
                 columnspacing=1.0, handletextpad=0.4)

    # B: kidney, two experiments that removed the same share
    model(ax[1], "kidney", "N", K["kidney"], "kidney-like, model")
    measured(ax[1], "kidney", K["kidney"], "o", "Muthukrishnan 2018",
             "mouse, CITED1 driver", dx=-0.035)
    measured(ax[1], "kidney", K["ink"], "D", "Cebrian 2014",
             "mouse, Gdnf driver", dx=0.035)
    frame(ax[1], "Kidney")
    ax[1].set_ylabel("Share of the units that results")
    ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=1)

    # C: how far one organ sits from its pair, model and people
    seeds = np.asarray(V["seedwise"]["kidney"]["mass"], float)
    dm = 100 * (np.log(seeds) - np.log(seeds).mean())
    groups = [("constructs\ngrown alike", dm, K["kidney"], True)]
    if V.get("pairs"):
        dp = 100 * np.asarray(V["pairs"]["deviation"], float)
        groups.append(("a person's\ntwo kidneys", dp, K["grey"], False))
    rng = np.random.default_rng(3)
    for i, (lab, d, col, filled) in enumerate(groups):
        x = i + rng.uniform(-0.16, 0.16, len(d))
        ax[2].plot(x, d, "o", ms=2.8, mew=0.7, mec=col, mfc=col if filled else "white",
                   alpha=0.85, zorder=2)
        sd = d.std(ddof=1)
        ax[2].plot([i + 0.30, i + 0.30], [-sd, sd], "-", lw=1.4, color=K["ink"], zorder=3)
        ax[2].text(i + 0.34, sd, f"{sd:.1f} %", ha="left", va="center",
                   fontsize=F.SZ["annot"], color=K["ink"])
        rows.extend(["kidney", lab.replace("\n", " "), "deviation from its pair, %",
                     j, v, "", ""] for j, v in enumerate(d))
    ax[2].axhline(0, lw=0.6, color=K["grey"], zorder=1)
    ax[2].set_xticks(range(len(groups)))
    ax[2].set_xticklabels([g[0] for g in groups])
    ax[2].set_xlim(-0.5, len(groups) - 0.15)
    ax[2].set_ylabel("Deviation of one organ (%)")
    ax[2].set_title("Two organs grown alike", fontsize=F.SZ["label"])

    for a in ax:
        F.tidy(a)
    F.panel_grid(fig, [list(ax)], ["a", "b", "c"])
    F.audit(fig, "Fig10_validation")
    fig.savefig(str(OUT / "Fig10_validation.pdf"))
    fig.savefig(str(OUT / "Fig10_validation.png"), dpi=300)
    plt.close(fig)
    with open(SD / "Fig10_validation.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["regime", "series", "quantity", "x", "y", "lo", "hi"])
        wr.writerows(rows)
    print("wrote Fig10_validation")
