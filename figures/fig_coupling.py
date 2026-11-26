"""Figure 7 of the manuscript: the extension in which progenitors make their tip branch.

Lines are the model, ablated over control; open symbols and grey bands are
what two mouse experiments measured. Model days are not embryonic days, so the
measurements are drawn at the right-hand edge and not along the time axis.
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

V = json.load(open(RES / "kidney_coupling.json")) if (RES / "kidney_coupling.json").exists() else None
if V is None:
    print("no results/kidney_coupling.json; run scripts/37_kidney_coupling.py first")
else:
    ta, T = V["ablate_day"], V["window"]
    names = list(V["nucleation"])
    low = names[-1]                       # where perfusion limits nucleation least
    m0 = V["measured"]["no response (Cebrian 2014)"]
    m1 = V["measured"]["response (Muthukrishnan 2018)"]
    t = np.arange(ta, T + 1) - ta
    sl = slice(ta - 1, T)
    fig, ax = plt.subplots(1, 3, figsize=(F.WIDTH_FULL, 60 * F.MM))
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.38, top=0.88, wspace=0.40)
    rows = []

    # A: what the ablation does to the pool, the cap and the sites
    c = V["conditions"][low]["0"]
    for key, lab, col, ls in (("pool_by_day", "progenitors", K["cells"], "-"),
                              ("cap_by_day", "per site", K["ink"], "-"),
                              ("sites_by_day", "branching sites", K["kidney"], "-")):
        y = np.asarray(c[key])[sl]
        ax[0].plot(t, y, ls, lw=1.1, color=col, label=lab)
        rows += [["a", lab, int(d), float(v)] for d, v in zip(t, y)]
    xr = t[-1] + 4
    ax[0].plot([xr] * 3, m0["cap"], "o", ms=4, mfc="white", mec=K["ink"], mew=1.0,
               clip_on=False, label="cells per tip, measured")
    ax[0].plot([xr + 5] * 3, m0["tips"], "D", ms=3.6, mfc="white", mec=K["kidney"],
               mew=1.0, clip_on=False, label="tips, measured")
    ax[0].axhline(1, lw=0.6, color=K["grey"], zorder=0)
    ax[0].set_xlim(0, xr + 9); ax[0].set_ylim(0.3, 1.5)
    ax[0].set_xlabel("Days since the ablation"); ax[0].set_ylabel("Ablated over control")
    ax[0].set_title("No proliferative response", fontsize=F.SZ["label"])
    ax[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2,
                 columnspacing=1.0, handletextpad=0.4)

    # B: units, with and without the response, where perfusion limits least
    for resp, lab, ls in (("0", "no response", "-"), ("1", "division doubled", (0, (3, 1.5)))):
        y = np.asarray(V["conditions"][low][resp]["units_by_day"])[sl]
        ax[1].plot(t, y, ls=ls, lw=1.1, color=K["kidney"], label=lab)
        rows += [["b", lab, int(d), float(v)] for d, v in zip(t, y)]
    ax[1].axhspan(min(m0["units"]), max(m0["units"]), color=K["pale"], lw=0, zorder=0)
    ax[1].axhspan(min(m1["units"]), max(m1["units"]), color=K["pale"], lw=0, zorder=0)
    ax[1].text(t[-1], np.mean(m0["units"]), "measured,\nno response", ha="right",
               va="center", fontsize=F.SZ["small"], color=K["grey"])
    ax[1].text(1.5, np.mean(m1["units"]), "measured,\nresponse", ha="left",
               va="center", fontsize=F.SZ["small"], color=K["grey"])
    ax[1].axhline(1, lw=0.6, color=K["grey"], zorder=0)
    ax[1].set_xlim(0, t[-1]); ax[1].set_ylim(0.3, 1.25)
    ax[1].set_xlabel("Days since the ablation"); ax[1].set_ylabel("Units (ratio)")
    ax[1].set_title("Units", fontsize=F.SZ["label"])
    ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=1)

    # C: the end of the window, every setting of the nucleation law
    for j, n in enumerate(names):
        for resp, dx, filled in (("0", -0.13, True), ("1", 0.13, False)):
            r = V["conditions"][n][resp]["units"]
            ax[2].errorbar([j + dx], [r["value"]],
                           yerr=[[r["value"] - r["lo"]], [r["hi"] - r["value"]]],
                           fmt="o", ms=3.6, mew=1.0, mec=K["kidney"],
                           mfc=K["kidney"] if filled else "white", ecolor=K["kidney"],
                           elinewidth=0.8, capsize=0)
            rows.append(["c", f"{n}, response {resp}", r["value"], f"{r['lo']}:{r['hi']}"])
    ax[2].axhspan(min(m0["units"]), max(m0["units"]), color=K["pale"], lw=0, zorder=0)
    ax[2].axhspan(min(m1["units"]), max(m1["units"]), color=K["pale"], lw=0, zorder=0)
    ax[2].axhline(1, lw=0.6, color=K["grey"], zorder=0)
    ax[2].set_xticks(range(len(names)))
    ax[2].set_xticklabels([f"{V['nucleation'][n]['q_set']:g}\n{V['nucleation'][n]['nu']:g}"
                           for n in names])
    ax[2].set_xlabel("Threshold (top), rate of the law")
    ax[2].set_xlim(-0.5, len(names) - 0.5); ax[2].set_ylim(0.3, 1.25)
    ax[2].set_ylabel("Units at the end (ratio)")
    ax[2].set_title("Where perfusion limits less", fontsize=F.SZ["label"])
    ax[2].plot([], [], "o", ms=3.6, mec=K["kidney"], mfc=K["kidney"], label="no response")
    ax[2].plot([], [], "o", ms=3.6, mec=K["kidney"], mfc="white", mew=1.0,
               label="division doubled")
    ax[2].legend(loc="upper center", bbox_to_anchor=(0.5, -0.40), ncol=2,
                 columnspacing=1.0, handletextpad=0.4)

    for a in ax:
        F.tidy(a)
    F.panel_grid(fig, [list(ax)], ["a", "b", "c"])
    F.audit(fig, "Fig11_coupling")
    fig.savefig(str(OUT / "Fig11_coupling.pdf"))
    fig.savefig(str(OUT / "Fig11_coupling.png"), dpi=300)
    plt.close(fig)
    with open(SD / "Fig11_coupling.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["panel", "series", "x", "y"])
        wr.writerows(rows)
    print("wrote Fig11_coupling")
