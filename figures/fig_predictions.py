"""Figure 6 of the second manuscript: predictions for ablations not yet done.

Filled symbols joined by a line are the model at the doses already compared
with a measurement, open squares are those measurements, and stars are the
model at doses and times nobody has reported.
"""
import csv
import json

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

need = [RES / "validation.json", RES / "predictions.json"]
if not all(p.exists() for p in need):
    print("needs results/validation.json and results/predictions.json; run "
          "scripts/36_validation.py and scripts/38_predictions.py first")
else:
    V, Pr = (json.load(open(p)) for p in need)
    fig, ax = plt.subplots(1, 3, figsize=(F.WIDTH_FULL, 60 * F.MM))
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.36, top=0.88, wspace=0.42)
    rows = []

    def bar(a, x, r, col, marker, ms, filled=True, label=None):
        a.errorbar([x], [r["value"]], yerr=[[r["value"] - r["lo"]], [r["hi"] - r["value"]]],
                   fmt=marker, ms=ms, mew=1.0, mec=col, mfc=col if filled else "white",
                   ecolor=col, elinewidth=0.8, capsize=0, zorder=4, label=label)

    # A: dose, pancreas-like and liver-like
    for reg, col in (("pancreas", K["pancreas"]), ("liver", K["liver"])):
        known = sorted((float(k), v["mass"]) for k, v in V["ablation"][reg].items())
        new = sorted((float(k), v) for k, v in Pr["dose"][reg].items())
        both = sorted(known + new) + [(1.0, dict(value=1.0, lo=1.0, hi=1.0))]
        ax[0].plot([p[0] for p in both], [p[1]["value"] for p in both], "-", lw=0.9,
                   color=col, zorder=2)
        for i, (x, r) in enumerate(known):
            bar(ax[0], x, r, col, "o", 3.0, label=f"{reg}-like, compared" if i == 0 else None)
            rows.append(["a", reg, "model, compared", x, r["value"], r["lo"], r["hi"]])
        for i, (x, r) in enumerate(new):
            bar(ax[0], x, r, col, "*", 7.5, label=f"{reg}-like, predicted" if i == 0 else None)
            rows.append(["a", reg, "model, predicted", x, r["value"], r["lo"], r["hi"]])
        for r in V["published"][reg]:
            ax[0].plot([r["kept"]], [r["left"]], "s", ms=4.4, mfc="white", mec=col,
                       mew=1.0, zorder=3)
            rows.append(["a", reg, "measured", r["kept"], r["left"], r.get("lo", ""), r.get("hi", "")])
    ax[0].plot([0, 1], [0, 1], ls=(0, (2, 2)), lw=0.8, color=K["grey"], zorder=1)
    ax[0].plot([], [], "s", ms=4.4, mfc="white", mec=K["ink"], mew=1.0, label="measured")
    ax[0].set_xlim(0.08, 1.05); ax[0].set_ylim(0, 1.12)
    ax[0].set_xlabel("Share of the progenitor pool left")
    ax[0].set_ylabel("Share of the organ that results")
    ax[0].set_title("Dose", fontsize=F.SZ["label"])
    ax[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2,
                 columnspacing=0.8, handletextpad=0.3)

    # B: depth of a kidney ablation, with and without a proliferative response
    for resp, filled, lab in (("0", True, "no response"), ("1", False, "division doubled")):
        pts = sorted((float(k.split("|")[1]), v) for k, v in Pr["kidney"][resp].items()
                     if k.startswith("45|"))
        ax[1].plot([1 - p[0] for p in pts[::-1]] + [1.0],
                   [p[1]["value"] for p in pts[::-1]] + [1.0],
                   ls="-" if filled else (0, (3, 1.5)), lw=0.9, color=K["kidney"], zorder=2)
        for share, r in pts:
            done = abs(share - 0.4) < 1e-9
            bar(ax[1], 1 - share, r, K["kidney"], "o" if done else "*", 3.2 if done else 7.5,
                filled=filled)
            rows.append(["b", f"response {resp}", "compared" if done else "predicted",
                         1 - share, r["value"], r["lo"], r["hi"]])
    for r in V["published"]["kidney"]:
        ax[1].plot([r["kept"]], [r["left"]], "s", ms=4.4, mfc="white", mec=K["ink"],
                   mew=1.0, zorder=3)
    ax[1].plot([0, 1], [0, 1], ls=(0, (2, 2)), lw=0.8, color=K["grey"], zorder=1)
    ax[1].plot([], [], "-", color=K["kidney"], lw=0.9, label="no response")
    ax[1].plot([], [], ls=(0, (3, 1.5)), color=K["kidney"], lw=0.9, label="division doubled")
    ax[1].plot([], [], "*", ms=7.5, mfc=K["kidney"], mec=K["kidney"], label="predicted")
    ax[1].plot([], [], "s", ms=4.4, mfc="white", mec=K["ink"], mew=1.0, label="measured")
    ax[1].set_xlim(0.3, 1.05); ax[1].set_ylim(0, 1.12)
    ax[1].set_xlabel("Share of the progenitors left")
    ax[1].set_ylabel("Share of the units that results")
    ax[1].set_title("Depth of a kidney ablation", fontsize=F.SZ["label"])
    ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2,
                 columnspacing=0.8, handletextpad=0.3)

    # C: the same ablation early and late
    days = sorted({int(k.split("|")[0]) for k in Pr["kidney"]["0"] if k.endswith("|0.4")})
    for resp, filled, dx in (("0", True, -0.12), ("1", False, 0.12)):
        for j, d in enumerate(days):
            r = Pr["kidney"][resp][f"{d}|0.4"]
            done = d == 45
            bar(ax[2], j + dx, r, K["kidney"], "o" if done else "*", 3.2 if done else 7.5,
                filled=filled)
            rows.append(["c", f"response {resp}", f"day {d}", d, r["value"], r["lo"], r["hi"]])
    ax[2].axhline(1, lw=0.6, color=K["grey"], zorder=0)
    ax[2].set_xticks(range(len(days)))
    ax[2].set_xticklabels([f"day {d}" for d in days])
    ax[2].set_xlim(-0.5, len(days) - 0.5); ax[2].set_ylim(0, 1.12)
    ax[2].set_xlabel("Day of the ablation (window of 120)")
    ax[2].set_ylabel("Share of the units that results")
    ax[2].set_title("Timing of a 40 % ablation", fontsize=F.SZ["label"])

    for a in ax:
        F.tidy(a)
    F.panel_grid(fig, [list(ax)], ["a", "b", "c"])
    F.audit(fig, "Fig12_predictions")
    fig.savefig(str(OUT / "Fig12_predictions.pdf"))
    fig.savefig(str(OUT / "Fig12_predictions.png"), dpi=300)
    plt.close(fig)
    with open(SD / "Fig12_predictions.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["panel", "series", "kind", "x", "y", "lo", "hi"])
        wr.writerows(rows)
    print("wrote Fig12_predictions")
