"""Figure: what a developing kidney does to its pool, against what the law says."""
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
K = dict(ink="#2b2b2b", grey="#9a9a9a", kidney="#2e6da4", warm="#b8860b",
         model="#c0392b", liver="#4a7c3f")

T = json.load(open(RES / "transcriptome.json"))
w = np.array(T["kidney"]["week"], float)
npc = np.array(T["kidney"]["npc"], float)
ind = np.array(T["kidney"]["induced"], float)
neph = np.array(T["kidney"]["nephron"], float)
lin = np.array(T["kidney"]["lineage"], float)
share, per = npc / lin, ind / npc

fig = plt.figure(figsize=(F.WIDTH_FULL, 92 * F.MM))
gs = fig.add_gridspec(2, 6, left=0.075, right=0.985, bottom=0.09, top=0.93,
                      wspace=1.45, hspace=0.62)
ax = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]),
      fig.add_subplot(gs[0, 4:6]), fig.add_subplot(gs[1, 0:3]),
      fig.add_subplot(gs[1, 3:6])]

# A: the lineage, week by week
bot = np.zeros(len(w))
for v, lab, col in ((npc, "progenitor", K["kidney"]),
                    (ind, "being induced", K["warm"]),
                    (neph, "unit formed", K["grey"])):
    ax[0].bar(np.arange(len(w)), v / lin, bottom=bot, width=0.62, color=col, lw=0)
    bot += v / lin
for i, n in enumerate(lin):
    ax[0].text(i, 1.02, f"{int(n)}", ha="center", fontsize=5, color=K["ink"])
ax[0].set_xticks(range(len(w))); ax[0].set_xticklabels([f"wk {int(x)}" for x in w], fontsize=6)
ax[0].set_ylabel("Share of the nephron lineage"); ax[0].set_ylim(0, 1.1)
ax[0].set_title("Human fetal kidney", fontsize=7, color=K["kidney"])
ax[0].legend(handles=[Patch(facecolor=c, label=l) for l, c in
                      (("progenitor", K["kidney"]), ("being induced", K["warm"]),
                       ("unit formed", K["grey"]))],
             frameon=False, fontsize=5.5, loc="upper center", handlelength=1.1,
             handleheight=0.8, borderpad=0.1, ncol=3, columnspacing=0.8,
             bbox_to_anchor=(0.5, -0.17))

# B: one sample per stage, and neither ratio orders with stage
for y, col in ((share, K["kidney"]), (per, K["warm"])):
    ax[1].plot(w, y, "o", ms=3.4, color=col, mew=0)
k2 = T.get("kidney_second_study")
if k2:
    s2 = np.array(k2["npc"], float) / np.maximum(np.array(k2["lineage"], float), 1)
    ax[1].plot(np.array(k2["week"], float), s2, "s", ms=3.4, color=K["kidney"],
               mfc="none", mew=0.8)
    ax[1].annotate("another study", (float(k2["week"][0]), float(s2[0])),
                   textcoords="offset points", xytext=(4, 4), fontsize=5,
                   color=K["kidney"])
sp = T["spearman"]
ax[1].set_yscale("log"); ax[1].set_xlabel("Gestational week")
ax[1].set_ylabel("Ratio"); ax[1].set_xlim(7.5, 20.5)
ax[1].set_title(f"$\\rho$ = {sp['npc_share'][0]:+.2f}, P = {sp['npc_share'][1]:.2f}"
                f"   |   $\\rho$ = {sp['induced_per_npc'][0]:+.2f}, "
                f"P = {sp['induced_per_npc'][1]:.2f}", fontsize=6, color=K["ink"])
ax[1].legend(handles=[Patch(facecolor=c, label=l) for l, c in
                      (("progenitor share", K["kidney"]),
                       ("induced per progenitor", K["warm"]))],
             frameon=False, fontsize=5.5, loc="lower left", handlelength=1.1,
             handleheight=0.8, borderpad=0.1)

# C: where the model's pool goes, which is not what "spent" suggests
# Only the budget organs. The liver renews, so its pool is replaced many times
# over and a share of what it started with is not a meaningful denominator.
regs = ("kidney", "pancreas")
rows_c = []
_rp = RES / "paper" / "regime_reps.json"
RP = json.load(open(_rp)) if _rp.exists() else {}
if not RP:
    print("fig_transcriptome: no regime_reps.json, panel c falls back to one run "
          "per regime; run `python scripts/17_paper_runs.py reps` for the means")


def fate(r):
    """Mean over the replicates, with the one-run file as the fallback."""
    v = RP.get(r)
    if v is None:
        fa = np.load(RES / "paper" / f"regime_{r}.npz")["fate"]
        return float(fa[0]), float(fa[1]), float(fa[3])
    f = v["fate"]
    return f["consumed"]["mean"], f["lost"]["mean"], f["left"]["mean"]


for i, r in enumerate(regs):
    cons, lost, leftv = fate(r)
    tot = max(cons + lost + leftv, 1e-9)
    bot = 0.0
    for v, col in ((cons, K["kidney"]), (lost, K["model"]), (leftv, K["grey"])):
        ax[2].bar(i, v / tot, bottom=bot, width=0.6, color=col, lw=0)
        bot += v / tot
    rows_c.append([r, cons, lost, leftv])
ax[2].set_xticks(range(len(regs)))
ax[2].set_xticklabels(regs, fontsize=6)
ax[2].set_ylabel("Share of the initial pool"); ax[2].set_ylim(0, 1.0)
_n = RP.get("kidney", {}).get("seeds")
ax[2].set_title(f"Where the pool goes ({_n} replicates)" if _n
                else "Where the pool goes, in the model", fontsize=7, color=K["ink"])
rows_c.append(["liver", *fate("liver")])
ax[2].legend(handles=[Patch(facecolor=K["kidney"], label="consumed into units"),
                      Patch(facecolor=K["model"], label="lost to attrition"),
                      Patch(facecolor=K["grey"], label="still held")],
             frameon=False, fontsize=5.5, loc="upper center", handlelength=1.1,
             handleheight=0.8, borderpad=0.1, ncol=3, columnspacing=0.8,
             bbox_to_anchor=(0.5, -0.17))

# D: the two mouse organs, progenitor share against stage
pan, lv = T["pancreas"], T["liver"]
for d, key, lo_k, hi_k, lab, col in (
        (pan, "frac_strict", "s_lo", "s_hi", "pancreas, multipotent", K["warm"]),
        (lv, "frac_prog", "lo", "hi", "liver, hepatoblast", K["liver"])):
    x = np.array([float(v[1:]) for v in d["stage"]], float)
    y = np.array(d[key], float)
    e = np.vstack([y - np.array(d[lo_k], float), np.array(d[hi_k], float) - y])
    ax[3].errorbar(x, np.maximum(y, 3e-5), yerr=np.abs(e), fmt="o-", ms=2.6,
                   lw=0.9, color=col, elinewidth=0.7, capsize=1.2, mew=0)
ax[3].set_yscale("log"); ax[3].set_ylim(2e-5, 2)
ax[3].set_xlabel("Embryonic day"); ax[3].set_ylabel("Progenitor share\nof the lineage")
ax[3].set_title("Mouse pancreas and liver", fontsize=7, color=K["ink"])
ax[3].legend(handles=[Patch(facecolor=c, label=l) for l, c in
                      (("pancreas, multipotent", K["warm"]),
                       ("liver, hepatoblast", K["liver"]))],
             frameon=False, fontsize=5.5, loc="lower left", handlelength=1.1,
             handleheight=0.8, borderpad=0.1)

# E: when each pool goes, against how much of the organ's window has passed
rows_e = []
for name, stages, frac, col in (
        ("kidney", w, npc / lin, K["kidney"]),
        ("pancreas", [float(v[1:]) for v in pan["stage"]], np.array(pan["frac_strict"], float),
         K["warm"]),
        ("liver", [float(v[1:]) for v in lv["stage"]], np.array(lv["frac_prog"], float),
         K["liver"])):
    win = T["windows"][name]
    x = (np.array(stages, float) - win["start"]) / (win["end"] - win["start"])
    y = np.array(frac, float) / max(float(np.array(frac, float)[0]), 1e-12)
    ax[4].plot(x, np.maximum(y, 3e-5), "o-", ms=2.6, lw=0.9, color=col, mew=0)
    rows_e += [[name, float(a), float(b)] for a, b in zip(x, y)]
ax[4].axvline(1.0, color=K["ink"], lw=0.6, ls=(0, (2, 2)))
ax[4].text(1.02, 2e-4, "window ends", fontsize=5.5, color=K["ink"], rotation=90,
           va="bottom")
ax[4].set_yscale("log"); ax[4].set_ylim(2e-5, 2); ax[4].set_xlim(-0.05, 1.25)
ax[4].set_xlabel("Fraction of the unit-forming window elapsed")
ax[4].set_ylabel("Progenitor share,\nrelative to the first stage")
ax[4].set_title("When each pool goes", fontsize=7, color=K["ink"])
ax[4].legend(handles=[Patch(facecolor=c, label=l) for l, c in
                      (("kidney", K["kidney"]), ("pancreas", K["warm"]),
                       ("liver", K["liver"]))],
             frameon=False, fontsize=5.5, loc="lower left", handlelength=1.1,
             handleheight=0.8, borderpad=0.1, ncol=3, columnspacing=0.8)

F.panel_grid(fig, [ax[:3], ax[3:]], ["a", "b", "c", "d", "e"])
fig.savefig(str(OUT / "Fig9_transcriptome.pdf"))
fig.savefig(str(OUT / "Fig9_transcriptome.png"), dpi=300)
plt.close(fig)
with open(SD / "Fig9_windows.csv", "w", newline="") as f:
    wr = csv.writer(f)
    wr.writerow(["organ", "window_fraction_elapsed", "progenitor_share_relative"])
    wr.writerows(rows_e)
with open(SD / "Fig9_pool_fate.csv", "w", newline="") as f:
    wr = csv.writer(f)
    wr.writerow(["regime", "consumed", "lost", "still_held"])
    wr.writerows(rows_c)
with open(SD / "Fig9_kidney_trends.csv", "w", newline="") as f:
    wr = csv.writer(f)
    wr.writerow(["week", "lineage", "progenitor", "induced", "nephron",
                 "progenitor_share", "induced_per_progenitor"])
    wr.writerows(zip(w.astype(int), lin.astype(int), npc.astype(int),
                     ind.astype(int), neph.astype(int), share, per))
print("wrote Fig9_transcriptome")
