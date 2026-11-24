"""Figures 1 and 2 of the manuscript, each a schematic and its data together.

The manuscript should be followable from the figures alone, so the picture of
the mechanism sits in the same figure as the measurement it explains rather
than in one of its own. Figure 1 states the law and then shows a construct
grown under it and where each pool ends up. Figure 2 asks what the binding
resource is worth: what seeding buys, and what relieving each bound releases.

Colour in the schematic panels marks the three resources and not the three
regimes: supply purple, competent cells teal, units and the room between them
ink. The data panels use the regime colours of the rest of the paper.
"""
import csv
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyBboxPatch, Patch, Wedge

from organgrow import figstyle as F
from organgrow.config import ROOT

F.apply()
RES, OUT, SD = ROOT / "results", ROOT / "figures" / "out", ROOT / "source_data"
OUT.mkdir(parents=True, exist_ok=True)
SD.mkdir(parents=True, exist_ok=True)
K = F.C
SUP, CEL, ROO = K["supply"], K["cells"], K["room"]      # the three resources
SM, TN = F.SZ["small"], 5.4
REG = (("liver", "renewing", K["liver"]),
       ("kidney", "spent,\nwith loss", K["kidney"]),
       ("pancreas", "spent,\nno loss", K["pancreas"]))


def aspect(fig, ax):
    b = ax.get_position()
    return (b.height * fig.get_figheight()) / (b.width * fig.get_figwidth())


def circ(fig, ax, c, r, **kw):
    ax.add_patch(Ellipse(c, 2 * r * aspect(fig, ax), 2 * r, **kw))


def blank(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")


def load(name):
    p = RES / name
    return json.load(open(p)) if p.exists() else None


# ==================================================== Figure 1: the model
snap = RES / "paper" / "snapshot_liver.npz"
reps = RES / "paper" / "regime_reps.json"
if snap.exists() and reps.exists():
    RP = json.load(open(reps))
    fig = plt.figure(figsize=(F.WIDTH_FULL, 112 * F.MM))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.86], width_ratios=[1.0, 1.0, 0.95],
                          left=0.015, right=0.985, top=0.945, bottom=0.045,
                          wspace=0.34, hspace=0.26)
    axA, axB, axC = (fig.add_subplot(gs[0, i]) for i in range(3))
    axD, axE, axF = (fig.add_subplot(gs[1, i]) for i in range(3))
    for a in (axA, axB, axC, axF):
        blank(a)
    fig.canvas.draw()

    # ------------------------------------------------------------ A: the law
    axA.text(0.0, 1.0, "A tip holds a supply and a pool",
             fontsize=F.SZ["label"], va="top", ha="left")
    for p0, p1, lw, col in (((0.02, 0.58), (0.20, 0.58), 2.4, ROO),
                            ((0.20, 0.58), (0.34, 0.78), 1.7, ROO),
                            ((0.20, 0.58), (0.34, 0.38), 1.7, ROO),
                            ((0.34, 0.38), (0.48, 0.31), 1.1, ROO),
                            ((0.34, 0.78), (0.48, 0.84), 1.1, ROO),
                            ((0.34, 0.78), (0.54, 0.60), 1.4, SUP)):
        axA.plot([p0[0], p1[0]], [p0[1], p1[1]], lw=lw, color=col,
                 solid_capstyle="round", zorder=1)
    tip = (0.46, 0.60)
    circ(fig, axA, tip, 0.235, fc=K["pale"], ec=K["grey"], lw=0.6,
         ls=(0, (2.5, 2)), alpha=0.5, zorder=0)
    axA.text(tip[0] + 0.06, tip[1] - 0.285, "the tissue this tip serves",
             fontsize=TN, color=K["grey"], ha="center", va="top")
    circ(fig, axA, tip, 0.028, fc=SUP, ec="none", zorder=4)
    axA.annotate("", xy=(0.425, 0.655), xytext=(0.325, 0.760),
                 arrowprops=dict(arrowstyle="-|>", lw=1.3, color=SUP,
                                 mutation_scale=7))
    axA.text(0.015, 0.90, "$q_i$ supply per unit", fontsize=TN, color=SUP,
             ha="left", va="center")
    rng = np.random.default_rng(3)
    ka = aspect(fig, axA)
    for _ in range(9):
        rr, th = 0.10 * np.sqrt(rng.random()), rng.random() * 2 * np.pi
        circ(fig, axA, (tip[0] + rr * np.cos(th) * ka, tip[1] + rr * np.sin(th)),
             0.016, fc=CEL, ec="none", alpha=0.9, zorder=5)
    axA.annotate("", xy=(tip[0] - 0.07, tip[1] - 0.03), xytext=(0.18, 0.33),
                 arrowprops=dict(arrowstyle="-", lw=0.5, color=CEL))
    axA.text(0.17, 0.325, "$P_i$ competent\ncells", fontsize=TN, color=CEL,
             ha="right", va="center", linespacing=1.3)
    u1, u2 = (tip[0] + 0.125, tip[1] + 0.125), (tip[0] + 0.105, tip[1] - 0.135)
    for u in (u1, u2):
        circ(fig, axA, u, 0.044, fc="white", ec=ROO, lw=1.0, zorder=5)
    axA.annotate("", xy=u2, xytext=u1,
                 arrowprops=dict(arrowstyle="<|-|>", lw=0.7, color=ROO,
                                 shrinkA=5, shrinkB=5, mutation_scale=6))
    axA.text(u1[0] + 0.085, (u1[1] + u2[1]) / 2, "$\\lambda_{\\min}$",
             fontsize=SM, color=ROO, ha="left", va="center")
    axA.text(0.99, 0.80, "units built,\nat least $\\lambda_{\\min}$ apart",
             fontsize=TN, color=ROO, ha="right", va="top", linespacing=1.35)
    axA.add_patch(FancyBboxPatch((0.02, 0.015), 0.70, 0.095,
                                 boxstyle="round,pad=0.012", fc="white",
                                 ec=ROO, lw=0.7, zorder=6))
    axA.text(0.37, 0.063, r"$\lambda_i = \nu\,P_i\,(q_i - q_{\mathrm{set}})_+$",
             fontsize=F.SZ["label"], ha="center", va="center", zorder=7)

    # ------------------------------------------------- B: what can run out
    axB.text(0.0, 1.0, "Three things can run out, each local",
             fontsize=F.SZ["label"], va="top", ha="left")
    ROWS = [(0.760, SUP, "Perfusion",
             "$q_i$ falls to $q_{\\mathrm{set}}$. Cells are there,\nthe supply to use them is not."),
            (0.455, CEL, "Competent cells",
             "$P_i = 0$. The supply is there,\nnothing is left to build with."),
            (0.150, ROO, "Room",
             "no site lies $\\lambda_{\\min}$ from its\nneighbours. Both are there.")]
    for y, col, head, body in ROWS:
        axB.add_patch(FancyBboxPatch((0.0, y - 0.125), 0.99, 0.255,
                                     boxstyle="round,pad=0.008", fc="white",
                                     ec=col, lw=0.8))
        axB.text(0.05, y + 0.058, head, fontsize=SM, color=col, va="center",
                 ha="left", fontweight="bold")
        axB.text(0.05, y - 0.045, body, fontsize=TN, color=K["ink"],
                 va="center", ha="left", linespacing=1.45)
    kb = aspect(fig, axB)
    g = (0.855, 0.760)
    axB.add_patch(Wedge(g, 0.072, 180, 360, width=0.024, fc=K["pale"], ec="none"))
    axB.add_patch(Wedge(g, 0.072, 180, 197, width=0.024, fc=SUP, ec="none"))
    axB.plot([g[0] - 0.075, g[0] - 0.047], [g[1], g[1]], lw=1.0, color=ROO)
    axB.text(g[0], g[1] - 0.055, "$q_{\\mathrm{set}}$", fontsize=TN,
             ha="center", va="top")
    circ(fig, axB, (0.855, 0.455), 0.075, fc="white", ec=CEL, lw=0.9,
         ls=(0, (2, 1.6)))
    axB.text(0.855, 0.455, "0", fontsize=SM, color=CEL, ha="center", va="center")
    for dx in (-0.042 * kb, 0.042 * kb):
        circ(fig, axB, (0.855 + dx, 0.165), 0.042, fc="white", ec=ROO, lw=0.9)
    axB.plot([0.855, 0.855], [0.113, 0.217], lw=0.7, color=ROO, ls=(0, (1.5, 1.2)))

    # ------------------------------------------- C: the single property
    axC.text(0.0, 1.0, "One property is varied", fontsize=F.SZ["label"],
             va="top", ha="left")
    axC.text(0.0, 0.875, "whether the cell that forms a unit is\nspent in forming it",
             fontsize=TN, color=K["grey"], va="top", ha="left", linespacing=1.4)
    kc = aspect(fig, axC)
    for x0, lab, spent in ((0.0, "spent   $w = 0$", True),
                           (0.52, "renewing   $w = 1$", False)):
        axC.text(x0 + 0.23, 0.675, lab, fontsize=SM, color=ROO, ha="center")
        for i in range(4):
            circ(fig, axC, (x0 + 0.085 + 0.062 * kc * i, 0.545), 0.026,
                 fc=CEL, ec="none")
        axC.annotate("", xy=(x0 + 0.115, 0.325), xytext=(x0 + 0.115, 0.475),
                     arrowprops=dict(arrowstyle="-|>", lw=0.8, color=ROO,
                                     mutation_scale=7))
        axC.text(x0 + 0.155, 0.40, "a cell\nforms a unit", fontsize=TN,
                 color=K["grey"], ha="left", va="center", linespacing=1.35)
        for i in range(3 if spent else 4):
            solid = spent or i < 3
            circ(fig, axC, (x0 + 0.085 + 0.062 * kc * i, 0.225), 0.026,
                 fc=CEL if solid else "white", ec="none" if solid else CEL,
                 lw=0.9, ls="-" if solid else (0, (1.6, 1.2)))
        circ(fig, axC, (x0 + 0.085 + 0.062 * kc * 4.6, 0.225), 0.036,
             fc="white", ec=ROO, lw=1.0)
        axC.text(x0 + 0.23, 0.065,
                 "the pool is drawn down" if spent else "the pool is restored",
                 fontsize=TN, color=ROO, ha="center", va="center")

    # --------------------------------------- D: a construct grown under it
    d = np.load(snap)
    scale = 7.4
    st, en, r = d["starts"] * scale, d["ends"] * scale, d["radius"]
    wv = 0.25 + 2.2 * (r - r.min()) / max(np.ptp(r), 1e-9)
    for a, b, lw in zip(st, en, wv):
        axD.plot([a[0], b[0]], [a[1], b[1]], color=ROO, lw=lw, alpha=0.55,
                 solid_capstyle="round")
    U = d["U"]
    axD.scatter(U[:, 0], U[:, 1], s=2.5, color=K["liver"], alpha=0.35, lw=0)
    tip_xy = d["ends"][d["tips"]] * scale
    strand = (d["q_tip"] <= d["q_set"]) & (d["P_tip"] > 0)
    axD.scatter(tip_xy[strand, 0], tip_xy[strand, 1], s=7, color=SUP, lw=0)
    axD.set_aspect("equal"); axD.axis("off")
    axD.set_title("A construct at arrest", fontsize=F.SZ["label"])
    axD.legend(handles=[Patch(facecolor=K["liver"], label="units"),
                        Patch(facecolor=SUP, label="tips holding cells below $q_{\\mathrm{set}}$")],
               loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=1)

    # ------------------------------------------------ E: where each pool is
    PARTS = (("spent", "spent", ROO), ("supply_held", "held by perfusion", SUP),
             ("space_held", "held by room", K["grey"]), ("active", "still active", K["pale"]))
    rows_e = []
    for i, (reg, lab, col) in enumerate(REG):
        bot = 0.0
        for key, _, c in PARTS:
            v = RP[reg]["closure"][key]["mean"]
            axE.bar(i, v, 0.62, bottom=bot, color=c, lw=0)
            bot += v
        fa = RP[reg]["fate"]
        rows_e.append([reg, *[RP[reg]["closure"][k]["mean"] for k, _, _ in PARTS],
                       fa["consumed"]["mean"], fa["lost"]["mean"]])
    kid = RP["kidney"]["fate"]
    axE.annotate(f"of what it spent,\n{kid['consumed']['mean']:.2f} built units and\n"
                 f"{kid['lost']['mean']:.2f} was lost",
                 xy=(1, 0.70), xytext=(1.0, 1.16), ha="center", fontsize=TN,
                 color=K["ink"], linespacing=1.4,
                 arrowprops=dict(arrowstyle="-|>", lw=0.6, color=K["ink"],
                                 mutation_scale=5))
    axE.set_xticks(range(3))
    axE.set_xticklabels([r[1] for r in REG], fontsize=TN)
    axE.set_ylabel("Share of the initial pool"); axE.set_ylim(0, 1.45)
    axE.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    axE.set_title("Where each pool ends up", fontsize=F.SZ["label"])
    F.tidy(axE)
    axE.legend(handles=[Patch(facecolor=c, label=l) for _, l, c in PARTS],
               loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=2,
               columnspacing=1.0)

    # ---------------------------------------------------- F: the three rows
    axF.text(0.0, 1.0, "and three regimes follow", fontsize=F.SZ["label"],
             va="top", ha="left")
    for x, h in ((0.46, "What holds it"), (0.86, "Restart?")):
        axF.text(x, 0.86, h, fontsize=TN, color=K["grey"], ha="center", va="top")
    FR = [(0.62, K["liver"], "renewing", "perfusion and room", "yes", K["good"]),
          (0.40, K["kidney"], "spent,\nwith loss", "nothing is left", "no", K["bad"]),
          (0.18, K["pancreas"], "spent,\nno loss", "perfusion at poor tips",
           "in part", K["warn"])]
    for y, col, name, binds, mark, mcol in FR:
        axF.plot([0.0, 0.99], [y + 0.11, y + 0.11], lw=0.5, color=K["pale"])
        axF.text(0.0, y, name, fontsize=TN, color=col, va="center", ha="left",
                 linespacing=1.45)
        axF.text(0.46, y, binds, fontsize=TN, color=K["ink"], va="center",
                 ha="center")
        axF.text(0.86, y, mark, fontsize=TN, color=mcol, va="center",
                 ha="center", fontweight="bold")
    axF.plot([0.0, 0.99], [0.07, 0.07], lw=0.5, color=K["pale"])
    axF.text(0.0, 0.015, "Panels A to C and F are schematic.", fontsize=TN,
             color=K["grey"], ha="left", va="top")

    F.panel_grid(fig, [[axA, axB, axC], [axD, axE, axF]],
                 ["a", "b", "c", "d", "e", "f"])
    F.audit(fig, "Fig1_model")
    fig.savefig(str(OUT / "Fig1_model.pdf"))
    fig.savefig(str(OUT / "Fig1_model.png"), dpi=300)
    plt.close(fig)
    with open(SD / "Fig1_pool_at_arrest.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["regime", "spent", "held_by_supply", "held_by_room",
                     "active", "consumed_into_units", "lost_to_attrition"])
        wr.writerows(rows_e)
    print("wrote Fig1_model")
else:
    print("Fig1_model skipped: needs results/paper/snapshot_liver.npz and regime_reps.json")

# ============================= Figure 2: what the binding resource is worth
E = load("elasticity.json")
Rb = load("robustness.json")
I = load("interventions.json")
if E and Rb and I:
    QS = ["ref", "nu_lo", "nu_hi", "q_lo", "q_hi"]
    fig, ax = plt.subplots(2, 3, figsize=(F.WIDTH_FULL, 112 * F.MM))
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.085, top=0.93,
                        wspace=0.42, hspace=0.62)
    (a0, a1, a2), (a3, a4, a5) = ax

    # A and B: the two elasticities, five settings each
    for a, key, lab, lims in ((a0, "P0", "$\\varepsilon_{P_0}$, seeded pool", (-0.15, 1.2)),
                              (a1, "q_set", "$\\varepsilon_{q}$, perfusion setpoint", (-3.3, 0.4))):
        for i, (reg, _, col) in enumerate(REG):
            v = [E[f"{reg}_{st}"][key] for st in QS if f"{reg}_{st}" in E and key in E[f"{reg}_{st}"]]
            xs = i + np.linspace(-0.22, 0.22, max(len(v), 1))
            a.errorbar(xs, [u["mean"] for u in v], yerr=[u["se"] for u in v],
                       fmt="o", ms=3.2, color=col, lw=0.8, capsize=0)
        a.axhline(0, color=K["grey"], lw=0.6, ls=(0, (2, 2)))
        a.set_xticks(range(3)); a.set_xticklabels([r[1] for r in REG], fontsize=TN)
        a.set_ylabel(lab); a.set_ylim(*lims); F.tidy(a)
    a0.set_title("What seeding more cells buys", fontsize=F.SZ["label"])
    a1.set_title("What lowering the threshold buys", fontsize=F.SZ["label"])

    # C: the ordering under 13 variants of the auxiliary rules
    order = ["reference", "root_surface", "root_random", "no_deficit", "no_load",
             "no_crowding", "no_mesenchyme", "territory", "binomial",
             "binomial_half_day", "binomial_quarter_day", "competing_fates",
             "no_self_renewal"]
    labels = ["reference", "root at surface", "root at random", "no deficit sprouting",
              "no load sprouting", "no crowding sprouting", "no mesenchymal signal",
              "territory assignment", "binomial nucleation", "binomial, half-day",
              "binomial, quarter-day", "competing fates", "pool consumable"]
    rows_c = []
    for j, v in enumerate(order):
        if v not in Rb:
            continue
        L, Pn = Rb[v]["liver"], Rb[v]["pancreas"]
        c = Rb[v].get("contrast")
        dd = c["mean"] if c else Pn["eps_pool"] - L["eps_pool"]
        se = c["se"] if c else float(np.hypot(Pn["eps_se"], L["eps_se"]))
        a2.errorbar(dd, j, xerr=se, fmt="o", ms=3.0, color=K["ink"], lw=0.8, capsize=0)
        a2.text(1.93, j, f"z = {dd / se:.1f}", va="center", ha="right",
                fontsize=TN, color=K["grey"])
        rows_c.append([v, dd, se, dd / se])
    a2.axvline(0, color=K["grey"], lw=0.6, ls=(0, (2, 2)))
    a2.set_yticks(range(len(order)))
    a2.set_yticklabels(labels, fontsize=TN)
    a2.invert_yaxis(); a2.set_xlim(-0.6, 2.0)
    a2.set_xlabel("spent-no-loss minus renewing")
    a2.set_title("The ordering, 13 variants", fontsize=F.SZ["label"])
    F.tidy(a2)

    # D, E, F: what relieving each bound does
    COND = (("none", "unchanged"), ("supply", "supply relieved"), ("room", "room relieved"))
    SHADE = (0.45, 0.75, 1.0)

    def tint(hexcol, f):
        c = np.array([int(hexcol[i:i + 2], 16) for i in (1, 3, 5)], float) / 255
        return tuple(1 - f * (1 - c))

    x = np.arange(3); w = 0.26
    rows_d = []
    for j, (cond, lab) in enumerate(COND):
        v = [100 * I[f"{r[0]}_{cond}"]["units_added"] for r in REG]
        e = [100 * I[f"{r[0]}_{cond}"]["units_added_se"] for r in REG]
        a3.bar(x + (j - 1) * w, v, w, yerr=e, lw=0, capsize=1.5,
               error_kw=dict(elinewidth=0.6),
               color=[tint(r[2], SHADE[j]) for r in REG])
        rows_d += [[r[0], cond, a, b] for r, a, b in zip(REG, v, e)]
        rel = [I[f"{r[0]}_{cond}"]["released"] for r in REG]
        rse = [I[f"{r[0]}_{cond}"]["se"] for r in REG]
        a4.bar(x + (j - 1) * w, [0 if u != u else u for u in rel], w,
               yerr=[0 if u != u else u for u in rse], lw=0, capsize=1.5,
               error_kw=dict(elinewidth=0.6),
               color=[tint(r[2], SHADE[j]) for r in REG])
    a3.set_ylabel("Units added by day 180 (%)")
    a3.set_title("What each relief builds", fontsize=F.SZ["label"])
    a4.set_ylabel("Share of the day-140 pool released")
    a4.set_title("What each relief releases", fontsize=F.SZ["label"])
    a4.set_ylim(0, None)
    for i, (reg, _, _) in enumerate(REG):
        if not I[f"{reg}_none"]["defined"]:
            a4.text(i, a4.get_ylim()[1] * 0.06, "nothing held", ha="center",
                    va="bottom", fontsize=TN, color=K["ink"], style="italic")

    pool = np.array([I[f"{r[0]}_none"]["pool"] for r in REG], float)
    lo = np.array([I[f"{r[0]}_none"]["pool_lo"] for r in REG], float)
    hi = np.array([I[f"{r[0]}_none"]["pool_hi"] for r in REG], float)
    flo = np.array([I[f"{r[0]}_none"]["pool_floor"] for r in REG], float)
    base = 0.05
    a5.bar(x, np.maximum(pool, base), 0.5, lw=0,
           color=[tint(r[2], 0.85) for r in REG])
    a5.errorbar(x, np.maximum(pool, base),
                yerr=[np.maximum(pool, base) - np.maximum(lo, base),
                      np.maximum(hi, base) - np.maximum(pool, base)],
                fmt="none", ecolor=K["ink"], elinewidth=0.6, capsize=1.5)
    for i, f in enumerate(flo):
        a5.plot([i - 0.3, i + 0.3], [f, f], ls=(0, (2, 1.5)), lw=0.8, color=K["ink"])
    for i, v in enumerate(pool):
        if v < base:
            a5.text(i, base * 1.25, "0", ha="center", va="bottom", fontsize=TN,
                    color=K["ink"])
    a5.set_yscale("log"); a5.set_ylim(base, 1e4)
    a5.set_ylabel("Progenitors held on day 140")
    a5.set_title("Pool left to release", fontsize=F.SZ["label"])

    for a in (a3, a4, a5):
        a.set_xticks(x); a.set_xticklabels([r[1] for r in REG], fontsize=TN)
        F.tidy(a)
    fig.legend(handles=[Patch(facecolor=tint(K["grey"], SHADE[j]), label=lab)
                        for j, (_, lab) in enumerate(COND)]
               + [plt.Line2D([0], [0], ls=(0, (2, 1.5)), lw=0.8, color=K["ink"],
                             label="pool a share needs")],
               loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=4,
               columnspacing=1.6)
    F.panel_grid(fig, [[a0, a1, a2], [a3, a4, a5]],
                 ["a", "b", "c", "d", "e", "f"])
    F.audit(fig, "Fig2_bounds")
    fig.savefig(str(OUT / "Fig2_bounds.pdf"))
    fig.savefig(str(OUT / "Fig2_bounds.png"), dpi=300)
    plt.close(fig)
    with open(SD / "Fig2_ordering.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["variant", "paired_difference", "se", "z"])
        wr.writerows(rows_c)
    with open(SD / "Fig2_interventions.csv", "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["regime", "intervention", "units_added_pct", "se_pct"])
        wr.writerows(rows_d)
    print("wrote Fig2_bounds")
else:
    print("Fig2_bounds skipped: needs elasticity.json, robustness.json and "
          "interventions.json")
