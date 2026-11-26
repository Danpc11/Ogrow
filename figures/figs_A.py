"""The six figures and the graphical abstract of MANUSCRIPT_A.md.

Each figure makes one point, and its first panel is a drawing of what the
other panels measure. Throughout: filled symbols and solid lines are the
model, open symbols are measurements from the literature, and a colour means
one thing: blue, vermillion and green are the liver-like, kidney-like and
pancreas-like sources, and purple, orange and black are perfusion,
progenitors and room.

Usage:  python figures/figs_A.py            every figure whose results exist
        python figures/figs_A.py 3 4 S1     only those

B1 is not a figure of this manuscript: it is the panel set on recovery after
injury, kept here until the engineering manuscript has its own script.
"""
import csv
import json
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

from organgrow import figstyle as F
from organgrow.config import ROOT

F.apply()
RES, OUT, SD = ROOT / "results", ROOT / "figures" / "out", ROOT / "source_data"
OUT.mkdir(parents=True, exist_ok=True)
SD.mkdir(parents=True, exist_ok=True)
K, MM, W = F.C, F.MM, F.WIDTH_FULL
REG = ("liver", "kidney", "pancreas")
NAME = {"liver": "liver-like", "kidney": "kidney-like", "pancreas": "pancreas-like"}


# ----------------------------------------------------------------- helpers
def load(name):
    p = RES / name
    return json.load(open(p)) if p.exists() else None


def letter(fig, ax, s, dx=-0.01, dy=0.012):
    """A panel letter against the axes' tight box, in figure coordinates."""
    fig.canvas.draw()
    bb = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(bb.x0 + dx, bb.y1 + dy, s, fontsize=F.SZ["panel"] + 1, fontweight="bold",
             ha="left", va="bottom")


def blank(ax, xlim=(0, 1), ylim=(0, 1), equal=True):
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    if equal:
        ax.set_aspect("equal", adjustable="box")
    ax.axis("off")


def arrow(ax, a, b, col=None, lw=0.9, style="-|>", scale=6, ls="-", rad=0.0, z=5):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle=style, mutation_scale=scale, lw=lw,
                                 color=col or K["ink"], linestyle=ls, zorder=z,
                                 connectionstyle=f"arc3,rad={rad}", shrinkA=0, shrinkB=0))


def note(ax, x, y, s, size="small", col=None, ha="center", va="center", **kw):
    ax.text(x, y, s, fontsize=F.SZ[size], color=col or K["ink"], ha=ha, va=va, **kw)


def tree(ax, base, angle, length, depth, spread=0.55, shrink=0.74, lw=2.2, col=None,
         tips=None, rng=None, z=2):
    """A branching tree drawn by recursion; returns the tips through `tips`."""
    rng = rng or np.random.default_rng(1)
    end = (base[0] + length * np.cos(angle), base[1] + length * np.sin(angle))
    ax.plot([base[0], end[0]], [base[1], end[1]], "-", lw=lw, color=col or K["supply"],
            solid_capstyle="round", zorder=z)
    if depth == 0:
        if tips is not None:
            tips.append(end)
        return
    for sgn in (-1, 1):
        tree(ax, end, angle + sgn * spread * (0.8 + 0.4 * rng.random()), length * shrink,
             depth - 1, spread, shrink, max(lw * 0.72, 0.5), col, tips, rng, z)


def cap(ax, c, n, r=0.035, dot=0.008, col=None, rng=None, alpha=1.0, z=6):
    """Progenitors around a tip: n dots scattered in a disc."""
    rng = rng or np.random.default_rng(2)
    for _ in range(n):
        a, d = rng.uniform(0, 2 * np.pi), r * np.sqrt(rng.random())
        ax.add_patch(Circle((c[0] + d * np.cos(a), c[1] + d * np.sin(a)), dot,
                            fc=col or K["cells"], ec="none", alpha=alpha, zorder=z))


def units(ax, c, n, r=0.06, size=0.013, rng=None, fc=None, z=4):
    """Functional units near a tip: open circles."""
    rng = rng or np.random.default_rng(3)
    for _ in range(n):
        a, d = rng.uniform(0, 2 * np.pi), r * (0.45 + 0.55 * rng.random())
        ax.add_patch(Circle((c[0] + d * np.cos(a), c[1] + d * np.sin(a)), size,
                            fc=fc or "white", ec=K["ink"], lw=0.5, zorder=z))


def bar(ax, x, r, col, marker="o", ms=3.4, filled=True, **kw):
    ax.errorbar([x], [r["value"]], yerr=[[r["value"] - r["lo"]], [r["hi"] - r["value"]]],
                fmt=marker, ms=ms, mew=0.9, mec=col, mfc=col if filled else "white",
                ecolor=col, elinewidth=0.8, capsize=0, zorder=4, **kw)


def save(fig, name, rows=None, header=None):
    F.audit(fig, name, quiet=True)
    fig.savefig(str(OUT / f"{name}.pdf"))
    fig.savefig(str(OUT / f"{name}.png"), dpi=300)
    plt.close(fig)
    if rows:
        with open(SD / f"{name}.csv", "w", newline="") as f:
            wr = csv.writer(f)
            wr.writerow(header)
            wr.writerows(rows)
    print(f"wrote {name}")


# ------------------------------------------------------------------ Figure 1
def fig1():
    """One local law, and the three ways it stops."""
    reps = load("paper/regime_reps.json")
    snap = RES / "paper" / "snapshot_liver.npz"
    if reps is None or not snap.exists():
        return print("Figure 1 needs results/paper (scripts/17_paper_runs.py)")
    fig = plt.figure(figsize=(W, 128 * MM))
    gs = GridSpec(2, 24, figure=fig, left=0.06, right=0.99, top=0.95, bottom=0.12,
                  hspace=0.40, wspace=1.0, height_ratios=[1.05, 1.0])
    rows = []

    # A: the tissue and the rule at one tip
    a = fig.add_subplot(gs[0, 0:15]); blank(a, (0, 1.75), (0, 1))
    a.add_patch(FancyBboxPatch((0.03, 0.06), 0.74, 0.88, boxstyle="round,pad=0,rounding_size=0.2",
                               fc=K["box"], ec=K["pale"], lw=0.8, zorder=0))
    tips, rng = [], np.random.default_rng(4)
    a.plot([0.40, 0.40], [-0.02, 0.14], "-", lw=3.0, color=K["supply"], solid_capstyle="butt", zorder=2)
    tree(a, (0.40, 0.12), np.pi / 2, 0.20, 3, spread=0.62, lw=2.6, tips=tips, rng=rng)
    for i, t in enumerate(tips):
        units(a, t, 3 + (i % 3), r=0.055, size=0.011, rng=rng)
        cap(a, t, 7, r=0.026, dot=0.006, rng=rng)
    note(a, 0.40, 0.985, "growing tissue", "label")
    note(a, 0.45, 0.085, "perfusion", "small", K["supply"], ha="left")
    # the zoom
    t0 = max(tips, key=lambda p: p[0])
    a.add_patch(Circle(t0, 0.075, fc="none", ec=K["ink"], lw=0.6, ls=(0, (2, 1.5)), zorder=7))
    zc = (1.28, 0.60)
    a.add_patch(Circle(zc, 0.33, fc="white", ec=K["ink"], lw=0.7, zorder=1))
    a.plot([t0[0] + 0.05, zc[0] - 0.30], [t0[1] + 0.055, zc[1] + 0.14], lw=0.5, color=K["grey"], zorder=0)
    a.plot([t0[0] + 0.05, zc[0] - 0.30], [t0[1] - 0.055, zc[1] - 0.14], lw=0.5, color=K["grey"], zorder=0)
    a.plot([zc[0] - 0.30, zc[0] - 0.03], [zc[1] - 0.12, zc[1]], "-", lw=5, color=K["supply"],
           solid_capstyle="round", zorder=2)
    tip = (zc[0] - 0.03, zc[1])
    cap(a, tip, 16, r=0.085, dot=0.013, rng=rng)
    for dx, dy in ((0.17, 0.13), (0.21, -0.03), (0.12, -0.17), (0.02, 0.19)):
        a.add_patch(Circle((tip[0] + dx, tip[1] + dy), 0.036, fc="white", ec=K["ink"], lw=0.7, zorder=4))
    a.add_patch(Circle((tip[0] + 0.24, tip[1] + 0.12), 0.036, fc="none", ec=K["ink"], lw=0.7,
                       ls=(0, (1.5, 1.2)), zorder=4))
    arrow(a, (tip[0] + 0.07, tip[1] + 0.05), (tip[0] + 0.205, tip[1] + 0.105), K["cells"], lw=0.8)
    note(a, zc[0] - 0.19, zc[1] - 0.20, "tip", "small", K["supply"])
    note(a, tip[0] - 0.09, tip[1] + 0.125, "progenitors, $P$", "small", K["cells"])
    note(a, tip[0] + 0.20, tip[1] - 0.125, "units", "small")
    note(a, zc[0], 0.975, "one tip", "label")
    # the law
    note(a, zc[0], 0.20, r"units formed per day  $=\ \nu\,P\,(q-q^{*})_{+}$", "label")
    note(a, zc[0], 0.125, "until one of three runs out at the tip", "small", K["grey"])
    for x, s_, c_ in ((zc[0] - 0.33, "perfusion\n$q \\rightarrow q^{*}$", K["supply"]),
                      (zc[0], "progenitors\n$P \\rightarrow 0$", K["cells"]),
                      (zc[0] + 0.33, "room\nno site left", K["room"])):
        note(a, x, 0.035, s_, "small", c_, va="center", fontweight="bold")

    # B: a construct grown under the law
    b = fig.add_subplot(gs[0, 16:24])
    z = np.load(snap)
    sc = 7.4                               # tree coordinates are in units of lam
    st, en, rad, U = z["starts"] * sc, z["ends"] * sc, z["radius"], z["U"]
    for s0, e0, r0 in zip(st, en, rad):
        b.plot([s0[0], e0[0]], [s0[1], e0[1]], "-", color=K["supply"],
               lw=0.3 + 2.4 * r0 / rad.max(), solid_capstyle="round", zorder=2)
    b.scatter(U[:, 0], U[:, 1], s=9, fc="white", ec=K["ink"], lw=0.45, zorder=3)
    lo, hi = U[:, :2].min(0), U[:, :2].max(0)
    pad = 0.06 * (hi - lo).max()
    b.set_xlim(lo[0] - pad, hi[0] + pad); b.set_ylim(lo[1] - pad, hi[1] + pad)
    b.set_aspect("equal", adjustable="box"); b.axis("off")
    b.set_title(f"One liver-like construct at arrest:\n{len(U)} units on {len(st)} vessel segments",
                fontsize=F.SZ["label"])

    # C, D: the three sources in time
    c = fig.add_subplot(gs[1, 0:6]); d = fig.add_subplot(gs[1, 8:14])
    for reg in REG:
        runs = sorted((RES / "paper").glob(f"regime_{reg}*.npz"))
        N = np.array([np.load(p)["N"] for p in runs], float)
        P = np.array([np.load(p)["P"] for p in runs], float)
        t = np.arange(1, N.shape[1] + 1)
        for ax_, Y in ((c, N), (d, P / P[:, :1])):
            m = Y.mean(0)
            ax_.fill_between(t, np.percentile(Y, 10, 0), np.percentile(Y, 90, 0), color=K[reg],
                             alpha=0.18, lw=0)
            ax_.plot(t, m, "-", color=K[reg], lw=1.2, label=NAME[reg])
        rows += [[reg, "units", int(x), float(y)] for x, y in zip(t, N.mean(0))]
        rows += [[reg, "pool over seeded pool", int(x), float(y)] for x, y in zip(t, (P / P[:, :1]).mean(0))]
    c.set_xlabel("Day"); c.set_ylabel("Units"); c.set_ylim(0, None)
    c.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0))
    c.set_title("Arrest without a global size variable", fontsize=F.SZ["label"])
    d.set_xlabel("Day"); d.set_ylabel("Progenitors / seeded"); d.set_yscale("symlog", linthresh=0.05)
    d.set_ylim(0, 3); d.set_yticks([0, 0.1, 1]); d.set_yticklabels(["0", "0.1", "1"])
    d.set_title("The source renews, decays, or is spent", fontsize=F.SZ["label"])

    # E: what holds the pool at arrest
    e = fig.add_subplot(gs[1, 18:24])
    parts = (("spent", "spent or lost", K["pale"]), ("supply_held", "held: no perfusion", K["supply"]),
             ("space_held", "held: no room", K["room"]))
    for i, reg in enumerate(REG):
        left = 0.0
        tot = sum(reps[reg]["closure"][k]["mean"] for k, _, _ in parts)
        for k, lab, col in parts:
            v = reps[reg]["closure"][k]["mean"] / tot
            e.barh(i, v, left=left, color=col, edgecolor="white", lw=0.6, height=0.62,
                   label=lab if i == 0 else None)
            if v > 0.12:
                e.text(left + v / 2, i, f"{v:.2f}", ha="center", va="center", fontsize=F.SZ["small"],
                       color="white" if col != K["pale"] else K["ink"])
            rows.append([reg, lab, "", float(v)])
            left += v
    e.set_yticks(range(3)); e.set_yticklabels([NAME[r] for r in REG]); e.invert_yaxis()
    for lab_, reg in zip(e.get_yticklabels(), REG):
        lab_.set_color(K[reg])
    e.set_xlim(0, 1); e.set_xlabel("Share of the progenitor pool at arrest")
    e.legend(loc="upper center", bbox_to_anchor=(0.45, -0.24), ncol=3, columnspacing=0.8,
             handletextpad=0.4, handlelength=0.9)
    e.set_title("Each stops against a different bound", fontsize=F.SZ["label"])
    for ax_ in (c, d, e):
        F.tidy(ax_)
    for ax_, s in ((a, "A"), (b, "B"), (c, "C"), (d, "D"), (e, "E")):
        letter(fig, ax_, s)
    save(fig, "FigA1_law", rows, ["regime", "quantity", "day", "value"])


# ------------------------------------------------------------------ Figure 2
def fig2():
    """Which property of the cell source decides the bound, and what follows."""
    R, El, Iv = load("regimes.json"), load("elasticity.json"), load("interventions.json")
    if R is None or El is None or Iv is None or "_elasticity_map" not in R:
        return print("Figure 2 needs results/regimes.json (with the elasticity map), "
                     "elasticity.json and interventions.json")
    meta = R["_meta"]
    fig = plt.figure(figsize=(W, 150 * MM))
    gs = GridSpec(2, 24, figure=fig, left=0.10, right=0.99, top=0.93, bottom=0.08,
                  wspace=1.0, hspace=0.62, height_ratios=[1.0, 0.85])
    rows = []

    # A: the factorial, as a map of which bound ends growth
    a = fig.add_subplot(gs[0, 0:10])
    pools, turn = meta["pool"], [t["name"] for t in meta["turnover"]]
    named = {(0.0, 60.0, "static"): "pancreas", (0.0, 300.0, "decaying"): "kidney",
             (1.0, 2600.0, "balanced"): "liver"}
    ylab, y = [], 0
    for w in meta["renewal"]:
        for tn in turn:
            for x, P0 in enumerate(pools):
                c = R[f"w{w:.1f}_P{P0:.0f}_{tn}"]
                shares = {"spent": c["spent"], "supply": c["supply_held"], "room": c["space_held"]}
                hold = max(shares, key=shares.get)
                col = {"spent": K["pale"], "supply": K["supply"], "room": K["room"]}[hold]
                a.add_patch(Rectangle((x - 0.46, y - 0.44), 0.92, 0.88, fc=col, ec="white", lw=0.8))
                a.text(x, y, f"{c['N']:.0f}", ha="center", va="center", fontsize=F.SZ["small"],
                       color=K["ink"] if hold == "spent" else "white")
                reg = named.get((w, P0, tn))
                if reg:
                    a.add_patch(Rectangle((x - 0.46, y - 0.44), 0.92, 0.88, fc="none", ec=K[reg],
                                          lw=1.8, zorder=5))
                rows.append(["a", w, P0, tn, c["N"], hold])
            ylab.append(tn); y += 1
        y += 0.5
    a.set_xlim(-0.6, 2.6); a.set_ylim(y - 0.9, -0.6)
    a.set_xticks(range(3)); a.set_xticklabels([f"{p:.0f}" for p in pools])
    a.set_yticks([0, 1, 2, 3.5, 4.5, 5.5]); a.set_yticklabels(ylab)
    a.set_xlabel("Progenitors seeded")
    a.tick_params(length=0)
    for sp in a.spines.values():
        sp.set_visible(False)
    a.text(-1.55, 1, "cell is spent", rotation=90, ha="center", va="center", fontsize=F.SZ["label"])
    a.text(-1.55, 4.5, "cell renews", rotation=90, ha="center", va="center", fontsize=F.SZ["label"])
    a.text(-1.25, -0.95, "turnover", ha="left", va="center", fontsize=F.SZ["small"], color=K["grey"])
    a.set_title("What ends growth (number: units)", fontsize=F.SZ["label"], pad=10)
    for lab_, col in (("progenitors gone", K["pale"]), ("perfusion", K["supply"])):
        a.plot([], [], "s", ms=6, color=col, label=lab_)
    for reg in REG:
        a.plot([], [], "s", ms=6, mfc="none", mec=K[reg], mew=1.5, label=NAME[reg])
    a.legend(loc="upper center", bbox_to_anchor=(0.42, -0.14), ncol=3, columnspacing=0.8,
             handletextpad=0.3)

    # B: how many settings each axis flips
    b = fig.add_subplot(gs[0, 12:17])
    flips = [("turnover", meta["turnover_flips"], meta["other_cells"]),
             ("renewal", meta["renewal_flips"], meta["cells"]),
             ("pool size", meta["pool_flips"], meta["other_cells"])]
    for i, (nm, k, n) in enumerate(flips):
        b.bar(i, k / n, color=K["ink"] if i == 0 else K["grey"], width=0.62)
        b.text(i, k / n + 0.03, f"{k} of {n}", ha="center", va="bottom", fontsize=F.SZ["small"])
        rows.append(["b", nm, k, n, k / n, ""])
    b.set_xticks(range(3)); b.set_xticklabels([f[0] for f in flips], rotation=25, ha="right")
    b.set_ylim(0, 1.05); b.set_ylabel("Settings in which changing it\nchanges the bound")
    b.set_title("Turnover decides most", fontsize=F.SZ["label"], pad=10)

    # C: what the unit count answers to
    c = fig.add_subplot(gs[0, 19:24])
    for i, reg in enumerate(REG):
        e = El[f"{reg}_ref"]
        for j, (key, mk, fill) in enumerate((("P0", "o", True), ("q_set", "s", False))):
            v, se = abs(e[key]["mean"]), e[key]["se"]
            c.errorbar([v], [i + (j - 0.5) * 0.28], xerr=[se], fmt=mk, ms=4.2, mew=0.9, mec=K[reg],
                       mfc=K[reg] if fill else "white", ecolor=K[reg], elinewidth=0.8, capsize=0)
            rows.append(["c", reg, key, "", e[key]["mean"], se])
    c.plot([], [], "o", ms=4.2, color=K["ink"], label="to the pool seeded")
    c.plot([], [], "s", ms=4.2, mfc="white", mec=K["ink"], label="to the perfusion threshold")
    c.axvline(1, lw=0.6, ls=(0, (2, 2)), color=K["grey"])
    c.set_yticks(range(3)); c.set_yticklabels([NAME[r] for r in REG]); c.invert_yaxis()
    for lab_, reg in zip(c.get_yticklabels(), REG):
        lab_.set_color(K[reg])
    c.set_xlim(0, 2.3); c.set_xlabel(r"Elasticity, $|\varepsilon|$")
    c.legend(loc="upper center", bbox_to_anchor=(0.30, -0.22), ncol=1)
    c.set_title("Three sources,\nthree sensitivities", fontsize=F.SZ["label"], pad=6)
    # D: the elasticity against nearness to capacity; the domain moves one and not the other
    d = fig.add_subplot(gs[1, 1:11])
    Cp = load("capacity.json")
    col = {"spent": K["pale"], "supply": K["supply"], "room": K["room"], "still active": "white"}
    cap0 = Cp["capacity"]["0.6"] if Cp else None
    xof = (lambda n_: n_ / cap0) if Cp else (lambda n_: n_)
    for w in meta["renewal"]:
        for P0 in pools:
            for tn in turn:
                c_ = R[f"w{w:.1f}_P{P0:.0f}_{tn}"]
                reg = named.get((w, P0, tn))
                d.errorbar([xof(c_["N"])], [c_["eps_pool"]], yerr=[c_["eps_pool_se"]], fmt="o",
                           ms=5.6 if reg else 4.6, mfc=col[c_["holds"]], mec=K[reg] if reg else K["grey"],
                           mew=1.6 if reg else 0.6, ecolor=K["grey"], elinewidth=0.6, capsize=0, zorder=3)
                rows.append(["d", w, P0, tn, c_["eps_pool"], xof(c_["N"])])
    for lab_, k_ in (("progenitors gone at arrest", "spent"), ("held by perfusion", "supply")):
        d.plot([], [], "o", ms=4.6, mfc=col[k_], mec=K["grey"], mew=0.6, label=lab_)
    if Cp:
        for reg, dm in Cp["domain"].items():
            pts = sorted((float(k_), v) for k_, v in dm.items() if k_ != "change")
            xs, ys = [v["fill"] for _, v in pts], [v["eps"] for _, v in pts]
            d.plot(xs, ys, "-", lw=1.1, color=K[reg], zorder=4)
            d.plot(xs, ys, "D", ms=4.4, mfc="white", mec=K[reg], mew=1.2, zorder=5,
                   label="same pool, domain changed" if reg == list(Cp["domain"])[0] else None)
            arrow(d, (xs[0], ys[0]), (xs[0] + 0.55 * (xs[1] - xs[0]), ys[0] + 0.55 * (ys[1] - ys[0])),
                  K[reg], lw=1.1, scale=7, z=6)
            for (sc_, v_) in pts:
                rows.append(["d", f"{reg}, domain {sc_:g}", v_["units"], "", v_["eps"], v_["fill"]])
        d.text(0.97, 0.95, f"Spearman {Cp['all']['rho_fill']:+.2f} with units/capacity,\n"
               f"{Cp['all']['rho_units']:+.2f} with units alone", transform=d.transAxes,
               ha="right", va="top", fontsize=F.SZ["small"])
        d.set_xlabel("Units at arrest / capacity of the domain")
        d.set_title("Nearness to capacity, not size, weakens\nthe dependence on progenitors",
                    fontsize=F.SZ["label"])
    else:
        d.set_xlabel("Units at arrest")
        d.set_title("Size answers the pool less as the tissue nears\nwhat perfusion supports",
                    fontsize=F.SZ["label"])
    d.set_ylabel(r"Elasticity of unit count to the pool, $\varepsilon_P$")
    d.set_ylim(0, 1.0); d.set_xlim(0, None)
    d.legend(loc="lower left")

    # E: what can restart an arrested construct
    e2 = fig.add_subplot(gs[1, 14:24])
    for i_, reg in enumerate(REG):
        for j_, (kind, colr) in enumerate((("none", K["pale"]), ("supply", K["supply"]), ("room", K["room"]))):
            r_ = Iv[f"{reg}_{kind}"]
            if r_["defined"]:
                e2.bar(i_ + (j_ - 1) * 0.26, 100 * r_["units_added"], yerr=100 * r_["units_added_se"],
                       width=0.24, color=colr, error_kw=dict(lw=0.7, ecolor=K["grey"]),
                       label={"none": "left alone", "supply": "perfusion threshold lowered",
                              "room": "spacing lowered"}[kind] if i_ == 0 else None)
            rows.append(["e", reg, kind, "", r_["units_added"] if r_["defined"] else "", ""])
    e2.text(1, 1.5, "no progenitors left:\nnothing to restart", ha="center", va="bottom",
            fontsize=F.SZ["small"], color=K["kidney"])
    e2.set_xticks(range(3)); e2.set_xticklabels([NAME[r_] for r_ in REG])
    for lab_, reg in zip(e2.get_xticklabels(), REG):
        lab_.set_color(K[reg])
    e2.set_ylabel("Units gained in 40 days (%)")
    e2.set_title("Growth restarts only where\nprogenitors remain", fontsize=F.SZ["label"])
    e2.legend(loc="upper right")
    for ax_ in (b, c, d, e2):
        F.tidy(ax_)
    for ax_, s_ in ((a, "A"), (b, "B"), (c, "C"), (d, "D"), (e2, "E")):
        letter(fig, ax_, s_)
    save(fig, "FigA2_axes", rows, ["panel", "a", "b", "c", "value", "extra"])


# ------------------------------------------------------------------ Figure 3
def fig3():
    """The model beside three published perturbations and a human cohort."""
    V = load("validation.json")
    if V is None:
        return print("Figure 3 needs results/validation.json (scripts/36_validation.py)")
    fig = plt.figure(figsize=(W, 122 * MM))
    gs = GridSpec(2, 24, figure=fig, left=0.07, right=0.99, top=0.96, bottom=0.14,
                  hspace=0.34, wspace=1.0, height_ratios=[0.62, 1.0])
    rows = []

    # A: the experiment and its counterpart in the model
    a = fig.add_subplot(gs[0, :]); blank(a, (0, 3.6), (0, 1))
    rng = np.random.default_rng(7)

    def bud(x, n, dead=0):
        a.add_patch(Circle((x, 0.52), 0.21, fc=K["box"], ec=K["pale"], lw=0.8))
        pts = [(x + 0.14 * np.sqrt(rng.random()) * np.cos(t), 0.52 + 0.14 * np.sqrt(rng.random()) * np.sin(t))
               for t in rng.uniform(0, 2 * np.pi, n)]
        for i, p_ in enumerate(pts):
            if i < dead:
                a.plot([p_[0] - 0.014, p_[0] + 0.014], [p_[1] - 0.014, p_[1] + 0.014], lw=0.7, color=K["grey"])
                a.plot([p_[0] - 0.014, p_[0] + 0.014], [p_[1] + 0.014, p_[1] - 0.014], lw=0.7, color=K["grey"])
            else:
                a.add_patch(Circle(p_, 0.016, fc=K["cells"], ec="none"))

    def organ(x, r, col, lab):
        a.add_patch(Circle((x, 0.52), r, fc=col, ec="none", alpha=0.85))
        note(a, x, 0.13, lab, "small")

    note(a, 0.08, 0.93, "Experiment", "label", ha="left", fontweight="bold")
    bud(0.30, 22); note(a, 0.30, 0.13, "organ bud with\nits progenitors", "small")
    arrow(a, (0.58, 0.52), (0.94, 0.52)); note(a, 0.76, 0.70, "toxin, in a\ntime window", "small", K["grey"])
    bud(1.22, 22, dead=13); note(a, 1.22, 0.13, "a known share\nis killed", "small")
    arrow(a, (1.50, 0.52), (1.80, 0.52)); note(a, 1.65, 0.64, "growth", "small", K["grey"])
    organ(2.00, 0.105, K["pancreas"], "pancreas:\nsmall for life")
    organ(2.50, 0.19, K["liver"], "liver:\nnormal in days")
    a.plot([2.80, 2.80], [0.08, 0.92], lw=0.6, color=K["pale"])
    note(a, 2.88, 0.93, "Model", "label", ha="left", fontweight="bold")
    note(a, 3.22, 0.62, "seed the same\nshare of the pool,\nchange nothing else,\nread the size", "small")
    note(a, 3.22, 0.14, "parameters\nnot optimized\nto these data", "small", K["grey"])

    def model(ax_, reg, key, col, label):
        pts = sorted((float(k), v[key]) for k, v in V["ablation"][reg].items())
        pts.append((1.0, dict(value=1.0, lo=1.0, hi=1.0)))
        ax_.plot([p_[0] for p_ in pts], [p_[1]["value"] for p_ in pts], "-", lw=1.0, color=col, zorder=2)
        for i, (x, r) in enumerate(pts):
            bar(ax_, x, r, col, label=label if i == 0 else None)
            rows.append([reg, "model", key, x, r["value"], r["lo"], r["hi"]])

    def measured(ax_, reg, col, marker, study, label, dx=0.0):
        pts = [r for r in V["published"][reg] if r["study"] == study]
        for i, r in enumerate(pts):
            err = [[r["left"] - r["lo"]], [r["hi"] - r["left"]]] if "lo" in r else None
            ax_.errorbar([r["kept"] + dx], [r["left"]], yerr=err, fmt=marker, ms=5.2, mfc="white",
                         mec=col, mew=1.1, ecolor=col, elinewidth=0.8, capsize=1.8, zorder=5,
                         label=label if i == 0 else None)
            rows.append([reg, study, r["at"], r["kept"], r["left"], r.get("lo", ""), r.get("hi", "")])

    def frame(ax_, title, ylab):
        ax_.plot([0, 1], [0, 1], ls=(0, (2, 2)), lw=0.8, color=K["grey"], zorder=1)
        ax_.set_xlim(0.15, 1.05); ax_.set_ylim(0, 1.12)
        ax_.set_xlabel("Share of progenitors left"); ax_.set_ylabel(ylab)
        ax_.set_title(title, fontsize=F.SZ["label"])

    b = fig.add_subplot(gs[1, 1:11])
    model(b, "pancreas", "mass", K["pancreas"], "pancreas-like, model")
    model(b, "liver", "mass", K["liver"], "liver-like, model")
    measured(b, "pancreas", K["pancreas"], "s", "Stanger 2007", "mouse pancreas")
    measured(b, "liver", K["liver"], "s", "Stanger 2007", "mouse liver")
    frame(b, "Pancreas and liver", "Share of the organ that results")
    b.text(0.93, 0.80, "in proportion", rotation=38, fontsize=F.SZ["small"], color=K["grey"],
           ha="right", va="bottom")
    b.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=2, columnspacing=0.8,
             handletextpad=0.3)

    c = fig.add_subplot(gs[1, 14:24])
    model(c, "kidney", "N", K["kidney"], "kidney-like, model")
    measured(c, "kidney", K["kidney"], "o", "Muthukrishnan 2018", "mouse, survivors divide faster", dx=-0.03)
    measured(c, "kidney", K["ink"], "D", "Cebrian 2014", "mouse, they do not", dx=0.03)
    frame(c, "Kidney: one model, two experiments", "Share of the nephrons that results")
    c.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=1)

    for ax_ in (b, c):
        F.tidy(ax_)
    for ax_, s_ in ((a, "A"), (b, "B"), (c, "C")):
        letter(fig, ax_, s_)
    save(fig, "FigA3_measurements", rows, ["regime", "series", "quantity", "x", "y", "lo", "hi"])


# --------------------------------------------------------------- Figure S1
def figS1():
    """Two organs grown alike: the model's seeds beside a person's two kidneys."""
    V = load("validation.json")
    if V is None:
        return print("Figure S1 needs results/validation.json")
    fig = plt.figure(figsize=(88 * MM, 62 * MM))
    d = fig.add_axes([0.26, 0.20, 0.70, 0.64])
    rows, rng = [], np.random.default_rng(7)
    seeds = np.asarray(V["seedwise"]["kidney"]["mass"], float)
    groups = [("constructs\ngrown alike", 100 * (np.log(seeds) - np.log(seeds).mean()), K["kidney"], True)]
    if V.get("pairs"):
        groups.append(("a person's\ntwo kidneys", 100 * np.asarray(V["pairs"]["deviation"], float),
                       K["grey"], False))
    for i, (lab, dv, col, filled) in enumerate(groups):
        x = i + rng.uniform(-0.16, 0.16, len(dv))
        d.plot(x, dv, "o", ms=3, mew=0.7, mec=col, mfc=col if filled else "white", alpha=0.9, zorder=2)
        sd = dv.std(ddof=1)
        d.plot([i + 0.30, i + 0.30], [-sd, sd], "-", lw=1.5, color=K["ink"], zorder=3)
        d.text(i, 1.01, f"s.d. {sd:.1f} %", ha="center", va="bottom", fontsize=F.SZ["small"],
               transform=d.get_xaxis_transform())
        rows += [["pairs", lab.replace("\n", " "), "deviation, %", j, float(v_), "", ""] for j, v_ in enumerate(dv)]
    d.axhline(0, lw=0.6, color=K["grey"], zorder=1)
    d.set_xticks(range(len(groups))); d.set_xticklabels([g[0] for g in groups])
    d.set_xlim(-0.5, len(groups) - 0.25); d.set_ylabel("Deviation of one organ\nfrom its pair (%)")
    d.set_title("Two organs grown alike", fontsize=F.SZ["label"], pad=14)
    F.tidy(d)
    save(fig, "FigAS1_pairs", rows, ["regime", "series", "quantity", "x", "y", "lo", "hi"])


# ------------------------------------------------------------------ Figure 4
def fig4():
    """What the kidney needs that the reference model lacks."""
    V = load("kidney_coupling.json")
    if V is None:
        return print("Figure 4 needs results/kidney_coupling.json (scripts/37_kidney_coupling.py)")
    ta, T = V["ablate_day"], V["window"]
    names = list(V["nucleation"]); low = names[-1]
    m0 = V["measured"]["no response (Cebrian 2014)"]
    m1 = V["measured"]["response (Muthukrishnan 2018)"]
    t = np.arange(ta, T + 1) - ta; sl = slice(ta - 1, T)
    fig = plt.figure(figsize=(W, 126 * MM))
    gs = GridSpec(2, 24, figure=fig, left=0.07, right=0.99, top=0.96, bottom=0.14,
                  hspace=0.36, wspace=1.0, height_ratios=[0.72, 1.0])
    rows, rng = [], np.random.default_rng(11)

    # A: reciprocal induction, before and after an ablation
    a = fig.add_subplot(gs[0, :]); blank(a, (0, 3.6), (0, 1))

    def stage(x, n_cap, n_tips, title, sub_):
        a.plot([x, x], [0.10, 0.36], "-", lw=3.2, color=K["supply"], solid_capstyle="round")
        ends = [(x, 0.36)] if n_tips == 1 else [(x - 0.11, 0.56), (x + 0.11, 0.56)]
        if n_tips == 2:
            for e_ in ends:
                a.plot([x, e_[0]], [0.36, e_[1]], "-", lw=2.4, color=K["supply"], solid_capstyle="round")
        for e_ in ends:
            cap(a, (e_[0], e_[1] + 0.03), n_cap, r=0.085, dot=0.014, rng=rng)
        note(a, x, 0.93, title, "small", fontweight="bold")
        note(a, x, 0.02, sub_, "small", K["grey"])

    note(a, 0.02, 0.50, "intact", "label", rotation=90, fontweight="bold")
    stage(0.30, 9, 1, "one cap", "a tip and its\nprogenitors")
    arrow(a, (0.47, 0.42), (0.62, 0.42)); note(a, 0.545, 0.50, "divide", "small", K["grey"])
    stage(0.80, 18, 1, "two caps", "the cap doubles")
    arrow(a, (0.97, 0.42), (1.12, 0.42)); note(a, 1.045, 0.50, "branch", "small", K["grey"])
    stage(1.34, 9, 2, "the tip divides", "each branch\nkeeps one cap")
    a.plot([1.66, 1.66], [0.05, 0.95], lw=0.6, color=K["pale"])
    note(a, 1.72, 0.50, "ablated", "label", rotation=90, fontweight="bold")
    stage(1.98, 5, 1, "cap reduced", "40 % of the\nprogenitors killed")
    arrow(a, (2.15, 0.42), (2.34, 0.42)); note(a, 2.245, 0.52, "divide,\nno branching", "small", K["grey"])
    stage(2.54, 9, 1, "cap restored", "the tip waited")
    arrow(a, (2.71, 0.42), (2.90, 0.42)); note(a, 2.805, 0.50, "then branch", "small", K["grey"])
    stage(3.14, 9, 2, "fewer tips, for good", "the delay is\nnever made up")

    # B: pool, cap and sites after the ablation
    b = fig.add_subplot(gs[1, 0:6])
    c0 = V["conditions"][low]["0"]
    for key, lab, col in (("pool_by_day", "progenitors, total", K["cells"]),
                          ("cap_by_day", "progenitors per tip", K["ink"]),
                          ("sites_by_day", "tips", K["supply"])):
        y = np.asarray(c0[key])[sl]
        b.plot(t, y, "-", lw=1.1, color=col, label=lab)
        rows += [["b", lab, int(d_), float(v_)] for d_, v_ in zip(t, y)]
    xr = t[-1] + 5
    b.plot([xr] * 3, m0["cap"], "o", ms=4.2, mfc="white", mec=K["ink"], mew=1.0, clip_on=False,
           label="per tip, mouse")
    b.plot([xr + 6] * 3, m0["tips"], "D", ms=3.8, mfc="white", mec=K["supply"], mew=1.0,
           clip_on=False, label="tips, mouse")
    b.axhline(1, lw=0.6, color=K["grey"], zorder=0)
    b.set_xlim(0, xr + 10); b.set_ylim(0.3, 1.5)
    b.set_xlabel("Model days since the ablation"); b.set_ylabel("Ablated / control")
    b.text(xr + 3, 1.19, "mouse,\ndays 1 to 3", ha="center", va="bottom", fontsize=F.SZ["small"],
           color=K["grey"])
    b.set_title("The cap recovers; the tree does not", fontsize=F.SZ["label"])
    b.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=2, columnspacing=0.8,
             handletextpad=0.4)

    # C: units, with and without a proliferative response
    c = fig.add_subplot(gs[1, 8:14])
    for resp, lab, ls in (("0", "survivors divide as before", "-"),
                          ("1", "survivors divide twice as fast", (0, (3, 1.5)))):
        y = np.asarray(V["conditions"][low][resp]["units_by_day"])[sl]
        c.plot(t, y, ls=ls, lw=1.2, color=K["kidney"], label=lab)
        rows += [["c", lab, int(d_), float(v_)] for d_, v_ in zip(t, y)]
    for m_, lab in ((m0, "measured:\nno response"), (m1, "measured:\nresponse")):
        c.axhspan(min(m_["units"]), max(m_["units"]), color=K["pale"], lw=0, zorder=0)
        c.text(t[-1] - 1, np.mean(m_["units"]), lab, ha="right", va="center",
               fontsize=F.SZ["small"], color=K["grey"])
    c.axhline(1, lw=0.6, color=K["grey"], zorder=0)
    c.set_xlim(0, t[-1]); c.set_ylim(0.3, 1.25)
    c.set_xlabel("Model days since the ablation"); c.set_ylabel("Nephron-like units, ablated / control")
    c.set_title("One switch separates the experiments", fontsize=F.SZ["label"])
    c.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=1)

    # D: the deficit against how far perfusion limits nucleation
    d = fig.add_subplot(gs[1, 16:24])
    for j, n in enumerate(names):
        for resp, dx, filled in (("0", -0.13, True), ("1", 0.13, False)):
            bar(d, j + dx, V["conditions"][n][resp]["units"], K["kidney"], filled=filled, ms=3.8)
            r = V["conditions"][n][resp]["units"]
            rows.append(["d", f"{n}, response {resp}", j, r["value"]])
    for m_ in (m0, m1):
        d.axhspan(min(m_["units"]), max(m_["units"]), color=K["pale"], lw=0, zorder=0)
    d.axhline(1, lw=0.6, color=K["grey"], zorder=0)
    d.set_xticks(range(len(names)))
    d.set_xticklabels([f"$q^{{*}}$ {V['nucleation'][n]['q_set']:g}\n$\\nu$ {V['nucleation'][n]['nu']:g}"
                       for n in names])
    d.set_xlabel("Setting of the nucleation law\n" + r"perfusion limits less $\rightarrow$")
    d.set_xlim(-0.5, len(names) - 0.5); d.set_ylim(0.3, 1.25)
    d.set_ylabel("Units at the end, ablated / control")
    d.set_title("The deficit deepens where\nperfusion limits less", fontsize=F.SZ["label"])
    d.plot([], [], "o", ms=3.8, color=K["kidney"], label="as before")
    d.plot([], [], "o", ms=3.8, mfc="white", mec=K["kidney"], mew=0.9, label="twice as fast")
    d.legend(loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=2, columnspacing=0.8,
             handletextpad=0.3)
    for ax_ in (b, c, d):
        F.tidy(ax_)
    for ax_, s_ in ((a, "A"), (b, "B"), (c, "C"), (d, "D")):
        letter(fig, ax_, s_)
    save(fig, "FigA4_kidney", rows, ["panel", "series", "x", "y"])


# ------------------------------------------------------------------ Figure 5
def figB1():
    """For the engineering manuscript: restart, and what sets the size regained."""
    Iv, Ij, Fl = load("interventions.json"), load("injury.json"), load("inflow_law.json")
    if Iv is None or Ij is None or Fl is None:
        return print("Figure 5 needs results/interventions.json, injury.json and inflow_law.json")
    fig = plt.figure(figsize=(W, 122 * MM))
    gs = GridSpec(2, 24, figure=fig, left=0.07, right=0.99, top=0.95, bottom=0.10,
                  hspace=0.50, wspace=1.0)
    rows = []

    # A: restart
    a = fig.add_subplot(gs[0, 0:7])
    for i, reg in enumerate(REG):
        for j, (kind, col) in enumerate((("none", K["pale"]), ("supply", K["supply"]), ("room", K["room"]))):
            r = Iv[f"{reg}_{kind}"]
            x = i + (j - 1) * 0.26
            if r["defined"]:
                a.bar(x, 100 * r["units_added"], yerr=100 * r["units_added_se"], width=0.24, color=col,
                      error_kw=dict(lw=0.7, ecolor=K["grey"]),
                      label={"none": "left alone", "supply": "perfusion relieved",
                             "room": "room relieved"}[kind] if i == 0 else None)
            rows.append(["a", reg, kind, r["units_added"] if r["defined"] else "", ""])
    a.text(1, 1.5, "pool spent:\nnothing to restart", ha="center", va="bottom",
           fontsize=F.SZ["small"], color=K["kidney"])
    a.set_xticks(range(3)); a.set_xticklabels([NAME[r] for r in REG])
    for lab_, reg in zip(a.get_xticklabels(), REG):
        lab_.set_color(K[reg])
    a.set_ylabel("Units gained in 40 days (%)")
    a.set_title("Growth restarts only where cells remain", fontsize=F.SZ["label"])
    a.legend(loc="upper right")

    # B: the two ways to drive the circuit
    b = fig.add_subplot(gs[0, 9:15]); blank(b, (0, 1.5), (0, 1))
    for x, title, sym, txt, col in ((0.36, "pressure held", "$\\Delta P$", "flow rises with\nevery unit added", K["grey"]),
                                    (1.14, "flow held", "$Q$", "the tissue shares\na fixed delivery", K["supply"])):
        b.add_patch(Circle((x, 0.70), 0.13, fc="white", ec=col, lw=1.2))
        note(b, x, 0.70, sym, "label", col)
        b.plot([x, x], [0.57, 0.40], "-", lw=2.4, color=K["supply"], solid_capstyle="round")
        for dx_ in (-0.12, 0.12):
            b.plot([x, x + dx_], [0.40, 0.28], "-", lw=1.6, color=K["supply"], solid_capstyle="round")
            b.add_patch(Circle((x + dx_ * 1.25, 0.22), 0.05, fc="white", ec=K["ink"], lw=0.6))
        note(b, x, 0.95, title, "small", fontweight="bold")
        note(b, x, 0.05, txt, "small", K["grey"])
    b.set_title("Two ways to perfuse", fontsize=F.SZ["label"])

    # C: size regained after a 30 % loss
    c = fig.add_subplot(gs[0, 17:24])
    conds = (("as grown", "left alone", "o", False), ("reactivated", "pressure held", "s", True),
             ("reactivated, fixed inflow", "flow held", "D", True))
    rng = np.random.default_rng(5)
    for i, reg in enumerate(("liver", "kidney")):
        for j, (key, lab, mk, fill) in enumerate(conds):
            v = np.asarray(Ij[f"{reg}_{key}"]["mass"]["seeds"], float)
            x = i + (j - 1) * 0.27
            c.plot(x + rng.uniform(-0.05, 0.05, len(v)), v, mk, ms=2.6, mew=0.6, mec=K[reg],
                   mfc=K[reg] if fill else "white", alpha=0.75, zorder=2)
            c.plot([x - 0.09, x + 0.09], [v.mean()] * 2, "-", lw=1.4, color=K["ink"], zorder=3)
            rows += [["c", reg, lab, float(q), ""] for q in v]
    for (key, lab, mk, fill) in conds:
        c.plot([], [], mk, ms=3.4, mec=K["ink"], mfc=K["ink"] if fill else "white", mew=0.6, label=lab)
    c.axhline(1, lw=0.6, ls=(0, (2, 2)), color=K["grey"], zorder=1)
    c.set_xticks(range(2)); c.set_xticklabels([NAME["liver"], NAME["kidney"]])
    for lab_, reg in zip(c.get_xticklabels(), ("liver", "kidney")):
        lab_.set_color(K[reg])
    c.set_xlim(-0.55, 1.55); c.set_ylabel("Mass / mass before the loss")
    c.set_title("After a 30 % loss", fontsize=F.SZ["label"])
    c.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, columnspacing=0.6,
             handletextpad=0.2)

    # D: the size regained follows the flow delivered
    d = fig.add_subplot(gs[1, 0:7])
    for reg in ("liver", "kidney"):
        x = np.array([r["scale"] for r in Fl[reg]["rows"]])
        y = np.array([r["mass_rel"] for r in Fl[reg]["rows"]])
        e_ = np.array([r["mass_rel_sd"] for r in Fl[reg]["rows"]])
        d.errorbar(x, y, yerr=e_, fmt="o", ms=3.6, color=K[reg], elinewidth=0.8, capsize=0,
                   label=f"{NAME[reg]}, slope {Fl[reg]['slope']:.2f}")
        d.plot(x, y[2] * (x / x[2]) ** Fl[reg]["slope"], "-", lw=0.9, color=K[reg])
        rows += [["d", reg, float(a_), float(b_), float(c_)] for a_, b_, c_ in zip(x, y, e_)]
    d.plot([0.55, 1.5], [0.55, 1.5], ls=(0, (2, 2)), lw=0.8, color=K["grey"])
    d.set_xscale("log"); d.set_yscale("log")
    d.set_xticks([0.6, 0.8, 1.0, 1.4]); d.set_xticklabels(["0.6", "0.8", "1", "1.4"])
    d.set_yticks([0.6, 0.8, 1.0, 1.4]); d.set_yticklabels(["0.6", "0.8", "1", "1.4"])
    d.minorticks_off()
    d.set_xlabel("Flow delivered / flow before the loss"); d.set_ylabel("Mass regained / mass before")
    d.set_title("The flow delivered sets the size", fontsize=F.SZ["label"])
    d.legend(loc="upper left")

    # E: mass returns, units do not
    e = fig.add_subplot(gs[1, 9:15])
    for i, reg in enumerate(("liver", "kidney")):
        r = Ij[f"{reg}_reactivated, fixed inflow"]
        for j, (key, lab, col) in enumerate((("units", "units", K["pale"]), ("mass", "mass", K[reg]))):
            v = np.asarray(r[key]["seeds"], float)
            e.bar(i + (j - 0.5) * 0.34, v.mean(), width=0.32, color=col,
                  yerr=v.std(ddof=1), error_kw=dict(lw=0.7, ecolor=K["grey"]))
            e.text(i + (j - 0.5) * 0.34, 0.04, lab, rotation=90, ha="center", va="bottom",
                   fontsize=F.SZ["small"], color=K["ink"] if j == 0 else "white")
            rows.append(["e", reg, key, float(v.mean()), float(v.std(ddof=1))])
    e.axhline(1, lw=0.6, ls=(0, (2, 2)), color=K["grey"])
    e.set_xticks(range(2)); e.set_xticklabels([NAME["liver"], NAME["kidney"]])
    for lab_, reg in zip(e.get_xticklabels(), ("liver", "kidney")):
        lab_.set_color(K[reg])
    e.set_ylim(0, 1.2); e.set_ylabel("Relative to before the loss")
    e.set_title("Mass returns, units do not", fontsize=F.SZ["label"])

    # F: compensation by the units that remain
    f = fig.add_subplot(gs[1, 17:24])
    kept = np.linspace(0.45, 1.0, 50)
    f.plot(kept, 1 / kept, "-", lw=0.9, color=K["grey"], label="full compensation")
    for loss in ("0.3", "0.5"):
        r = Ij[f"kidney_compensation_{loss}"]
        u, m = np.asarray(r["units"]["seeds"], float), np.asarray(r["mass_per_unit"]["seeds"], float)
        f.plot(u, m, "o", ms=3, mec=K["kidney"], mfc=K["kidney"], alpha=0.6, mew=0)
        f.plot([u.mean()], [m.mean()], "o", ms=5.5, mfc="white", mec=K["kidney"], mew=1.2,
               label="kidney-like, model" if loss == "0.3" else None)
        rows.append(["f", "kidney", loss, float(u.mean()), float(m.mean())])
    f.plot([1.0], [1.0], "o", ms=3, color=K["ink"])
    f.set_xlim(0.42, 1.04); f.set_ylim(0.9, 2.3)
    f.set_xlabel("Units kept"); f.set_ylabel("Mass per unit / before the loss")
    f.set_title("Remaining units enlarge", fontsize=F.SZ["label"])
    f.legend(loc="upper right")
    for ax_ in (a, c, d, e, f):
        F.tidy(ax_)
    for ax_, s_ in ((a, "A"), (b, "B"), (c, "C"), (d, "D"), (e, "E"), (f, "F")):
        letter(fig, ax_, s_)
    save(fig, "FigB1_after_arrest", rows, ["panel", "a", "b", "c", "d"])


# ------------------------------------------------------------------ Figure 6
def fig5():
    """Predictions for ablations not yet done, and the reading the paper proposes."""
    V, Pr, Kc = load("validation.json"), load("predictions.json"), load("kidney_coupling.json")
    if V is None or Pr is None or Kc is None:
        return print("Figure 5 needs validation.json, predictions.json and kidney_coupling.json")
    fig = plt.figure(figsize=(W, 128 * MM))
    gs = GridSpec(2, 24, figure=fig, left=0.07, right=0.99, top=0.95, bottom=0.07,
                  hspace=0.62, wspace=1.0, height_ratios=[1.0, 0.62])
    rows = []

    # A: dose
    a = fig.add_subplot(gs[0, 0:7])
    for reg, col in (("pancreas", K["pancreas"]), ("liver", K["liver"])):
        known = sorted((float(k), v["mass"]) for k, v in V["ablation"][reg].items())
        new = sorted((float(k), v) for k, v in Pr["dose"][reg].items())
        both = sorted(known + new) + [(1.0, dict(value=1.0, lo=1.0, hi=1.0))]
        a.plot([p_[0] for p_ in both], [p_[1]["value"] for p_ in both], "-", lw=0.9, color=col, zorder=2)
        for x, r in known:
            bar(a, x, r, col, ms=2.8)
            rows.append(["a", reg, "model, compared", x, r["value"], r["lo"], r["hi"]])
        for x, r in new:
            alt = Pr.get("dose_alt", {}).get(reg, {}).get(f"{x:g}")
            if alt:      # the same dose under four other settings of the law
                a.plot([x + 0.022] * 2, [min(alt.values()), max(alt.values())], "-", lw=2.4,
                       color=K["pale"], solid_capstyle="butt", zorder=1)
            bar(a, x, r, col, marker="*", ms=8)
            rows.append(["a", reg, "model, predicted", x, r["value"], r["lo"], r["hi"]])
        for r in V["published"][reg]:
            a.plot([r["kept"]], [r["left"]], "s", ms=4.6, mfc="white", mec=col, mew=1.0, zorder=3)
    a.plot([0, 1], [0, 1], ls=(0, (2, 2)), lw=0.8, color=K["grey"], zorder=1)
    a.text(0.62, 1.01, "liver", color=K["liver"], fontsize=F.SZ["small"], ha="center")
    a.text(0.70, 0.52, "pancreas", color=K["pancreas"], fontsize=F.SZ["small"], ha="left")
    a.set_xlim(0.08, 1.05); a.set_ylim(0, 1.12)
    a.set_xlabel("Share of progenitors left"); a.set_ylabel("Share of the organ that results")
    a.set_title("Dose", fontsize=F.SZ["label"])

    # B: depth of a kidney ablation
    b = fig.add_subplot(gs[0, 9:16])
    for resp, filled in (("0", True), ("1", False)):
        pts = sorted((1 - float(k.split("|")[1]), v) for k, v in Pr["kidney"][resp].items()
                     if k.startswith("45|"))
        b.plot([p_[0] for p_ in pts] + [1.0], [p_[1]["value"] for p_ in pts] + [1.0],
               ls="-" if filled else (0, (3, 1.5)), lw=0.9, color=K["kidney"], zorder=2)
        for x, r in pts:
            done = abs(x - 0.6) < 1e-9
            bar(b, x, r, K["kidney"], marker="o" if done else "*", ms=3.2 if done else 8, filled=filled)
            rows.append(["b", f"response {resp}", "compared" if done else "predicted", x,
                         r["value"], r["lo"], r["hi"]])
    for r in V["published"]["kidney"]:
        b.plot([r["kept"]], [r["left"]], "s", ms=4.6, mfc="white", mec=K["ink"], mew=1.0, zorder=3)
    b.plot([0, 1], [0, 1], ls=(0, (2, 2)), lw=0.8, color=K["grey"], zorder=1)
    b.set_xlim(0.3, 1.05); b.set_ylim(0, 1.12)
    b.set_xlabel("Share of progenitors left"); b.set_ylabel("Share of the nephrons that results")
    b.set_title("Depth of a kidney ablation", fontsize=F.SZ["label"])

    # C: timing
    c = fig.add_subplot(gs[0, 18:24])
    days = sorted({int(k.split("|")[0]) for k in Pr["kidney"]["0"] if k.endswith("|0.4")})
    for resp, filled, dx in (("0", True, -0.12), ("1", False, 0.12)):
        for j, d_ in enumerate(days):
            r = Pr["kidney"][resp][f"{d_}|0.4"]
            done = d_ == 45
            bar(c, j + dx, r, K["kidney"], marker="o" if done else "*", ms=3.2 if done else 8, filled=filled)
            rows.append(["c", f"response {resp}", f"day {d_}", d_, r["value"], r["lo"], r["hi"]])
    c.axhline(1, lw=0.6, color=K["grey"], zorder=0)
    c.set_xticks(range(len(days))); c.set_xticklabels(["early", "as\npublished", "late"])
    c.set_xlim(-0.5, len(days) - 0.5); c.set_ylim(0, 1.12)
    c.set_xlabel("When 40 % are killed"); c.set_ylabel("Share of the nephrons that results")
    c.set_title("Timing", fontsize=F.SZ["label"])
    h = [plt.Line2D([], [], marker="*", ms=8, ls="", color=K["ink"], label="predicted"),
         plt.Line2D([], [], marker="o", ms=3.2, ls="", color=K["ink"], label="model, already compared"),
         plt.Line2D([], [], marker="s", ms=4.6, ls="", mfc="white", mec=K["ink"], label="measured"),
         plt.Line2D([], [], ls="-", lw=0.9, color=K["kidney"], label="survivors divide as before"),
         plt.Line2D([], [], ls=(0, (3, 1.5)), lw=0.9, color=K["kidney"], label="twice as fast")]
    fig.legend(handles=h, loc="upper center", bbox_to_anchor=(0.53, 0.435), ncol=5,
               columnspacing=1.2, handletextpad=0.4)

    # D: the reading: how size answers progenitor loss names the limiting resource
    d = fig.add_subplot(gs[1, 1:23])
    el = lambda left, kept: float(np.log(left) / np.log(kept))
    pub = V["published"]
    meas = [("liver", el(pub["liver"][0]["left"], pub["liver"][0]["kept"]), K["liver"]),
            ("kidney,\nsurvivors divide faster", np.mean([el(r["left"], r["kept"]) for r in pub["kidney"]
                                                           if r["study"].startswith("Muth")]), K["kidney"]),
            ("pancreas", el(pub["pancreas"][1]["left"], pub["pancreas"][1]["kept"]), K["pancreas"]),
            ("kidney,\nthey do not", np.mean([el(r["left"], r["kept"]) for r in pub["kidney"]
                                              if r["study"].startswith("Ceb")]), K["kidney"])]
    low = list(Kc["nucleation"])[-1]
    mod = [("liver-like", el(V["ablation"]["liver"]["0.35"]["mass"]["value"], 0.35), K["liver"]),
           ("coupled,\nfaster", el(Kc["conditions"][low]["1"]["units"]["value"], 0.6), K["kidney"]),
           ("pancreas-like", el(V["ablation"]["pancreas"]["0.39"]["mass"]["value"], 0.39), K["pancreas"]),
           ("coupled,\nas before", el(Kc["conditions"][low]["0"]["units"]["value"], 0.6), K["kidney"])]
    d.axhline(0, lw=0.8, color=K["ink"], zorder=1)
    d.axvspan(-0.08, 0.5, color=K["supply"], alpha=0.10, lw=0, zorder=0)
    d.axvspan(0.5, 1.55, color=K["cells"], alpha=0.12, lw=0, zorder=0)
    d.text(0.21, 1.42, "near perfusion-set size:\nsize hardly follows", ha="center", va="top",
           fontsize=F.SZ["small"], color=K["supply"], fontweight="bold")
    d.text(1.02, 1.42, "far from it:\nsize follows progenitor number", ha="center", va="top",
           fontsize=F.SZ["small"], color=K["ink"], fontweight="bold")
    for (lab, x, col), up in zip(meas, (1, 1, 1, 1)):
        d.plot([x], [0.28], "s", ms=5.5, mfc="white", mec=col, mew=1.2, zorder=4, clip_on=False)
        d.text(x, 0.42, lab, ha="center", va="bottom", fontsize=F.SZ["small"], color=col)
        rows.append(["d", "measured", lab.replace("\n", " "), x, "", "", ""])
    for i, (lab, x, col) in enumerate(mod):
        d.plot([x], [-0.28], "o", ms=5, color=col, zorder=4, clip_on=False)
        d.text(x + (-0.02 if i % 2 == 0 else 0.02), -0.42, lab, ha="right" if i % 2 == 0 else "left",
               va="top", fontsize=F.SZ["small"], color=col)
        rows.append(["d", "model", lab.replace("\n", " "), x, "", "", ""])
    d.text(-0.10, 0.28, "measured", ha="right", va="center", fontsize=F.SZ["small"])
    d.text(-0.10, -0.28, "model", ha="right", va="center", fontsize=F.SZ["small"])
    d.set_xlim(-0.08, 1.55); d.set_ylim(-1.0, 1.5)
    d.set_yticks([]); d.spines["left"].set_visible(False)
    d.set_xticks([0, 0.5, 1.0, 1.5])
    d.set_xlabel("Effective elasticity of organ size to progenitor number over one finite ablation:  "
                 r"$\log$(share of organ) / $\log$(share of progenitors left)")
    d.set_title("The organs on one axis: how strongly size follows progenitor number",
                fontsize=F.SZ["label"])
    for ax_ in (a, b, c):
        F.tidy(ax_)
    for ax_, s_ in ((a, "A"), (b, "B"), (c, "C"), (d, "D")):
        letter(fig, ax_, s_)
    save(fig, "FigA5_predictions", rows, ["panel", "series", "kind", "x", "y", "lo", "hi"])


# -------------------------------------------------------- graphical abstract
def abstract():
    """A square summary: one rule, three resources, three answers to a loss."""
    fig = plt.figure(figsize=(110 * MM, 110 * MM))
    a = fig.add_axes([0.02, 0.02, 0.96, 0.96]); blank(a)
    rng = np.random.default_rng(21)
    note(a, 0.5, 0.965, "One growth rule at every tip", "panel", fontweight="bold")
    tips = []
    a.plot([0.5, 0.5], [0.62, 0.69], "-", lw=3, color=K["supply"], solid_capstyle="round")
    tree(a, (0.5, 0.68), np.pi / 2, 0.075, 2, spread=0.75, lw=2.4, tips=tips, rng=rng)
    for t_ in tips:
        cap(a, t_, 6, r=0.018, dot=0.005, rng=rng)
        units(a, (t_[0], t_[1] + 0.012), 2, r=0.032, size=0.008, rng=rng)
    note(a, 0.5, 0.585, r"units form at a rate  $\nu\,P\,(q-q^{*})_{+}$  under constraints that act together",
         "small")
    cols = ((0.18, "perfusion", K["supply"], "liver-like", K["liver"], 0.92, "barely shrinks"),
            (0.50, "progenitors", K["cells"], "pancreas-like", K["pancreas"], 0.44, "shrinks in proportion"),
            (0.82, "progenitors that\ndrive branching", K["cells"], "kidney-like", K["kidney"], 0.62,
             "fewer tips,\nfewer nephrons"))
    for x0_, x1_ in ((0.325, 0.355), (0.645, 0.675)):      # the three act on one another
        arrow(a, (x0_, 0.445), (x1_, 0.445), K["ink"], lw=0.9, style="<|-|>", scale=6)
    for x, res, rc, reg, gc, size, txt in cols:
        arrow(a, (0.5 + (x - 0.5) * 0.25, 0.555), (x, 0.50), K["grey"], lw=0.7)
        a.add_patch(FancyBboxPatch((x - 0.135, 0.405), 0.27, 0.08, boxstyle="round,pad=0,rounding_size=0.02",
                                   fc="white", ec=rc, lw=1.1))
        note(a, x, 0.445, res, "small", rc, fontweight="bold")
        note(a, x, 0.375, "limits most", "small", K["grey"])
        note(a, x, 0.325, reg, "small", gc, fontweight="bold")
        a.add_patch(Circle((x - 0.055, 0.215), 0.048, fc=gc, ec="none", alpha=0.9))
        arrow(a, (x - 0.002, 0.215), (x + 0.03, 0.215), K["ink"], lw=0.7, scale=5)
        a.add_patch(Circle((x + 0.075, 0.215), 0.048 * np.sqrt(size), fc=gc, ec="none", alpha=0.9))
        note(a, x, 0.115, txt, "small")
    note(a, 0.5, 0.275, "after removing most of the progenitors", "small", K["grey"])
    note(a, 0.5, 0.035, "Coupled local constraints set where growth arrests\nand how size answers a loss of progenitors",
         "label", fontweight="bold")
    fig.savefig(str(OUT / "FigA0_graphical_abstract.pdf"))
    fig.savefig(str(OUT / "FigA0_graphical_abstract.png"), dpi=300)
    plt.close(fig)
    print("wrote FigA0_graphical_abstract")


FIGS = {"0": abstract, "1": fig1, "2": fig2, "3": fig3, "4": fig4, "5": fig5, "S1": figS1, "B1": figB1}

def manifest():
    """Record what the figures were drawn from: version, commit, platform, result files.

    A figure and the text can part company when one is drawn on another machine
    or from older results. The manifest says which results stood behind the
    files in `figures/out/`, so that a mismatch can be traced without redrawing.
    """
    import hashlib
    import platform
    import subprocess
    import organgrow
    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                                capture_output=True, text=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "status", "--porcelain", "--", "src", "scripts", "config"],
                                    cwd=ROOT, capture_output=True, text=True).stdout.strip())
    except OSError:
        commit, dirty = "unknown", False
    src = {}
    for name in ("paper/regime_reps.json", "regimes.json", "elasticity.json", "interventions.json",
                 "validation.json", "kidney_coupling.json", "predictions.json", "capacity.json"):
        f = RES / name
        if f.exists():
            src[name] = hashlib.sha256(f.read_bytes()).hexdigest()[:16]
    key = {}
    for name, path, pick in (("elasticity against units reached, factorial", "regimes.json",
                              lambda r: r.get("_elasticity_map", {}).get("rho_units")),
                             ("elasticity against units over capacity", "capacity.json",
                              lambda r: r.get("all", {}).get("rho_fill"))):
        r = load(path)
        if r:
            key[name] = pick(r)
    m = dict(version=organgrow.__version__, commit=commit + ("+uncommitted" if dirty else ""),
             platform=f"{platform.system()} {platform.machine()}", python=platform.python_version(),
             numpy=np.__version__, results=src, numbers_drawn=key)
    json.dump(m, open(OUT / "FigA_manifest.json", "w"), indent=1)
    print(f"wrote FigA_manifest.json: version {m['version']}, commit {m['commit']}, {m['platform']}")


if __name__ == "__main__":
    for k in (sys.argv[1:] or sorted(FIGS)):
        FIGS[k]()
    manifest()
