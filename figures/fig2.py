"""Figure 2 — the same engine generates liver, kidney and pancreas."""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from organgrow import grow, sdf

from organgrow import figstyle as F
from organgrow.config import ROOT

F.apply()
RES = ROOT / "results"
OUTDIR = ROOT / "figures" / "out"
OUTDIR.mkdir(parents=True, exist_ok=True)
K = F.C
MM = F.MM


ORG = [("liver", "Liver", K["liver"], 190),
       ("kidney", "Kidney", K["kidney"], 150),
       ("lung", "Lung", K["lung"], 175),
       ("pancreas", "Pancreas", K["pancreas"], 165)]
SNAP = (10, 25, 45, 70)
NREP = 10

runs, reps = {}, {}
for kind, _, _, _ in ORG:
    runs[kind] = grow(kind, seed=0, days=70, snap_at=SNAP)
    reps[kind] = [grow(kind, seed=s, days=70, snap_at=(),
                              sigma_init=0.45, sigma_param=0.04)[4]
                  for s in range(1, NREP + 1)]

fig = plt.figure(figsize=(F.WIDTH_FULL, 180 * F.MM))
gs = fig.add_gridspec(5, 4, height_ratios=[1, 1, 1, 1, 1.22],
                      hspace=0.26, wspace=0.52,
                      left=0.075, right=0.985, top=0.972, bottom=0.062)


def proj(P, rot=0.75, tilt=0.34):
    ca, sa, cb, sb = np.cos(rot), np.sin(rot), np.cos(tilt), np.sin(tilt)
    X = P[:, 0] * ca - P[:, 2] * sa
    Z = P[:, 0] * sa + P[:, 2] * ca
    Y = P[:, 1] * cb - Z * sb
    D = P[:, 1] * sb + Z * cb
    f = 600 / (600 + D)
    return X * f, Y * f, D, f


def silhouette(ax, kind, R):
    g = np.linspace(-R * 1.5, R * 1.5, 300)
    X, Y = np.meshgrid(g, g)
    ax.contour(X, Y, sdf(kind, X, Y, np.zeros_like(X)), levels=[0],
               colors=[K["pale"]], linewidths=0.7, zorder=0)


row_axes = []
for r, (kind, name, col, R) in enumerate(ORG):
    row_axes.append([])
    tr = runs[kind][4]
    for c, day in enumerate(SNAP):
        ax = fig.add_subplot(gs[r, c]); row_axes[r].append(ax)
        U, m, _ = tr.snapshots[day]
        keep = U[:, 2] <= 0
        Uk, mk = U[keep], m[keep]
        silhouette(ax, kind, R)
        if len(Uk):
            x, y, d, f = proj(Uk)
            o = np.argsort(-d)
            ax.scatter(x[o], y[o], s=np.clip(mk[o] * 170 * f[o] ** 2, 0.5, 9),
                       c=col, lw=0, alpha=0.70, zorder=2)
        g = np.linspace(-R * 1.5, R * 1.5, 300)
        XX, YY = np.meshgrid(g, g)
        inside = sdf(kind, XX, YY, np.zeros_like(XX)) < 0
        xs, ys = XX[inside], YY[inside]
        pad = 0.10 * max(np.ptp(xs), np.ptp(ys))
        ax.set_xlim(xs.min() - pad, xs.max() + pad)
        ax.set_ylim(ys.min() - pad, ys.max() + pad)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        if r == 0:
            ax.text(0.5, 1.02, "day %d" % day, transform=ax.transAxes,
                    ha="center", va="bottom", fontsize=6.2, color=K["ink"])
        ax.text(0.98, 0.03, "$N$ = %d" % U.shape[0], transform=ax.transAxes,
                ha="right", va="bottom", fontsize=5.9, color=K["grey"])
        if c == 0:
            ax.text(0.02, 0.97, name, transform=ax.transAxes, fontsize=6.6,
                    color=col, fontweight="bold", va="top", ha="left")

axd = fig.add_subplot(gs[4, 0])
for kind, name, col, _ in ORG:
    t = np.asarray(runs[kind][4].t)
    M = np.vstack([np.asarray(r.N) for r in reps[kind]])
    axd.fill_between(t, np.percentile(M, 10, axis=0), np.percentile(M, 90, axis=0),
                     color=col, alpha=0.18, lw=0)
    axd.plot(t, np.asarray(runs[kind][4].N), color=col, lw=1.1, label=name)
axd.set_xlabel("Day"); axd.set_ylabel("Units, $N$")
axd.legend(frameon=False, handlelength=1.4, borderpad=0.15, labelspacing=0.25,
           handletextpad=0.45, loc="upper left")

axe = fig.add_subplot(gs[4, 1])
for kind, name, col, _ in ORG:
    tr = runs[kind][4]
    t = np.asarray(tr.t)
    P = np.asarray(tr.P, dtype=float)
    if np.all(np.isnan(P)):
        cu = np.asarray(tr.cu)
        axe.plot(t, cu / cu[0], color=col, lw=1.1, ls=(0, (3, 1.6)))
    else:
        axe.plot(t, P / P[0], color=col, lw=1.1)
axe.set_xlabel("Day"); axe.set_ylabel("Fabrication capacity\n(fraction of initial)")
axe.set_ylim(-0.03, 1.08)
axe.plot([], [], color=K["grey"], lw=1.1, label="pool, $P/P_0$")
axe.plot([], [], color=K["grey"], lw=1.1, ls=(0, (3, 1.6)), label="spare, $c_u/c_u^0$")
axe.legend(frameon=False, handlelength=1.6, borderpad=0.15, labelspacing=0.25,
           handletextpad=0.45, loc="upper right")

axf = fig.add_subplot(gs[4, 2])
data = [runs[k][1] / np.median(runs[k][1]) for k, _, _, _ in ORG]
vp = axf.violinplot(data, positions=list(range(1, len(ORG) + 1)), widths=0.74,
                    showextrema=False, showmedians=False)
for b, (_, _, col, _) in zip(vp["bodies"], ORG):
    b.set_facecolor(col); b.set_alpha(0.30); b.set_edgecolor(col); b.set_linewidth(0.7)
rng = np.random.default_rng(0)
for i, (d, (_, _, col, _)) in enumerate(zip(data, ORG), start=1):
    sub = d[rng.choice(len(d), size=min(len(d), 200), replace=False)]
    axf.scatter(i + rng.normal(0, 0.055, len(sub)), sub, s=0.7, c=col, lw=0, alpha=0.35)
    axf.plot([i - 0.26, i + 0.26], [np.median(d)] * 2, color=K["ink"], lw=1.0)
axf.set_xticks(range(1, len(ORG) + 1))
axf.set_xticklabels([o[1] for o in ORG], fontsize=F.SZ["small"], rotation=18, ha="right")
axf.set_ylabel("Unit mass / median")

axg = fig.add_subplot(gs[4, 3])
for kind, name, col, _ in ORG:
    t = np.asarray(runs[kind][4].t)
    H = np.vstack([np.asarray(r.hypo, dtype=float) for r in reps[kind]])
    axg.fill_between(t, np.percentile(H, 10, axis=0), np.percentile(H, 90, axis=0),
                     color=col, alpha=0.18, lw=0)
    axg.plot(t, np.median(H, axis=0), color=col, lw=1.1)
axg.set_xlabel("Day")
axg.set_ylabel("Under-perfused units\n(fraction, $c<0.5\\,c_u$)")
axg.set_ylim(bottom=0)

F.panel_grid(fig, row_axes + [[axd, axe, axf, axg]],
             ["a","","","", "b","","","", "c","","","", "d","","","",
              "e","f","g","h"])
fig.savefig(str(OUTDIR / "Fig2.pdf"), dpi=600, bbox_inches="tight")
fig.savefig(str(OUTDIR / "Fig2.png"), dpi=400, bbox_inches="tight")
print("ok")
