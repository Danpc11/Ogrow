"""Closed-model figures of the paper: regimes (2), calibration (4), vasculature (5)."""
import json

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from organgrow import figstyle as F
from organgrow.config import ROOT, CFG

F.apply()
RES, OUT = ROOT / "results" / "paper", ROOT / "figures" / "out"
OUT.mkdir(parents=True, exist_ok=True)
K = F.C
if not (RES / "regime_liver.npz").exists():
    raise SystemExit("closed-model figures skipped: run scripts/17_paper_runs.py")
L = lambda t: np.load(RES / f"{t}.npz")
REG = [("liver", "Liver, supply-bound", K["liver"]),
       ("kidney", "Kidney, pool-bound with loss", K["kidney"]),
       ("pancreas", "Pancreas, pool-bound budget", K["pancreas"])]
numbers = {}

# ---------------------------------------------------------------- Figure 2
fig, ax = plt.subplots(1, 3, figsize=(F.WIDTH_FULL, 58 * F.MM))
fig.subplots_adjust(left=0.07, right=0.985, bottom=0.2, top=0.84, wspace=0.42)
for reg, name, col in REG:
    d = L(f"regime_{reg}"); t = np.arange(1, len(d["N"]) + 1)
    ax[0].plot(t, d["N"], color=col, label=name)
    ax[1].plot(t, d["q"] / 0.6, color=col)
    ax[2].plot(t, np.maximum(d["P"] / d["P"][0], 1e-3), color=col)
    numbers[f"regime_{reg}"] = dict(N=float(d["N"][-1]), q_over_set=float(d["q"][-1] / 0.6),
                                    pool_left=float(d["P"][-1] / d["P"][0]))
ax[0].set_ylabel("Units, $N$"); ax[0].legend(loc="lower right", fontsize=F.SZ["small"])
ax[1].axhline(1, color=K["grey"], lw=0.6, ls=(0, (2, 2)))
ax[1].set_ylabel("Supply per unit / setpoint")
ax[2].set_yscale("log"); ax[2].set_ylabel("Progenitor pool / initial")
for a in ax:
    a.set_xlabel("Day")
F.panel_grid(fig, [list(ax)], ["a", "b", "c"])
fig.savefig(str(OUT / "Fig2_closed.pdf")); fig.savefig(str(OUT / "Fig2_closed.png"), dpi=300)
plt.close(fig)

# ---------------------------------------------------------------- Figure 4
coh = np.vstack([L(f"cohort_k{k}")["N"] for k in range(8)])
rate = [L(f"rate_k{k}") for k in range(8)]
Nr = np.array([r["N"][-1] for r in rate]); nur = np.array([float(r["nu"]) for r in rate])
n = coh.mean(0); k0 = int(np.argmin(np.abs(n - n[-1] / 1.82)))
cv = lambda a: a.std(ddof=1) / a.mean()
rng = np.random.default_rng(5)
mat = []
for sd in (2, 4, 6, 8, 10):
    ci = np.mean([cv(coh[np.arange(8), np.clip(np.round(rng.normal(0, sd, 8)).astype(int) + k0,
                                                 0, coh.shape[1] - 1)]) for _ in range(300)])
    mat.append((sd, ci, cv(coh[:, -1])))
R = np.vstack([r["N"] for r in rate]); nR = R.mean(0); kR = int(np.argmin(np.abs(nR - nR[-1] / 1.82)))
rate_pt = (cv(R[:, kR]), cv(R[:, -1]))
rho, p = stats.spearmanr(np.log(nur), Nr)
numbers["calibration"] = dict(maturity=[dict(offset_sd=a, cv_start=b, cv_final=c) for a, b, c in mat],
                              rate=dict(cv_start=rate_pt[0], cv_final=rate_pt[1]),
                              rate_vs_final_spearman=float(rho), rate_vs_final_p=float(p),
                              mouse=dict(cv_start=0.229, cv_final=0.074))

fig, ax = plt.subplots(1, 2, figsize=(F.WIDTH_FULL * 0.72, 62 * F.MM))
fig.subplots_adjust(left=0.1, right=0.98, bottom=0.2, top=0.86, wspace=0.4)
lx = np.log10(nur)
ax[0].scatter(lx, Nr, s=14, color=K["liver"], lw=0)
g = np.polyfit(lx, Nr, 1); xx = np.linspace(lx.min(), lx.max(), 10)
ax[0].plot(xx, np.polyval(g, xx), color=K["ink"], lw=0.8)
ax[0].set_xlabel(r"Nucleation rate, $\log_{10}\nu$"); ax[0].set_ylabel("Final units, $N$")
ax[0].text(0.03, 0.95, f"Spearman {rho:.2f}, $P$ = {p:.3f}", transform=ax[0].transAxes,
           va="top", fontsize=F.SZ["small"])
ax[1].plot([m[1] for m in mat], [m[2] for m in mat], "o-", color=K["kidney"], ms=3,
           label="maturity spread")
ax[1].plot(*rate_pt, "s", color=K["liver"], ms=5, label="rate spread")
ax[1].plot(0.229, 0.074, "*", color=K["ink"], ms=9, label="mouse liver")
ax[1].set_xlabel("CV of $N$ at start of interval"); ax[1].set_ylabel("CV of $N$ when closed")
ax[1].legend(loc="upper left", fontsize=F.SZ["small"])
F.panel_grid(fig, [list(ax)], ["a", "b"])
fig.savefig(str(OUT / "Fig4_closed.pdf")); fig.savefig(str(OUT / "Fig4_closed.png"), dpi=300)
plt.close(fig)

# ---------------------------------------------------------------- Figure 5
b = CFG["cost_exponent_b"]; x = b + 2
fig, ax = plt.subplots(1, 3, figsize=(F.WIDTH_FULL, 58 * F.MM))
fig.subplots_adjust(left=0.07, right=0.985, bottom=0.2, top=0.84, wspace=0.45)
for i, (reg, name, col) in enumerate(REG[:2]):
    mr = L(f"regime_{reg}")["murray"]
    ax[0].scatter(i + rng.uniform(-0.18, 0.18, len(mr)), mr, s=5, color=col, alpha=0.5, lw=0)
    ax[0].plot([i - 0.25, i + 0.25], [np.median(mr)] * 2, color=K["ink"], lw=1)
    numbers[f"murray_{reg}"] = dict(n=int(len(mr)), median=float(np.median(mr)),
                                    q25=float(np.percentile(mr, 25)), q75=float(np.percentile(mr, 75)))
ax[0].axhline(1, color=K["grey"], lw=0.6, ls=(0, (2, 2)))
ax[0].set_xticks([0, 1]); ax[0].set_xticklabels(["Liver", "Kidney"])
ax[0].set_ylabel(rf"$r_p^{{{x:.1f}}}/\sum r_c^{{{x:.1f}}}$")
lr = 2.1; Rt = 600.0; r = np.logspace(-1.3, 0.8, 300)
tau = (1.0 / (lr / r ** 4 + Rt)) / r ** 3; rc = (lr / (3 * Rt)) ** 0.25
ax[1].loglog(r / rc, tau / tau.max(), color=K["ink"])
ax[1].axvline(1, color=K["grey"], lw=0.6, ls=(0, (2, 2)))
ax[1].set_xlabel("Radius / critical radius"); ax[1].set_ylabel("Wall shear / peak")
ax[1].text(0.3, 0.12, "closes", transform=ax[1].transAxes, fontsize=F.SZ["small"], ha="center")
ax[1].text(0.78, 0.12, "holds", transform=ax[1].transAxes, fontsize=F.SZ["small"], ha="center")
bb = np.linspace(0.5, 1.5, 100)
ax[2].plot(bb, 2 * 2 ** (-4 / (bb + 2)), color=K["ink"])
ax[2].axhline(1, color=K["grey"], lw=0.6, ls=(0, (2, 2)))
ax[2].axvline(b, color=K["liver"], lw=0.6)
ax[2].set_xlabel("Wall cost exponent, $b$"); ax[2].set_ylabel("Daughters' conductance / parent")
F.panel_grid(fig, [list(ax)], ["a", "b", "c"])
fig.savefig(str(OUT / "Fig5_closed.pdf")); fig.savefig(str(OUT / "Fig5_closed.png"), dpi=300)
plt.close(fig)

json.dump(numbers, open(ROOT / "results" / "paper_numbers.json", "w"), indent=1)
print(json.dumps(numbers, indent=1))
