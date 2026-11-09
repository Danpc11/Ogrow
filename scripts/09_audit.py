#!/usr/bin/env python
"""Audit: every number a figure states must come from config or from results.

Reports the value each figure displays beside the value the pipeline
produces, so a stale caption cannot survive a rerun unnoticed.
"""
import json
import re
from itertools import combinations
from pathlib import Path
import numpy as np
from organgrow.config import CFG, ROOT

RES = ROOT / "results"
FIGS = ROOT / "figures"
bad = []


def check(label, shown, actual, tol=0.02):
    ok = (shown is None) or (abs(shown - actual) <= tol * max(abs(actual), 1e-9))
    mark = "ok  " if ok else "MISMATCH"
    print(f"  {mark} {label:<46} figure={shown}  pipeline={actual}")
    if not ok:
        bad.append(label)


def literal(fig, pattern):
    """Pull a number out of a figure script, so stale text is caught."""
    s = (FIGS / fig).read_text()
    m = re.search(pattern, s)
    return float(m.group(1)) if m else None


print("Fig 1  engine cost")
cost = json.load(open(RES / "cost.json"))
print(f"  tree traversal is O(E): {cost[0]['tree']*1e3:.3f} to "
      f"{cost[-1]['tree']*1e3:.1f} ms over {cost[0]['E']} to {cost[-1]['E']} edges")
ratio = [c["sparse"] / c["tree"] for c in cost if c["sparse"] == c["sparse"]]
print(f"  sparse/tree ratio: {min(ratio):.1f} to {max(ratio):.1f}x")

print("\nFig 3  closure mode and variance")
M = np.load(RES / "f3_main.npz")
for k in M.files:
    A = M[k]
    print(f"  {k:<9} sd(N) day 1 = {A[:,0].std(ddof=1):6.2f}   "
          f"day 70 = {A[:,-1].std(ddof=1):6.2f}   "
          f"ratio = {A[:,-1].std(ddof=1)/A[:,0].std(ddof=1):.3f}")
check("fig3 sigma_init stated in panel a",
      literal("fig3.py", r"sigma_\{init\}=([\d.]+)"), CFG["disorder"]["sigma_init"])

print("\nFig 4  self-averaging exponents")
A0 = json.load(open(RES / "scaling_nodisorder.json"))
N0 = np.array([r["N"] for r in A0])
for key in ("cv_N", "cv_mean_unit", "cv_mass", "vessel_identity",
            "density_discord"):
    y = np.array([r[key] for r in A0])
    a = -np.polyfit(np.log(N0), np.log(y), 1)[0]
    print(f"  {key:<18} a = {a:+.3f}")
sp = literal("fig4.py", r"sigma_\{param\}=([\d.]+)")
check("fig4 sigma_param stated in panel d", sp,
      CFG["disorder"]["sigma_param"] if sp is not None else CFG["disorder"]["sigma_param"])

print("\nFig 5  identifiability")
I = json.load(open(RES / "ident.json"))
for lab, v in I["sets"].items():
    c = v["cond"]
    print(f"  {lab:<24} rank {v['rank']}/{len(I['params'])}  "
          f"cond {c:.0f}" if c else f"  {lab:<24} rank {v['rank']}")
c1 = I["sets"]["N trajectory"]["cond"]
c2 = I["sets"]["N + per-unit supply"]["cond"]
shown = re.search(r"condition ([\d.]+) to ([\d.]+)", (FIGS / "fig5.py").read_text())
check("fig5 condition before", float(shown.group(1)) if shown else None, c1, 0.01)
check("fig5 condition after", float(shown.group(2)) if shown else None, c2, 0.01)
nrep = I["sets"]["N + per-unit supply"]["nrep"]
worst = I["params"][int(max(range(len(nrep)), key=lambda i: nrep[i]))]
best = I["params"][int(min(range(len(nrep)), key=lambda i: nrep[i]))]
lam = nrep[I["params"].index(worst)]
beta = nrep[I["params"].index(best)]
shown = re.search(r"(\d+) replicates vs (\d+)", (FIGS / "fig5.py").read_text())
check("fig5 replicates, costliest parameter", float(shown.group(1)) if shown else None, lam, 0.01)
check("fig5 replicates, cheapest parameter", float(shown.group(2)) if shown else None, beta, 0.01)

print("\nFig 6  variance partition")
if not (RES / "partition.json").exists():
    print("  skipped: needs the human cohort (see data/README.md)")
else:
    P = json.load(open(RES / "partition.json"))
    print(f"  {P['kidneys']} kidneys, {P['donors']} donors")
    print(f"  within {P['intra']:.3f}  between {P['inter']:.3f}  "
          f"ratio {P['ratio']:.3f}  P = {P['p_permutation']:.4f}")
    k = len(P["per_invariant"])
    alpha = 0.05 / k                      # Bonferroni over the invariants
    print(f"  paired Wilcoxon, Bonferroni alpha = {alpha:.4f}")
    nsig = 0
    for nm, v in P["per_invariant"].items():
        p = v["p_wilcoxon"]
        mark = "*" if p < alpha else " "
        nsig += p < alpha
        print(f"  {nm:<14} between {100*v['between_frac']:5.1f}%   "
              f"L-R {v['lr_rel_diff']:+6.1f}%   P = {p:.4f} {mark}")
    print(f"  {nsig} of {k} survive correction")
    shown = re.search(r"(\d+) of (\d+) survive", (FIGS / "fig6.py").read_text())
    if shown:
        check("fig6 invariants surviving correction",
              float(shown.group(1)), nsig, 0.001)

print("\n" + "=" * 60)
if bad:
    print("MISMATCHES:", ", ".join(bad))
    raise SystemExit(1)
print("every stated number matches the pipeline")
