#!/usr/bin/env python
"""Does the model recover from injury, and does function follow structure?

The organ is grown to closure (day 140), then a part of it is destroyed and
growth runs on for 80 more days. Injuring after closure, rather than during
growth, is what separates recovery from growth that would have happened
anyway; the control is the same run with no injury.

Three injuries, each a plane cut rather than a random removal, as a lobectomy
is. "parenchyma" removes the units past the cut and leaves the tree.
"vascular" occludes the vessels past the cut, which stops perfusion downstream
without removing tissue. "combined" does both.

Four readouts, two structural and two functional, all already in the model:
units, total unit mass, the mass whose supply meets the threshold, and the
flow through the root. Recovery of a quantity X is

    R_X(t) = (X(t) - X_post) / (X_pre - X_post)

so R = 0 just after the injury and R = 1 at the value before it.

Usage:  python scripts/24_injury.py        run (cached) and report
"""
import json

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow

OUT = ROOT / "results" / "injury"
OUT.mkdir(parents=True, exist_ok=True)
S, NU, CLOSE, END, FRAC = 0.6, 8.0, 140, 220, 0.30
SEEDS = (300, 301, 302)
REG = {"liver": dict(P0=2600.0, rho=0.10, delta=0.10, q_set=0.6, renewal=1.0),
       "kidney": dict(P0=300.0, rho=0.00, delta=0.05, q_set=0.6, renewal=0.0)}
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=6000)
KINDS = ("control", "parenchyma", "vascular", "combined")
# Recovery needs two things the growth model does not have on its own. Injury
# must restore a unit's ability to enlarge, which growth factors do and the
# decaying growth rate forbids; and the circulation must cap how far it goes,
# which a fixed inlet pressure does not, since the tree can grow its
# conductance and draw more. Each condition is run with and without.
CONDS = {"as grown": dict(),
         "reactivated": dict(injury_reactivates=1.0),
         "reactivated, fixed inflow": dict(injury_reactivates=1.0, inflow_limited=True)}
KEYS = ("N", "mass", "perfused", "Q")


def run(reg, kind, seed, cond="as grown", frac=FRAC):
    tag = f"{reg}_{kind}_{cond.replace(' ', '').replace(',', '')}_f{frac}_s{seed}"
    f = OUT / f"{tag}.json"
    if f.exists():
        return json.load(open(f))
    kw = {} if kind == "control" else dict(injure_day=CLOSE, injury=kind,
                                           injury_fraction=frac, **CONDS[cond])
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S), 190 * S,
                        days=END, rng=np.random.default_rng(seed), nu=NU,
                        **FIX, **REG[reg], **kw)
    series = {k: [float(x) for x in getattr(tr, k)] for k in KEYS}
    r = dict(regime=reg, kind=kind, cond=cond, fraction=frac, seed=seed,
             series=series, injury=tr.injury, closure=tr.closure)
    json.dump(r, open(f, "w"))
    return r


def recovery(rows, key):
    """Recovery of `key` at the end, averaged over seeds, with its spread.

    Defined only where the injury actually lowered the quantity, by more than
    2 % of its value before: an injury that does not touch a quantity has
    nothing to recover, and the ratio would divide by noise.
    """
    out = []
    for r in rows:
        s = r["series"][key]
        pre, post, end = s[CLOSE - 2], s[CLOSE - 1], s[-1]
        out.append((end - post) / (pre - post) if pre - post > 0.02 * abs(pre) else np.nan)
    if all(np.isnan(out)):
        return np.nan, np.nan
    return float(np.nanmean(out)), float(np.nanstd(out, ddof=1) / np.sqrt(len(out)))


def versus_control(rows, ctrl, key):
    """Where the quantity ends, against the same run with no injury."""
    v = [r["series"][key][-1] for r in rows]
    c = [r["series"][key][-1] for r in ctrl]
    return float(np.mean(v) / max(np.mean(c), 1e-9))


def levels(rows, key):
    """Where a quantity ends, relative to the day before the injury."""
    v = [r["series"][key][-1] / max(r["series"][key][CLOSE - 2], 1e-9) for r in rows]
    return float(np.mean(v)), float(np.std(v, ddof=1) / np.sqrt(len(v)))


if __name__ == "__main__":
    res = {}
    for reg in REG:
        print(f"\n=== {reg}: {int(100*FRAC)} % lost on day {CLOSE}, read on day {END}; "
              f"all relative to the day before ===")
        print(f"{'condition':>26} {'units':>14} {'mass':>14} {'mass per unit':>14} "
              f"{'perfused mass':>14}")
        for cond in CONDS:
            rows = [run(reg, "parenchyma", s, cond) for s in SEEDS]
            n, _ = levels(rows, "N"); mm, _ = levels(rows, "mass")
            pf, _ = levels(rows, "perfused")
            per = mm / max(n, 1e-9)
            print(f"{cond:>26} {n:>14.2f} {mm:>14.2f} {per:>14.2f} {pf:>14.2f}")
            res[f"{reg}_{cond}"] = dict(units=n, mass=mm, mass_per_unit=per, perfused=pf)
    # how much of the loss the kidney makes up by enlarging what it keeps
    print("\n=== kidney, reactivated and inflow limited: compensation against loss ===")
    print(f"{'lost':>8} {'nephrons kept':>15} {'mass per nephron':>18} "
          f"{'full compensation':>19} {'mass':>8}")
    for frac in (0.3, 0.5):
        rows = [run("kidney", "parenchyma", s, "reactivated, fixed inflow", frac)
                for s in SEEDS]
        n, _ = levels(rows, "N"); mm, _ = levels(rows, "mass")
        print(f"{100*frac:>7.0f} % {n:>15.2f} {mm / max(n, 1e-9):>18.2f} "
              f"{1/(1-frac):>19.2f} {mm:>8.2f}")
        res[f"kidney_compensation_{frac}"] = dict(units=n, mass=mm,
                                                  mass_per_unit=mm / max(n, 1e-9),
                                                  full=1 / (1 - frac))
    json.dump(res, open(ROOT / "results" / "injury.json", "w"), indent=1)
