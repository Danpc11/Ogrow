#!/usr/bin/env python
"""Does the model recover from injury, and does function follow structure?

The organ is grown to closure (day 140), then a part of it is destroyed and
growth runs on for 80 more days. Injuring after closure, rather than during
growth, is what separates recovery from growth that would have happened
anyway; every quantity is read against the day before the injury in the same
run, so no separate control is needed.

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
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow

OUT = ROOT / "results" / "injury"
OUT.mkdir(parents=True, exist_ok=True)
S, NU, CLOSE, END, FRAC = 0.6, 8.0, 140, 220, 0.30
SEEDS = tuple(range(300, 324))     # 24, as elsewhere in this work
REG = {"liver": dict(P0=2600.0, rho=0.10, delta=0.10, q_set=0.6, renewal=1.0),
       "kidney": dict(P0=300.0, rho=0.00, delta=0.05, q_set=0.6, renewal=0.0)}
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=6000)
# Recovery needs two things the growth model does not have on its own. Injury
# must restore a unit's ability to enlarge, which growth factors do and the
# decaying growth rate forbids; and the circulation must cap how far it goes,
# which a fixed inlet pressure does not, since the tree can grow its
# conductance and draw more. Each condition is run with and without.
SERIES = {"N": "N", "mass": "mass"}     # injury record key -> series key

CONDS = {"as grown": dict(),
         "reactivated": dict(injury_reactivates=1.0),
         "reactivated, fixed inflow": dict(injury_reactivates=1.0, inflow_limited=True)}
KEYS = ("N", "mass", "perfused", "Q")


def run(reg, kind, seed, cond="as grown", frac=FRAC):
    tag = f"{reg}_{kind}_{cond.replace(' ', '').replace(',', '')}_f{frac}_s{seed}"
    f = OUT / f"{tag}.json"
    if f.exists():
        return json.load(open(f))
    kw = dict(injure_day=CLOSE, injury=kind, injury_fraction=frac, **CONDS[cond])
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S), 190 * S,
                        days=END, rng=np.random.default_rng(seed), nu=NU,
                        **FIX, **REG[reg], **kw)
    series = {k: [float(x) for x in getattr(tr, k)] for k in KEYS}
    r = dict(regime=reg, kind=kind, cond=cond, fraction=frac, seed=seed,
             series=series, injury=tr.injury, closure=tr.closure)
    json.dump(r, open(f, "w"))
    return r


def levels(rows, key):
    """Where a quantity ends in each seed, relative to the day before the injury.

    The seeds are returned rather than summarised: a mean alone hides how wide
    the draw is, and the spread turns out to be the size of the difference
    between two machines running the same seed, so a reader needs to see it to
    judge the digits.
    """
    return np.array([r["series"][key][-1] / max(r["series"][key][CLOSE - 2], 1e-9)
                     for r in rows], float)


def made_up(rows, key):
    """The share of what the injury took that is back by the end, per seed.

    This is the quantity that separates rebuilding from restoring, so it is
    computed here rather than left to the prose. It is defined only where the
    injury actually took something: below a fiftieth of the pre-injury value
    the denominator is noise, and the seed is dropped.
    """
    out = []
    for r in rows:
        pre = r["injury"]["pre"][key]
        post = r["injury"]["post"][key]
        if pre - post > 0.02 * pre:
            out.append((r["series"][SERIES[key]][-1] - post) / (pre - post))
    return np.array(out, float)


def spread(v):
    """Mean and per-seed values, for the table and for the record."""
    return dict(mean=float(np.mean(v)), seeds=[float(x) for x in v])


def _warm(a):
    """Populate the cache for one run. Each tag is written by one worker only."""
    run(*a)


if __name__ == "__main__":
    grid = [(reg, "parenchyma", s, cond, FRAC)
            for reg in REG for cond in CONDS for s in SEEDS]
    grid += [("kidney", "parenchyma", s, "reactivated, fixed inflow", frac)
             for frac in (0.3, 0.5) for s in SEEDS]
    grid = [g for g in grid
            if not (OUT / (f"{g[0]}_{g[1]}_{g[3].replace(' ', '').replace(',', '')}"
                           f"_f{g[4]}_s{g[2]}.json")).exists()]
    if grid:
        n = max(1, min(int(os.environ.get("OGROW_JOBS", 0)) or (os.cpu_count() or 1),
                       len(grid)))
        print(f"running {len(grid)} simulations on {n} worker(s)")
        if n > 1:
            with ProcessPoolExecutor(max_workers=n) as pool:
                for _ in pool.map(_warm, grid):
                    pass
        else:
            for g in grid:
                _warm(g)

    res = {}
    for reg in REG:
        print(f"\n=== {reg}: {int(100*FRAC)} % lost on day {CLOSE}, read on day {END}; "
              f"all relative to the day before ===")
        print(f"{'condition':>26} {'units':>14} {'mass':>14} {'mass per unit':>14} "
              f"{'perfused mass':>14}")
        for cond in CONDS:
            rows = [run(reg, "parenchyma", s, cond) for s in SEEDS]
            n, mm = levels(rows, "N"), levels(rows, "mass")
            pf = levels(rows, "perfused")
            per = mm / np.maximum(n, 1e-9)      # per seed, not a ratio of means
            print(f"{cond:>26} {n.mean():>14.2f} {mm.mean():>14.2f} "
                  f"{per.mean():>14.2f} {pf.mean():>14.2f}")
            res[f"{reg}_{cond}"] = dict(units=spread(n), mass=spread(mm),
                                        mass_per_unit=spread(per), perfused=spread(pf),
                                        units_made_up=spread(made_up(rows, "N")),
                                        mass_made_up=spread(made_up(rows, "mass")))
    # how much of the loss the kidney makes up by enlarging what it keeps
    print("\n=== kidney, reactivated and inflow limited: compensation against loss ===")
    print(f"{'lost':>8} {'nephrons kept':>15} {'mass per nephron':>18} "
          f"{'full compensation':>19} {'mass':>8}")
    for frac in (0.3, 0.5):
        rows = [run("kidney", "parenchyma", s, "reactivated, fixed inflow", frac)
                for s in SEEDS]
        n, mm = levels(rows, "N"), levels(rows, "mass")
        per = mm / np.maximum(n, 1e-9)
        print(f"{100*frac:>7.0f} % {n.mean():>15.2f} {per.mean():>18.2f} "
              f"{1/(1-frac):>19.2f} {mm.mean():>8.2f}   seeds "
              f"{' '.join(f'{x:.2f}' for x in per)}")
        res[f"kidney_compensation_{frac}"] = dict(units=spread(n), mass=spread(mm),
                                                  mass_per_unit=spread(per),
                                                  full=1 / (1 - frac))
    json.dump(res, open(ROOT / "results" / "injury.json", "w"), indent=1)
