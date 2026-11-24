#!/usr/bin/env python
"""Is the scaling exponent a property of the model, or of one nucleation rate?

Unit count rises with domain volume as a power well below 1, so unit density
falls as the organ grows (`scripts/23_scaling.py`). The reading in THEORY.md
is that what binds a larger organ is how fast the vascular tree spreads and
not how much room its container has: a unit is placed within a fixed distance
of a tip, so progenitors that come up behind the front find no tip to sit at.

If that reading is right the exponent is not a constant of the model. It is
set by the nucleation rate against the speed of the front, and lowering nu
should let the tree get ahead and carry the exponent towards 1. This sweep
asks that question directly: the geometry sweep of 23_scaling.py, run at
several nu, with the exponent fitted at each.

A rising exponent would say that the rate at which units are seeded should be
matched to the speed of vascularisation rather than maximised. A flat one
would say the loss of density is structural and the reading above is wrong.
Either answer is worth having; the second is worth having more.

This is a diagnostic, not part of the manuscript pipeline, so `make paper`
does not run it. OGROW_SEEDS sets the replicates (default 3, for scoping;
use 24 to match the rest of this work) and OGROW_JOBS the workers.

Usage:  python scripts/25_front.py
"""
import json
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow

OUT = ROOT / "results" / "front"
OUT.mkdir(parents=True, exist_ok=True)
LAM, R0, D0 = 7.4, 190.0, 140
SCALES = (0.4, 0.6, 0.8, 1.0)
NUS = (2.0, 4.0, 8.0, 16.0)       # 8 is the reference rate used everywhere else
NSEED = int(os.environ.get("OGROW_SEEDS", 3))
SEEDS = tuple(range(300, 300 + NSEED))
LIVER = dict(rho=0.10, delta=0.10, q_set=0.6, renewal=1.0, P0=2600.0)
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=30000)
PACK = 0.74                       # close packing of equal spheres


def run(nu, S, seed):
    """One liver at one domain scale and one nucleation rate. Cached."""
    f = OUT / f"nu{nu}_S{S}_s{seed}.json"
    if f.exists():
        return json.load(open(f))
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S), R0 * S,
                        days=D0, rng=np.random.default_rng(seed), nu=nu,
                        **FIX, **LIVER)
    ceiling = (R0 * S / LAM) ** 3 * PACK
    r = dict(nu=nu, S=S, seed=seed, N=len(U), vessels=len(vt),
             fill=len(U) / ceiling, **tr.closure)
    json.dump(r, open(f, "w"))
    return r


def _warm(a):
    """Populate the cache for one run; each file is written by one worker only."""
    run(*a)


if __name__ == "__main__":
    grid = [g for g in [(nu, S, sd) for nu in NUS for S in SCALES for sd in SEEDS]
            if not (OUT / f"nu{g[0]}_S{g[1]}_s{g[2]}.json").exists()]
    if grid:
        n = max(1, min(int(os.environ.get("OGROW_JOBS", 0)) or (os.cpu_count() or 1),
                       len(grid)))
        print(f"running {len(grid)} simulations on {n} worker(s)", flush=True)
        if n > 1:
            with ProcessPoolExecutor(max_workers=n) as pool:
                for _ in pool.map(_warm, grid):
                    pass
        else:
            for g in grid:
                _warm(g)

    res = {}
    print(f"\n{NSEED} seeds per point, liver regime, domain scaled alone\n")
    print(f"{'nu':>5} {'S':>5} {'N':>16} {'supply':>7} {'room':>6} {'fill':>7}")
    for nu in NUS:
        rows = [run(nu, S, s) for S in SCALES for s in SEEDS]
        for S in SCALES:
            v = [x for x in rows if x["S"] == S]
            N = np.array([x["N"] for x in v], float)
            mm = (lambda k, v=v: float(np.mean([x[k] for x in v])))
            print(f"{nu:>5.1f} {S:>5} {N.mean():>7.0f} [{N.min():>4.0f},{N.max():>4.0f}] "
                  f"{mm('supply_held'):>7.2f} {mm('space_held'):>6.2f} "
                  f"{100 * mm('fill'):>6.2f} %")
        V = np.log([x["S"] ** 3 for x in rows])
        b, _ = np.polyfit(V, np.log([x["N"] for x in rows]), 1)
        res[f"nu{nu}"] = dict(nu=nu, exponent=float(b), seeds=NSEED, rows=rows)
        print(f"{'':>5} exponent V^{b:.2f}\n")

    e = [(res[f'nu{nu}']['nu'], res[f'nu{nu}']['exponent']) for nu in NUS]
    print("nucleation rate against the exponent: "
          + ", ".join(f"nu={a:.0f} -> {b:.2f}" for a, b in e))
    lo, hi = e[0][1], e[-1][1]
    print("a slower rate raises the exponent, so the loss of unit density with"
          "\nsize is set by the rate against the front and not by the geometry"
          if lo > hi + 0.05 else
          "the exponent does not rise as the rate falls, so the loss of unit"
          "\ndensity with size is not a matter of outrunning the front")
    json.dump(res, open(ROOT / "results" / "front.json", "w"), indent=1)
