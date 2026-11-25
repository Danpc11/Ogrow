#!/usr/bin/env python
"""Why a budget cannot be spent to the last cell, and what the two repairs cost.

A progenitor is tied to the tip it sits at. Where supply there falls short it
strands, so a pool that is a budget does not empty by being used. The claim in
KNOWN_LIMITATIONS.md is that this is structural rather than a matter of
calibration, and the three sweeps here are what it rests on.

    budget    the kidney regime over initial pool and loss rate. With no loss
              the pool always strands, and a larger pool strands more.
    floor     a residual drive below the supply threshold, so a stranded
              progenitor is not held for ever. It empties a budget and wrecks
              the renewing organ.
    elasticity  whether the ordering of the three organs depends on the loss
              term, measured on the kidney with and without it.

This is a diagnostic, not part of the manuscript pipeline, so `make paper`
does not run it.

Usage:  python scripts/27_pool_fate.py [budget|floor|elasticity ...]
"""
import sys

import numpy as np

from organgrow import sdf
from organgrow.closed import grow

S, DAYS, NU, H = 0.6, 140, 8.0, 0.3
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=6000)
# The three published regimes, as the other scripts define them. The kidney's
# loss rate read 0.0 here, so job_floor measured a kidney that is not the
# paper's; the sweeps that vary the loss rate pass their own and are
# unaffected. check_consistency.py now holds the canonical table and fails on
# any disagreement between the scripts that copy it.
BASE = {"liver": dict(P0=2600.0, rho=0.10, delta=0.10, renewal=1.0),
        "kidney": dict(P0=300.0, rho=0.00, delta=0.05, renewal=0.0),
        "pancreas": dict(P0=60.0, rho=0.00, delta=0.00, renewal=0.0)}
SEEDS = (300, 301, 302)


def run(seed, **p):
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S), 190 * S,
                        days=DAYS, rng=np.random.default_rng(seed), nu=NU,
                        q_set=0.6, **FIX, **p)
    return tr


def job_budget():
    print("\nthe kidney regime over initial pool and loss rate, three seeds\n")
    print(f"{'P0':>5} {'delta':>6} {'N':>6} {'consumed':>9} {'lost':>7} {'held':>7} "
          f"{'by supply':>10}")
    for P0 in (80, 120, 200, 300):
        for delta in (0.0, 0.005, 0.02, 0.05):
            r = [run(s, P0=float(P0), rho=0.0, delta=delta, renewal=0.0) for s in SEEDS]
            f = [t.fate for t in r]
            print(f"{P0:5d} {delta:6.3f} {np.mean([t.N[-1] for t in r]):6.0f} "
                  f"{np.mean([x['consumed'] for x in f]):8.0%} "
                  f"{np.mean([x['lost'] for x in f]):6.0%} "
                  f"{np.mean([x['left'] for x in f]):6.0%} "
                  f"{np.mean([t.closure['supply_held'] for t in r]):10.2f}")
    print("with no loss the pool always strands, and a larger pool strands more")


def job_floor():
    print("\na residual drive below the supply threshold, three seeds\n")
    print(f"{'regime':9s} {'q_floor':>8} {'N':>6} {'consumed':>9} {'held':>7} "
          f"{'perfused share':>15}")
    for reg, p in BASE.items():
        for fl in (0.0, 0.01, 0.03, 0.08):
            r = [run(s, q_floor=fl, **p) for s in SEEDS]
            f = [t.fate for t in r]
            print(f"{reg:9s} {fl:8.2f} {np.mean([t.N[-1] for t in r]):6.0f} "
                  f"{np.mean([x['consumed'] for x in f]):8.0%} "
                  f"{np.mean([x['left'] for x in f]):6.0%} "
                  f"{np.mean([t.perfused[-1] / max(t.mass[-1], 1e-9) for t in r]):15.2f}")
    print("it empties a budget and fills the renewing organ with units the"
          "\ncirculation does not reach")


def job_elasticity():
    print("\nelasticity of the kidney's final count to its initial pool, "
          "six paired replicates\n")
    for delta in (0.05, 0.0):
        e = []
        for seed in range(300, 306):
            hi = run(seed, P0=300 * np.exp(H), rho=0.0, delta=delta, renewal=0.0).N[-1]
            lo = run(seed, P0=300 * np.exp(-H), rho=0.0, delta=delta, renewal=0.0).N[-1]
            e.append((np.log(hi) - np.log(lo)) / (2 * H))
        e = np.array(e)
        print(f"  delta = {delta:.2f}: {e.mean():+.2f} ± "
              f"{e.std(ddof=1) / np.sqrt(len(e)):.2f}   seeds {np.round(e, 2)}")
    print("the ordering of the three organs does not rest on the loss term")


JOBS = dict(budget=job_budget, floor=job_floor, elasticity=job_elasticity)

if __name__ == "__main__":
    for name in (sys.argv[1:] or list(JOBS)):
        JOBS[name]()
