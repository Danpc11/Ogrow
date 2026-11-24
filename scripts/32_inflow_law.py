#!/usr/bin/env python
"""Does the mass a construct lands on follow the delivery it is given?

`scripts/24_injury.py` reports that holding the root flow at its pre-injury
value brings the mass back to the pre-injury value and stops it. That number
is not a prediction. The ceiling on a unit is

    m_j* = inflow_target * s_j / phi_set,        sum_j s_j = 1,

so the equilibrium total is inflow_target / phi_set, and 24_injury.py sets
phi_set to the intact organ's own delivery per unit mass. With the inflow held
at its pre-injury value that forces the total back to the pre-injury mass by
arithmetic, for every seed, which is also why the spread of the mass ratio
contracts.

This script asks the question that is not circular. It measures the delivery
per unit mass once, on the intact ensemble, and holds that one number fixed
across every seed and every condition. Then it varies the delivery the
circulation provides. If the model is saying something, the mass a construct
lands on should follow

    M_final  ~  Q_delivered / phi_demand

with the pre-injury mass no longer built into either side, and the return to
the pre-injury mass at scale 1 becomes a result rather than a definition.

A construct fed at constant pressure has no such ceiling and is the control.

Usage:  python scripts/32_inflow_law.py [REGIME ...]     default: both

    OGROW_SEEDS   replicates, default 24
    OGROW_JOBS    workers
"""
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow

OUT = ROOT / "results" / "inflow_law"
OUT.mkdir(parents=True, exist_ok=True)
S, NU, CLOSE, END = 0.6, 8.0, 140, 220
SCALES = (0.6, 0.8, 1.0, 1.2, 1.4)
NSEED = int(os.environ.get("OGROW_SEEDS", 24))
SEEDS = tuple(range(300, 300 + NSEED))
REG = {"liver": dict(P0=2600.0, rho=0.10, delta=0.10, q_set=0.6, renewal=1.0),
       "kidney": dict(P0=300.0, rho=0.00, delta=0.05, q_set=0.6, renewal=0.0)}
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=6000)
LOSS = 0.3


def run(reg, seed, scale, phi):
    """One construct: grown, injured on day CLOSE, followed to END."""
    tag = f"{reg}_s{seed}_x{scale}" + ("" if phi is None else "_phi")
    f = OUT / f"{tag}.json"
    if f.exists():
        return json.load(open(f))
    kw = dict(injure_day=CLOSE, injury="parenchyma", injury_fraction=LOSS,
              injury_reactivates=1.0, inflow_limited=True, inflow_scale=scale)
    if phi is not None:
        kw["phi_demand"] = float(phi)
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S), 190 * S,
                        days=END, rng=np.random.default_rng(seed), nu=NU,
                        **FIX, **REG[reg], **kw)
    pre = tr.injury["pre"]
    r = dict(regime=reg, seed=seed, scale=scale, phi=phi,
             pre_mass=float(pre["mass"]), pre_Q=float(tr.Q[CLOSE - 1]),
             mass=float(np.asarray(tr.mass)[-1]), N=len(U),
             Q=float(np.asarray(tr.Q)[-1]))
    json.dump(r, open(f, "w"))
    return r


def intact(reg, seed):
    """The intact organ at arrest, which is where phi_demand is measured."""
    f = OUT / f"intact_{reg}_s{seed}.json"
    if f.exists():
        return json.load(open(f))
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S), 190 * S,
                        days=CLOSE, rng=np.random.default_rng(seed), nu=NU,
                        **FIX, **REG[reg])
    r = dict(regime=reg, seed=seed, mass=float(np.asarray(tr.mass)[-1]),
             Q=float(np.asarray(tr.Q)[-1]), N=len(U))
    r["phi"] = r["Q"] / max(r["mass"], 1e-12)
    json.dump(r, open(f, "w"))
    return r


def _warm(a):
    (intact(*a[1:]) if a[0] == "intact" else run(*a[1:]))


def jobs(n, grid):
    if not grid:
        return
    k = max(1, min(int(os.environ.get("OGROW_JOBS", 0)) or (os.cpu_count() or 1),
                   len(grid)))
    print(f"running {len(grid)} simulations on {k} worker(s)", flush=True)
    if k > 1:
        with ProcessPoolExecutor(max_workers=k) as pool:
            for _ in pool.map(_warm, grid):
                pass
    else:
        for g in grid:
            _warm(g)


if __name__ == "__main__":
    regimes = [a for a in sys.argv[1:] if a in REG] or list(REG)
    res = {}
    for reg in regimes:
        jobs(reg, [("intact", reg, s) for s in SEEDS
                   if not (OUT / f"intact_{reg}_s{s}.json").exists()])
        ref = [intact(reg, s) for s in SEEDS]
        phi = float(np.mean([r["phi"] for r in ref]))
        phi_sd = float(np.std([r["phi"] for r in ref], ddof=1))
        print(f"\n{reg}: delivery per unit mass of the intact ensemble, "
              f"phi = {phi:.5f} ± {phi_sd:.5f} over {len(ref)} seeds")

        jobs(reg, [("run", reg, s, x, phi) for x in SCALES for s in SEEDS
                   if not (OUT / f"{reg}_s{s}_x{x}_phi.json").exists()])

        print(f"\n{'inflow':>8} {'mass / pre-injury':>22} {'mass / (Q/phi)':>18} "
              f"{'units / pre-injury':>19}")
        rows = []
        for x in SCALES:
            v = [run(reg, s, x, phi) for s in SEEDS]
            rel = np.array([r["mass"] / max(r["pre_mass"], 1e-12) for r in v])
            law = np.array([r["mass"] / max(x * r["pre_Q"] / phi, 1e-12) for r in v])
            print(f"{x:>8.1f} {rel.mean():>14.3f} [{rel.min():.2f}, {rel.max():.2f}] "
                  f"{law.mean():>12.3f} [{law.min():.2f}, {law.max():.2f}]")
            rows.append(dict(scale=x, mass_rel=float(rel.mean()),
                             mass_rel_sd=float(rel.std(ddof=1)),
                             mass_over_law=float(law.mean()),
                             mass_over_law_sd=float(law.std(ddof=1)),
                             seeds=[float(a) for a in rel]))
        res[reg] = dict(phi=phi, phi_sd=phi_sd, seeds=len(SEEDS), rows=rows)

        a = np.log([r["scale"] for r in rows])
        b = np.log([r["mass_rel"] for r in rows])
        slope = float(np.polyfit(a, b, 1)[0])
        res[reg]["slope"] = slope
        print(f"\n  mass scales as (delivered flow)^{slope:.2f}")
        print("  An exponent of 1 is what the ceiling gives by construction, so "
              "this is not\n  a test of the budget. What it measures is how "
              "completely a construct\n  reaches the budget in a finite window: "
              "the departure from 1, and the\n  overshoot and shortfall in the "
              "last column, are the result.")
    json.dump(res, open(ROOT / "results" / "inflow_law.json", "w"), indent=1)
