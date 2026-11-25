#!/usr/bin/env python
"""What changes when a unit that loses its supply can disappear.

`scripts/32_inflow_law.py` reports that a construct held at a fixed inflow
settles in total while its units do not: the residual sum |m - ceiling| over
sum ceiling stays around a third, and most of the units end at the mass floor.
The reason is a mechanism the reference model does not have. A unit whose share
of the delivery collapses has no way to disappear, so it sits at the floor
instead of involuting, and the construct ends as a few large units among many
vestigial ones.

This is the sensitivity that asks what the result looks like when that
mechanism is supplied. `involution` is a rate per day at which a unit whose
supply has fallen below `involution_q` of the setpoint loses mass; below half
the mass a unit is born with it is removed. The reference model is the rate of
zero, so the first row of each table is the published result.

The rate is not calibrated against any measurement. What makes the sensitivity
worth reporting is that the answer barely moves across a tenfold range of it,
so the conclusion rests on the mechanism existing and not on its speed.

Everything is at unit flow, where the published result sits. The demand phi is
measured once on an intact ensemble of a disjoint block of seeds, exactly as in
32, and does not depend on involution because it is read before the injury.

Usage:  python scripts/33_involution.py [REGIME ...]     default: both

    OGROW_SEEDS         replicates, default 24
    OGROW_JOBS          workers
    OGROW_ADOPT_CACHE   see organgrow.runcache
"""
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow, _p
from organgrow.runcache import fingerprint, cached

OUT = ROOT / "results" / "involution"
OUT.mkdir(parents=True, exist_ok=True)
S, NU, CLOSE, END = 0.6, 8.0, 140, 220
RATES = (0.0, 0.02, 0.05, 0.10, 0.20)
NSEED = int(os.environ.get("OGROW_SEEDS", 24))
SEEDS = tuple(range(300, 300 + NSEED))          # where the sweep is evaluated
CAL = tuple(range(400, 400 + NSEED))            # where phi is calibrated
REG = {"liver": dict(P0=2600.0, rho=0.10, delta=0.10, q_set=0.6, renewal=1.0),
       "kidney": dict(P0=300.0, rho=0.00, delta=0.05, q_set=0.6, renewal=0.0)}
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=6000)
LOSS, SCALE = 0.3, 1.0
M0 = float(_p("m0", None))


def _kw(reg, rate, phi):
    return dict(injure_day=CLOSE, injury="parenchyma", injury_fraction=LOSS,
                injury_reactivates=1.0, inflow_limited=True,
                inflow_scale=SCALE, phi_demand=float(phi), involution=float(rate))


def intact(reg, seed):
    """The intact organ at arrest, which is where phi is measured."""
    f = OUT / f"intact_{reg}_s{seed}.json"
    fp = fingerprint(regime=reg, seed=seed, S=S, NU=NU, days=CLOSE,
                     params={**FIX, **REG[reg]})
    hit = cached(f, fp)
    if hit is not None:
        return hit
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S), 190 * S,
                        days=CLOSE, rng=np.random.default_rng(seed), nu=NU,
                        **FIX, **REG[reg])
    r = dict(regime=reg, seed=seed, mass=float(np.asarray(tr.mass)[-1]),
             Q=float(np.asarray(tr.Q)[-1]), N=len(U), fingerprint=fp)
    r["phi"] = r["Q"] / max(r["mass"], 1e-12)
    json.dump(r, open(f, "w"))
    return r


def run(reg, seed, rate, phi):
    """One construct, injured and then followed with involution at `rate`."""
    f = OUT / f"{reg}_s{seed}_i{rate}.json"
    kw = _kw(reg, rate, phi)
    fp = fingerprint(regime=reg, seed=seed, S=S, NU=NU, days=END,
                     params={**FIX, **REG[reg], **kw})
    hit = cached(f, fp)
    if hit is not None:
        return hit
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S), 190 * S,
                        days=END, rng=np.random.default_rng(seed), nu=NU,
                        **FIX, **REG[reg], **kw)
    pre = tr.injury["pre"]
    M = np.asarray(tr.mass, float)
    c = tr.final.get("ceiling")
    c = None if c is None else np.asarray(c, float)
    r = dict(regime=reg, seed=seed, rate=rate, phi=phi,
             pre_mass=float(pre["mass"]), pre_N=int(pre["N"]),
             pre_Q=float(tr.Q[CLOSE - 1]), mass=float(M[-1]), N=len(U),
             residual=(float("nan") if c is None
                       else float(np.abs(m - c).sum() / max(c.sum(), 1e-12))),
             floor_share=float((m <= 1.05 * M0).mean()) if len(m) else float("nan"),
             drift_per_day=float(abs(M[-1] - M[-11]) / max(M[-1], 1e-12) / 10.0),
             balance=float(tr.fate["balance"]), fingerprint=fp)
    json.dump(r, open(f, "w"))
    return r


def _warm(a):
    (intact(*a[1:]) if a[0] == "intact" else run(*a[1:]))


def jobs(grid):
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
        jobs([("intact", reg, s) for s in CAL
              if cached(OUT / f"intact_{reg}_s{s}.json",
                        fingerprint(regime=reg, seed=s, S=S, NU=NU, days=CLOSE,
                                    params={**FIX, **REG[reg]})) is None])
        ref = [intact(reg, s) for s in CAL]
        phi = float(np.mean([r["phi"] for r in ref]))
        print(f"\n{reg}: a 30 % loss at unit flow, demand phi = {phi:.5f} "
              f"calibrated on seeds {CAL[0]} to {CAL[-1]}, evaluated on seeds "
              f"{SEEDS[0]} to {SEEDS[-1]}")

        jobs([("run", reg, s, x, phi) for x in RATES for s in SEEDS
              if cached(OUT / f"{reg}_s{s}_i{x}.json",
                        fingerprint(regime=reg, seed=s, S=S, NU=NU, days=END,
                                    params={**FIX, **REG[reg],
                                            **_kw(reg, x, phi)})) is None])

        print(f"\n{'involution':>11} {'mass / pre-injury':>22} {'residual':>9} "
              f"{'at floor':>9} {'units / pre-injury':>19} {'balance':>8}")
        rows = []
        for x in RATES:
            v = [run(reg, s, x, phi) for s in SEEDS]
            rel = np.array([r["mass"] / max(r["pre_mass"], 1e-12) for r in v])
            un = np.array([r["N"] / max(r["pre_N"], 1) for r in v])
            eps = np.array([r["residual"] for r in v])
            fl = np.array([r["floor_share"] for r in v])
            print(f"{x:>11.2f} {rel.mean():>14.3f} [{rel.min():.2f}, {rel.max():.2f}] "
                  f"{np.nanmean(eps):>9.3f} {np.nanmean(fl):>9.2f} "
                  f"{un.mean():>13.2f} [{un.min():.2f}, {un.max():.2f}] "
                  f"{max(abs(r['balance']) for r in v):>8.1g}", flush=True)
            rows.append(dict(rate=x, mass_rel=float(rel.mean()),
                             mass_rel_sd=float(rel.std(ddof=1)),
                             mass_rel_lo=float(rel.min()), mass_rel_hi=float(rel.max()),
                             units_rel=float(un.mean()),
                             units_rel_lo=float(un.min()), units_rel_hi=float(un.max()),
                             residual=float(np.nanmean(eps)),
                             floor_share=float(np.nanmean(fl)),
                             drift=float(max(r["drift_per_day"] for r in v))))
        on = [r for r in rows if r["rate"] > 0]
        res[reg] = dict(phi=phi, seeds=len(SEEDS), rows=rows,
                        calibration_seeds=[CAL[0], CAL[-1]],
                        evaluation_seeds=[SEEDS[0], SEEDS[-1]],
                        units_rel_on=[min(r["units_rel"] for r in on),
                                      max(r["units_rel"] for r in on)],
                        floor_on=[min(r["floor_share"] for r in on),
                                  max(r["floor_share"] for r in on)],
                        residual_on=[min(r["residual"] for r in on),
                                     max(r["residual"] for r in on)])
        d = res[reg]
        print(f"  Across a tenfold range of the rate the construct keeps "
              f"{d['units_rel_on'][0]:.2f} to {d['units_rel_on'][1]:.2f} of its "
              f"units\n  against {rows[0]['units_rel']:.2f} with the mechanism "
              f"absent, the share at the mass floor\n  falls from "
              f"{rows[0]['floor_share']:.2f} to between {d['floor_on'][0]:.2f} and "
              f"{d['floor_on'][1]:.2f}, and the residual from "
              f"{rows[0]['residual']:.2f}\n  to between {d['residual_on'][0]:.2f} "
              f"and {d['residual_on'][1]:.2f}. The rate sets how fast, not where "
              f"it lands.", flush=True)
    json.dump(res, open(ROOT / "results" / "involution.json", "w"), indent=1)
