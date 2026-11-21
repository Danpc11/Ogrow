#!/usr/bin/env python
"""Exploratory, from the development of the closed model; its parameters are
from earlier versions and the paper does not use it (see PAPER.md).

Does unit number self-average? A size sweep that checks its own controls.

Self-averaging is a statement about how the relative spread R = sd/mean of an
observable behaves as the system grows with everything else held fixed:

    R ~ N^(-1/2)   self-averaging, as for a sum of independent parts
    R ~ const      non-self-averaging, as under quenched disorder
    R grows        reinforcement, as in a Polya urn above its threshold

Three earlier attempts failed for the same reason: the knob that changed N
also changed the mechanism. Varying P0 or q_set moves the bound that binds;
varying lam changes units per vessel from 1.3 to 0.4 and pins the tree at its
vessel cap. None of those measure size.

Here the domain is scaled by s and every extensive quantity with it (pool,
design size, vessel cap, all as s^3), while every intensive one (lam, nu,
q_set, b, eta) stays fixed. Then four quantities that should not depend on
size are checked before the exponent is reported: units per vessel, median
radius, final supply per unit, and the fraction of the pool left. If any of
them drifts by more than TOLERANCE across the sweep the script says so and
reports no exponent, because whatever it measured, it was not size.
"""
import json
import time

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow

OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

R_REF = 190.0
import sys
SCALES = ([float(a) for a in sys.argv[1:]] if len(sys.argv) > 1
          else [0.6, 0.8, 1.0, 1.3, 1.7])  # N spans ~ s^3
NREP = 8
DAYS = 140
TOLERANCE = 0.25                          # max relative drift of a control

REF = dict(P0=600.0, max_vessels=1200)
# intensive latents, held fixed across sizes at their reference values
LATENT = dict(R_unit=650.0, units_per_branch=650.0 / 600.0)
FIXED = dict(nu=8.0, q_set=0.6, rho=0.1, delta=0.1)


def run_size(s, nrep=NREP, seed0=300):
    """Replicate organs at domain scale s, extensive quantities scaled as s^3."""
    vol = s ** 3
    kw = dict(FIXED, **LATENT, P0=REF["P0"] * vol,
              max_vessels=int(round(REF["max_vessels"] * vol)))
    R = R_REF * s
    out = []
    for k in range(nrep):
        rng = np.random.default_rng(seed0 + k)
        U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / s, y / s, z / s),
                            R, days=DAYS, rng=rng, **kw)
        r = np.asarray(vt.radius, float)
        P = np.asarray(tr.P, float)
        out.append(dict(N=len(U), vessels=len(vt),
                        units_per_vessel=len(U) / max(len(vt), 1),
                        r_median=float(np.median(r)),
                        q_final=float(tr.q[-1]),
                        pool_left=float(P[-1] / max(P[0], 1e-9)),
                        capped=bool(len(vt) >= kw["max_vessels"])))
    return out


def summarise(rows):
    N = np.array([r["N"] for r in rows], float)
    agg = dict(N_mean=float(N.mean()), N_sd=float(N.std(ddof=1)),
               R=float(N.std(ddof=1) / max(N.mean(), 1e-9)),
               any_capped=any(r["capped"] for r in rows))
    for key in ("units_per_vessel", "r_median", "q_final", "pool_left"):
        agg[key] = float(np.mean([r[key] for r in rows]))
    return agg


def main():
    table = []
    print(f"{'s':>5} {'N':>7} {'R':>8} {'N/ves':>7} {'r_med':>7} "
          f"{'q_f':>7} {'pool':>7} {'cap':>5} {'t':>6}")
    for s in SCALES:
        t0 = time.time()
        a = summarise(run_size(s))
        a["s"] = s
        table.append(a)
        print(f"{s:>5} {a['N_mean']:>7.0f} {a['R']:>8.4f} "
              f"{a['units_per_vessel']:>7.2f} {a['r_median']:>7.3f} "
              f"{a['q_final']:>7.3f} {a['pool_left']:>7.3f} "
              f"{str(a['any_capped']):>5} {time.time()-t0:>5.0f}s")

    # ---- the controls, checked before anything is concluded --------------
    print()
    ok = True
    for key in ("units_per_vessel", "r_median", "q_final", "pool_left"):
        v = np.array([a[key] for a in table])
        drift = float((v.max() - v.min()) / max(abs(v.mean()), 1e-9))
        flag = "ok" if drift <= TOLERANCE else "DRIFTS"
        ok &= drift <= TOLERANCE
        print(f"  control {key:<18} drift {100*drift:5.1f} %   {flag}")
    if any(a["any_capped"] for a in table):
        print("  control vessel cap            reached   DRIFTS")
        ok = False

    res = dict(scales=SCALES, days=DAYS, nrep=NREP, tolerance=TOLERANCE,
               reference=REF, fixed=FIXED, table=table, controls_ok=ok)
    if ok:
        N = np.array([a["N_mean"] for a in table])
        R = np.array([a["R"] for a in table])
        g, cov = np.polyfit(np.log(N), np.log(R), 1, cov=True)
        res["exponent"] = float(g[0])
        res["exponent_se"] = float(np.sqrt(cov[0, 0]))
        print(f"\n  R ~ N^{g[0]:.3f} +/- {np.sqrt(cov[0,0]):.3f} "
              f"over N = {N.min():.0f} to {N.max():.0f}")
        print("  -0.5 self-averaging, 0 quenched, > 0 reinforcement")
    else:
        print("\n  controls drift: no exponent is reported, because the sweep")
        print("  did not change size alone")

    json.dump(res, open(OUT / "size_sweep.json", "w"), indent=1)
    print("wrote results/size_sweep.json")


if __name__ == "__main__":
    main()
