#!/usr/bin/env python
"""How the closure mode and the unit count change with organ size.

Unit size is conserved far better than organ size across mammals, so the
minimum spacing lam stays fixed while the domain radius grows: a larger organ
holds more units of the same size, and the number that fit is V / lam^3.

Two sweeps. "geometry" changes the domain alone. "allometric" also carries the
standard scalings with body mass M: the competent pool with organ volume
(S^3), the developmental window as M^(1/4) (S^(3/4)) and cellular rates as
M^(-1/4) (S^(-3/4)), so that the product of rate and window is scale
invariant. The readouts are the unit count against domain volume, the share of
the pool held by each bound, and how much of the domain's close packing the
units occupy.

Usage:  python scripts/23_scaling.py          run (cached) and report
"""
import json

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow

OUT = ROOT / "results" / "scaling"
OUT.mkdir(parents=True, exist_ok=True)
LAM, R0, S0, D0, NU0, P00 = 7.4, 190.0, 0.6, 140, 8.0, 2600.0
SCALES, SEEDS = (0.4, 0.6, 0.8, 1.0), (300, 301)
LIVER = dict(rho=0.10, delta=0.10, q_set=0.6, renewal=1.0)
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=30000)
PACK = 0.74                      # close packing of equal spheres


def run(mode, S, seed):
    f = OUT / f"{mode}_S{S}_s{seed}.json"
    if f.exists():
        return json.load(open(f))
    g = S / S0
    days = D0 if mode == "geometry" else int(round(D0 * g ** 0.75))
    nu = NU0 if mode == "geometry" else NU0 * g ** -0.75
    P0 = P00 if mode == "geometry" else P00 * g ** 3
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S), R0 * S,
                        days=days, rng=np.random.default_rng(seed), nu=nu, P0=P0,
                        **FIX, **LIVER)
    ceiling = (R0 * S / LAM) ** 3 * PACK
    r = dict(mode=mode, S=S, seed=seed, days=days, nu=nu, P0=P0, N=len(U),
             vessels=len(vt), ceiling=ceiling, fill=len(U) / ceiling, **tr.closure)
    json.dump(r, open(f, "w"))
    return r


if __name__ == "__main__":
    res = {}
    print(f"{'sweep':>11} {'S':>4} {'days':>5} {'N':>7} {'supply':>7} {'room':>6} "
          f"{'spent':>6} {'fill':>7}")
    for mode in ("geometry", "allometric"):
        rows = [run(mode, S, s) for S in SCALES for s in SEEDS]
        for S in SCALES:
            v = [x for x in rows if x["S"] == S]
            m = lambda k: float(np.mean([x[k] for x in v]))
            print(f"{mode:>11} {S:>4} {v[0]['days']:>5} {m('N'):>7.1f} {m('supply_held'):>7.2f} "
                  f"{m('space_held'):>6.2f} {m('spent'):>6.2f} {100*m('fill'):>6.2f} %")
        V = np.log([x["S"] ** 3 for x in rows]); N = np.log([x["N"] for x in rows])
        b = float(np.polyfit(V, N, 1)[0])
        print(f"{'':>11} unit count scales as V^{b:.2f} "
              f"(V^1 is what a conserved unit size requires)\n")
        res[mode] = dict(exponent=b, rows=rows)
    json.dump(res, open(ROOT / "results" / "scaling.json", "w"), indent=1)
