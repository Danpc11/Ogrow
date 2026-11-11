#!/usr/bin/env python
"""Which bound binds first, and what that does to the variance.

In the closed model organ identity is not a branch in the code. One law
governs nucleation at every tip,

    lambda_i = nu * P_i * (q_i - q_set)_+ ,   dP_i/dt = (rho - delta) P_i - k_i

so a growing organ runs into one of two ceilings: it runs out of supply per
unit, or it runs out of progenitors. Which one it meets is set by two
numbers, nu and rho - delta, and that is what this sweep maps.

The prediction the sweep tests: variance in unit number contracts only where
the supply bound binds, because only there does the organ read a quantity
that its own growth lowers. Where the progenitor pool runs out first, the
spread is carried to the end untouched. Liver, kidney and pancreas are then
positions on this plane rather than three separate mechanisms.
"""
import json

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow

OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

KIND = "liver"          # the domain; identity here is (nu, rho - delta)
R = 190
DAYS = 140
NREP = 10
P0 = 600.0
Q_SET = 0.6

NU = [1.0, 2.0, 4.0, 8.0, 16.0]
NET = [0.10, 0.05, 0.0, -0.05, -0.10, -0.20]    # rho - delta


def cohort(nu, net, nrep=NREP, seed0=300):
    """Replicate organs at one point of the plane."""
    N, P, Q = [], [], []
    rho, delta = (net, 0.0) if net >= 0 else (0.0, -net)
    for s in range(nrep):
        rng = np.random.default_rng(seed0 + s)
        U, m, vt, tr = grow(lambda x, y, z: sdf(KIND, x, y, z), R,
                            days=DAYS, rng=rng, nu=nu, q_set=Q_SET,
                            rho=rho, delta=delta, P0=P0)
        N.append(np.asarray(tr.N, float))
        P.append(np.asarray(tr.P, float))
        Q.append(np.asarray(tr.q, float))
    return np.vstack(N), np.vstack(P), np.vstack(Q)


def classify(N, P, Q):
    """Which ceiling the organ met, and what happened to the spread."""
    n, p, q = N.mean(0), P.mean(0), Q.mean(0)
    # the pool is spent if it has fallen to a small fraction of its start
    pool_spent = bool(p[-1] < 0.05 * p[0])
    # supply binds if per-unit supply has come down onto the setpoint
    supply_met = bool(q[-1] <= Q_SET * 1.1)
    if supply_met and not pool_spent:
        bound = "supply"
    elif pool_spent and not supply_met:
        bound = "progenitor"
    elif pool_spent and supply_met:
        bound = "both"
    else:
        bound = "neither"
    sd_i, sd_f = N[:, 0].std(ddof=1), N[:, -1].std(ddof=1)
    return dict(bound=bound, N_final=float(n[-1]),
                sd_initial=float(sd_i), sd_final=float(sd_f),
                sd_ratio=float(sd_f / max(sd_i, 1e-9)),
                cv_final=float(sd_f / max(n[-1], 1e-9)),
                pool_left=float(p[-1] / max(p[0], 1e-9)),
                q_final=float(q[-1]))


def main():
    rows = []
    print(f"{'nu':>6} {'rho-delta':>10} {'N_f':>7} {'bound':>12} "
          f"{'sd_f/sd_i':>10} {'pool left':>10}")
    for nu in NU:
        for net in NET:
            N, P, Q = cohort(nu, net)
            r = classify(N, P, Q)
            r.update(nu=nu, net=net)
            rows.append(r)
            print(f"{nu:>6} {net:>10} {r['N_final']:>7.0f} "
                  f"{r['bound']:>12} {r['sd_ratio']:>10.3f} "
                  f"{r['pool_left']:>10.3f}")

    by = {}
    for r in rows:
        by.setdefault(r["bound"], []).append(r["sd_ratio"])
    print()
    for b, v in sorted(by.items()):
        print(f"  {b:>12}: sd_f/sd_i = {np.mean(v):.3f} "
              f"(n = {len(v)}, range {min(v):.3f} to {max(v):.3f})")

    json.dump(dict(nu=NU, net=NET, days=DAYS, q_set=Q_SET, P0=P0,
                   nrep=NREP, grid=rows),
              open(OUT / "phase_diagram.json", "w"), indent=1)
    print("\nwrote results/phase_diagram.json")


if __name__ == "__main__":
    main()
