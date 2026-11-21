#!/usr/bin/env python
"""Exploratory, from the development of the closed model; its parameters are
from earlier versions and the paper does not use it (see PAPER.md).

Calibrate the closed model's liver against the measured mouse contraction.

Target (Azizoglu et al. 2025, mouse liver, postnatal): lobule number grows
1.82-fold while its coefficient of variation across animals contracts
3.09-fold, and the window closes.

The organ is grown in the liver regime: the progenitor pool is large enough
not to bind, so the supply bound is what stops it. For each (nu, q_set) a
cohort of replicates is grown. The last day is the closed organ; the earlier
day at which the mean count was 1/1.82 of the final one is the start of the
measured interval; the contraction is the ratio of the two CVs.

Usage:
    python scripts/12_calibrate_closed.py NU QSET K0 K1   replicates K0..K1-1
    python scripts/12_calibrate_closed.py NU QSET         one point, saved
    python scripts/12_calibrate_closed.py                 collect saved points
"""
import json
import sys

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow

OUT = ROOT / "results" / "closed_calibration"
OUT.mkdir(parents=True, exist_ok=True)

GROWTH, CONTRACTION = 1.82, 3.09
S = 0.6                 # domain scale of the calibration organ
P0 = 2600.0             # large enough that the progenitor bound does not bind
DAYS, NREP, SEED0 = 140, 8, 300


REPS = OUT / "reps"
REPS.mkdir(exist_ok=True)


def run_rep(nu, q_set, k):
    """One replicate, saved, so a long point can be built over several calls."""
    f = REPS / f"nu{nu}_q{q_set}_k{k}.npz"
    if f.exists():
        return
    rng = np.random.default_rng(SEED0 + k)
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S),
                        190 * S, days=DAYS, rng=rng, nu=nu, q_set=q_set,
                        rho=0.1, delta=0.1, P0=P0, R_unit=650.0,
                        units_per_branch=650 / 600.0, max_vessels=6000)
    np.savez(f, N=np.asarray(tr.N, float), q=np.asarray(tr.q, float))


def cohort(nu, q_set):
    fs = [REPS / f"nu{nu}_q{q_set}_k{k}.npz" for k in range(NREP)]
    have = [f for f in fs if f.exists()]
    if len(have) < NREP:
        return None
    d = [np.load(f) for f in have]
    return np.vstack([x["N"] for x in d]), np.vstack([x["q"] for x in d])


def score(N, Q, q_set):
    n = N.mean(0)
    k = int(np.argmin(np.abs(n - n[-1] / GROWTH)))
    cv = lambda d: N[:, d].std(ddof=1) / max(N[:, d].mean(), 1e-9)
    tail = (n[-1] - n[int(0.8 * len(n))]) / max(n[-1], 1.0)
    return dict(day_start=k + 1, contraction=float(cv(k) / max(cv(-1), 1e-9)),
                cv_start=float(cv(k)), cv_final=float(cv(-1)),
                N_final=float(n[-1]), q_final=float(Q.mean(0)[-1]),
                closed=bool(tail < 0.05 and Q.mean(0)[-1] <= 1.1 * q_set),
                growth_last_fifth=float(tail),
                N_final_each=[float(x) for x in N[:, -1]])


if __name__ == "__main__":
    if len(sys.argv) == 5:                       # run replicates k0..k1-1
        nu, q_set, k0, k1 = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
        for k in range(k0, min(k1, NREP)):
            run_rep(nu, q_set, k)
        print(f"nu={nu} q_set={q_set}: replicates {k0} to {min(k1, NREP) - 1} done")
    elif len(sys.argv) == 3:
        nu, q_set = float(sys.argv[1]), float(sys.argv[2])
        c = cohort(nu, q_set)
        if c is None:
            for k in range(NREP):
                run_rep(nu, q_set, k)
            c = cohort(nu, q_set)
        r = score(*c, q_set)
        r.update(nu=nu, q_set=q_set)
        json.dump(r, open(OUT / f"nu{nu}_q{q_set}.json", "w"), indent=1)
        print(f"nu={nu:<5} q_set={q_set:<5} N_f={r['N_final']:6.0f} "
              f"contraction={r['contraction']:5.2f} CV_f={r['cv_final']:.4f} "
              f"q_f={r['q_final']:.3f} closed={r['closed']} "
              f"(last fifth +{100*r['growth_last_fifth']:.1f} %)")
    else:
        rows = [json.load(open(f)) for f in sorted(OUT.glob("*.json"))]
        ok = [r for r in rows if r["closed"]]
        for r in sorted(rows, key=lambda r: (r["nu"], r["q_set"])):
            print(f"nu={r['nu']:<5} q_set={r['q_set']:<5} contraction="
                  f"{r['contraction']:5.2f} closed={r['closed']}")
        if ok:
            best = min(ok, key=lambda r: abs(np.log(r["contraction"] / CONTRACTION)))
            print(f"\nbest closed: nu={best['nu']}, q_set={best['q_set']}, "
                  f"contraction {best['contraction']:.2f} against {CONTRACTION}")
        else:
            print("\nno point closes the window")
