#!/usr/bin/env python
"""Calibrate the closed engine against the measured mouse liver contraction.

Three conditions, asked for together rather than one after another, because
in the open engine they fought each other:

  1. unit number grows 1.82-fold          (80 -> 146 lobules, Azizoglu 2025)
  2. the CV of unit number contracts 3.09-fold over that growth
  3. the window closes inside the run: supply per unit reaches q_set and the
     unit count stops rising

Only three numbers are free: nu, q_set and the length of the run. Everything
else is derived (the setpoint from the vessel-tissue coupling point, the
critical radius from the shear condition, the scales from the
nondimensionalisation) or fixed by the regime (rho = delta for a liver, where
capacity binds before the progenitor pool).

With Poisson nucleation the terminal floor is demographic, which is what the
animal has: counting 146 lobules carries a floor of 1/sqrt(146) = 0.083,
just above the 0.074 reported. If the closed engine reproduces the
contraction without a large parameter spread, the tension the open engine had
between closing and contracting was an artefact of putting all the terminal
variance into sigma_param.
"""
import json
import itertools
import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow

OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

TARGET_GROWTH = 1.82          # N_f / N_i
TARGET_CONTRACTION = 3.09     # CV_i / CV_f at that growth
NREP = 14
KIND = "liver"
R = 190


def run(nu, q_set, days, P0=600.0, nrep=NREP, seed0=300):
    """One cohort: unit-count trajectories for `nrep` replicate organs."""
    traces = []
    for s in range(nrep):
        rng = np.random.default_rng(seed0 + s)
        U, m, vt, tr = grow(lambda x, y, z: sdf(KIND, x, y, z), R,
                            days=days, rng=rng, nu=nu, q_set=q_set,
                            rho=0.1, delta=0.1, P0=P0)
        traces.append((np.asarray(tr.N, float), np.asarray(tr.q, float)))
    N = np.vstack([t[0] for t in traces])
    q = np.vstack([t[1] for t in traces])
    return N, q


def score(N, q, q_set):
    """How far a cohort sits from the three conditions."""
    mean = N.mean(0)
    if mean[0] <= 0 or mean[-1] / mean[0] < TARGET_GROWTH:
        return None
    k = int(np.argmin(np.abs(mean / mean[0] - TARGET_GROWTH)))
    cv_i = N[:, 0].std(ddof=1) / mean[0]
    cv_k = N[:, k].std(ddof=1) / mean[k]
    if cv_k <= 0:
        return None
    contraction = cv_i / cv_k
    # the window has closed if the count has stopped rising and supply per
    # unit has come down to the setpoint
    tail = (mean[-1] - mean[int(0.8 * len(mean))]) / max(mean[-1], 1.0)
    closed = bool(tail < 0.05 and q.mean(0)[-1] <= q_set * 1.1)
    return dict(day_at_growth=k + 1, contraction=float(contraction),
                closed=closed, N_final=float(mean[-1]),
                cv_initial=float(cv_i), cv_at_growth=float(cv_k),
                cv_final=float(N[:, -1].std(ddof=1) / mean[-1]),
                q_final=float(q.mean(0)[-1]),
                error=float(abs(np.log(contraction / TARGET_CONTRACTION))))


def main():
    # P0 is large enough that the progenitor pool does not run out: the
    # liver is the organ where capacity binds first, and a small pool puts it
    # in the kidney's regime instead, which is what an earlier sweep hit.
    grid = dict(nu=[0.5, 1.0, 2.0, 4.0, 8.0],
                q_set=[0.2, 0.4, 0.6, 0.8],
                days=[70])
    rows, best = [], None
    print(f"{'nu':>6} {'q_set':>7} {'days':>6} {'N_f':>7} {'contr':>8} "
          f"{'CV_f':>8} {'closed':>7} {'err':>6}")
    for nu, q_set, days in itertools.product(grid["nu"], grid["q_set"],
                                             grid["days"]):
        N, q = run(nu, q_set, days)
        s = score(N, q, q_set)
        if s is None:
            continue
        s.update(nu=nu, q_set=q_set, days=days)
        rows.append(s)
        print(f"{nu:>6} {q_set:>7} {days:>6} {s['N_final']:>7.0f} "
              f"{s['contraction']:>8.2f} {s['cv_final']:>8.4f} "
              f"{str(s['closed']):>7} {s['error']:>6.2f}")
        # a run that does not close is not a candidate, however well its
        # contraction matches: the window closing is the condition the open
        # engine could never satisfy at the same time
        if s["closed"] and (best is None or s["error"] < best["error"]):
            best = s

    res = dict(target_growth=TARGET_GROWTH,
               target_contraction=TARGET_CONTRACTION,
               poisson_floor_146=float(1 / np.sqrt(146)),
               observed_cv_final=0.074, grid=rows, best=best)
    json.dump(res, open(OUT / "closed_calibration.json", "w"), indent=1)

    print()
    if best is None:
        print("no setting closes the window and matches the contraction")
        print("that is a statement about the model, not about the search")
    else:
        print(f"best: nu = {best['nu']}, q_set = {best['q_set']}, "
              f"days = {best['days']}")
        print(f"  contraction {best['contraction']:.2f} "
              f"against {TARGET_CONTRACTION} measured")
        print(f"  terminal CV {best['cv_final']:.4f} against 0.074 measured, "
              f"Poisson floor {1/np.sqrt(146):.3f}")
    print("wrote results/closed_calibration.json")


if __name__ == "__main__":
    main()
