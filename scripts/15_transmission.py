#!/usr/bin/env python
"""Which individual spread reaches the final unit count, in each regime.

The open engine answered "which closure mode destroys variance" by
construction: exhaustion and limiting pool read nothing, so they conserved the
spread exactly. In the closed model all organs run one law and differ only in
which bound binds first, so the answer can be measured.

Each replicate gets its own progenitor pool, P0 * exp(s_P z1), or its own
supply setpoint, q_set * exp(s_q z2), with the same z for the same replicate
across conditions. The theory predicts a diagonal transmission:

    supply-bound (liver)        q_set spread reaches N, P0 spread does not
    progenitor-bound (kidney)   P0 spread reaches N, q_set spread does not

Transmission is sqrt(CV_N^2 - CV_baseline^2) / input CV, so demographic noise,
present in every condition, is removed.

Usage:
    python scripts/15_transmission.py REGIME CONDITION K0 K1   run replicates
    python scripts/15_transmission.py                          report
"""
import json
import sys

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow

OUT = ROOT / "results" / "transmission"
OUT.mkdir(parents=True, exist_ok=True)
import os
S, DAYS, SEED0 = 0.6, 140, 300
NREP = int(os.environ.get('OG_NREP', 8))
SIG_P, SIG_Q = 0.30, 0.30   # equal input spreads, so the two transmissions are comparable
REGIMES = {
    "liver":  dict(P0=2600.0, rho=0.10, delta=0.10, q_set=0.6),   # supply binds
    "kidney": dict(P0=300.0,  rho=0.00, delta=0.05, q_set=0.6),   # pool binds, with loss
    "pancreas": dict(P0=60.0, rho=0.00, delta=0.00, q_set=0.6),   # a budget, no loss
}
CONDITIONS = ("baseline", "pool", "setpoint")
Z = np.random.default_rng(2024).normal(size=(max(NREP, 64), 2))


def run_rep(regime, cond, k):
    f = OUT / f"{regime}_{cond}_k{k}.json"
    if f.exists():
        return
    p = dict(REGIMES[regime])
    if cond == "pool":
        p["P0"] *= float(np.exp(SIG_P * Z[k, 0]))
    if cond == "setpoint":
        p["q_set"] *= float(np.exp(SIG_Q * Z[k, 1]))
    rng = np.random.default_rng(SEED0 + k)
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S),
                        190 * S, days=DAYS, rng=rng, nu=8.0, R_unit=650.0,
                        units_per_branch=650 / 600.0, max_vessels=6000, **p)
    P = np.asarray(tr.P, float)
    json.dump(dict(N=len(U), q_final=float(tr.q[-1]),
                   pool_final=float(P[-1]), **p), open(f, "w"))


def cv(regime, cond):
    fs = [OUT / f"{regime}_{cond}_k{k}.json" for k in range(NREP)]
    if not all(f.exists() for f in fs):
        return None, None
    rows = [json.load(open(f)) for f in fs]
    N = np.array([r["N"] for r in rows], float)
    return float(N.std(ddof=1) / N.mean()), rows


if __name__ == "__main__":
    if len(sys.argv) == 5:
        reg, cond, k0, k1 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
        for k in range(k0, min(k1, NREP)):
            run_rep(reg, cond, k)
        print(f"{reg} {cond}: replicates {k0} to {min(k1, NREP) - 1} done")
        raise SystemExit(0)
    res = {}
    print(f"{'regime':>8} {'input':>9} {'CV of N':>8} {'transmission':>13} {'bound'}")
    for reg in REGIMES:
        base, _ = cv(reg, "baseline")
        for cond, sig in (("baseline", 0.0), ("pool", SIG_P), ("setpoint", SIG_Q)):
            c, rows = cv(reg, cond)
            if c is None:
                print(f"{reg:>8} {cond:>9}  (incomplete)")
                continue
            tr = (np.sqrt(max(c ** 2 - base ** 2, 0.0)) / sig) if sig else float("nan")
            lo = hi = float("nan")
            if sig:
                _, brows = cv(reg, "baseline")
                Nb = np.array([r["N"] for r in brows], float)
                Nc = np.array([r["N"] for r in rows], float)
                rng = np.random.default_rng(1)
                bs = []
                for _ in range(2000):
                    a = rng.choice(Nc, len(Nc)); b = rng.choice(Nb, len(Nb))
                    ca, cb = a.std(ddof=1) / a.mean(), b.std(ddof=1) / b.mean()
                    bs.append(np.sqrt(max(ca ** 2 - cb ** 2, 0.0)) / sig)
                lo, hi = np.percentile(bs, [2.5, 97.5])
            # what holds an organ's pool is read from its distribution over
            # tips (closed.closure_state, scripts/21_robustness.py), not from
            # a mean supply, so no label is given here
            bound = ""
            res[f"{reg}_{cond}"] = dict(cv=c, transmission=tr, ci95=[float(lo), float(hi)], bound=bound,
                                        N=[r["N"] for r in rows])
            print(f"{reg:>8} {cond:>9} {c:>8.3f} {tr:>13.2f}  [{lo:.2f}, {hi:.2f}]  {bound}")
    json.dump(res, open(ROOT / "results" / "transmission.json", "w"), indent=1)
