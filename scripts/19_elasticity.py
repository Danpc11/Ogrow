#!/usr/bin/env python
"""Elasticity of the final unit count to the pool and to the supply setpoint.

For each replicate the same seed is run with one input multiplied by exp(+h)
and by exp(-h); the elasticity is (ln N+ - ln N-) / 2h. Common random numbers
remove most of the demographic noise that a comparison of spreads between
cohorts has to average away, and the elasticity is what a mean-field
derivation predicts. For small spreads it equals the transmission of
scripts/15_transmission.py.

Usage:
    python scripts/19_elasticity.py REGIME SETTING K0 K1    run replicates
    python scripts/19_elasticity.py                         report
"""
import json
import os
from concurrent.futures import ProcessPoolExecutor
import sys

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow
from organgrow.runcache import fingerprint, stale, save_npz

OUT = ROOT / "results" / "elasticity"
OUT.mkdir(parents=True, exist_ok=True)
S, DAYS, SEED0, H = 0.6, 140, 300, 0.3
# How many replicates the report reads. 24 is what the rest of this work uses;
# the loop used to read a hard-coded 64, which silently ignored anything past
# it and disagreed with what the Makefile writes.
NSEED = int(os.environ.get("OGROW_SEEDS", 24))
BASE = {
    "liver":    dict(P0=2600.0, rho=0.10, delta=0.10, q_set=0.6, renewal=1.0, nu=8.0),
    "kidney":   dict(P0=300.0,  rho=0.00, delta=0.05, q_set=0.6, renewal=0.0, nu=8.0),
    "pancreas": dict(P0=60.0,   rho=0.00, delta=0.00, q_set=0.6, renewal=0.0, nu=8.0),
}
SETTINGS = {"ref": {}, "nu_lo": dict(nu=4.0), "nu_hi": dict(nu=16.0),
            "q_lo": dict(q_set=0.5), "q_hi": dict(q_set=0.7)}
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=6000)


def spec(reg, setting, k, inp, sign):
    """The file, the growth parameters and the fingerprint of one perturbation."""
    p = dict(BASE[reg]); p.update(SETTINGS[setting])
    p[inp] = p[inp] * float(np.exp(sign * H))
    f = OUT / f"{reg}_{setting}_{inp}{'+' if sign > 0 else '-'}_k{k}.npz"
    return f, p, fingerprint(seed=SEED0 + k, S=S, days=DAYS, params={**FIX, **p})


def run(reg, setting, k, inp, sign):
    f, p, fp = spec(reg, setting, k, inp, sign)
    if not stale(f, fp):
        return
    rng = np.random.default_rng(SEED0 + k)
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S),
                        190 * S, days=DAYS, rng=rng, **FIX, **p)
    fin = tr.final
    save_npz(f, fp, N=np.asarray(tr.N, float), C=np.asarray(tr.C, float),
             Q=np.asarray(tr.Q, float), q=np.asarray(tr.q, float),
             P=np.asarray(tr.P, float),
             q_tip=fin["q"][fin["tips"]], P_tip=fin["P"][fin["tips"]],
             q_set=p["q_set"])


def _warm(a):
    """Populate the cache for one run; each file is written by one worker only."""
    run(*a)


def elasticities(reg, setting):
    out = {}
    for inp in ("P0", "q_set"):
        e = []
        for k in range(NSEED):
            a, _, fa = spec(reg, setting, k, inp, +1)
            b, _, fb = spec(reg, setting, k, inp, -1)
            if stale(a, fa) or stale(b, fb):
                continue
            with np.load(a) as za, np.load(b) as zb:
                e.append((np.log(za["N"][-1]) - np.log(zb["N"][-1])) / (2 * H))
        if e:
            e = np.array(e)
            out[inp] = dict(mean=float(e.mean()), se=float(e.std(ddof=1) / np.sqrt(len(e)))
                            if len(e) > 1 else float("nan"), n=len(e))
    return out


if __name__ == "__main__":
    if len(sys.argv) == 5:
        reg, setting, k0, k1 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
        grid = [(reg, setting, k, inp, sign) for k in range(k0, k1)
                for inp in ("P0", "q_set") for sign in (+1, -1)]
        n = max(1, min(int(os.environ.get("OGROW_JOBS", 0)) or (os.cpu_count() or 1),
                       len(grid)))
        if n > 1:
            with ProcessPoolExecutor(max_workers=n) as pool:
                for _ in pool.map(_warm, grid):
                    pass
        else:
            for g in grid:
                _warm(g)
        print(f"{reg} {setting}: replicates {k0} to {k1 - 1} done")
        raise SystemExit(0)
    res = {}
    print(f"{'regime':>9} {'setting':>8} {'eps_pool':>16} {'eps_setpoint':>16}")
    for reg in BASE:
        for st in SETTINGS:
            e = elasticities(reg, st)
            if not e:
                continue
            res[f"{reg}_{st}"] = e
            fmt = lambda d: f"{d['mean']:+.2f} ± {d['se']:.2f} (n={d['n']})" if d else "-"
            print(f"{reg:>9} {st:>8} {fmt(e.get('P0')):>16} {fmt(e.get('q_set')):>16}")
    json.dump(res, open(ROOT / "results" / "elasticity.json", "w"), indent=1)
