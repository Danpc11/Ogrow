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
import sys

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow

OUT = ROOT / "results" / "elasticity"
OUT.mkdir(parents=True, exist_ok=True)
S, DAYS, SEED0, H = 0.6, 140, 300, 0.3
BASE = {
    "liver":    dict(P0=2600.0, rho=0.10, delta=0.10, q_set=0.6, renewal=1.0, nu=8.0),
    "kidney":   dict(P0=300.0,  rho=0.00, delta=0.05, q_set=0.6, renewal=0.0, nu=8.0),
    "pancreas": dict(P0=60.0,   rho=0.00, delta=0.00, q_set=0.6, renewal=0.0, nu=8.0),
}
SETTINGS = {"ref": {}, "nu_lo": dict(nu=4.0), "nu_hi": dict(nu=16.0),
            "q_lo": dict(q_set=0.5), "q_hi": dict(q_set=0.7)}
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=6000)


def run(reg, setting, k, inp, sign):
    f = OUT / f"{reg}_{setting}_{inp}{'+' if sign > 0 else '-'}_k{k}.npz"
    if f.exists():
        return
    p = dict(BASE[reg]); p.update(SETTINGS[setting])
    p[inp] = p[inp] * float(np.exp(sign * H))
    rng = np.random.default_rng(SEED0 + k)
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S),
                        190 * S, days=DAYS, rng=rng, **FIX, **p)
    fin = tr.final
    np.savez(f, N=np.asarray(tr.N, float), C=np.asarray(tr.C, float), Q=np.asarray(tr.Q, float),
             q=np.asarray(tr.q, float), P=np.asarray(tr.P, float),
             q_tip=fin["q"][fin["tips"]], P_tip=fin["P"][fin["tips"]],
             q_set=p["q_set"])


def elasticities(reg, setting):
    out = {}
    for inp in ("P0", "q_set"):
        e = []
        for k in range(64):
            a = OUT / f"{reg}_{setting}_{inp}+_k{k}.npz"
            b = OUT / f"{reg}_{setting}_{inp}-_k{k}.npz"
            if not (a.exists() and b.exists()):
                continue
            e.append((np.log(np.load(a)["N"][-1]) - np.log(np.load(b)["N"][-1])) / (2 * H))
        if e:
            e = np.array(e)
            out[inp] = dict(mean=float(e.mean()), se=float(e.std(ddof=1) / np.sqrt(len(e)))
                            if len(e) > 1 else float("nan"), n=len(e))
    return out


if __name__ == "__main__":
    if len(sys.argv) == 5:
        reg, setting, k0, k1 = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
        for k in range(k0, k1):
            for inp in ("P0", "q_set"):
                for sign in (+1, -1):
                    run(reg, setting, k, inp, sign)
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
