#!/usr/bin/env python
"""Does the model recover from injury, and does function follow structure?

The organ is grown to closure (day 140), then a part of it is destroyed and
growth runs on for 80 more days. Injuring after closure, rather than during
growth, is what separates recovery from growth that would have happened
anyway; the control is the same run with no injury.

Three injuries, each a plane cut rather than a random removal, as a lobectomy
is. "parenchyma" removes the units past the cut and leaves the tree.
"vascular" occludes the vessels past the cut, which stops perfusion downstream
without removing tissue. "combined" does both.

Four readouts, two structural and two functional, all already in the model:
units, total unit mass, the mass whose supply meets the threshold, and the
flow through the root. Recovery of a quantity X is

    R_X(t) = (X(t) - X_post) / (X_pre - X_post)

so R = 0 just after the injury and R = 1 at the value before it.

Usage:  python scripts/24_injury.py        run (cached) and report
"""
import json

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow

OUT = ROOT / "results" / "injury"
OUT.mkdir(parents=True, exist_ok=True)
S, NU, CLOSE, END, FRAC = 0.6, 8.0, 140, 220, 0.30
SEEDS = (300, 301, 302)
REG = {"liver": dict(P0=2600.0, rho=0.10, delta=0.10, q_set=0.6, renewal=1.0),
       "kidney": dict(P0=300.0, rho=0.00, delta=0.05, q_set=0.6, renewal=0.0)}
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=6000)
KINDS = ("control", "parenchyma", "vascular", "combined")
KEYS = ("N", "mass", "perfused", "Q")


def run(reg, kind, seed):
    f = OUT / f"{reg}_{kind}_s{seed}.json"
    if f.exists():
        return json.load(open(f))
    kw = {} if kind == "control" else dict(injure_day=CLOSE, injury=kind,
                                           injury_fraction=FRAC)
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S), 190 * S,
                        days=END, rng=np.random.default_rng(seed), nu=NU,
                        **FIX, **REG[reg], **kw)
    series = {k: [float(x) for x in getattr(tr, k)] for k in KEYS}
    r = dict(regime=reg, kind=kind, seed=seed, series=series,
             injury=tr.injury, closure=tr.closure)
    json.dump(r, open(f, "w"))
    return r


def recovery(rows, key):
    """Recovery of `key` at the end, averaged over seeds, with its spread.

    Defined only where the injury actually lowered the quantity, by more than
    2 % of its value before: an injury that does not touch a quantity has
    nothing to recover, and the ratio would divide by noise.
    """
    out = []
    for r in rows:
        s = r["series"][key]
        pre, post, end = s[CLOSE - 2], s[CLOSE - 1], s[-1]
        out.append((end - post) / (pre - post) if pre - post > 0.02 * abs(pre) else np.nan)
    if all(np.isnan(out)):
        return np.nan, np.nan
    return float(np.nanmean(out)), float(np.nanstd(out, ddof=1) / np.sqrt(len(out)))


def versus_control(rows, ctrl, key):
    """Where the quantity ends, against the same run with no injury."""
    v = [r["series"][key][-1] for r in rows]
    c = [r["series"][key][-1] for r in ctrl]
    return float(np.mean(v) / max(np.mean(c), 1e-9))


if __name__ == "__main__":
    res = {}
    for reg in REG:
        print(f"\n=== {reg}: injured on day {CLOSE}, {int(100*FRAC)} % of the organ, "
              f"read on day {END} ===")
        print(f"{'injury':>11} {'lost':>22} | {'recovery of':>12} "
              f"{'units':>14} {'mass':>14} {'perfused mass':>14} {'root flow':>14}")
        for kind in KINDS:
            rows = [run(reg, kind, s) for s in SEEDS]
            if kind == "control":
                control_rows = rows
                drift = {k: float(np.mean([r["series"][k][-1] / max(r["series"][k][CLOSE-2], 1e-9)
                                           for r in rows])) for k in KEYS}
                print(f"{'control':>11} {'nothing':>22} | {'drift x':>12} "
                      + "".join(f"{drift[k]:>14.3f}" for k in KEYS))
                res[f"{reg}_control"] = drift
                continue
            inj = rows[0]["injury"]
            lost = f"{inj['pre']['N']}->{inj['post']['N']} units"
            cells = {k: recovery(rows, k) for k in KEYS}
            rel = {k: versus_control(rows, control_rows, k) for k in KEYS}
            fmt = lambda k: ("      none" if np.isnan(cells[k][0])
                             else f"{cells[k][0]:>6.2f}±{cells[k][1]:.2f}")
            print(f"{kind:>11} {lost:>22} | {'recovered':>12} "
                  + "".join(f"{fmt(k):>14}" for k in KEYS))
            print(f"{'':>11} {'':>22} | {'of control':>12} "
                  + "".join(f"{rel[k]:>14.2f}" for k in KEYS))
            res[f"{reg}_{kind}"] = {k: dict(recovered=cells[k][0], se=cells[k][1],
                                            of_control=rel[k]) for k in KEYS}
    json.dump(res, open(ROOT / "results" / "injury.json", "w"), indent=1)
