#!/usr/bin/env python
"""Is closure irreversible, and which constraint holds the remaining pool?

On day 140, when growth has stopped, the state is changed and growth runs on
for 40 days: "supply" lowers the threshold by 20 %, which raises the hydraulic
availability of every tip relative to it without touching the pool; "room"
shrinks the least spacing between units by 20 %, which frees space without
touching perfusion or the pool; "none" lets growth continue unchanged. The
readout is the share of the pool held on day 140 that becomes units by day
180. If supply holds a share s of the pool and room a share r, relieving
supply should release about s, relieving room about r, and nothing should
release almost nothing.

Usage:  python scripts/22_interventions.py         run (cached) and report
"""
import json

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow

OUT = ROOT / "results" / "interventions"
OUT.mkdir(parents=True, exist_ok=True)
S, NU, DAY0, DAY1, SEEDS = 0.6, 8.0, 140, 180, range(300, 306)
REG = {"liver": dict(P0=2600.0, rho=0.10, delta=0.10, q_set=0.6),
       "kidney": dict(P0=300.0, rho=0.00, delta=0.05, q_set=0.6),
       "pancreas": dict(P0=60.0, rho=0.00, delta=0.00, q_set=0.6)}
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=6000)


def run(reg, cond, seed):
    f = OUT / f"{reg}_{cond}_s{seed}.json"
    if f.exists():
        return json.load(open(f))
    kw = {} if cond == "none" else dict(intervene_day=DAY0 + 1, intervene=cond, intervene_factor=0.8)
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S), 190 * S, days=DAY1,
                        rng=np.random.default_rng(seed), nu=NU, **FIX, **REG[reg], **kw)
    N, P = np.asarray(tr.N, float), np.asarray(tr.P, float)
    r = dict(N0=float(N[DAY0 - 1]), N1=float(N[-1]), P0=float(P[DAY0 - 1]),
             released=float((N[-1] - N[DAY0 - 1]) / max(P[DAY0 - 1], 1e-9)))
    json.dump(r, open(f, "w"))
    return r


if __name__ == "__main__":
    res = {}
    print(f"{'regime':>9} {'intervention':>13} {'pool on day 140':>16} {'share released by day 180':>27}")
    for reg in REG:
        for cond in ("none", "supply", "room"):
            rows = [run(reg, cond, s) for s in SEEDS]
            rel = np.array([r["released"] for r in rows])
            add = np.array([(r["N1"] - r["N0"]) / max(r["N0"], 1e-9) for r in rows])
            res[f"{reg}_{cond}"] = dict(pool=float(np.mean([r["P0"] for r in rows])),
                                        released=float(rel.mean()),
                                        se=float(rel.std(ddof=1) / np.sqrt(len(rel))),
                                        each=[float(x) for x in rel],
                                        units_added=float(add.mean()),
                                        units_added_se=float(add.std(ddof=1) / np.sqrt(len(add))))
            print(f"{reg:>9} {cond:>13} {res[f'{reg}_{cond}']['pool']:>16.0f} "
                  f"{rel.mean():>20.3f} ± {res[f'{reg}_{cond}']['se']:.3f}")
    json.dump(res, open(ROOT / "results" / "interventions.json", "w"), indent=1)
