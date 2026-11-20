#!/usr/bin/env python
"""Does the separation of regimes depend on the auxiliary rules?

Each variant changes one rule that is not part of the nucleation law, and the
same two readouts are measured in the three regimes:

  closure     where the progenitor pool stands when growth stops: spent, held
              at tips short of supply, held at tips with no room, active
  pool        elasticity of the final unit count to the progenitor pool,
              from the same seed at P0 * exp(+-0.3)

Variants: the vascular root at the deepest point (reference), at the surface,
or at random; each sprouting rule removed in turn (deficit, load, crowding,
mesenchymal signal); units assigned by vascular territory instead of by
nearest tip.

Usage:
    python scripts/21_robustness.py VARIANT [VARIANT ...]   run
    python scripts/21_robustness.py                         report
"""
import json
import sys

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow

OUT = ROOT / "results" / "robustness"
OUT.mkdir(parents=True, exist_ok=True)
S, DAYS, NU, SEEDS, H = 0.6, 140, 8.0, range(300, 306), 0.3
REG = {"liver": dict(P0=2600.0, rho=0.10, delta=0.10, q_set=0.6),
       "kidney": dict(P0=300.0, rho=0.00, delta=0.05, q_set=0.6),
       "pancreas": dict(P0=60.0, rho=0.00, delta=0.00, q_set=0.6)}
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=6000)
VARIANTS = {
    "reference": {},
    "root_surface": dict(root_site="surface"),
    "root_random": dict(root_site="random"),
    "no_deficit": dict(sprout=0.0),
    "no_load": dict(sprout_load=0.0),
    "no_crowding": dict(sprout_room=0.0),
    "no_mesenchyme": dict(attractor_density=0.0),
    "territory": dict(assign_mode="segment"),
}


def run(variant, reg, seed, sign):
    f = OUT / f"{variant}_{reg}_s{seed}_{sign}.json"
    if f.exists():
        return json.load(open(f))
    p = dict(REG[reg]); p["P0"] *= float(np.exp({"0": 0.0, "+": H, "-": -H}[sign]))
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S), 190 * S,
                        days=DAYS, rng=np.random.default_rng(seed), nu=NU,
                        **FIX, **p, **VARIANTS[variant])
    r = dict(N=len(U), vessels=len(vt), territory_diff=tr.final["territory_diff"],
             **tr.closure)
    json.dump(r, open(f, "w"))
    return r


def summary(variant):
    out = {}
    for reg in REG:
        base = [run(variant, reg, s, "0") for s in SEEDS]
        eps = [(np.log(max(run(variant, reg, s, "+")["N"], 1)) -
                np.log(max(run(variant, reg, s, "-")["N"], 1))) / (2 * H) for s in SEEDS]
        out[reg] = dict(N=float(np.mean([b["N"] for b in base])),
                        **{k: float(np.mean([b[k] for b in base]))
                           for k in ("spent", "supply_held", "space_held", "active", "territory_diff")},
                        eps_pool=float(np.mean(eps)), eps_se=float(np.std(eps, ddof=1) / np.sqrt(len(eps))),
                        eps_each=[float(e) for e in eps])
    return out


if __name__ == "__main__":
    if len(sys.argv) > 1:
        for v in sys.argv[1:]:
            summary(v); print(f"{v}: done", flush=True)
        raise SystemExit(0)
    res = {v: summary(v) for v in VARIANTS
           if all((OUT / f"{v}_{r}_s{s}_{g}.json").exists()
                  for r in REG for s in SEEDS for g in "0+-")}
    json.dump(res, open(ROOT / "results" / "robustness.json", "w"), indent=1)
    print(f"{'variant':>14} {'regime':>9} {'N':>5} {'spent':>6} {'supply':>7} {'room':>6} {'eps_pool':>14}  order")
    for v, R in res.items():
        order = R["liver"]["eps_pool"] < R["kidney"]["eps_pool"] < R["pancreas"]["eps_pool"]
        for i, reg in enumerate(REG):
            x = R[reg]
            print(f"{v if i == 0 else '':>14} {reg:>9} {x['N']:>5.0f} {x['spent']:>6.2f} {x['supply_held']:>7.2f} "
                  f"{x['space_held']:>6.2f} {x['eps_pool']:>+7.2f} ± {x['eps_se']:.2f}  "
                  f"{('holds' if order else 'BROKEN') if i == 0 else ''}")
