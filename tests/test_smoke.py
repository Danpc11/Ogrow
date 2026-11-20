#!/usr/bin/env python
"""Smoke tests: a few seconds, no data, no figures.

Run with `make test` or `python tests/test_smoke.py`. Each check is a property
that an earlier version of the code violated, so a failure here is a known
kind of regression rather than a style preference.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from organgrow import grow, sdf                      # noqa: E402
from organgrow.config import CFG                     # noqa: E402
from organgrow.hemodynamics import VascularTree      # noqa: E402
from organgrow import closed                         # noqa: E402

FAIL = []


def check(name, cond, detail=""):
    print(f"  {'ok  ' if cond else 'FAIL'} {name} {detail}")
    if not cond:
        FAIL.append(name)


print("engine")
for kind, ident in CFG["identity"].items():
    U, m, T, V, tr = grow(kind, seed=1, days=15, snap_at=())
    d = cKDTree(U).query(U, k=2)[0][:, 1]
    check(f"{kind}: grows", len(U) > 10, f"N={len(U)}")
    check(f"{kind}: spacing >= lam", d.min() >= ident["lam"] * 0.999,
          f"min nn={d.min():.2f} lam={ident['lam']}")
    check(f"{kind}: masses above floor", np.mean(m <= 0.0201) < 0.5,
          f"floor fraction={np.mean(m <= 0.0201):.2f}")
    check(f"{kind}: finite trace", np.all(np.isfinite(tr.cu)))

print("hemodynamics")
rng = np.random.default_rng(0)
vt = VascularTree([0, 0, 0], [0, 2, 0], 1.0)
for _ in range(60):
    p = int(rng.integers(len(vt)))
    b = vt.end[p]
    vt.add(p, b.copy(), b + rng.normal(0, 1.5, 3), vt.radius[p] * 0.79)
tips = vt.tips()
tc = {i: 0.01 for i in tips[: len(tips) // 2]}      # half the tips carry tissue
f, p_in, C = vt.perfuse(tc, dP=1.0)
check("flow conserved at every branch",
      all(abs(f[k] - sum(f[c] for c in range(len(vt)) if vt.parent[c] == k)) < 1e-9
          for k in range(len(vt)) if not vt.is_tip[k]))
check("dead-end tips carry no flow",
      all(f[i] == 0.0 for i in tips[len(tips) // 2:]))
check("conductance bounded by the tissue", C <= sum(tc.values()) * (1 + 1e-9),
      f"C={C:.4g} tissue={sum(tc.values()):.4g}")
check("pressures fall along the tree",
      np.all(vt.pressure_out <= vt.pressure_in + 1e-12) and np.all(vt.pressure_out >= -1e-12))
check("critical radius positive where tissue hangs",
      np.all(vt.r_crit[vt.below > 0] > 0))

print("closed model")
rng = np.random.default_rng(1)
U, m, vt, tr = closed.grow(lambda x, y, z: sdf("liver", x, y, z), 190, days=30, rng=rng,
                           nu=1.0, q_set=0.5, rho=0.0, delta=0.0, P0=300.0)
leak = CFG["closed"]["m_leak"] ** CFG["closed"]["xi"] / CFG["closed"]["N_scale"]
tissue = (m ** CFG["closed"]["xi"]).sum() / CFG["closed"]["N_scale"] + leak * len(vt.tips())
check("grows", tr.N[-1] > 50, f"N={tr.N[-1]}")
check("tree conductance bounded by tissue plus leak", tr.C[-1] <= tissue * (1 + 1e-9),
      f"C={tr.C[-1]:.4g} bound={tissue:.4g}")
check("no runaway calibre", max(vt.radius) < 20 * 2.0, f"r max={max(vt.radius):.2f}")
check("supply stays in [0, 1]", 0.0 <= min(tr.q) and max(tr.q) <= 1.0 + 1e-9)
d = cKDTree(U).query(U, k=2)[0][:, 1]
check("unit spacing >= lam", d.min() >= 7.4 * 0.999, f"min nn={d.min():.2f}")

print("derivations")
# the critical radius is where shear actually peaks: one segment of length l
# in series with tissue of conductance G, shear read from perfuse() itself
from organgrow.hemodynamics import VascularTree
l, G = 2.0, 0.05
rs = np.logspace(-1.0, 0.6, 400); sh = []
for r in rs:
    t1 = VascularTree(np.zeros(3), np.array([0.0, l, 0.0]), r)
    t1.perfuse({0: G}); sh.append(float(t1.shear[0]))
r_peak = rs[int(np.argmax(sh))]
check("critical radius is the shear peak", abs(r_peak / float(t1.r_crit[0]) - 1) < 0.02,
      f"peak {r_peak:.3f}, r_crit {float(t1.r_crit[0]):.3f}")
# flow is conserved at a node that both branches and carries tissue
t2 = VascularTree(np.zeros(3), np.array([0.0, 1.0, 0.0]), 1.0)
t2.add(0, np.array([0.0, 1.0, 0.0]), np.array([0.0, 2.0, 0.0]), 0.8)
t2.add(0, np.array([0.0, 1.0, 0.0]), np.array([1.0, 1.0, 0.0]), 0.8)
Gi = 0.3
t2.perfuse({0: Gi, 1: 0.5, 2: 0.2})
local = t2.pressure_out[0] * Gi
check("flow conserved where an internal segment carries tissue",
      abs(t2.flow[0] - (t2.flow[1] + t2.flow[2] + local)) < 1e-9 * t2.flow[0],
      f"in {t2.flow[0]:.6f}, out {t2.flow[1] + t2.flow[2] + local:.6f}")
# progenitors are conserved exactly in every realisation: with no division,
# what was there at the start is now units, progenitors left, or lost
from organgrow import sdf as _sdf
from organgrow.closed import grow as _grow, closure_state
bal = []
for k in range(3):
    U_, m_, vt_, tr_ = _grow(lambda x, y, z: _sdf("liver", x / 0.6, y / 0.6, z / 0.6), 114,
                             days=60, rng=np.random.default_rng(10 + k), nu=8.0, q_set=0.6,
                             rho=0.0, delta=0.03, P0=40.0, R_unit=650.0,
                             units_per_branch=650 / 600.0, max_vessels=3000)
    bal.append(40 - (len(U_) + int(round(tr_.P[-1])) + tr_.lost[-1]))
check("progenitors conserved exactly in each realisation", all(b == 0 for b in bal),
      f"P0 - (N + P + lost) = {bal}")
# the closure state tells a tip held by lack of room from one held by supply
st = closure_state({0: 5, 1: 5, 2: 5, 3: 5}, np.array([0.5, 0.9, 0.9, 0.9]), None,
                   0.6, np.array([3, 4, 0, 2]), np.array([0, 4, 0, 1]), 40.0)
check("closure separates supply, room and active",
      abs(st["supply_held"] - 5 / 40) < 1e-12 and abs(st["space_held"] - 5 / 40) < 1e-12
      and abs(st["active"] - 10 / 40) < 1e-12, str(st))
# an unreachable shear setpoint is an error, not a silently substituted radius
try:
    _grow(lambda x, y, z: _sdf("liver", x, y, z), 114, days=2, rng=np.random.default_rng(0),
          nu=8.0, q_set=0.6, rho=0.0, delta=0.0, P0=10.0, tau_set=1e9)
    check("unreachable setpoint raises", False, "no error")
except ValueError:
    check("unreachable setpoint raises", True)

# units assigned by vascular territory still hang from tips, and an unknown
# root site is an error rather than a silent default
U_, m_, vt_, tr_ = _grow(lambda x, y, z: _sdf("liver", x / 0.6, y / 0.6, z / 0.6), 114,
                         days=40, rng=np.random.default_rng(5), nu=8.0, q_set=0.6, rho=0.0,
                         delta=0.0, P0=40.0, R_unit=650.0, units_per_branch=650 / 600.0,
                         max_vessels=3000, assign_mode="segment", root_site="surface")
from organgrow.closed import _assign
own = _assign(vt_, U_, vt_.tips(), 7.4, "segment")
check("territory assignment returns tips", set(own.tolist()) <= set(vt_.tips()))
try:
    _grow(lambda x, y, z: _sdf("liver", x, y, z), 114, days=2, rng=np.random.default_rng(0),
          nu=8.0, q_set=0.6, rho=0.0, delta=0.0, P0=10.0, root_site="nowhere")
    check("unknown root site raises", False, "no error")
except ValueError:
    check("unknown root site raises", True)

print()
if FAIL:
    print("FAILED:", ", ".join(FAIL))
    sys.exit(1)
print("all smoke tests passed")
