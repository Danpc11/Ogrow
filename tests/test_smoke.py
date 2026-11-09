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

print()
if FAIL:
    print("FAILED:", ", ".join(FAIL))
    sys.exit(1)
print("all smoke tests passed")
