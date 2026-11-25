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
st = closure_state({0: 5, 1: 5, 2: 5, 3: 5}, np.array([0.5, 0.9, 0.9, 0.9]),
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

# the critical radius satisfies tau'(r_c) = 0 and tau''(r_c) < 0 analytically
_l, _G = 2.0, 0.5
_rc = (_l * _G / 3.0) ** 0.25
_tau = lambda r: (1.0 / (_l / r ** 4 + 1.0 / _G)) / r ** 3
_h = 1e-4 * _rc
_d1 = (_tau(_rc + _h) - _tau(_rc - _h)) / (2 * _h)
_d2 = (_tau(_rc + _h) - 2 * _tau(_rc) + _tau(_rc - _h)) / _h ** 2
check("shear has zero slope and negative curvature at r_c",
      abs(_d1) < 1e-6 * _tau(_rc) / _rc and _d2 < 0, f"tau' {_d1:.2e}, tau'' {_d2:.2e}")
# hydraulic consistency: the flow the root draws is the flow the tissue takes
t3 = VascularTree(np.zeros(3), np.array([0.0, 1.0, 0.0]), 1.0)
t3.add(0, np.array([0.0, 1.0, 0.0]), np.array([0.0, 2.0, 0.0]), 0.7)
t3.add(0, np.array([0.0, 1.0, 0.0]), np.array([1.0, 1.0, 0.0]), 0.6)
t3.add(1, np.array([0.0, 2.0, 0.0]), np.array([0.0, 3.0, 0.0]), 0.5)
loads = {0: 0.1, 2: 0.4, 3: 0.3}
t3.perfuse(loads)
taken = sum(t3.pressure_out[k] * g for k, g in loads.items())
check("root flow equals the flow the tissue takes", abs(t3.flow[0] - taken) < 1e-9 * t3.flow[0],
      f"root {t3.flow[0]:.6f}, tissue {taken:.6f}")
# progenitor accounting with division: P0 + born = units + left + lost
acc = []
for k in range(2):
    U_, m_, vt_, tr_ = _grow(lambda x, y, z: _sdf("liver", x / 0.6, y / 0.6, z / 0.6), 114,
                             days=50, rng=np.random.default_rng(20 + k), nu=8.0, q_set=0.6,
                             rho=0.05, delta=0.03, P0=40.0, R_unit=650.0,
                             units_per_branch=650 / 600.0, max_vessels=3000)
    acc.append(40 + tr_.born[-1] - (len(U_) + int(round(tr_.P[-1])) + tr_.lost[-1]))
check("progenitors balance with division and loss", all(a == 0 for a in acc), f"{acc}")

# binomial nucleation with a split day conserves progenitors just as exactly
accb = []
for k in range(2):
    U_, m_, vt_, tr_ = _grow(lambda x, y, z: _sdf("liver", x / 0.6, y / 0.6, z / 0.6), 114,
                             days=40, rng=np.random.default_rng(40 + k), nu=8.0, q_set=0.6,
                             rho=0.05, delta=0.03, P0=40.0, R_unit=650.0,
                             units_per_branch=650 / 600.0, max_vessels=3000,
                             nucleation="binomial", substeps=2)
    accb.append(40 + tr_.born[-1] - (len(U_) + int(round(tr_.P[-1])) + tr_.lost[-1]))
check("binomial nucleation with substeps balances progenitors", all(a == 0 for a in accb), f"{accb}")

# competing fates conserve progenitors exactly as well
accc = []
for k in range(2):
    U_, m_, vt_, tr_ = _grow(lambda x, y, z: _sdf("liver", x / 0.6, y / 0.6, z / 0.6), 114,
                             days=40, rng=np.random.default_rng(60 + k), nu=8.0, q_set=0.6,
                             rho=0.05, delta=0.03, P0=40.0, R_unit=650.0,
                             units_per_branch=650 / 600.0, max_vessels=3000, turnover="competing")
    accc.append(40 + tr_.born[-1] - (len(U_) + int(round(tr_.P[-1])) + tr_.lost[-1]))
check("competing fates balance progenitors", all(a == 0 for a in accc), f"{accc}")

# where the differentiated cell is its own progenitor (renewal = 1, the liver),
# the pool is not drawn down by unit formation and the books still balance
accr = []
for k in range(2):
    U_, m_, vt_, tr_ = _grow(lambda x, y, z: _sdf("liver", x / 0.6, y / 0.6, z / 0.6), 114,
                             days=40, rng=np.random.default_rng(80 + k), nu=8.0, q_set=0.6,
                             rho=0.05, delta=0.03, P0=40.0, R_unit=650.0,
                             units_per_branch=650 / 600.0, max_vessels=3000, renewal=1.0)
    accr.append(40 + tr_.born[-1] + tr_.renewed[-1]
                - (len(U_) + int(round(tr_.P[-1])) + tr_.lost[-1]))
check("self-renewal balances progenitors", all(a == 0 for a in accr), f"{accr}")
# and a self-renewing pool is never spent by making units
U_, m_, vt_, tr_ = _grow(lambda x, y, z: _sdf("liver", x / 0.6, y / 0.6, z / 0.6), 114,
                         days=40, rng=np.random.default_rng(90), nu=8.0, q_set=0.6,
                         rho=0.0, delta=0.0, P0=40.0, R_unit=650.0,
                         units_per_branch=650 / 600.0, max_vessels=3000, renewal=1.0)
check("a self-renewing pool is not spent", int(round(tr_.P[-1])) == 40,
      f"P0 40, final {tr_.P[-1]:.0f}, units {len(U_)}")

# an occluded segment carries no flow, and neither does what hangs below it
t4 = VascularTree(np.zeros(3), np.array([0.0, 1.0, 0.0]), 1.0)
t4.add(0, np.array([0.0, 1.0, 0.0]), np.array([0.0, 2.0, 0.0]), 0.8)
t4.add(0, np.array([0.0, 1.0, 0.0]), np.array([1.0, 1.0, 0.0]), 0.8)
t4.add(1, np.array([0.0, 2.0, 0.0]), np.array([0.0, 3.0, 0.0]), 0.6)
t4.perfuse({2: 0.4, 3: 0.3})
open_flow = float(t4.flow[2])
t4.occluded[1] = True
t4.perfuse({2: 0.4, 3: 0.3})
check("an occluded segment and its subtree carry no flow",
      abs(float(t4.flow[1])) < 1e-12 and abs(float(t4.flow[3])) < 1e-12
      and float(t4.flow[2]) > 0, f"occluded {t4.flow[1]:.2e}, below {t4.flow[3]:.2e}")
check("occlusion leaves the other branch conducting", float(t4.flow[2]) > 0.5 * open_flow,
      f"{t4.flow[2]:.5f} against {open_flow:.5f} before")
# a plane cut removes the fraction of units it was asked for
U_, m_, vt_, tr_ = _grow(lambda x, y, z: _sdf("liver", x / 0.6, y / 0.6, z / 0.6), 114,
                         days=45, rng=np.random.default_rng(7), nu=8.0, q_set=0.6,
                         rho=0.0, delta=0.0, P0=400.0, R_unit=650.0,
                         units_per_branch=650 / 600.0, max_vessels=3000,
                         injure_day=40, injury="parenchyma", injury_fraction=0.3)
i_ = tr_.injury
got = 1 - i_["post"]["N"] / max(i_["pre"]["N"], 1)
check("a plane cut removes the fraction asked for", abs(got - 0.3) < 0.03,
      f"asked 0.30, removed {got:.3f} ({i_['pre']['N']} to {i_['post']['N']})")
check("the injury is recorded with its day and kind",
      i_["day"] == 40 and i_["kind"] == "parenchyma", str({k: i_[k] for k in ("day", "kind")}))
try:
    _grow(lambda x, y, z: _sdf("liver", x, y, z), 114, days=3, rng=np.random.default_rng(0),
          nu=8.0, q_set=0.6, rho=0.0, delta=0.0, P0=10.0, injure_day=2, injury="meteorite")
    check("unknown injury raises", False, "no error")
except ValueError:
    check("unknown injury raises", True)

# growth reactivated after the injury rebuilds; a fixed inflow stops the rebuild
def _after_injury(**kw):
    U2, m2, vt2, tr2 = _grow(lambda x, y, z: _sdf("liver", x / 0.6, y / 0.6, z / 0.6), 114,
                             days=70, rng=np.random.default_rng(7), nu=8.0, q_set=0.6,
                             rho=0.0, delta=0.0, P0=400.0, R_unit=650.0,
                             units_per_branch=650 / 600.0, max_vessels=3000,
                             injure_day=40, injury="parenchyma", injury_fraction=0.3, **kw)
    return float(m2.sum()) / tr2.injury["pre"]["mass"]

# the declared boundary condition has to hold exactly, not only in the ceiling
_U3, _m3, _vt3, _tr3 = _grow(lambda x, y, z: _sdf("liver", x / 0.6, y / 0.6, z / 0.6), 114,
                             days=70, rng=np.random.default_rng(7), nu=8.0, q_set=0.6,
                             rho=0.0, delta=0.0, P0=400.0, R_unit=650.0,
                             units_per_branch=650 / 600.0, max_vessels=3000,
                             injure_day=40, injury="parenchyma", injury_fraction=0.3,
                             injury_reactivates=1.0, inflow_limited=True)
_pre_Q = _tr3.Q[38]
_drift = max(abs(x - _pre_Q) for x in _tr3.Q[40:]) / max(_pre_Q, 1e-12)
check("a fixed inflow holds the root flow at its pre-injury value",
      _drift < 1e-9, f"largest departure {_drift:.2e} of {_pre_Q:.6f}")

_quiet = _after_injury()
_react = _after_injury(injury_reactivates=1.0)
_held = _after_injury(injury_reactivates=1.0, inflow_limited=True)
check("without reactivation the organ stays below its old mass", _quiet < 0.95,
      f"mass {_quiet:.2f} of pre-injury")
check("reactivating growth rebuilds the mass", _react > _quiet + 0.2 and _react > 0.95,
      f"{_quiet:.2f} without, {_react:.2f} with")
check("at a fixed inlet pressure the rebuild overshoots", _react > 1.3,
      f"mass {_react:.2f} of pre-injury")
check("a fixed inflow stops the rebuild at the old size", 0.9 < _held < 1.15,
      f"mass {_held:.2f} of pre-injury")

# The size it stops at is the delivered flow over what the tissue needs per
# unit mass. Left to itself the second is read off the organ's own pre-injury
# state, which forces it back to exactly that mass; these two check that the
# dial is real and that leaving it alone changes nothing.
_more = _after_injury(injury_reactivates=1.0, inflow_limited=True, inflow_scale=1.4)
_less = _after_injury(injury_reactivates=1.0, inflow_limited=True, inflow_scale=0.6)
check("a larger delivered flow lands the organ on a larger mass",
      _more > _held + 0.1 and _less < _held - 0.05,
      f"0.6x {_less:.2f}, 1.0x {_held:.2f}, 1.4x {_more:.2f} of pre-injury")
_same = _after_injury(injury_reactivates=1.0, inflow_limited=True,
                      inflow_scale=1.0, phi_demand=None)
check("leaving the demand unset leaves the reference behaviour untouched",
      abs(_same - _held) < 1e-12, f"{_same:.6f} against {_held:.6f}")

# A plateau in the total is not the same as every unit reaching its ceiling.
# The ceilings are reported so the two can be told apart; this checks that they
# are there, one per unit, and that the total mass can match the budget while
# the units are individually still far from theirs.
_Uk, _mk, _vtk, _trk = _grow(lambda x, y, z: _sdf("liver", x / 0.6, y / 0.6, z / 0.6), 114,
                             days=70, rng=np.random.default_rng(7), nu=8.0, q_set=0.6,
                             rho=0.0, delta=0.0, P0=400.0, R_unit=650.0,
                             units_per_branch=650 / 600.0, max_vessels=3000,
                             injure_day=40, injury="parenchyma", injury_fraction=0.3,
                             injury_reactivates=1.0, inflow_limited=True)
_cl = _trk.final["ceiling"]
check("the ceiling is reported, one per unit",
      _cl is not None and len(_cl) == len(_mk), f"units {len(_mk)}, ceilings "
      f"{'none' if _cl is None else len(_cl)}")
_eps = float(np.abs(_mk - _cl).sum() / max(_cl.sum(), 1e-12))
check("the total can sit on the budget while the units do not",
      abs(float(_mk.sum()) / max(float(_cl.sum()), 1e-12) - 1.0) < 0.25 and _eps > 0.02,
      f"total {float(_mk.sum()) / max(float(_cl.sum()), 1e-12):.3f} of the "
      f"budget, residual {_eps:.3f}")

# the final per-segment state covers every segment, so tips index into it
from organgrow.closed import _q_full as _qf
_fin = tr_.final
check("final q and P cover every segment",
      len(_fin["q"]) == len(vt_.parent) and len(_fin["P"]) == len(vt_.parent),
      f"segments {len(vt_.parent)}, q {len(_fin['q'])}, P {len(_fin['P'])}")
try:
    _ = _fin["q"][_fin["tips"]], _fin["P"][_fin["tips"]]
    check("tips index into the final q and P", True)
except IndexError as e:
    check("tips index into the final q and P", False, str(e))
check("a segment born after the solve takes its parent's supply",
      list(_qf(np.array([0.4, 0.7]), np.array([-1, 0, 1, 0]), 4)) == [0.4, 0.7, 0.7, 0.4],
      str(list(_qf(np.array([0.4, 0.7]), np.array([-1, 0, 1, 0]), 4))))

# where the pool went has to add up, and a pool that empties is not the same
# as a pool that was used
for _lab, _kw in (("renewing", dict(rho=0.10, delta=0.10, renewal=1.0)),
                  ("budget", dict(rho=0.0, delta=0.05, renewal=0.0))):
    _U4, _m4, _vt4, _tr4 = _grow(lambda x, y, z: _sdf("liver", x / 0.6, y / 0.6, z / 0.6), 114,
                                 days=50, rng=np.random.default_rng(11), nu=8.0, q_set=0.6,
                                 P0=300.0, R_unit=650.0, units_per_branch=650 / 600.0,
                                 max_vessels=3000, **_kw)
    _f = _tr4.fate
    check(f"the pool balances in a {_lab} organ", abs(_f["balance"]) < 1e-6,
          f"P0 + born - (consumed + left + lost + culled) = {_f['balance']:.3g}")
    check(f"consumed and lost are reported apart in a {_lab} organ",
          _f["consumed"] >= 0 and _f["lost"] >= 0,
          f"consumed {_f['consumed']:.2f}, lost {_f['lost']:.2f}")
    # The books have to close when units are destroyed too. A unit removed by
    # an injury was still made from a cell, so consumed counts formations and
    # not survivors, and pool cells destroyed with their niche are culled, not
    # lost to attrition.
    for _ik in ("parenchyma", "vascular"):
        _U5, _m5, _vt5, _tr5 = _grow(lambda x, y, z: _sdf("liver", x / 0.6, y / 0.6, z / 0.6), 114,
                                     days=70, rng=np.random.default_rng(11), nu=8.0,
                                     q_set=0.6, P0=300.0, R_unit=650.0,
                                     units_per_branch=650 / 600.0, max_vessels=3000,
                                     injure_day=40, injury=_ik, injury_fraction=0.3,
                                     **_kw)
        _g5 = _tr5.fate
        check(f"the pool balances after a {_ik} injury in a {_lab} organ",
              abs(_g5["balance"]) < 1e-6,
              f"balance {_g5['balance']:.3g}, formed {_g5['units_formed']}, "
              f"left {_g5['units_left']}, culled {_g5['culled']:.2f}")

# q_floor defaults to zero, so the reference model is unchanged by its presence
_Ua, _ma, _vta, _tra = _grow(lambda x, y, z: _sdf("liver", x / 0.6, y / 0.6, z / 0.6), 114,
                             days=40, rng=np.random.default_rng(5), nu=8.0, q_set=0.6,
                             rho=0.0, delta=0.0, P0=120.0, R_unit=650.0,
                             units_per_branch=650 / 600.0, max_vessels=3000)
_Ub, _mb, _vtb, _trb = _grow(lambda x, y, z: _sdf("liver", x / 0.6, y / 0.6, z / 0.6), 114,
                             days=40, rng=np.random.default_rng(5), nu=8.0, q_set=0.6,
                             rho=0.0, delta=0.0, P0=120.0, R_unit=650.0, q_floor=0.0,
                             units_per_branch=650 / 600.0, max_vessels=3000)
check("q_floor at zero leaves the reference model untouched",
      _tra.N[-1] == _trb.N[-1] and len(_Ua) == len(_Ub),
      f"{_tra.N[-1]} against {_trb.N[-1]}")
_Uc, _mc, _vtc, _trc = _grow(lambda x, y, z: _sdf("liver", x / 0.6, y / 0.6, z / 0.6), 114,
                             days=40, rng=np.random.default_rng(5), nu=8.0, q_set=0.6,
                             rho=0.0, delta=0.0, P0=120.0, R_unit=650.0, q_floor=0.05,
                             units_per_branch=650 / 600.0, max_vessels=3000)
check("a residual drive below threshold spends more of the pool",
      _trc.fate["consumed"] > _tra.fate["consumed"] + 0.05,
      f"consumed {_tra.fate['consumed']:.2f} without it, {_trc.fate['consumed']:.2f} with it")

# involution defaults to zero, so the reference model is unchanged by it
_Ud, _md, _vtd, _trd = _grow(lambda x, y, z: _sdf("liver", x / 0.6, y / 0.6, z / 0.6), 114,
                             days=40, rng=np.random.default_rng(5), nu=8.0, q_set=0.6,
                             rho=0.0, delta=0.0, P0=120.0, R_unit=650.0, involution=0.0,
                             units_per_branch=650 / 600.0, max_vessels=3000)
check("involution at zero leaves the reference model untouched",
      len(_Ud) == len(_Ua) and float(_md.sum()) == float(_ma.sum()),
      f"{len(_Ua)} units and mass {float(_ma.sum()):.6f} against {len(_Ud)} "
      f"and {float(_md.sum()):.6f}")


# A unit that loses its supply is removed rather than left at the mass floor.
# Checked where the floor population is largest: a renewing construct held at a
# fixed inflow after injury, which is where the share rule strands units on
# tips that deliver almost nothing.
def _inv(rate):
    U6, m6, vt6, tr6 = _grow(lambda x, y, z: _sdf("liver", x / 0.6, y / 0.6, z / 0.6), 114,
                             days=70, rng=np.random.default_rng(7), nu=8.0, q_set=0.6,
                             rho=0.10, delta=0.10, P0=400.0, R_unit=650.0,
                             units_per_branch=650 / 600.0, max_vessels=3000,
                             injure_day=40, injury="parenchyma", injury_fraction=0.3,
                             injury_reactivates=1.0, inflow_limited=True,
                             involution=rate)
    q6 = np.asarray(tr6.final["q"], float)[np.asarray(tr6.final["tip_of"], int)]
    return (float((m6 <= 1.05 * 0.02).mean()), len(m6), tr6.fate,
            float(q6.min() / 0.6) if len(q6) else float("nan"))


_f0, _n0, _fa0, _q0 = _inv(0.0)
_f1, _n1, _fa1, _q1 = _inv(0.05)
check("involution clears units stranded at the mass floor", _f1 < _f0 - 0.1,
      f"{_f0:.2f} of units at the floor without it, {_f1:.2f} with it")
check("involution removes units rather than renaming the count",
      _n1 < _n0 and _fa1["units_formed"] > _fa1["units_left"],
      f"{_n0} units without it, {_n1} with it, of {_fa1['units_formed']} ever formed")
check("the pool still balances when units are removed",
      abs(_fa1["balance"]) < 1e-6, f"balance {_fa1['balance']:.3g}")
check("what involution leaves behind is supplied", _q1 > _q0,
      f"worst supply {_q0:.3f} of the setpoint without it, {_q1:.3f} with it")

print()
if FAIL:
    print("FAILED:", ", ".join(FAIL))
    sys.exit(1)
print("all smoke tests passed")
