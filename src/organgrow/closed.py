"""A closed growth model: one nucleation law, organ identity as a regime.

Units. Three scales are set and everything else is a ratio:

    pressure   dP = 1          the perfusion pressure across the organ
    length     lam = 1         the minimum spacing between units
    radius     r_root = 1      so the root segment's conductance is 1/l_root

With those, conductance is in units of r_root^4 / lam and flow in units of
dP * that. A unit of mass 1 is given terminal conductance 1 / R0, and R0 is
set to the unit count at which the organ should become supply-limited, so
that q = Q / N is of order 1 when the tree and the tissue are matched. q_set
is then a fraction of 1 and nu is per unit of progenitor per unit of excess
supply. The shear setpoint follows from the same scales: it is the shear the
root carries at the design point, dP / (R0 * r_root^3), so it is derived
rather than chosen.

The earlier engine carried three closure modes as separate branches, one of
them a feedback loop and two of them decay curves. Asking which mode destroys
variance then had the answer built in: an open loop conserves the spread
because it reads nothing.

Here there is one law. At each tip i, with supply per unit q_i = Q_i / n_i,

    lambda_i = nu * P_i * (q_i - q_set)_+        units nucleated per unit time
    mu_i     = sprout * (q_set - q_i)_+          vessels sprouted per unit time
    k_i      ~ Poisson(lambda_i dt)              demographic noise
    dP_i/dt  = (rho - delta) * P_i - lambda_i    the progenitor pool

One signal, q_i - q_set, drives both responses with opposite signs: excess
supply makes tissue, a shortfall makes vessels. Neither runs at a rate of its
own, so the tree and the tissue settle against each other.

and the three regimes fall out of two numbers:

    rho ~ delta                  pool holds; capacity binds first  (liver)
    delta > rho                  pool runs out first               (kidney, lung)
    rho = 0, P_i(0) fixed        a budget, spent once              (pancreas)

Identity is therefore a position on the axis (rho - delta, P_0), and the
question of which bound binds first becomes empirical rather than assumed.
Flow is pressure-driven, so a branch carries what its calibre allows, and
the body grows the domain rather than the supply.
"""
from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from .hemodynamics import VascularTree


@dataclass
class Trace:
    t: list = field(default_factory=list)
    N: list = field(default_factory=list)
    P: list = field(default_factory=list)
    q: list = field(default_factory=list)       # supply per unit
    C: list = field(default_factory=list)       # tree conductance
    mass: list = field(default_factory=list)
    bound: list = field(default_factory=list)   # "supply" or "progenitor"


def grow(sdf, R, *, days, rng, nu, q_set, rho, delta, P0, dP=1.0,
         eta=0.12, tau_set=None, pressure_exp=0.5, xi=1.0, kappa=1.0,
         g0=0.30, tau_g=12.0, N_scale=650.0, lam=7.4, s_rate=0.0,
         root_radius=1.0, shrink=0.25, vascular_fraction=0.08, crowding=1.0,
         sprout=4.0,
         max_units=4000, max_vessels=1200):
    """Grow one organ. `sdf(x, y, z)` is negative inside the domain."""
    box = np.array([R * 1.35, R * 1.25, R * 0.85])
    # terminal conductance per unit of mass, set so that q = Q / N is of
    # order one when the tissue has grown to N_scale units
    R0 = float(N_scale)
    # q is measured against the flow one mature unit would draw at the design
    # point: a whole organ of N_scale units has terminal conductance 1, so
    # each unit draws dP / N_scale. Dividing by that puts q of order 1 and
    # makes q_set a fraction of a unit's design supply.
    q_scale = dP / R0
    # room the vessels may take before they compete with the parenchyma,
    # in tree units (lengths in lam, radii in r_root)
    vol_budget = vascular_fraction * (2.0 * np.prod(box / lam))
    # the shear setpoint is not free either: it is the shear the root would
    # carry at the target organ size, so the tree is neither starved nor
    # oversized at its design point
    # Shear on a segment in series with its tissue is not monotone in the
    # radius: tau = Q/r^3 rises while the segment's own resistance dominates
    # and falls once the tissue's does, with a maximum at
    #     r_c = (3 l / R_t)^(1/4)
    # Below r_c the feedback is positive and the vessel closes; above it the
    # feedback is negative and the vessel holds. A setpoint must therefore sit
    # below the peak, and segments must start on its stable side.
    l_root = 2.1
    # A segment does not see the whole tissue resistance: with N_scale units
    # spread over the tips of a tree, each terminal branch carries the
    # resistance of the handful of units below it. Using the single-segment
    # value would put the critical radius orders of magnitude off.
    units_per_tip = max(N_scale / max(max_vessels * 0.5, 1.0), 1.0)
    R_seg = R0 / units_per_tip
    r_crit = (3.0 * l_root / R_seg) ** 0.25
    tau_peak = (dP / (l_root / r_crit ** 4 + R_seg)) / r_crit ** 3
    if tau_set is None:
        tau_set = 0.3 * tau_peak
    # the radius at which a segment in series with its tissue sits at the
    # setpoint: segments are born here, on the stable side, so that
    # adaptation moves them rather than rescuing them
    from scipy.optimize import brentq
    _tau = lambda x: (dP / (l_root / x ** 4 + R_seg)) / x ** 3
    try:
        r_eq = brentq(lambda x: _tau(x) - tau_set, r_crit * 1.001, 100.0)
    except ValueError:
        r_eq = 2.0 * r_crit

    def inside(P, s=1.0):
        return sdf(P[:, 0] / s, P[:, 1] / s, P[:, 2] / s) < 0

    seed = rng.uniform(-1, 1, (20000, 3)) * box
    s0 = seed[inside(seed)][0]
    # the tree lives in units of lam: lengths divided by lam, radii by r_root
    vt = VascularTree(s0 / lam, s0 / lam + np.array([0.0, l_root, 0.0]),
                      r_eq)
    scale = lam                      # tree coordinates -> domain coordinates

    U = np.zeros((0, 3))            # unit positions
    m = np.zeros(0)                 # unit masses
    tip_of = np.zeros(0, int)       # which tip feeds each unit
    P = {0: float(P0)}              # progenitor pool per tip
    tr = Trace()

    for day in range(1, days + 1):
        s = (1.0 + s_rate * day)    # the body grows the space, not the supply

        # ---- perfusion: the tree decides what each branch carries ----------
        if len(U):
            # tissue calls flow toward itself through its own conductance
            tc = {}
            for t_, mj in zip(tip_of, m):
                tc[int(t_)] = tc.get(int(t_), 0.0) + float(mj) ** xi / R0
        else:
            tc = {k: 1.0 / R0 for k in vt.tips()}
        f, p, C = vt.perfuse(tc, dP=dP)

        # ---- calibres chase their shear setpoint ---------------------------
        vt.relax(tau_set, rate=eta, pressure_exp=pressure_exp, p_in=dP,
                 pressures=p, r_min=r_crit, r_max=20.0 * r_eq,
                 shrink=shrink,
                 volume_budget=vol_budget, crowding=crowding)
        # nothing is allowed below the critical radius: there the feedback
        # is positive and the segment would close
        vt.radius = list(np.maximum(np.asarray(vt.radius), r_crit))

        # ---- nucleation: one law, demographic noise ------------------------
        tips = vt.tips()
        n_at = {}
        for t_ in tip_of:
            n_at[int(t_)] = n_at.get(int(t_), 0) + 1
        new_pos, new_tip = [], []
        for i in tips:
            Pi = P.get(i, 0.0)
            if Pi <= 0:
                continue
            q_i = f[i] / (q_scale * max(n_at.get(i, 0), 1))
            drive = max(q_i - q_set, 0.0)
            lam_i = nu * Pi * drive
            k = int(rng.poisson(lam_i)) if lam_i > 0 else 0
            for _ in range(k):
                if len(U) + len(new_pos) >= max_units:
                    break
                c = vt.end[i] * scale + rng.normal(0, lam, 3)
                if sdf(c[0] / s, c[1] / s, c[2] / s) >= 0:
                    continue
                if len(U) and np.min(((U - c) ** 2).sum(-1)) < lam * lam:
                    continue
                if new_pos and np.min(((np.asarray(new_pos) - c) ** 2).sum(-1)) < lam * lam:
                    continue
                new_pos.append(c); new_tip.append(i)
            P[i] = max(Pi + (rho - delta) * Pi - k, 0.0)
        if new_pos:
            U = np.vstack([U, np.asarray(new_pos)])
            m = np.concatenate([m, np.full(len(new_pos), 0.02)])
            tip_of = np.concatenate([tip_of, np.asarray(new_tip, int)])

        # ---- unit growth: the ceiling is the flow the unit receives --------
        if len(U):
            Qj = np.array([f[int(t_)] / (q_scale * max(n_at.get(int(t_), 1), 1))
                           for t_ in tip_of])
            g = g0 * np.exp(-day / tau_g)
            m = m + g * m * (1.0 - m / np.maximum(kappa * Qj, 1e-6))
            m = np.maximum(m, 0.02)

        # ---- the tree sprouts where the tissue is short of supply ---------
        # One signal drives both responses, with opposite signs. Where supply
        # per unit exceeds the setpoint there is room for another unit; where
        # it falls short the tissue calls a vessel. Branching is therefore
        # governed by the same quantity as nucleation instead of running at a
        # rate of its own, which is what lets the tree and the tissue settle
        # against each other.
        if len(vt) < max_vessels and len(U):
            for i in list(tips):
                n_i = n_at.get(i, 0)
                if n_i == 0:
                    continue
                q_i = f[i] / (q_scale * n_i)
                deficit = max(q_set - q_i, 0.0)
                k_v = int(rng.poisson(sprout * deficit)) if deficit > 0 else 0
                for _ in range(k_v):
                    if len(vt) >= max_vessels:
                        break
                    sel = np.flatnonzero(tip_of == i)
                    if not len(sel):
                        break
                    u = U[sel[rng.integers(len(sel))]]
                    b = vt.end[i]
                    nk = vt.add(i, b.copy(), b + (u / scale - b) * 0.55,
                                max(vt.radius[i] * 0.79, r_crit * 1.05))
                    P[nk] = P.get(i, 0.0) * 0.5
                    P[i] = P.get(i, 0.0) * 0.5

        Ptot = float(sum(P.values()))
        qbar = float(np.mean([f[i] / (q_scale * max(n_at.get(i, 1), 1))
                              for i in tips]))
        tr.t.append(day); tr.N.append(len(U)); tr.P.append(Ptot)
        tr.q.append(qbar); tr.C.append(C); tr.mass.append(float(m.sum()))
        tr.bound.append("progenitor" if Ptot < 1e-3 else
                        ("supply" if qbar <= q_set * 1.05 else "growing"))

    return U, m, vt, tr
