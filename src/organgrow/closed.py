"""A closed growth model: one nucleation law, organ identity as a regime.

Units. Three scales are set and everything else is a ratio:

    pressure   dP = 1          the perfusion pressure across the organ
    length     lam = 1         the minimum spacing between units
    radius     r_root = 1      so the root segment's conductance is 1/l_root

With those, conductance is in units of r_root^4 / lam and flow in units of
dP times that. A mature unit (mass 1) has terminal conductance 1 / N_scale,
so an organ of N_scale mature units has tissue conductance 1 and the tree
and the tissue are matched at that size. The shear setpoint follows from the
same scales: a fraction of the peak shear a segment can carry in series with
its tissue at the design point.

The earlier engine carried three closure modes as separate branches, one of
them a feedback loop and two of them decay curves. Asking which mode destroys
variance then had the answer built in: an open loop conserves the spread
because it reads nothing.

Here there is one law. At each tip i the tree delivers an outlet pressure
p_i. The supply a unit can draw there is what a mature unit would draw at
that pressure, in units of the design-point draw:

    q_i = p_i / dP

This is read from the tree, not from the tissue already present, which is
what keeps the loop closed without being circular: the flow a branch carries
is the tissue's demand, but the pressure it delivers is the tree's doing.

    lambda_i = nu * P_i * (q_i - q_set)_+        units nucleated per unit time
    mu_i     = sprout * (q_set - q_i)_+          vessels sprouted per unit time
             + sprout_room * (placements refused for lack of space)
    k_i      ~ Poisson(lambda_i dt)              demographic noise
    dP_i/dt  = (rho - delta) * P_i - lambda_i    the progenitor pool
    dm_j/dt  = g(t) m_j (1 - m_j / (kappa q_i))  the unit grows to its supply

One signal, q_i - q_set, drives both responses with opposite signs: excess
supply makes tissue, a shortfall makes vessels. The three regimes fall out
of two numbers:

    rho ~ delta                  pool holds; capacity binds first  (liver)
    delta > rho                  pool runs out first               (kidney, lung)
    rho = 0, P_i(0) fixed        a budget, spent once              (pancreas)

Identity is a position on the axis (rho - delta, P_0), and which bound binds
first is an outcome. Flow is pressure-driven, so a branch carries what its
calibre allows, and the body grows the domain rather than the supply.

A critical radius falls out of the physics. Shear on a segment in series
with the conductance below it is not monotone in the radius: tau = Q / r^3
rises while the segment's own resistance dominates and falls once the
tissue's does, with a maximum at

    r_c = (3 l / R_below)^(1/4)

Below r_c the calibre feedback is positive and the vessel closes; above it
the vessel holds. It is computed per segment at every step from the real
conductance below it, and it is the floor a segment may not cross.
"""
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq
from scipy.spatial import cKDTree

from .config import CFG
from .hemodynamics import VascularTree

_CL = CFG.get("closed", {})


@dataclass
class Trace:
    t: list = field(default_factory=list)
    N: list = field(default_factory=list)
    P: list = field(default_factory=list)
    q: list = field(default_factory=list)       # mean supply per tip
    C: list = field(default_factory=list)       # tree conductance
    mass: list = field(default_factory=list)
    vessels: list = field(default_factory=list)
    bound: list = field(default_factory=list)   # "growing", "supply" or "progenitor"


def _p(name, value):
    """A keyword given explicitly wins; otherwise the config value."""
    return _CL[name] if value is None else value


def grow(sdf, R, *, days, rng, nu, q_set, rho, delta, P0, dP=1.0,
         lam=7.4, s_rate=0.0, tau_set=None,
         N_scale=None, m0=None, m_leak=None, xi=None, kappa=None, g0=None,
         tau_g=None, eta=None, pressure_exp=None, shrink=None,
         tau_set_fraction=None, daughter_ratio=None, sprout_reach=None,
         sprout=None, sprout_room=None, vascular_fraction=None, crowding=None, l_root=None,
         l_sprout=None, nucleation_spread=None,
         max_units=None, max_vessels=None):
    """Grow one organ. `sdf(x, y, z)` is negative inside the domain.

    The identity parameters (nu, q_set, rho, delta, P0) are required; every
    other parameter defaults to config/params.yaml, section `closed`.
    """
    N_scale = float(_p("N_scale", N_scale)); m0 = _p("m0", m0)
    m_leak = _p("m_leak", m_leak); xi = _p("xi", xi); kappa = _p("kappa", kappa)
    g0 = _p("g0", g0); tau_g = _p("tau_g", tau_g); eta = _p("eta", eta)
    pressure_exp = _p("pressure_exp", pressure_exp); shrink = _p("shrink", shrink)
    tau_set_fraction = _p("tau_set_fraction", tau_set_fraction)
    daughter_ratio = _p("daughter_ratio", daughter_ratio)
    sprout_reach = _p("sprout_reach", sprout_reach); sprout = _p("sprout", sprout)
    sprout_room = _p("sprout_room", sprout_room)
    vascular_fraction = _p("vascular_fraction", vascular_fraction)
    crowding = _p("crowding", crowding); l_root = _p("l_root", l_root)
    l_sprout = _p("l_sprout", l_sprout)
    nucleation_spread = _p("nucleation_spread", nucleation_spread)
    max_units = int(_p("max_units", max_units))
    max_vessels = int(_p("max_vessels", max_vessels))

    box = np.array([R * 1.35, R * 1.25, R * 0.85])
    R0 = N_scale                       # resistance of one mature unit
    # room the vessels may take before they compete with the parenchyma,
    # in tree units (lengths in lam, radii in r_root)
    vol_budget = vascular_fraction * (2.0 * np.prod(box / lam))

    # Design point: a terminal branch carries the handful of units below it.
    # Its shear in series with that tissue peaks at r_c; the setpoint sits
    # below the peak and segments are born on the stable side of it.
    units_per_tip = max(N_scale / max(max_vessels * 0.5, 1.0), 1.0)
    R_seg = R0 / units_per_tip
    r_crit0 = (3.0 * l_root / R_seg) ** 0.25
    _tau = lambda x: (dP / (l_root / x ** 4 + R_seg)) / x ** 3
    tau_peak = _tau(r_crit0)
    if tau_set is None:
        tau_set = tau_set_fraction * tau_peak
    try:
        r_eq = brentq(lambda x: _tau(x) - tau_set, r_crit0 * 1.001, 100.0)
    except ValueError:
        r_eq = 2.0 * r_crit0

    def inside(X, s=1.0):
        return sdf(X[:, 0] / s, X[:, 1] / s, X[:, 2] / s) < 0

    seed = rng.uniform(-1, 1, (20000, 3)) * box
    s0 = seed[inside(seed)][0]
    scale = lam                       # tree coordinates -> domain coordinates
    vt = VascularTree(s0 / scale, s0 / scale + np.array([0.0, l_root, 0.0]), r_eq)

    U = np.zeros((0, 3))              # unit positions
    m = np.zeros(0)                   # unit masses
    tip_of = np.zeros(0, int)         # which tip feeds each unit
    P = {0: float(P0)}                # progenitor pool per tip
    tr = Trace()
    leak = m_leak ** xi / R0          # what a tip carries before its first unit

    for day in range(1, days + 1):
        s = 1.0 + s_rate * day        # the body grows the space, not the supply
        tips = vt.tips()
        n_at = np.bincount(tip_of, minlength=len(vt)) if len(tip_of) else np.zeros(len(vt), int)
        # tissue hangs along a vessel, not only at its end: a segment that has
        # branched keeps the units and the progenitors it already carried, so
        # every segment with tissue or pool is a site of nucleation and sprouting
        active = [k for k in range(len(vt)) if n_at[k] > 0 or P.get(k, 0.0) > 0]

        # ---- perfusion: the tree decides what each branch carries ----------
        tc = {i: leak for i in tips}
        for t_, mj in zip(tip_of, m):
            tc[int(t_)] = tc.get(int(t_), 0.0) + float(mj) ** xi / R0
        f, p_in, C = vt.perfuse(tc, dP=dP)
        p_out = vt.pressure_out
        q = np.maximum(p_out, 0.0) / dP            # supply per unit, by tip

        # ---- calibres chase their shear setpoint ---------------------------
        # the transmural pressure of a segment is the mean of inlet and outlet
        vt.relax(tau_set, rate=eta, pressure_exp=pressure_exp, p_in=dP,
                 pressures=0.5 * (p_in + p_out),
                 r_min=vt.r_crit, r_max=20.0 * r_eq, shrink=shrink,
                 volume_budget=vol_budget, crowding=crowding)

        # ---- nucleation: one law, demographic noise ------------------------
        new_pos, new_tip = [], []
        no_room = np.zeros(len(vt), int)   # placements refused for lack of space
        for i in active:
            Pi = P.get(i, 0.0)
            if Pi <= 0:
                continue
            drive = max(q[i] - q_set, 0.0)
            lam_i = nu * Pi * drive
            k = int(rng.poisson(lam_i)) if lam_i > 0 else 0
            placed = 0
            for _ in range(k):
                if len(U) + len(new_pos) >= max_units:
                    break
                c = vt.end[i] * scale + rng.normal(0, nucleation_spread * lam, 3)
                if sdf(c[0] / s, c[1] / s, c[2] / s) >= 0:
                    no_room[i] += 1
                    continue
                if len(U) and np.min(((U - c) ** 2).sum(-1)) < lam * lam:
                    no_room[i] += 1
                    continue
                if new_pos and np.min(((np.asarray(new_pos) - c) ** 2).sum(-1)) < lam * lam:
                    no_room[i] += 1
                    continue
                new_pos.append(c); new_tip.append(i)
                placed += 1
            # progenitors are spent by the units made, not by the ones refused
            P[i] = max(Pi + (rho - delta) * Pi - placed, 0.0)
        if new_pos:
            U = np.vstack([U, np.asarray(new_pos)])
            m = np.concatenate([m, np.full(len(new_pos), m0)])
            tip_of = np.concatenate([tip_of, np.asarray(new_tip, int)])
            n_at = np.bincount(tip_of, minlength=len(vt))

        # ---- unit growth: the ceiling is the supply at the unit's tip ------
        if len(U):
            g = g0 * np.exp(-day / tau_g)
            ceiling = np.maximum(kappa * q[tip_of], 1e-6)
            m = m + g * m * (1.0 - m / ceiling)
            m = np.maximum(m, m0)

        # ---- the tree sprouts where the tissue is short of supply ---------
        # The same signal drives both responses with opposite signs: where
        # supply per unit exceeds the setpoint there is room for another
        # unit; where it falls short the tissue calls a vessel. A sprout
        # grows toward one of the units that called it and takes with it the
        # units now closer to its end, and the progenitor pool splits in the
        # same proportion.
        if len(vt) < max_vessels and len(U):
            utree = cKDTree(U)
            for i in active:
                # a second, geometric reason to branch: the tissue at this
                # segment had supply to make units but no room to put them, so
                # the tree must reach into space that is still empty. The
                # sprout points away from the units around the segment's end
                # and reaches past the cloud they occupy.
                if no_room[i] > 0:
                    k_r = int(rng.poisson(sprout_room * no_room[i]))
                    for _ in range(k_r):
                        if len(vt) >= max_vessels:
                            break
                        b = vt.end[i]
                        near = utree.query_ball_point(b * scale, r=3.0 * lam)
                        away = (b - U[near].mean(0) / scale) if near else rng.normal(size=3)
                        away = away + 0.5 * rng.normal(size=3)
                        away /= max(np.linalg.norm(away), 1e-9)
                        end = b + away * l_sprout
                        if sdf(end[0] * scale / s, end[1] * scale / s, end[2] * scale / s) >= 0:
                            continue
                        nk = vt.add(i, b.copy(), end,
                                    max(vt.radius[i] * daughter_ratio, vt.r_crit[i] * 1.05))
                        # the new tip inherits half the pool: it is where the
                        # next units will go
                        P[nk] = P.get(i, 0.0) * 0.5
                        P[i] = P.get(i, 0.0) * 0.5
                if n_at[i] == 0:
                    continue
                deficit = max(q_set - q[i], 0.0)
                k_v = int(rng.poisson(sprout * deficit)) if deficit > 0 else 0
                for _ in range(k_v):
                    if len(vt) >= max_vessels:
                        break
                    sel = np.flatnonzero(tip_of == i)
                    if not len(sel):
                        break
                    u = U[sel[rng.integers(len(sel))]]
                    b = vt.end[i]
                    end = b + (u / scale - b) * sprout_reach
                    nk = vt.add(i, b.copy(), end,
                                max(vt.radius[i] * daughter_ratio, vt.r_crit[i] * 1.05))
                    d_old = ((U[sel] - b * scale) ** 2).sum(-1)
                    d_new = ((U[sel] - end * scale) ** 2).sum(-1)
                    moved = sel[d_new < d_old]
                    tip_of[moved] = nk
                    frac = len(moved) / len(sel)
                    P[nk] = P.get(i, 0.0) * frac
                    P[i] = P.get(i, 0.0) * (1.0 - frac)
                    sel = np.flatnonzero(tip_of == i)
                    if not len(sel):
                        break
            n_at = np.bincount(tip_of, minlength=len(vt))

        Ptot = float(sum(P.values()))
        # q was read at the start of the day, before any sprout was added
        live = [i for i in active if n_at[i] > 0] or tips
        qbar = float(np.mean(q[live]))
        tr.t.append(day); tr.N.append(len(U)); tr.P.append(Ptot)
        tr.q.append(qbar); tr.C.append(C); tr.mass.append(float(m.sum()))
        tr.vessels.append(len(vt))
        tr.bound.append("progenitor" if Ptot < 1e-3 else
                        ("supply" if qbar <= q_set * 1.05 else "growing"))

    return U, m, vt, tr
