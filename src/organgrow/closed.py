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
    k_i      = min(Poisson(lambda_i dt), P_i)    units made; or Binomial(P_i,
                                                 1 - exp(-nu (q_i - q_set)_+ dt))
    P_i is a whole number of progenitors: each unit made spends one, and the
    survivors divide with probability rho dt and are lost with probability
    delta dt (independently, or as competing fates). P0 + born = units + left
    + lost holds exactly in every realisation.
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

    r_c = (l / (3 R_below))^(1/4)

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
    lost: list = field(default_factory=list)   # progenitors lost, cumulative
    born: list = field(default_factory=list)   # progenitors born, cumulative
    renewed: list = field(default_factory=list)  # progenitors returned, cumulative
    closure: dict = field(default_factory=dict)  # closure_by_trials at the end
    Q: list = field(default_factory=list)      # total flow through the root
    bound: list = field(default_factory=list)   # "growing", "supply" or "progenitor"


def _p(name, value):
    """A keyword given explicitly wins; otherwise the config value."""
    return _CL[name] if value is None else value


def _assign(vt, U, tips, scale, mode):
    """The tip that feeds each unit.

    "tip": the terminal branch whose end is nearest in space.
    "segment": the vessel segment nearest to the unit, then the nearest tip
    within that segment's subtree, so that a unit stays in the territory of
    the vessel it lies against rather than jumping to a tip of another
    branch that happens to end close by.
    """
    tips = np.asarray(tips, int)
    ends = np.asarray([vt.end[k] for k in range(len(vt.parent))], float) * scale
    if mode == "tip":
        return tips[cKDTree(ends[tips]).query(U)[1]]
    if mode != "segment":
        raise ValueError(f"unknown assign_mode {mode!r}")
    par = np.asarray(vt.parent)
    starts = np.where(par[:, None] >= 0, ends[np.maximum(par, 0)], ends)
    # exact distance from each unit to each straight segment: project, clip
    # the projection to the segment, measure
    ab = ends - starts
    L2 = np.maximum((ab ** 2).sum(1), 1e-12)
    seg = np.empty(len(U), int)
    for a0 in range(0, len(U), 256):
        u = U[a0:a0 + 256][:, None, :]
        t = np.clip(((u - starts[None]) * ab[None]).sum(-1) / L2[None], 0.0, 1.0)
        d2 = ((u - (starts[None] + t[..., None] * ab[None])) ** 2).sum(-1)
        seg[a0:a0 + 256] = np.argmin(d2, axis=1)
    kids = [[] for _ in range(len(par))]
    for k in range(1, len(par)):
        kids[par[k]].append(k)
    is_tip = np.zeros(len(par), bool); is_tip[tips] = True
    sub = {}
    def subtree_tips(k):
        if k not in sub:
            out, stack = [], [k]
            while stack:
                j = stack.pop()
                if is_tip[j]:
                    out.append(j)
                stack.extend(kids[j])
            sub[k] = np.asarray(out, int)
        return sub[k]
    out = np.empty(len(U), int)
    for s_ in np.unique(seg):
        idx = np.flatnonzero(seg == s_)
        cand = subtree_tips(int(s_))
        if len(cand) == 0:
            cand = tips
        out[idx] = cand[cKDTree(ends[cand]).query(U[idx])[1]]
    return out


def closure_state(P, q, vt, q_set, tried, no_room, P_initial):
    """Where the progenitor pool stands, read from its distribution over tips.

    Returns the share of the initial pool that is spent, held at tips whose
    supply is at or below the threshold (supply holds it), held at tips above
    the threshold whose placements were all refused (space holds it), and
    held at tips above the threshold that can still nucleate. A mean supply
    per unit hides which of these applies; the shares do not.
    """
    n = min(len(q), len(tried))
    held = np.array([P.get(k, 0.0) for k in range(n)], float)
    qq = np.asarray(q, float)[:n]
    t = np.asarray(tried, float)[:n]; r = np.asarray(no_room, float)[:n]
    below = held[qq <= q_set].sum()
    above = qq > q_set
    blocked = above & (t > 0) & (r >= t)
    tot0 = max(P_initial, 1e-9)
    held_tot = held.sum()
    return dict(spent=float(max(1.0 - held_tot / tot0, 0.0)),
                supply_held=float(below / tot0),
                space_held=float(held[blocked].sum() / tot0),
                active=float(held[above & ~blocked].sum() / tot0))


def closure_by_trials(P, q, vt, q_set, U, sdf_scaled, scale, spread, lam, nu,
                      P_initial, n_trials=40, seed=12345):
    """Where the pool stands at closure, with room tested rather than inferred.

    For every tip that holds progenitors and has supply above the threshold,
    n_trials placements are drawn with an independent generator, around the
    tip as nucleation draws them, and the share that would land inside the
    organ and at least lam from every unit estimates its placement
    probability p. Its pool counts as held by room in proportion 1 - p and as
    active in proportion p. Pool at tips at or below the threshold is held by
    supply. The residual nucleation capacity is
    sum nu P (q - q_set)+ p. This does not depend on what one stochastic day
    happened to try.
    """
    rng = np.random.default_rng(seed)
    par = np.asarray(vt.parent)
    n = len(par)
    held = np.array([P.get(k, 0.0) for k in range(n)], float)
    # tips born on the last day have no supply yet: they take their parent's
    qq = np.empty(n)
    qq[:len(q)] = np.asarray(q, float)[:n]
    for k in range(len(q), n):
        qq[k] = qq[par[k]]
    tree = cKDTree(U) if len(U) else None
    supply = held[qq <= q_set].sum(); room = active = lam_eff = 0.0
    for i in np.flatnonzero((qq > q_set) & (held > 0)):
        def _trial(n):
            c = vt.end[i] * scale + rng.normal(0, spread, (n, 3))
            ok = np.array([sdf_scaled(*x) < 0 for x in c])
            if tree is not None and ok.any():
                ok[ok] = tree.query(c[ok])[0] >= lam
            return ok
        ok = _trial(n_trials)
        if not ok.any():
            # no success in n_trials: draw nine times more before reading zero,
            # so that a true probability of 0.02 is read as zero with chance
            # 0.98^400 rather than 0.98^40
            ok = np.concatenate([ok, _trial(9 * n_trials)])
        p = float(ok.mean())
        room += held[i] * (1 - p); active += held[i] * p
        lam_eff += nu * held[i] * (qq[i] - q_set) * p
    tot0 = max(P_initial, 1e-9)
    return dict(spent=float(max(1.0 - held.sum() / tot0, 0.0)),
                supply_held=float(supply / tot0), space_held=float(room / tot0),
                active=float(active / tot0), residual_rate=float(lam_eff))


def grow(sdf, R, *, days, rng, nu, q_set, rho, delta, P0, dP=1.0,
         lam=7.4, s_rate=0.0, tau_set=None,
         N_scale=None, R_unit=None, units_per_branch=None, m0=None, m_leak=None, xi=None, kappa=None, g0=None,
         tau_g=None, eta=None, b_cost=None, shrink=None,
         tau_set_fraction=None, sprout_reach=None,
         sprout=None, sprout_room=None, sprout_load=None, junction_cost=None, vascular_fraction=None, crowding=None, l_root=None,
         l_sprout=None, nucleation_spread=None,
         max_units=None, max_vessels=None,
         attractor_density=None, attractor_reach=None, attractor_kill=None,
         root_site="deepest", assign_mode="tip",
         nucleation="poisson_capped", substeps=1, turnover="independent",
         renewal=0.0,
         intervene_day=None, intervene=None, intervene_factor=0.8):
    """Grow one organ. `sdf(x, y, z)` is negative inside the domain.

    The identity parameters (nu, q_set, rho, delta, P0) are required; every
    other parameter defaults to config/params.yaml, section `closed`.
    """
    N_scale = float(_p("N_scale", N_scale)); m0 = _p("m0", m0)
    m_leak = _p("m_leak", m_leak); xi = _p("xi", xi); kappa = _p("kappa", kappa)
    g0 = _p("g0", g0); tau_g = _p("tau_g", tau_g); eta = _p("eta", eta)
    # One number sets three exponents. b is the metabolic cost exponent of
    # the vessel wall; Murray's b = 1 is the special case where cost is blood
    # volume, and measured branching exponents of 2.7 to 3.0 correspond to
    # b between 0.7 and 1.0, with b ~ 0.9 where flow is pulsatile.
    #     branching   r_p^x = sum r_c^x    with  x = b + 2
    #     pressure    tau_set ~ p^(b/2)
    #     daughter    r_c / r_p = N^(-1/x) for an N-way split
    b_cost = CFG["cost_exponent_b"] if b_cost is None else b_cost
    branch_exp = b_cost + 2.0
    pressure_exp = b_cost / 2.0
    shrink = _p("shrink", shrink)
    tau_set_fraction = _p("tau_set_fraction", tau_set_fraction)
    # a symmetric bifurcation under the same cost law
    daughter_ratio = 2.0 ** (-1.0 / branch_exp)
    sprout_reach = _p("sprout_reach", sprout_reach); sprout = _p("sprout", sprout)
    sprout_room = _p("sprout_room", sprout_room)
    sprout_load = _p("sprout_load", sprout_load)
    junction_cost = _p("junction_cost", junction_cost)
    attractor_density = _p("attractor_density", attractor_density)
    attractor_reach = _p("attractor_reach", attractor_reach)
    attractor_kill = _p("attractor_kill", attractor_kill)
    r_max_factor = _p("r_max_factor", None)     # widest calibre, in units of r_eq
    sprout_floor = _p("sprout_floor", None)     # a sprout starts this far above r_crit

    def _worth_splitting(vt, i, n_here):
        """Does reaching further buy more than the junction costs?

        A branch does not buy conductance. Two daughters at the cost-law
        ratio conduct 2 * 2^(-4/(b+2)) of the parent, which is less than one:
        splitting always costs conductance. What it buys is reach, the tissue
        a terminal branch can serve, and a tip already carrying n units is
        serving n / units_per_tip times what one branch is meant to.

        So the comparison is between the tissue that would be relieved and
        the tissue the junction takes. Both are volumes, which is what makes
        them comparable without a new constant: the junction occupies
        junction_cost * r^3, and relieving the overload frees the units above
        what one branch serves, each of which occupies about lam^3.
        """
        r = float(vt.radius[i])
        excess = n_here - units_per_tip
        if excess <= 0:
            return False
        relieved = excess * lam ** 3 / scale ** 3   # in tree units
        cost = junction_cost * r ** 3
        return bool(relieved > cost)
    vascular_fraction = _p("vascular_fraction", vascular_fraction)
    crowding = _p("crowding", crowding); l_root = _p("l_root", l_root)
    l_sprout = _p("l_sprout", l_sprout)
    nucleation_spread = _p("nucleation_spread", nucleation_spread)
    max_units = int(_p("max_units", max_units))
    max_vessels = int(_p("max_vessels", max_vessels))

    box = np.array([R * 1.35, R * 1.25, R * 0.85])
    # N_scale used to do two jobs: it set the resistance of one unit, which is
    # intensive, and, divided by the vessel cap, how many units a terminal
    # branch serves, also intensive but built from two extensive numbers.
    # Scaling the organ therefore changed the physics of a single unit. They
    # are now two latent parameters of their own. Their defaults reproduce
    # the calibrated reference organ exactly; what they mean physically is a
    # question for the fit, not something assumed here.
    R0 = float(R_unit) if R_unit is not None else N_scale
    # room the vessels may take before they compete with the parenchyma,
    # in tree units (lengths in lam, radii in r_root)
    # The organ's volume, from the fraction of a fixed grid over the sampling
    # box (-box, +box), whose volume is 8 bx by bz, that falls inside the
    # domain. The vascular budget is a fraction of it, in units of lam^3.
    # (An earlier version used 2 * prod(box / lam), a quarter of the box and
    # not the organ.) The same volume sets the attractor count.
    _g = np.stack(np.meshgrid(*[np.linspace(-1, 1, 25)] * 3, indexing="ij"), -1).reshape(-1, 3) * box
    _vol_dom = float(np.prod(2.0 * box)) * float(np.mean([sdf(*p) < 0 for p in _g]))
    vol_budget = vascular_fraction * _vol_dom / lam ** 3

    # Design point: a terminal branch carries the handful of units below it.
    # Its shear in series with that tissue peaks at r_c; the setpoint sits
    # below the peak and segments are born on the stable side of it.
    units_per_tip = (float(units_per_branch) if units_per_branch is not None
                     else max(N_scale / max(max_vessels * 0.5, 1.0), 1.0))
    R_seg = R0 / units_per_tip
    r_crit0 = (l_root / (3.0 * R_seg)) ** 0.25    # shear peak, see hemodynamics.perfuse
    _tau = lambda x: (dP / (l_root / x ** 4 + R_seg)) / x ** 3
    tau_peak = _tau(r_crit0)
    if tau_set is None:
        tau_set = tau_set_fraction * tau_peak
    try:
        r_eq = brentq(lambda x: _tau(x) - tau_set, r_crit0 * 1.001, 100.0)
    except ValueError:
        raise ValueError(
            f"no stable equilibrium radius: the shear setpoint {tau_set:.4g} is "
            f"not below the peak shear {tau_peak:.4g} at the critical radius "
            f"{r_crit0:.4g}; the setpoint, the root length or the unit "
            f"resistance is inconsistent") from None

    def inside(X, s=1.0):
        return sdf(X[:, 0] / s, X[:, 1] / s, X[:, 2] / s) < 0

    seed = rng.uniform(-1, 1, (20000, 3)) * box
    # The vascular root enters at a fixed anatomical site, the hilum, not at a
    # random point. Drawing it at random was the largest source of
    # between-replicate spread in the closed model: a root that lands on the
    # boundary leaves most nucleation candidates outside the organ, so too few
    # units fit for the first branch to pay for its junction, and the tree
    # stays at one vessel; a root deep inside branches explosively. Spearman
    # 0.72 between root depth and final vessel count over eight seeds.
    # The root is placed at the deepest point of the domain, found on a fixed
    # grid so that it is identical across replicates and scales with s. The
    # random draw above is kept so the generator stream is unchanged.
    _g = np.stack(np.meshgrid(*[np.linspace(-1, 1, 41)] * 3,
                              indexing="ij"), -1).reshape(-1, 3) * box
    _d = np.array([sdf(*p) for p in _g])
    s0 = _g[int(np.argmin(_d))]
    # Alternatives, for testing how much the results depend on this choice:
    # "surface" enters from the organ's surface, walking from the deepest
    # point along -y to the boundary and stepping one spacing back inside, as
    # a hilum does; "random" draws an interior point per replicate.
    if root_site == "surface":
        d = np.array([0.0, -1.0, 0.0]); p = s0.copy()
        while sdf(*(p + d * lam * 0.5)) < 0:
            p = p + d * lam * 0.5
        s0 = p - d * lam
    elif root_site == "random":
        s0 = seed[inside(seed)][0]
    elif root_site != "deepest":
        raise ValueError(f"unknown root_site {root_site!r}")
    scale = lam                       # tree coordinates -> domain coordinates
    vt = VascularTree(s0 / scale, s0 / scale + np.array([0.0, l_root, 0.0]), r_eq)

    U = np.zeros((0, 3))              # unit positions
    m = np.zeros(0)                   # unit masses
    tip_of = np.zeros(0, int)         # which tip feeds each unit
    # Progenitors are counted as whole cells. Each day a tip's progenitors
    # differentiate into units (one progenitor per unit placed), divide with
    # probability rho and are lost with probability delta, all as discrete
    # events; a branching tip hands each of its progenitors to the new branch
    # with probability one half. With rho = 0 the balance
    # P0 = units made + progenitors left + progenitors lost holds exactly in
    # every realisation, not only on average (tests/test_smoke.py).
    P = {0: int(round(P0))}           # progenitor count per tip
    Ptot0 = float(P[0])
    lost = 0                          # progenitors lost, cumulative
    born = 0                          # progenitors born, cumulative
    renewed = 0                       # progenitors returned after making a unit
    tr = Trace()
    leak = m_leak ** xi / R0          # what a tip carries before its first unit

    # ---- the mesenchymal signal: unserved tissue that calls the tree --------
    # Nothing else in the model drives growth into unoccupied space: units
    # appear only next to tips and tips sprout toward tissue that already
    # exists, so without this the organ grew as a colony around its root and
    # left 61 to 87 % of its domain empty. The signal is formalised as space
    # colonisation (Runions and colleagues, leaf venation): attractor points
    # spread through the domain at a fixed density stand for mesenchyme not
    # yet served, each pulls the nearest tip within reach, and each is removed
    # once a vessel or a unit comes close. Density, reach and kill distance
    # are intensive, so the signal scales with the organ.
    # _vol_dom: the organ volume computed with the vascular budget above
    n_attr = int(round(attractor_density * _vol_dom / lam ** 3))
    A = rng.uniform(-1, 1, (max(n_attr * 3, 10), 3)) * box
    A = A[np.array([sdf(*p) < 0 for p in A])][:n_attr]

    spacing = lam                     # least distance between units
    for day in range(1, days + 1):
        # A causal test of closure: on intervene_day the state is changed and
        # growth continues. "supply" lowers the threshold by intervene_factor,
        # which raises the hydraulic availability of every tip relative to it
        # without touching the pool; "room" shrinks the least spacing between
        # units by intervene_factor, which frees space without touching
        # perfusion or the pool.
        if intervene_day is not None and day == intervene_day:
            if intervene == "supply":
                q_set = q_set * intervene_factor
            elif intervene == "room":
                spacing = spacing * intervene_factor
            else:
                raise ValueError(f"unknown intervention {intervene!r}")
        s = 1.0 + s_rate * day        # the body grows the space, not the supply
        tips = vt.tips()
        # Each unit is fed by the terminal branch nearest to it. Units used to
        # stay with the tip they were born at, so when that tip branched to
        # relieve its load the load never moved: the parent kept every unit,
        # kept passing the junction test, and kept sprouting. The root ended
        # with 270 daughters and the tree was a star around it, with its whole
        # domain beyond reach. Reassigning by proximity is what makes a branch
        # actually relieve the tissue that justified it.
        if len(U) and len(tips):
            tip_of = _assign(vt, U, tips, scale, assign_mode)
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
                 r_min=vt.r_crit, r_max=r_max_factor * r_eq, shrink=shrink,
                 volume_budget=vol_budget, crowding=crowding, junction_cost=junction_cost)

        # ---- nucleation: one law, demographic noise ------------------------
        new_pos, new_tip = [], []
        no_room = np.zeros(len(vt), int)   # placements refused for lack of space
        tried = np.zeros(len(vt), int)     # placements attempted
        # the day can be split into substeps for the progenitor processes
        # (nucleation, division, loss), with supply held at its value for the
        # day; substeps = 1 is the reference
        dt = 1.0 / substeps
        for _sub in range(substeps):
            for i in active:
                Pi = P.get(i, 0.0)
                if Pi <= 0:
                    continue
                drive = max(q[i] - q_set, 0.0)
                lam_i = nu * Pi * drive * dt
                if nucleation == "binomial":
                    # each progenitor differentiates with hazard nu * drive over
                    # the step, independently; the cap is automatic
                    k = int(rng.binomial(int(Pi), 1.0 - np.exp(-nu * drive * dt))) if drive > 0 else 0
                elif nucleation == "poisson_capped":
                    k = int(rng.poisson(lam_i)) if lam_i > 0 else 0
                else:
                    raise ValueError(f"unknown nucleation {nucleation!r}")
                # A unit takes one cell from the competent pool, so a tip cannot
                # make more units than it holds. Whether that cell is spent
                # depends on the organ: see `renewal` below. The draw used to be uncapped: a tip with two
                # progenitors could place five units and clamp its pool at zero,
                # which made a pool-bound organ exceed its own budget by up to
                # 40 %, and by an amount that depended on the drive, that is on
                # the setpoint.
                k = min(k, int(Pi))
                tried[i] += k
                placed = 0
                for _ in range(k):
                    if len(U) + len(new_pos) >= max_units:
                        break
                    c = vt.end[i] * scale + rng.normal(0, nucleation_spread * lam, 3)
                    if sdf(c[0] / s, c[1] / s, c[2] / s) >= 0:
                        no_room[i] += 1
                        continue
                    if len(U) and np.min(((U - c) ** 2).sum(-1)) < spacing * spacing:
                        no_room[i] += 1
                        continue
                    if new_pos and np.min(((np.asarray(new_pos) - c) ** 2).sum(-1)) < spacing * spacing:
                        no_room[i] += 1
                        continue
                    new_pos.append(c); new_tip.append(i)
                    placed += 1
                # Whether making a unit spends a progenitor is the organ's
                # defining property. Where the differentiated cell is its own
                # progenitor, as mononuclear hepatocytes are in the liver, the
                # cell that forms a unit stays competent and the pool is not
                # drawn down: renewal = 1. Where the progenitor is a separate,
                # consumable population, as the cap mesenchyme is in the kidney
                # and the multipotent progenitors are in the pancreas, it is
                # spent: renewal = 0. The survivors then divide and die.
                returned = int(rng.binomial(placed, min(renewal, 1.0))) if renewal > 0 and placed else 0
                renewed += returned
                left = int(Pi) - placed + returned
                if turnover == "independent":
                    born_i = int(rng.binomial(left, min(rho * dt, 1.0))) if rho > 0 and left else 0
                    died = int(rng.binomial(left, min(delta * dt, 1.0))) if delta > 0 and left else 0
                elif turnover == "competing":
                    # one fate per progenitor per step: divide, be lost, or
                    # neither, with hazards rho and delta competing
                    if left and rho + delta > 0:
                        pe = 1.0 - np.exp(-(rho + delta) * dt)
                        born_i, died, _ = rng.multinomial(
                            left, [pe * rho / (rho + delta), pe * delta / (rho + delta),
                                   1.0 - pe])
                        born_i, died = int(born_i), int(died)
                    else:
                        born_i = died = 0
                else:
                    raise ValueError(f"unknown turnover {turnover!r}")
                lost += died
                born += born_i
                P[i] = left + born_i - died
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
                # A third reason to branch: a tip carrying far more tissue
                # than a terminal branch is meant to serve must divide,
                # whether or not that tissue is short of supply. Without it a
                # tip with ample supply and room to place units never
                # branches, so the organ can sit on one vessel forever.
                load = n_at[i] / units_per_tip if i < len(n_at) else 0.0
                if load > 1.0 and _worth_splitting(vt, i, n_at[i]):
                    for _ in range(int(rng.poisson(sprout_load * (load - 1.0)))):
                        if (len(vt) >= max_vessels or
                                vt.volume(junction_cost) >= vol_budget):
                            break
                        sel = np.flatnonzero(tip_of == i)
                        if not len(sel):
                            break
                        u = U[sel[rng.integers(len(sel))]]
                        b = vt.end[i]
                        nk = vt.add(i, b.copy(),
                                    b + (u / scale - b) * sprout_reach,
                                    max(vt.radius[i] * daughter_ratio,
                                        vt.r_crit[i] * sprout_floor))
                        P[nk] = int(rng.binomial(int(P.get(i, 0)), 0.5))
                        P[i] = int(P.get(i, 0)) - P[nk]
                        n_at = np.bincount(tip_of, minlength=len(vt))

                if no_room[i] > 0 and _worth_splitting(vt, i, n_at[i]):
                    # Respond to how crowded the tip is, which is the share of
                    # placements refused, not their number. The number grows
                    # with the progenitor pool, so a pool of a few hundred
                    # made hundreds of vessels in one day and every tree ran
                    # to the vessel cap.
                    k_r = int(rng.poisson(sprout_room * no_room[i]
                                          / max(tried[i], 1)))
                    for _ in range(k_r):
                        if (len(vt) >= max_vessels or
                                vt.volume(junction_cost) >= vol_budget):
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
                                    max(vt.radius[i] * daughter_ratio, vt.r_crit[i] * sprout_floor))
                        # the new tip inherits half the pool: it is where the
                        # next units will go
                        P[nk] = int(rng.binomial(int(P.get(i, 0)), 0.5))
                        P[i] = int(P.get(i, 0)) - P[nk]
                if n_at[i] == 0:
                    continue
                deficit = max(q_set - q[i], 0.0)
                # The junction test stays on this path too. Removing it was
                # tried on the reasoning that tissue short of supply should call
                # a vessel whatever the tip's load; once the mesenchymal signal
                # let the organ reach its supply bound, deficit was positive at
                # every tip and the ungated path ran the tree to its cap (6000
                # vessels against 773 with it off, at the same filling).
                k_v = (int(rng.poisson(sprout * deficit))
                       if deficit > 0 and _worth_splitting(vt, i, n_at[i])
                       else 0)
                for _ in range(k_v):
                    # the same stop as every other sprouting rule
                    if len(vt) >= max_vessels or vt.volume(junction_cost) >= vol_budget:
                        break
                    sel = np.flatnonzero(tip_of == i)
                    if not len(sel):
                        break
                    u = U[sel[rng.integers(len(sel))]]
                    b = vt.end[i]
                    end = b + (u / scale - b) * sprout_reach
                    nk = vt.add(i, b.copy(), end,
                                max(vt.radius[i] * daughter_ratio, vt.r_crit[i] * sprout_floor))
                    d_old = ((U[sel] - b * scale) ** 2).sum(-1)
                    d_new = ((U[sel] - end * scale) ** 2).sum(-1)
                    moved = sel[d_new < d_old]
                    tip_of[moved] = nk
                    frac = len(moved) / len(sel)
                    P[nk] = int(rng.binomial(int(P.get(i, 0)), frac))
                    P[i] = int(P.get(i, 0)) - P[nk]
                    sel = np.flatnonzero(tip_of == i)
                    if not len(sel):
                        break
            n_at = np.bincount(tip_of, minlength=len(vt))

        # ---- tips extend into unserved tissue (mesenchymal signal) --------
        if len(A):
            tips_now = vt.tips()
            ends = np.asarray([vt.end[k] for k in tips_now]) * scale
            d, near = cKDTree(ends).query(A)
            pulled = d < attractor_reach * lam
            if pulled.any():
                for t_idx in np.unique(near[pulled]):
                    if (len(vt) >= max_vessels or
                            vt.volume(junction_cost) >= vol_budget):
                        break
                    i = tips_now[int(t_idx)]
                    pts = A[pulled & (near == t_idx)]
                    v = pts - ends[t_idx]
                    v /= np.maximum(np.linalg.norm(v, axis=1, keepdims=True), 1e-9)
                    step = v.mean(0)
                    nrm = np.linalg.norm(step)
                    if nrm < 1e-6:
                        continue
                    b = vt.end[i]
                    end = b + step / nrm * l_sprout
                    if sdf(*(end * scale / s)) >= 0:
                        continue
                    _rc = np.atleast_1d(vt.r_crit)
                    rc_i = _rc[i] if i < len(_rc) else _rc[min(vt.parent[i], len(_rc) - 1)]
                    nk = vt.add(i, b.copy(), end,
                                max(vt.radius[i] * daughter_ratio, rc_i * sprout_floor))
                    # extending, not branching: the progenitor niche travels
                    # with the tip, so the front does not starve
                    P[nk] = int(P.get(i, 0))
                    P[i] = 0
            # tissue reached by a vessel or a unit is served
            served = cKDTree(np.asarray(vt.endpoints()) * scale).query(A)[0] < attractor_kill * lam
            if len(U):
                served |= cKDTree(U).query(A)[0] < attractor_kill * lam
            A = A[~served]
            n_at = np.bincount(tip_of, minlength=len(vt))

        Ptot = float(sum(P.values()))
        # q was read at the start of the day, before any sprout was added
        live = [i for i in active if n_at[i] > 0] or tips
        qbar = float(np.mean(q[live]))
        tr.t.append(day); tr.N.append(len(U)); tr.P.append(Ptot); tr.lost.append(lost); tr.born.append(born); tr.renewed.append(renewed)
        tr.q.append(qbar); tr.C.append(C); tr.mass.append(float(m.sum()))
        tr.vessels.append(len(vt))
        tr.Q.append(float(vt.flow[0]) if len(vt.flow) else 0.0)
        tr.bound.append(closure_state(P, q, vt, q_set, tried, no_room, Ptot0))

    # final state per node, for analyses that need where supply and pool sit
    tr.closure = closure_by_trials(P, q, vt, q_set, U, lambda x, y, z: sdf(x / s, y / s, z / s),
                                   scale, nucleation_spread * lam, spacing, nu, Ptot0)
    _tips = vt.tips()
    _diff = (float(np.mean(_assign(vt, U, _tips, scale, "tip") != _assign(vt, U, _tips, scale, "segment")))
             if len(U) and len(_tips) else 0.0)
    tr.final = dict(territory_diff=_diff, q=np.asarray(q, float).copy(),
                    P=np.array([P.get(k, 0.0) for k in range(len(q))], float),
                    tips=np.asarray(vt.tips(), int))
    return U, m, vt, tr
