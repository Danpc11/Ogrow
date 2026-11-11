"""Flow, shear and calibre on a vascular tree, under a volume budget.

A tree has analytic resistance, so no linear system is needed: demand is
accumulated from the leaves to the root in one pass, and conductance is
combined from the leaves to the root in another. Both are O(E).

The three quantities that matter downstream:

    w_k   = r_k^4 / l_k                 Poiseuille conductance (Eq. 2)
    tau_k = f_k / r_k^3                 wall shear, up to a constant (Eq. 11)
    r_k   relaxes toward the radius at which tau_k equals its setpoint (Eq. 10)

Units are arbitrary and consistent; only ratios enter the model.
"""
import numpy as np


class VascularTree:
    """Segments with parent pointers, radii, flows and shear.

    Segment 0 is the root. Every other segment has exactly one parent, so
    children are implicit and both traversals are a single reversed loop over
    the segment index, which is created in topological order.
    """

    __slots__ = ("parent", "start", "end", "radius", "length",
                 "flow", "shear", "is_tip",
                 "pressure_in", "pressure_out", "conductance_total",
                 "r_crit", "below")

    def __init__(self, root_start, root_end, r0):
        self.parent = [-1]
        self.start = [np.asarray(root_start, float)]
        self.end = [np.asarray(root_end, float)]
        self.radius = [float(r0)]
        self.length = [max(float(np.linalg.norm(
            np.asarray(root_end, float) - np.asarray(root_start, float))), 1e-6)]
        self.flow = np.zeros(1)
        self.shear = np.zeros(1)
        self.is_tip = [True]
        self.pressure_in = np.zeros(1)
        self.pressure_out = np.zeros(1)
        self.conductance_total = 0.0
        self.r_crit = np.zeros(1)
        self.below = np.zeros(1)

    def __len__(self):
        return len(self.parent)

    def add(self, parent, start, end, radius):
        """Append a daughter segment and mark its parent as no longer a tip."""
        self.parent.append(int(parent))
        self.start.append(np.asarray(start, float))
        self.end.append(np.asarray(end, float))
        self.length.append(max(float(np.linalg.norm(
            np.asarray(end, float) - np.asarray(start, float))), 1e-6))
        self.radius.append(float(radius))
        self.flow = np.append(np.asarray(self.flow, float), 0.0)
        self.shear = np.append(np.asarray(self.shear, float), 0.0)
        self.is_tip.append(True)
        self.is_tip[int(parent)] = False
        return len(self.parent) - 1

    # ------------------------------------------------------------------ flow
    def pressures(self, p_in=1.0):
        """Pressure at each segment, from the root down.

        Shear-only adaptation is unstable: at a uniform setpoint r^3 is
        proportional to f, so conductance follows flow^(4/3), and with a fixed
        perfusion pressure the two reinforce each other without bound. The
        drop along the tree is what breaks the loop, because a segment's
        setpoint rises with the pressure it carries.
        """
        n = len(self.parent)
        f = np.asarray(self.flow)
        r = np.asarray(self.radius)
        l = np.asarray(self.length)
        w = r ** 4 / np.maximum(l, 1e-6)
        drop = f / np.maximum(w, 1e-12)
        p = np.zeros(n)
        par = np.asarray(self.parent)
        for k in range(1, n):                 # parents precede daughters
            p[k] = p[par[k]] - drop[k]
        # rescale so the root-to-tip span is [p_in, p_in/span]: the absolute
        # drop depends on units, only the profile along the tree is physical
        lo = p.min()
        if lo < p[0]:
            p = (p - lo) / (p[0] - lo)
        else:
            p = np.ones(n)
        return p_in * np.maximum(p, 1e-2)

    def solve(self, tip_demand):
        """Accumulate demand to the root, then set flow and shear per segment.

        `tip_demand` maps segment index to the demand drawn at that tip. The
        pass is leaves-to-root because segments are created in topological
        order, so a reversed loop visits every daughter before its parent.
        """
        n = len(self.parent)
        f = np.zeros(n)
        for k, d in tip_demand.items():
            f[k] += d
        par = np.asarray(self.parent)
        for k in range(n - 1, 0, -1):
            f[par[k]] += f[k]
        r = np.asarray(self.radius)
        self.flow = f
        # wall shear for Poiseuille flow: tau = 4 mu Q / (pi r^3), constants
        # folded away since only the ratio to the setpoint is used
        self.shear = f / np.maximum(r, 1e-6) ** 3
        return f

    # ----------------------------------------------------------- conductance
    def conductance(self):
        """Effective conductance from the root, by series-parallel reduction.

        Walking leaves to root, each segment's downstream conductance is its
        own in series with the parallel sum of its daughters'.
        """
        n = len(self.parent)
        r = np.asarray(self.radius)
        l = np.asarray(self.length)
        own = r ** 4 / np.maximum(l, 1e-6)
        below = np.zeros(n)                     # parallel sum of daughters
        eff = np.zeros(n)
        par = np.asarray(self.parent)
        for k in range(n - 1, -1, -1):
            eff[k] = own[k] if below[k] <= 0 else 1.0 / (1.0 / own[k] +
                                                         1.0 / below[k])
            if par[k] >= 0:
                below[par[k]] += eff[k]
        return float(eff[0])

    # ------------------------------------------------- pressure-driven flow
    def perfuse(self, terminal_conductance, dP=1.0):
        """Flow and pressure from a fixed perfusion pressure.

        Demand does not set the flow. The tree's geometry does: conductance is
        reduced leaves-to-root, the root draws Q = dP * C, and that flow is
        split among daughters in proportion to their conductances. A branch
        therefore carries what its calibre allows, which is what makes the
        shear law non-trivial and Murray's relation a result rather than a
        restatement of the rule.

        `terminal_conductance` maps a tip index to the conductance of the
        tissue hanging from it, so growing tissue calls flow toward itself
        without dictating it.
        """
        n = len(self.parent)
        r = np.asarray(self.radius)
        l = np.asarray(self.length)
        own = r ** 4 / np.maximum(l, 1e-6)
        par = np.asarray(self.parent)

        below = np.zeros(n)
        for k, c in terminal_conductance.items():
            below[int(k)] += float(c)
        eff = np.zeros(n)
        for k in range(n - 1, -1, -1):              # leaves to root
            # a tip with nothing below it is a dead end, not an open drain:
            # its conductance is zero, otherwise it would discharge to zero
            # pressure, its shear would rise with its radius and the calibre
            # loop would run away
            eff[k] = 0.0 if below[k] <= 0 else 1.0 / (1.0 / own[k] +
                                                      1.0 / below[k])
            if par[k] >= 0:
                below[par[k]] += eff[k]
        # the critical radius of each segment in series with everything
        # below it: shear peaks at (3 l / R_below)^(1/4); under it the
        # calibre feedback is positive, so it is the floor the segment may
        # not cross (Eq. in CLOSED_MODEL.md). Dead ends get the root's value.
        r_crit = np.where(below > 0, (3.0 * l * np.maximum(below, 1e-12)) ** 0.25, 0.0)
        r_crit[below <= 0] = 0.0
        self.r_crit = r_crit
        self.below = below

        kids = [[] for _ in range(n)]
        for k in range(1, n):
            kids[par[k]].append(k)

        f = np.zeros(n)
        p = np.zeros(n)
        pout = np.zeros(n)
        f[0] = dP * eff[0]
        p[0] = dP
        for k in range(n):                          # root to leaves
            drop = f[k] / max(own[k], 1e-12)
            pout[k] = p[k] - drop
            cs = kids[k]
            if cs:
                tot = sum(eff[c] for c in cs)
                for c in cs:
                    f[c] = f[k] * (eff[c] / max(tot, 1e-12))
                    p[c] = pout[k]
        self.flow = f
        self.shear = f / np.maximum(r, 1e-6) ** 3
        self.pressure_in = p
        self.pressure_out = pout             # what the tissue at a tip sees
        self.conductance_total = float(eff[0])
        return f, p, float(eff[0])

    # -------------------------------------------------------------- calibres
    def volume(self, junction_cost=0.0):
        """Tissue the vasculature takes, lumens plus junctions.

        Lumen volume alone is degenerate in the branching number: splitting a
        vessel in two costs nothing, so nothing stops a tip dividing. A
        junction occupies tissue of its own, of order r^3, and that steric
        term is what makes the network finite (Murray's cost function is flat
        in the branching number; a wall cost with a junction penalty is not).
        """
        r = np.asarray(self.radius)
        l = np.asarray(self.length)
        lumen = float(np.pi * np.sum(r ** 2 * l))
        if junction_cost:
            par = np.asarray(self.parent)
            n_kids = np.bincount(par[par >= 0], minlength=len(r))
            branching = n_kids >= 2
            lumen += junction_cost * float(np.sum(r[branching] ** 3))
        return lumen

    def relax(self, tau_set, rate=0.12, r_min=0.25, r_max=60.0,
              volume_budget=None, crowding=1.0, pressure_exp=0.0, p_in=1.0,
              pressures=None, shrink=0.0):
        """Adapt each radius toward its shear setpoint, stably.

        `shrink` is the constant regression tendency that makes the fixed
        point a balance rather than a runaway; `volume_budget` raises the
        setpoint as the vessels fill the space the parenchyma needs.
        """
        f = np.asarray(self.flow)
        r = np.asarray(self.radius)
        if pressure_exp and pressures is not None:
            # the setpoint rises with transmural pressure, from the cost of
            # mural mass at constant homeostatic wall stress. These are the
            # pressures the flow solve returned, not a rescaled profile.
            tau_set = tau_set * np.maximum(pressures / max(p_in, 1e-9),
                                           1e-3) ** pressure_exp
        if volume_budget:
            # mural mass costs space the parenchyma needs, so a crowded organ
            # runs its vessels at higher shear rather than widening them
            phi = min(self.volume() / volume_budget, 0.999)
            tau_set = tau_set / max(1.0 - phi, 1e-3) ** crowding
        ts = np.maximum(np.asarray(tau_set, float), 1e-9)
        # Shear-only adaptation with a power-law target is unstable: the
        # radius that satisfies the setpoint scales as f^(1/3) and f rises
        # with r^4, so the two reinforce each other. The published stable
        # form (Pries, Reglin & Secomb) responds to the logarithm of the
        # shear ratio, which bounds the step, and carries a shrinking
        # tendency so that a vessel carrying no flow regresses.
        tau = np.where(f > 0, f / np.maximum(r, 1e-6) ** 3, 0.0)
        # A segment whose flow stops falling would otherwise widen without
        # bound: less demand lowers its shear, which widens it, which lowers
        # the shear again. The shrinking tendency has to be able to win, so
        # the signal is capped on the widening side.
        signal = np.clip(np.log(np.maximum(tau, 1e-12) / ts), -4.0, 4.0) - shrink
        r_new = r * (1.0 + rate * np.tanh(signal))
        self.radius = list(np.clip(r_new, r_min, r_max))
        return np.asarray(self.radius)

    # ------------------------------------------------------------- geometry
    def endpoints(self):
        return np.asarray(self.end)

    def tips(self):
        return [k for k, t in enumerate(self.is_tip) if t]

    def radii(self):
        return np.asarray(self.radius)

    def r_eff(self):
        """Shear-weighted effective radius, <r^5>/<r^4> (Eq. 11).

        Reported rather than used: it is the quantity a homogenised treatment
        would need, and it depends on the whole calibre distribution.
        """
        r = np.asarray(self.radius)
        return float((r ** 5).mean() / max((r ** 4).mean(), 1e-12))
