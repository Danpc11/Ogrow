"""Flow, shear and calibre on a vascular tree.

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
                 "flow", "shear", "is_tip")

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

    # -------------------------------------------------------------- calibres
    def relax(self, tau_set, rate=0.12, r_min=0.25, r_max=60.0):
        """Move each radius toward the one that would put shear at setpoint.

        With tau proportional to f / r^3, the radius satisfying the setpoint is
        (f / tau_set)^(1/3); the step is a fraction of the gap, so a segment
        whose flow has just risen widens over several days rather than at once.
        """
        f = np.asarray(self.flow)
        r = np.asarray(self.radius)
        r_target = np.where(f > 0, (f / max(tau_set, 1e-9)) ** (1.0 / 3.0), r_min)
        r_new = r + rate * (r_target - r)
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
