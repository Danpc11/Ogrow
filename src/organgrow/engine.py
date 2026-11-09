"""Generative organ-growth engine.

Local rules only. Organ identity enters through one parameter: how the
fabrication window closes.

Closure modes
-------------
saturation   capacity saturation      dN/dt = nu * max(C/N - c_set, 0)
exhaustion   progenitor exhaustion    dN/dt = iota * P ;  dP/dt = -(k) P
pool         limiting pool            dN/dt drawn from a fixed budget
"""
from dataclasses import dataclass, field
import numpy as np
from scipy.spatial import cKDTree

# --------------------------------------------------------------- domains
def sdf(kind, x, y, z):
    """Signed distance to the tissue domain boundary; negative inside."""
    if kind == "liver":
        a = np.hypot(np.hypot((x - 30) / 1.45, (y + 5) / 0.88), z / 0.62) - 105
        b = np.hypot(np.hypot((x + 92) / 0.95, (y + 12) / 0.70), z / 0.55) - 62
        d = np.minimum(a, b)
        notch = np.hypot(np.hypot((x + 22) / 0.62, (y + 96) / 1.05), z / 1.6) - 54
        d = np.maximum(d, -notch)
        return np.maximum(d, y - 62)
    if kind == "kidney":
        e = np.hypot(np.hypot(x / 0.78, y / 1.30), z / 0.70) - 100
        cut = np.hypot(np.hypot((x + 92) / 0.80, y / 0.92), z / 1.30) - 72
        return np.maximum(e, -cut)
    if kind == "lung":
        # wedge with a concave hilum on the medial face and a domed base
        e = np.hypot(np.hypot(x / 0.92, (y + 10) / 1.45), z / 0.70) - 115
        hil = np.hypot(np.hypot((x + 86) / 0.70, (y + 20) / 0.85), z / 1.25) - 58
        return np.maximum(np.maximum(e, -hil), -(y + 150))
    # pancreas: tapering body curving from head to tail
    t = np.clip(x / 165.0, -1, 1)
    cy = 34 * np.sin(t * 1.35) - 6
    rad = 50 * (1 - 0.52 * (t * 0.5 + 0.5))
    return np.maximum(np.hypot(y - cy, z * 1.35) - rad, np.abs(x) - 165)


try:                                  # parameters live in config/params.yaml
    from .config import CFG
    IDENTITY = {k: dict(v) for k, v in CFG["identity"].items()}
    _ENG = CFG["engine"]
except Exception:                     # standalone use without the config file
    IDENTITY = {
        "liver":    dict(mode="saturation", R=190, lam=7.4, nu=400.0, c_set=0.055),
        "kidney":   dict(mode="exhaustion", R=150, lam=6.0, iota=0.165, decay=0.952),
        "pancreas": dict(mode="pool",       R=165, lam=6.6, draw=0.075),
        "lung":     dict(mode="exhaustion", R=175, lam=6.2, iota=0.110, decay=0.930,
                         vessels="bundle"),
    }
    _ENG = dict(days=70, n_attractors=9000, max_units=4200, beta=2.2,
                g0=0.30, tau_g=12.0, C0=14.0, alpha_allometric=1.155,
                tau_body=16.0)


@dataclass
class Trace:
    t: list = field(default_factory=list)
    N: list = field(default_factory=list)
    P: list = field(default_factory=list)
    cu: list = field(default_factory=list)
    hypo: list = field(default_factory=list)
    snapshots: dict = field(default_factory=dict)


def root_supply(day, C0=None, tau_b=None, alpha=None):
    """Root conductance grows with the body; exponent from measured allometry."""
    C0 = _ENG['C0'] if C0 is None else C0
    tau_b = _ENG['tau_body'] if tau_b is None else tau_b
    alpha = _ENG['alpha_allometric'] if alpha is None else alpha
    return C0 * (1.0 + day / tau_b) ** alpha


def grow(kind, seed=0, days=None, P0=620.0, C0=None, snap_at=(10, 25, 45, 70),
         beta=None, g0=None, tau_g=None, max_units=None, sigma_init=0.0, sigma_param=0.0, n_seed=12,
         n_att=9000, rules=('shear', 'exclusion', 'perunit')):
    """Two independent sources of variability, which the theory separates.

    sigma_init   spread of the INITIAL CONDITION (seeded units). This is what
                 capacity-limited nucleation can correct.
    sigma_param  spread of the PER-INDIVIDUAL PARAMETERS (perfusion setpoint,
                 progenitor pool). This is quenched disorder and sets the
                 floor CV_inf >= CV(c_set); it cannot be corrected.
    """
    """Run one replicate. Returns (units, primary tree, vessels, trace)."""
    days = _ENG['days'] if days is None else days
    C0 = _ENG['C0'] if C0 is None else C0
    beta = _ENG['beta'] if beta is None else beta
    g0 = _ENG['g0'] if g0 is None else g0
    tau_g = _ENG['tau_g'] if tau_g is None else tau_g
    max_units = _ENG['max_units'] if max_units is None else max_units
    n_att = _ENG['n_attractors'] if n_att == 9000 else n_att
    rng = np.random.default_rng(seed)
    cfg = dict(IDENTITY[kind])
    if sigma_param > 0:
        if "c_set" in cfg:
            cfg["c_set"] *= float(np.exp(rng.normal(0, sigma_param)))
        P0 = P0 * float(np.exp(rng.normal(0, sigma_param)))
        C0 = C0 * float(np.exp(rng.normal(0, sigma_param * 0.35)))
    n_seed_i = max(1, int(round(n_seed * np.exp(rng.normal(0, sigma_init)))))
    R, lam = cfg["R"], cfg["lam"]
    lam2 = lam * lam

    def inside(p):
        return sdf(kind, p[0], p[1], p[2]) < 0

    # ---- attractors: unoccupied territory to be colonised ----
    box = np.array([R * 1.35, R * 1.25, R * 0.85])

    def sample_attractors(n):
        """Batched rejection sampling: narrow domains have a low accept rate."""
        out = []
        got = 0
        while got < n:
            P = rng.uniform(-1, 1, (max(n, 4000), 3)) * box
            keep = sdf(kind, P[:, 0], P[:, 1], P[:, 2]) < 0
            if keep.any():
                out.append(P[keep]); got += int(keep.sum())
        return np.vstack(out)[:n]

    N_ATT = int(n_att)
    att = sample_attractors(N_ATT)

    # anatomical entry point: portal hilum (liver), renal hilum (kidney),
    # pancreatic head. Fall back to any interior point.
    SEED = {"liver": (-5.0, -10.0, 0.0),
            "kidney": (-30.0, 0.0, 0.0),
            "pancreas": (-120.0, -6.0, 0.0),
            "lung": (-70.0, -10.0, 0.0)}
    s0 = np.array(SEED[kind], dtype=float)
    if not inside(s0):
        P = rng.uniform(-1, 1, (20000, 3)) * box
        ok = sdf(kind, P[:, 0], P[:, 1], P[:, 2]) < 0
        if ok.any():
            s0 = P[ok][0]
    tips = [s0.copy()]
    U, m, T, V = [], [], [(s0, s0)], [(s0, s0)]
    for _ in range(n_seed_i):
        for _try in range(60):
            d = rng.normal(size=3); d /= np.linalg.norm(d)
            cand = s0 + d * rng.uniform(lam * 1.2, lam * 3.0)
            if not inside(cand):
                continue
            if U and np.min(((np.asarray(U) - cand) ** 2).sum(-1)) < lam * lam:
                continue
            U.append(cand); m.append(0.12)
            T.append((s0.copy(), cand.copy())); tips.append(cand.copy())
            break
    kill = lam * 1.25
    infl = lam * 16.0

    P = float(P0)
    tr = Trace()

    for day in range(1, days + 1):
        N = max(len(U), 1)
        Cd = beta * float(np.sum(m)) + 6.0 if m else 6.0
        C = 1.0 / (1.0 / root_supply(day, C0) + 1.0 / Cd)
        cu = C / N

        # ---- identity: budget of new units for today ----
        if cfg["mode"] == "saturation":
            # per-unit law (theory):  dN/dt = nu (C/N - c_set)
            # total-capacity law:     dN/dt = nu (C - c_set N)
            if "perunit" in rules:
                drive = cu - cfg["c_set"]
            else:
                drive = C - cfg["c_set"] * N
            k = int(min(max(cfg["nu"] * drive, 0.0), 90))
        elif cfg["mode"] == "exhaustion":
            k = int(min(cfg["iota"] * P, 90))
            P *= cfg["decay"]
        else:
            k = int(min(cfg["draw"] * P, 90))
            P = max(P - k * 1.6, 0.0)
        k = max(min(k, max_units - len(U)), 0)

        # ---- space colonisation: tips grow toward unclaimed attractors ----
        placed = 0
        guard = 0
        utree = cKDTree(np.asarray(U)) if U else None
        while placed < k and len(att) and guard < 40:
            guard += 1
            Tarr = np.asarray(tips)
            dmin, nearest = cKDTree(Tarr).query(att, k=1)
            claimed = dmin < infl
            if not claimed.any():
                claimed = dmin < infl * 3.0
                if not claimed.any():
                    break
            moved = []
            for ti in np.unique(nearest[claimed]):
                sel = claimed & (nearest == ti)
                dirs = att[sel] - tips[ti]
                v = dirs.sum(0)
                nv = np.linalg.norm(v)
                if nv < 1e-9:
                    continue
                cand = tips[ti] + v / nv * (lam * 2.1)
                if not inside(cand):
                    continue
                if ("exclusion" in rules and utree is not None
                        and utree.query(cand[None, :], k=1)[0][0] < lam):
                    continue
                T.append((tips[ti].copy(), cand.copy()))
                U.append(cand); m.append(0.09)
                moved.append(cand)
                placed += 1
                if placed >= k:
                    break
            if not moved:
                break
            tips.extend(moved)
            Uarr = np.asarray(U)
            utree = cKDTree(Uarr)
            att = att[utree.query(att, k=1)[0] > kill]
            if len(att) < N_ATT * 0.25:          # replenish: sampling must not
                fresh = sample_attractors(N_ATT // 2)   # cap growth
                far = utree.query(fresh, k=1)[0] > kill
                att = np.vstack([att, fresh[far]]) if far.any() else att

        if not U:
            tr.t.append(day); tr.N.append(0); tr.P.append(P); tr.cu.append(cu); tr.hypo.append(0.0)
            continue

        # ---- second tree ----
        # "demand": vessels grow toward poorly supplied units (liver, kidney,
        # pancreas, where the supply tree is also the exchange tree).
        # "bundle": vessels follow the primary tree segment by segment, as the
        # pulmonary artery follows the bronchus (lung).
        if len(V) < 1200:
            if cfg.get("vessels") == "bundle":
                for seg in T[len(V):len(V) + 5]:
                    a, b = np.asarray(seg[0]), np.asarray(seg[1])
                    off = rng.normal(0, lam * 0.35, 3)      # runs alongside
                    V.append((a + off, b + off))
            else:
                Varr = np.asarray([v[1] for v in V])
                for _ in range(5):
                    u = U[rng.integers(len(U))]
                    j = int(cKDTree(Varr).query(u[None, :], k=1)[1][0])
                    b = Varr[j]
                    V.append((b.copy(), b + (u - b) * 0.55))
                    Varr = np.asarray([v[1] for v in V])

        # ---- unit growth limited by local supply ----
        Varr = np.asarray([v[1] for v in V])
        Uarr = np.asarray(U)
        _dv = cKDTree(Varr).query(Uarr, k=1)[0]
        sup = cu * (np.exp(-_dv / 58.0) if 'shear' in rules else 1.0)
        marr = np.asarray(m)
        marr = marr + g0 * np.exp(-day / tau_g) * marr * (1 - marr / np.maximum(0.85 * sup, 1e-4))
        m = list(np.maximum(marr, 0.02))

        tr.t.append(day)
        tr.N.append(len(U))
        tr.P.append(P if cfg["mode"] != "saturation" else np.nan)
        tr.cu.append(cu)
        tr.hypo.append(float(np.mean(sup < 0.5 * cu)))
        if day in snap_at:
            tr.snapshots[day] = (np.asarray(U).copy(), np.asarray(m).copy(),
                                 [tuple(map(np.asarray, e)) for e in T])

    return np.asarray(U), np.asarray(m), T, V, tr
