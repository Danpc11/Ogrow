"""Invariant vector used to score a simulated organ."""
import numpy as np
from scipy.spatial import cKDTree

def invariants(U, m, T, V):
    U = np.asarray(U); m = np.asarray(m)
    Varr = np.asarray([v[1] for v in V])
    seg = np.array([np.linalg.norm(np.asarray(b) - np.asarray(a)) for a, b in T if np.any(np.asarray(b) != np.asarray(a))])
    d1 = cKDTree(U).query(U, k=2)[0][:, 1]              # nearest-neighbour spacing
    dv = cKDTree(Varr).query(U, k=1)[0]                 # unit-to-vessel distance
    hull = np.ptp(U, axis=0).prod()
    return np.array([
        np.log(len(m)),                      # unit number
        np.log(m.mean()),                    # mean unit mass
        m.std(ddof=1) / m.mean(),            # unit-mass dispersion
        np.log(m.sum() / max(hull, 1e-9)),   # mass per bounding volume
        np.log(d1.mean()),                   # mean unit spacing
        d1.std(ddof=1) / d1.mean(),          # spacing dispersion
        np.log(dv.mean()),                   # mean supply distance
        np.log(seg.mean()) if len(seg) else 0.0,   # mean segment length
    ])
