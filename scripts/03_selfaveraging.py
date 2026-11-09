#!/usr/bin/env python
"""Figure 4: which observables concentrate with organ size, and which do not.

Additive quantities should decay as N^-a with a > 0; properties of the
particular realisation should not decay at all.
"""
import json
import numpy as np
from scipy.spatial import cKDTree
from organgrow.config import CFG, ROOT, seeds
from organgrow import grow, sdf

OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
KIND = "liver"
DAYS = CFG["sweeps"]["scaling_days"]
d, e = CFG["disorder"], CFG["engine"]
R = CFG["identity"][KIND]["R"]


def probe_points(rng, n=1500):
    """Fixed probe grid inside the domain, shared by every replicate."""
    out = []
    while len(out) < n:
        P = rng.uniform(-1, 1, (4000, 3)) * np.array([R * 1.1, R, R * 0.6])
        keep = sdf(KIND, P[:, 0], P[:, 1], P[:, 2]) < 0
        out.extend(P[keep])
    return np.asarray(out[:n])


def observables(U, m, V):
    Varr = np.asarray([v[1] for v in V])
    sup = np.exp(-cKDTree(Varr).query(U, k=1)[0] / 58.0)
    return dict(N=float(len(m)), mass=float(np.sum(m)),
                mean_unit=float(np.mean(m)),
                dissipation=float(np.sum(np.asarray(m) / np.maximum(sup, 1e-6))))


def main(sigma_param):
    rng = np.random.default_rng(0)
    grid = probe_points(rng)
    rows = []
    for nd in DAYS:
        obs, labels, dens = [], [], []
        for s in seeds("fig4_scaling"):
            U, m, T, V, tr = grow(KIND, seed=s, days=nd, snap_at=(),
                                  sigma_init=d["sigma_init"], sigma_param=sigma_param)
            if len(U) < 8:
                continue
            obs.append(observables(U, m, V))
            Varr = np.asarray([v[1] for v in V])
            labels.append(cKDTree(Varr).query(grid, k=1)[1])      # which vessel
            dens.append(cKDTree(U).query(grid, k=min(8, len(U)))[0].mean(1))
        if len(obs) < 4:
            continue
        row = dict(days=nd, N=float(np.mean([o["N"] for o in obs])))
        for k in obs[0]:
            v = np.array([o[k] for o in obs])
            row["cv_" + k] = float(v.std(ddof=1) / v.mean())
        row["vessel_identity"] = float(np.mean(
            [np.mean(labels[i] != labels[j])
             for i in range(len(labels)) for j in range(i + 1, len(labels))]))
        row["density_discord"] = float(1.0 - np.mean(
            [np.corrcoef(dens[i], dens[j])[0, 1]
             for i in range(len(dens)) for j in range(i + 1, len(dens))]))
        rows.append(row)
        print(f"  N={row['N']:7.0f}  cv_N={row['cv_N']:.4f}  "
              f"vessel_identity={row['vessel_identity']:.3f}", flush=True)
    return rows


if __name__ == "__main__":
    for sp, tag in ((0.0, "nodisorder"), (d["sigma_param"], "disorder")):
        print(f"sigma_param = {sp}")
        json.dump(main(sp), open(OUT / f"scaling_{tag}.json", "w"), indent=1)
    print("wrote results/scaling_*.json")
