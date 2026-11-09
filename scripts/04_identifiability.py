#!/usr/bin/env python
"""Figure 5: which parameters are recoverable, and from which observables."""
import json
import numpy as np
from organgrow.config import CFG, ROOT, seeds
from organgrow import grow, IDENTITY

OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
KIND = "liver"
DAYS = [10, 20, 30, 45, 60, 70]
PARAMS = ["nu", "c_set", "lam", "beta", "C0"]
STEP = 0.06
NOISE = 0.10          # relative noise assumed on each observable
TARGET = 0.20         # parameter to within +/- 20 %

cfg, eng = CFG["identity"][KIND], CFG["engine"]
BASE = dict(nu=cfg["nu"], c_set=cfg["c_set"], lam=cfg["lam"],
            beta=eng["beta"], C0=eng["C0"])


def run(p, seed):
    keep = {k: IDENTITY[KIND][k] for k in ("nu", "c_set", "lam")}
    IDENTITY[KIND].update(nu=p["nu"], c_set=p["c_set"], lam=p["lam"])
    try:
        U, m, T, V, tr = grow(KIND, seed=seed, days=eng["days"], snap_at=(),
                              beta=p["beta"], C0=p["C0"])
    finally:
        IDENTITY[KIND].update(keep)
    i = [d - 1 for d in DAYS]
    # the second tree is extended at a fixed rate, so its size carries no
    # information about the parameters; it is not used as an observable.
    return dict(N=np.asarray(tr.N, float)[i],
                cu=np.asarray(tr.cu, float)[i])


def average(p):
    o = [run(p, s) for s in seeds("fig5_identifiability")]
    return {k: np.mean([x[k] for x in o], axis=0) for k in o[0]}


def vector(r, keys):
    return np.concatenate([np.log(np.maximum(r[k], 1e-9)) for k in keys])


base = average(BASE)
plus, minus = {}, {}
for name in PARAMS:
    for sign, store in ((+1, plus), (-1, minus)):
        p = dict(BASE); p[name] = BASE[name] * (1 + sign * STEP)
        store[name] = average(p)
    print("  sensitivity of", name, flush=True)

SETS = [("Final N only", ["N"], True),
        ("N trajectory", ["N"], False),
        ("N + per-unit supply", ["N", "cu"], False)]
out = {}
for label, keys, final_only in SETS:
    J = np.array([(vector(plus[n], keys) - vector(minus[n], keys)) / (2 * STEP)
                  for n in PARAMS]).T
    if final_only:
        J = J[[i * len(DAYS) + len(DAYS) - 1 for i in range(len(keys))], :]
    sv = np.linalg.svd(J, compute_uv=False)
    rank = int((sv > max(sv[0] * 1e-8, 1e-10)).sum())
    se = np.sqrt(np.diag(np.linalg.pinv(J.T @ J / NOISE ** 2)))
    out[label] = dict(sv=(sv / sv[0]).tolist(), rank=rank,
                      cond=float(sv[0] / sv[rank - 1]) if rank else None,
                      se=se.tolist(),
                      nrep=np.ceil((se / TARGET) ** 2).tolist())
    print(f"{label:>24} rank={rank}/{len(PARAMS)} cond={out[label]['cond']}")

Jfull = np.array([(vector(plus[n], ["N", "cu"]) -
                   vector(minus[n], ["N", "cu"])) / (2 * STEP)
                  for n in PARAMS])
np.savez(OUT / "sens.npz", J=Jfull, params=np.array(PARAMS),
         keys=np.array(["N", "cu"]), days=np.array(DAYS))
# the engine truncates N to an integer, so finite differences carry
# truncation noise; report how much the condition number moves with the seed
conds = []
for s_ in seeds("fig5_identifiability"):
    J1 = np.array([(vector(run(dict(BASE, **{n: BASE[n]*(1+STEP)}), s_), ["N","cu"]) -
                    vector(run(dict(BASE, **{n: BASE[n]*(1-STEP)}), s_), ["N","cu"]))
                   / (2*STEP) for n in PARAMS]).T
    sv_ = np.linalg.svd(J1, compute_uv=False)
    rk_ = int((sv_ > max(sv_[0]*1e-8, 1e-10)).sum())
    if rk_:
        conds.append(float(sv_[0]/sv_[rk_-1]))
out["N + per-unit supply"]["cond_per_seed"] = conds
print(f"  condition number across seeds: "
      f"{np.median(conds):.0f} (range {min(conds):.0f} to {max(conds):.0f})")
json.dump(dict(params=PARAMS, sets=out), open(OUT / "ident.json", "w"), indent=1)
print("wrote results/ident.json and results/sens.npz")
