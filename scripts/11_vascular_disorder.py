#!/usr/bin/env python
"""The vascular tree as a source of quenched disorder.

Three measurements, in the order they settle the question:

  gamma   does supply per unit rise or fall as units are added?
  eta     is the dispersion topological, or does it come from calibre relaxation?
  time    does the disorder diverge, or settle at a stationary level?
"""
import json
import numpy as np
from organgrow.config import CFG, ROOT, seeds
from organgrow import grow
import organgrow.engine as E

OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
d = CFG["disorder"]
DAYS = CFG["engine"]["days"]
RATES = CFG["sweeps"]["calibre_rate"]

res = {}

# -- gamma: the exponent in C ~ N^gamma, measured along each trajectory ------
print("exponent of supply against unit count")
res["gamma"] = {}
for kind in CFG["organs"]["growth_series"]:
    g = []
    for s in seeds("vascular_gamma"):
        U, m, T, V, tr = grow(kind, seed=s, days=DAYS, snap_at=())
        N = np.asarray(tr.N, float)
        C = np.asarray(tr.cu, float) * N
        ok = (N > 20) & (C > 0)
        if ok.sum() > 6:
            g.append(float(np.polyfit(np.log(N[ok]), np.log(C[ok]), 1)[0]))
    res["gamma"][kind] = dict(mean=float(np.mean(g)), sd=float(np.std(g, ddof=1)),
                              values=g)
    print(f"  {kind:>9} gamma = {np.mean(g):.3f} +/- {np.std(g, ddof=1):.3f}")

# -- eta: topology alone, or calibre relaxation? ----------------------------
print("\ndispersion of supply against calibre relaxation rate")
base = E._ENG["calibre_rate"]
res["eta_sweep"] = []
try:
    for rate in RATES:
        E._ENG["calibre_rate"] = rate
        C = []
        for s in seeds("vascular_eta"):
            # pancreas: unit count fixed by its budget, no seeded disorder, so
            # anything that disperses must come from the tree
            U, m, T, V, tr = grow("pancreas", seed=s, days=DAYS, snap_at=(),
                                  sigma_init=0.0, sigma_param=0.0)
            C.append(float(np.asarray(tr.cu, float)[-1] * np.asarray(tr.N, float)[-1]))
        C = np.array(C)
        # the dispersion itself is uncertain: a bootstrap over replicates
        rng_b = np.random.default_rng(0)
        boot = [np.std(x, ddof=1) / np.mean(x) for x in
                (rng_b.choice(C, len(C), replace=True) for _ in range(400))]
        res["eta_sweep"].append(dict(rate=float(rate),
                                     cv=float(C.std(ddof=1) / C.mean()),
                                     cv_lo=float(np.percentile(boot, 5)),
                                     cv_hi=float(np.percentile(boot, 95)),
                                     C=[float(x) for x in C]))
        print(f"  eta = {rate:<5} CV(C) = {C.std(ddof=1)/C.mean():.4f}  "
              f"[{np.percentile(boot, 5):.4f}, {np.percentile(boot, 95):.4f}]")
finally:
    E._ENG["calibre_rate"] = base      # whatever happens, the engine is restored

# -- time: does it diverge or settle? ---------------------------------------
print("\ndisorder over time")
R = [grow("liver", seed=s, days=100, snap_at=(),
          sigma_init=d["sigma_init"], sigma_param=0.0)[4] for s in seeds("vascular_time")]
res["over_time"] = []
for day in range(5, 101, 5):
    N = np.array([tr.N[day - 1] for tr in R], float)
    C = np.array([tr.cu[day - 1] * tr.N[day - 1] for tr in R], float)
    res["over_time"].append(dict(day=day, N=float(N.mean()),
                                 cv_N=float(N.std(ddof=1) / N.mean()),
                                 cv_C=float(C.std(ddof=1) / C.mean())))
for r in res["over_time"][::4]:
    print(f"  day {r['day']:>3}  N = {r['N']:>6.0f}  "
          f"CV(N) = {r['cv_N']:.4f}  CV(C) = {r['cv_C']:.4f}")

# -- Murray's law, measured not imposed -------------------------------------
from organgrow.hemodynamics import VascularTree
rng = np.random.default_rng(0)
vt = VascularTree([0, 0, 0], [0, 10, 0], CFG["engine"]["root_radius"])
for _ in range(200):
    j = int(rng.integers(len(vt)))
    b = vt.end[j]
    vt.add(j, b, b + rng.normal(0, 5, 3), vt.radius[j] * 0.79)
dem = {k: 1.0 for k in vt.tips()}
for _ in range(60):
    vt.solve(dem)
    vt.relax(CFG["engine"]["tau_set"], rate=CFG["engine"]["calibre_rate"])
r = vt.radii()
par = np.asarray(vt.parent)
kids = {}
for k in range(1, len(vt)):
    kids.setdefault(int(par[k]), []).append(k)
ratio = [float(r[p] ** 3 / sum(r[c] ** 3 for c in cs))
         for p, cs in kids.items() if len(cs) >= 2]
res["murray"] = dict(median=float(np.median(ratio)), n=len(ratio),
                     values=ratio)
print(f"\nMurray ratio r_p^3 / sum r_c^3 = {np.median(ratio):.2f} "
      f"over {len(ratio)} bifurcations (law predicts 1)")

json.dump(res, open(OUT / "vascular_disorder.json", "w"), indent=1)
print("\nwrote results/vascular_disorder.json")
