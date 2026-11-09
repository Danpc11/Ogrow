#!/usr/bin/env python
"""Figure 2: growth series and replicate ensembles for the three identities."""
import numpy as np
from organgrow.config import CFG, ROOT, seeds
from organgrow import grow

OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
SNAP = (10, 25, 45, 70)
d = CFG["disorder"]; e = CFG["engine"]

snaps, reps = {}, {}
for kind in CFG["organs"]["growth_series"]:
    U, m, T, V, tr = grow(kind, seed=CFG["seeds"]["fig2_representative"],
                          days=e["days"], snap_at=SNAP)
    snaps[kind] = tr.snapshots
    R = [grow(kind, seed=s, days=e["days"], snap_at=(),
              sigma_init=d["sigma_init"], sigma_param=d["sigma_param"])[4]
         for s in seeds("fig2_replicates")]
    reps[kind] = dict(N=np.array([r.N for r in R], float),
                      P=np.array([r.P for r in R], float),
                      hypo=np.array([r.hypo for r in R], float))
    print(f"{kind}: N={tr.N[-1]}  replicates={len(R)}", flush=True)

np.savez(OUT / "ens_main.npz",
         **{f"{k}_{f}": v[f] for k, v in reps.items() for f in v})
np.save(OUT / "ens_snapshots.npy", snaps, allow_pickle=True)
print("wrote results/ens_main.npz")

# reference cohort: invariant vectors of replicate organs, used by figures 4
# and by scripts/06_exchangeability.py
from organgrow.sim_invariants import invariants as _inv
ref = []
for s in seeds("fig2_replicates"):
    U, m, T, V, tr = grow("liver", seed=3000 + s, days=e["days"], snap_at=(),
                          sigma_init=d["sigma_init"], sigma_param=d["sigma_param"])
    if len(U) > 8:
        ref.append(_inv(U, m, T, V))
np.save(OUT / "ref_cohort.npy", np.asarray(ref))
print(f"reference cohort: {len(ref)} organs -> results/ref_cohort.npy")
