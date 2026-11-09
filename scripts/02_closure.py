#!/usr/bin/env python
"""Figure 3: does the closure mode destroy or merely dilute variance."""
import numpy as np
from organgrow.config import CFG, ROOT, seeds
from organgrow import grow

OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
d, e = CFG["disorder"], CFG["engine"]

main = {}
for kind in CFG["organs"]["closure_test"]:
    main[kind] = np.array([grow(kind, seed=s, days=e["days"], snap_at=(),
                                sigma_init=d["sigma_init"], sigma_param=0.0)[4].N
                           for s in seeds("fig3_main")], float)
    A = main[kind]
    print(f"{kind}: sd(N) day1={A[:,0].std(ddof=1):.2f} day70={A[:,-1].std(ddof=1):.2f}",
          flush=True)
np.savez(OUT / "f3_main.npz", **main)

sig = CFG["sweeps"]["sigma_param"]
for kind in ["liver", "kidney"]:
    A = np.zeros((len(sig), len(list(seeds("fig3_sweep"))), e["days"]))
    for i, sp in enumerate(sig):
        for j, s in enumerate(seeds("fig3_sweep")):
            A[i, j] = grow(kind, seed=s, days=e["days"], snap_at=(),
                           sigma_init=d["sigma_init"], sigma_param=sp)[4].N
        print(f"  {kind} sigma_param={sp} done", flush=True)
    np.savez(OUT / f"f3_sweep_{kind}.npz", sig=np.array(sig), N=A)
print("wrote results/f3_*.npz")
