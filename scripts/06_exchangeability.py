#!/usr/bin/env python
"""Is a simulated organ an exchangeable member of its cohort?

The expected pairwise distance in a self-normalised invariant space is
sqrt(2) for ANY exchangeable sample, so its value carries no information
about the cohort. What it does provide is a test: a simulation whose
distance to the cohort matches the cohort's own internal distance cannot
be told apart from a member of it.
"""
import json
from itertools import combinations
import numpy as np
from organgrow.config import ROOT

OUT = ROOT / "results"


def pairwise(V):
    sd = V.std(0, ddof=1); use = sd > 1e-12
    return np.array([np.sqrt(np.mean(((V[i][use] - V[j][use]) / sd[use]) ** 2))
                     for i, j in combinations(range(len(V)), 2)])


def null_reference(n, d, reps=400, seed=0):
    """Distance among exchangeable random vectors: the reference, not a result."""
    rng = np.random.default_rng(seed)
    return float(np.mean([pairwise(rng.normal(0, 1, (n, d))).mean()
                          for _ in range(reps)]))


if __name__ == "__main__":
    ref = OUT / "ref_cohort.npy"
    if not ref.exists():
        raise SystemExit("no reference cohort; run scripts/01_ensembles.py")
    V = np.load(ref)
    obs = pairwise(V)
    null = null_reference(*V.shape)
    res = dict(n=int(V.shape[0]), d=int(V.shape[1]),
               observed=float(obs.mean()),
               exchangeable_reference=null,
               sqrt2=float(np.sqrt(2)),
               note="the reference is a property of the normalisation, "
                    "not of the cohort; use it only to test exchangeability")
    json.dump(res, open(OUT / "exchangeability.json", "w"), indent=1)
    print(f"  observed {res['observed']:.3f}   "
          f"exchangeable reference {null:.3f}   sqrt(2) = {res['sqrt2']:.3f}")
