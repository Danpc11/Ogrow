#!/usr/bin/env python
"""Figure 6: how organ variance splits into systemic, positional and contingent.

Size is a property of the individual; shape is not. The left-right comparison
separates what the individual fixes from what the organ's surroundings fix.
"""
import json
from itertools import combinations
import numpy as np
from scipy import stats
from organgrow.config import CFG, ROOT

OUT = ROOT / "results"
src = OUT / "trusted_split.json"
if not src.exists():
    raise SystemExit("run scripts/07_trusted.py first")

D = json.load(open(src))
keys = CFG["human"]["invariants_ct"]
LOG = {"volume_cm3", "surface_cm2"}
ks = sorted(D)
V = np.array([[np.log(D[k][i]) if i in LOG else D[k][i] for i in keys] for k in ks])
sd = V.std(0, ddof=1); use = sd > 1e-12

idx = {}
for i, k in enumerate(ks):
    idx.setdefault(D[k]["pid"], []).append(i)
pairs = [(v[0], v[1]) for v in idx.values() if len(v) == 2]
dist = lambda i, j: float(np.sqrt(np.mean(((V[i][use] - V[j][use]) / sd[use]) ** 2)))
intra = np.array([dist(i, j) for i, j in pairs])
inter = np.array([dist(i, j) for i, j in combinations(range(len(ks)), 2)
                  if D[ks[i]]["pid"] != D[ks[j]]["pid"]])

# permutation null for the ratio: shuffle which kidney belongs to whom
rng = np.random.default_rng(0)
null = []
for _ in range(2000):
    p = rng.permutation(len(ks))
    a = np.mean([dist(p[i], p[j]) for i, j in pairs])
    b = np.mean([dist(p[i], p[j]) for i, j in
                 list(combinations(range(len(ks)), 2))[:600]])
    null.append(b / a)
null = np.array(null)
obs = inter.mean() / intra.mean()

res = dict(kidneys=len(ks), donors=len(idx),
           intra=float(intra.mean()), inter=float(inter.mean()), ratio=float(obs),
           null_mean=float(null.mean()), null_sd=float(null.std(ddof=1)),
           p_permutation=float(max((null >= obs).mean(), 1 / len(null))))

# per-invariant split, and the systematic left-right difference
res["per_invariant"] = {}
for a, nm in enumerate(keys):
    x = V[:, a]
    within = np.mean([(x[i] - x[j]) ** 2 / 2 for i, j in pairs])
    L = np.array([D[ks[i]][nm] for i, _ in pairs])
    R = np.array([D[ks[j]][nm] for _, j in pairs])
    res["per_invariant"][nm] = dict(
        between_frac=float(max(1 - within / x.var(ddof=1), 0)),
        lr_rel_diff=float(100 * (L.mean() - R.mean()) / ((L.mean() + R.mean()) / 2)),
        n_left_greater=int((L > R).sum()), n_pairs=len(L),
        p_sign=float(stats.binomtest(int((L > R).sum()), len(L), 0.5).pvalue),
        p_wilcoxon=float(stats.wilcoxon(L, R).pvalue))

json.dump(res, open(OUT / "partition.json", "w"), indent=1)
print(f"{res['kidneys']} kidneys, {res['donors']} donors")
print(f"  within {res['intra']:.3f}  between {res['inter']:.3f}  "
      f"ratio {res['ratio']:.2f}  P={res['p_permutation']:.4f}")
for nm, v in res["per_invariant"].items():
    print(f"  {nm:>14}: between {100*v['between_frac']:5.1f}%  "
          f"L-R {v['lr_rel_diff']:+6.1f}%  P={min(v['p_sign'], v['p_wilcoxon']):.4f}")
