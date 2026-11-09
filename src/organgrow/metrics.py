"""Distance in invariant space, and the irreducible floor."""
from itertools import combinations
import numpy as np


def _matrix(rows, keys, log_keys=()):
    return np.array([[np.log(r[k]) if k in log_keys else r[k] for k in keys]
                     for r in rows], float)


def epsilon(v, mu, sd):
    """Distance of one vector from a cohort, normalised per axis."""
    v, mu, sd = map(np.asarray, (v, mu, sd))
    use = sd > 1e-12
    return float(np.sqrt(np.mean(((v[use] - mu[use]) / sd[use]) ** 2)))


def floor(rows, keys, log_keys=()):
    """Mean pairwise distance within a cohort: the irreducible floor.

    Returns (mean, 5th pct, 95th pct, n_pairs). Axes with zero spread are
    dropped, since they carry no information about between-sample distance.
    """
    V = _matrix(rows, keys, log_keys)
    sd = V.std(0, ddof=1)
    use = sd > 1e-12
    e = np.array([np.sqrt(np.mean(((V[i][use] - V[j][use]) / sd[use]) ** 2))
                  for i, j in combinations(range(len(V)), 2)])
    return float(e.mean()), float(np.percentile(e, 5)), float(np.percentile(e, 95)), len(e)


def effective_axes(rows, keys, log_keys=()):
    """How many independent axes the invariant set really has."""
    V = _matrix(rows, keys, log_keys)
    sd = V.std(0, ddof=1)
    use = sd > 1e-12
    C = np.corrcoef(V[:, use].T)
    r = np.abs(C[np.triu_indices(use.sum(), 1)]).mean()
    return float(use.sum() / (1 + (use.sum() - 1) * r)), float(r)
