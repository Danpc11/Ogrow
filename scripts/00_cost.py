#!/usr/bin/env python
"""Figure 1e: cost of the flow solve, tree traversal vs sparse system."""
import json, time
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
from organgrow.config import ROOT

def _timeit(fn, *a):
    t = time.perf_counter(); fn(*a); return time.perf_counter() - t


def tree(depth):
    par, fr = [-1], [0]
    for _ in range(depth):
        nf = []
        for p in fr:
            for _ in range(2):
                par.append(p); nf.append(len(par) - 1)
        fr = nf
    return np.array(par), fr

def solve_tree(par, tips, Q=1.0):
    acc = np.zeros(len(par)); acc[tips] = Q / len(tips)
    for k in range(len(par) - 1, 0, -1):
        acc[par[k]] += acc[k]
    return acc[1:]

def solve_sparse(par, w, tips, Q=1.0):
    n = len(par); r, c, v = [], [], []
    for k in range(1, n):
        i, j, wij = par[k], k, w[k]
        r += [i, j, i, j]; c += [i, j, j, i]; v += [wij, wij, -wij, -wij]
    A = sp.coo_matrix((v, (r, c)), shape=(n, n)).tolil()
    b = np.zeros(n); b[0] = Q; b[tips] -= Q / len(tips)
    A[0, :] = 0; A[0, 0] = 1.0; b[0] = 0.0
    return spl.spsolve(A.tocsr(), b)

p0, t0 = tree(8); w0 = np.ones(len(p0))
for _ in range(3):                       # warm-up, discarded
    solve_tree(p0, t0); solve_sparse(p0, w0, t0)

out = []
for d in (8, 10, 12, 14, 16):
    par, tips = tree(d); w = np.ones(len(par))
    tt = np.median([_timeit(solve_tree, par, tips) for _ in range(7)])
    ss = (np.median([_timeit(solve_sparse, par, w, tips) for _ in range(3)])
          if len(par) <= 70000 else float("nan"))
    out.append(dict(E=len(par) - 1, tree=tt, sparse=ss))
    print(f"E={len(par)-1:>7}  tree={tt*1e3:8.3f} ms  sparse={ss*1e3:10.3f} ms")
json.dump(out, open(ROOT / "results" / "cost.json", "w"), indent=1)
