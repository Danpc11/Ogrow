#!/usr/bin/env python
"""Runs behind the closed-model figures of the paper (Figures 2, 4 and 5).

Every run is cached under results/paper/, so the script can be stopped and
resumed. The three regimes and their parameters are those of
scripts/15_transmission.py; nu = 8 throughout.

    regimes   one organ per regime, with daily traces           (Figure 2)
    reps      the same three regimes over NSEED seeds, for the
              replicate means of how each pool is spent          (Figure 2)
    cohort    eight livers with a common rate, daily traces     (Figure 4)
    rate      eight livers with an individual rate,
              sigma(log nu) = 1.0                               (Figure 4)
    murray    bifurcation ratios of the regime trees            (Figure 5)
    snapshot  geometry of one liver at closure                  (Figure 1)

Usage:  python scripts/17_paper_runs.py [JOB ...]     default: all
"""
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT, CFG
from organgrow.closed import grow

OUT = ROOT / "results" / "paper"
OUT.mkdir(parents=True, exist_ok=True)
S, DAYS, NU, NREP, SEED0 = 0.6, 140, 8.0, 8, 300
NSEED = 24      # replicates behind the regime means, as elsewhere in this work
REGIMES = {
    "liver":    dict(P0=2600.0, rho=0.10, delta=0.10, q_set=0.6, renewal=1.0),
    "kidney":   dict(P0=300.0,  rho=0.00, delta=0.05, q_set=0.6, renewal=0.0),
    "pancreas": dict(P0=60.0,   rho=0.00, delta=0.00, q_set=0.6, renewal=0.0),
}
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=6000)


def run(tag, seed, nu, **p):
    f = OUT / f"{tag}.npz"
    if f.exists():
        return
    rng = np.random.default_rng(seed)
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S),
                        190 * S, days=DAYS, rng=rng, nu=nu, **FIX, **p)
    r = np.asarray(vt.radius, float)
    par = np.asarray(vt.parent)
    x = CFG["cost_exponent_b"] + 2.0
    kids = {}
    for k in range(1, len(par)):
        kids.setdefault(int(par[k]), []).append(k)
    murray = [float(r[a] ** x / sum(r[c] ** x for c in cs))
              for a, cs in kids.items() if len(cs) >= 2]
    b = tr.closure
    fa = tr.fate
    np.savez(f, closure=np.array([b["spent"], b["supply_held"], b["space_held"], b["active"]]),
             fate=np.array([fa["consumed"], fa["lost"], fa["born"], fa["left"]]),
             N=np.asarray(tr.N, float), q=np.asarray(tr.q, float),
             P=np.asarray(tr.P, float), vessels=np.asarray(tr.vessels, float),
             murray=np.asarray(murray), nu=nu)
    print(f"  {tag}: N={len(U)}", flush=True)


def job_regimes():
    for reg, p in REGIMES.items():
        run(f"regime_{reg}", SEED0, NU, **p)


def _warm(a):
    """Populate the cache for one run; each file is written by one worker only."""
    tag, seed, p = a
    run(tag, seed, NU, **p)


def _tag(reg, k):
    """Seed 0 of every regime is the trace run, so it is not run twice."""
    return f"regime_{reg}" if k == 0 else f"regime_{reg}_k{k}"


def job_reps():
    """The three regimes over NSEED seeds, so Figure 2a is a mean and not one draw.

    How a pool is spent varies between replicates, and the panel used to show a
    single realisation. The per-run files are the same ones job_regimes writes,
    so seed 0 is shared and only the rest are new.
    """
    grid = [(_tag(reg, k), SEED0 + k, p)
            for reg, p in REGIMES.items() for k in range(NSEED)
            if not (OUT / f"{_tag(reg, k)}.npz").exists()]
    if grid:
        n = max(1, min(int(os.environ.get("OGROW_JOBS", 0)) or (os.cpu_count() or 1),
                       len(grid)))
        print(f"running {len(grid)} simulations on {n} worker(s)", flush=True)
        if n > 1:
            with ProcessPoolExecutor(max_workers=n) as pool:
                for _ in pool.map(_warm, grid):
                    pass
        else:
            for g in grid:
                _warm(g)

    keys = ("consumed", "lost", "born", "left")
    ckeys = ("spent", "supply_held", "space_held", "active")
    out = {}
    for reg in REGIMES:
        fa = np.array([np.load(OUT / f"{_tag(reg, k)}.npz")["fate"] for k in range(NSEED)])
        cl = np.array([np.load(OUT / f"{_tag(reg, k)}.npz")["closure"] for k in range(NSEED)])
        N = np.array([np.load(OUT / f"{_tag(reg, k)}.npz")["N"][-1] for k in range(NSEED)])

        def ms(A, names):
            return {n: dict(mean=float(A[:, i].mean()),
                            se=float(A[:, i].std(ddof=1) / np.sqrt(len(A))))
                    for i, n in enumerate(names)}

        out[reg] = dict(seeds=NSEED, fate=ms(fa, keys), closure=ms(cl, ckeys),
                        N=dict(mean=float(N.mean()),
                               se=float(N.std(ddof=1) / np.sqrt(len(N)))))
        f = out[reg]["fate"]
        print(f"  {reg:9s} N={out[reg]['N']['mean']:7.0f} "
              + "  ".join(f"{k} {f[k]['mean']:.3f}±{f[k]['se']:.3f}" for k in keys),
              flush=True)
    json.dump(out, open(OUT / "regime_reps.json", "w"), indent=1)


def job_cohort():
    for k in range(NREP):
        run(f"cohort_k{k}", SEED0 + k, NU, **REGIMES["liver"])


def job_rate():
    z = np.random.default_rng(99).normal(size=NREP)
    for k in range(NREP):
        run(f"rate_k{k}", SEED0 + k, NU * float(np.exp(1.0 * z[k])), **REGIMES["liver"])


def job_snapshot():
    """Geometry of one liver at closure, for the schematic of Figure 1."""
    f = OUT / "snapshot_liver.npz"
    if f.exists():
        return
    rng = np.random.default_rng(SEED0)
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S),
                        190 * S, days=DAYS, rng=rng, nu=NU, **FIX, **REGIMES["liver"])
    par = np.asarray(vt.parent)
    ends = np.asarray([vt.end[k] for k in range(len(par))], float)
    starts = np.asarray([vt.end[p] if p >= 0 else vt.end[k] for k, p in enumerate(par)], float)
    fin = tr.final
    tips = fin["tips"]
    np.savez(f, U=np.asarray(U, float), starts=starts, ends=ends,
             radius=np.asarray(vt.radius, float), tips=tips,
             q_tip=fin["q"][tips], P_tip=fin["P"][tips], q_set=REGIMES["liver"]["q_set"])
    print(f"  snapshot_liver: N={len(U)}, vessels={len(par)}", flush=True)


JOBS = dict(regimes=job_regimes, reps=job_reps, cohort=job_cohort, rate=job_rate,
            snapshot=job_snapshot)

if __name__ == "__main__":
    for j in (sys.argv[1:] or JOBS):
        if j in JOBS:
            JOBS[j]()
    print("done")
