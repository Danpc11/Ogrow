#!/usr/bin/env python
"""Runs behind the closed-model figures of the paper (Figures 2, 4 and 5).

Every run is cached under results/paper/, so the script can be stopped and
resumed. The three regimes and their parameters are those of
scripts/15_transmission.py; nu = 8 throughout.

    regimes   one organ per regime, with daily traces           (Figure 2)
    cohort    eight livers with a common rate, daily traces     (Figure 4)
    rate      eight livers with an individual rate,
              sigma(log nu) = 1.0                               (Figure 4)
    murray    bifurcation ratios of the regime trees            (Figure 5)
    snapshot  geometry of one liver at closure                  (Figure 1)

Usage:  python scripts/17_paper_runs.py [JOB ...]     default: all
"""
import sys

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT, CFG
from organgrow.closed import grow

OUT = ROOT / "results" / "paper"
OUT.mkdir(parents=True, exist_ok=True)
S, DAYS, NU, NREP, SEED0 = 0.6, 140, 8.0, 8, 300
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
    np.savez(f, closure=np.array([b["spent"], b["supply_held"], b["space_held"], b["active"]]),
             N=np.asarray(tr.N, float), q=np.asarray(tr.q, float),
             P=np.asarray(tr.P, float), vessels=np.asarray(tr.vessels, float),
             murray=np.asarray(murray), nu=nu)
    print(f"  {tag}: N={len(U)}", flush=True)


def job_regimes():
    for reg, p in REGIMES.items():
        run(f"regime_{reg}", SEED0, NU, **p)


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
    tips = fin["tips"][fin["tips"] < len(fin["q"])]
    np.savez(f, U=np.asarray(U, float), starts=starts, ends=ends,
             radius=np.asarray(vt.radius, float), tips=tips,
             q_tip=fin["q"][tips], P_tip=fin["P"][tips], q_set=REGIMES["liver"]["q_set"])
    print(f"  snapshot_liver: N={len(U)}, vessels={len(par)}", flush=True)


JOBS = dict(regimes=job_regimes, cohort=job_cohort, rate=job_rate, snapshot=job_snapshot)

if __name__ == "__main__":
    for j in (sys.argv[1:] or JOBS):
        if j in JOBS:
            JOBS[j]()
    print("done")
