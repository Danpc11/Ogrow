#!/usr/bin/env python
"""Does the binary property separate the regimes on its own?

The three regimes of this paper differ in several things at once: whether the
cell that forms a unit is spent (`renewal`), how large the seeded pool is
(`P0`), and how the pool turns over, which is a division rate and a loss rate
together (`rho`, `delta`).

    liver-like     renewal 1   P0 2600   rho 0.10  delta 0.10   balanced
    kidney-like    renewal 0   P0  300   rho 0     delta 0.05   decaying
    pancreas-like  renewal 0   P0   60   rho 0     delta 0      static

So the paper cannot credit the binary property with the separation without
showing what the others do. This varies each independently and reports, for
every cell, where the pool ends up: spent, held at tips whose supply is at or
below the setpoint (supply holds it), held at tips with supply to spare but no
room (space holds it), or still able to nucleate. Those four shares are what
"regime" means here, so a regime that survives a change of pool and turnover
at fixed renewal is a regime the binary property owns.

Turnover is one axis and not two. Division and loss were separated once, with
the division rate held at the liver's 0.10 while the loss rate swept down to
zero, and that was wrong twice over. It named cells after the kidney and the
pancreas that were not those regimes, since both divide at 0, and it paired a
dividing pool with no loss, which grows without bound: a refused nucleation
does not spend a progenitor, so at rho 0.10 and delta 0 the pool grows by a
factor of 620,000 over the window and the model attempts of the order of a
billion placements that are nearly all refused. Those cells cost half an hour
apiece and answered a question about no organ. Sweeping the published pairings
instead puts the three regimes at their own corners and leaves every cell
stationary in the pool.

Two readings matter and the script prints both. Within a column of the
factorial, does changing renewal alone move the closure? And across a row, does
changing the pool or the turnover at fixed renewal move it as much?

Usage:  python scripts/34_regimes.py

    OGROW_SEEDS   replicates, default 24
    OGROW_JOBS    workers
    OGROW_ADOPT_CACHE   see organgrow.runcache
"""
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow
from organgrow.runcache import fingerprint, cached

OUT = ROOT / "results" / "regimes"
OUT.mkdir(parents=True, exist_ok=True)
S, NU, DAYS = 0.6, 8.0, 140
NSEED = int(os.environ.get("OGROW_SEEDS", 24))
SEEDS = tuple(range(300, 300 + NSEED))
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=6000,
           q_set=0.6)
#: the values each axis takes, the published regimes among them
RENEWAL = (0.0, 1.0)
POOL = (60.0, 300.0, 2600.0)
#: division and loss travel together, as the three published sources pair them
TURNOVER = (("static", 0.00, 0.00),
            ("decaying", 0.00, 0.05),
            ("balanced", 0.10, 0.10))
TURN = {n: (r, d) for n, r, d in TURNOVER}
#: the elasticity of the final unit count to the seeded pool is the paper's
#: sharpest contrast between regimes, and it is measured at each regime's own
#: P0, so a renewing source at P0 2600 is compared with a spent one at 60. This
#: measures it again with the pool and the turnover matched, which is the only
#: way to say whether the binary property owns the contrast. H is the
#: perturbation scripts/19_elasticity.py uses.
H = 0.3
MATCHED = ((2600.0, "balanced"), (300.0, "decaying"), (60.0, "static"))

#: where the published regimes sit in that grid
NAMED = {(1.0, 2600.0, "balanced"): "liver-like",
         (0.0, 300.0, "decaying"): "kidney-like",
         (0.0, 60.0, "static"): "pancreas-like"}


def key(renewal, P0, turn, seed):
    return f"w{renewal}_P{P0:.0f}_{turn}_s{seed}"


def params(renewal, P0, turn):
    rho, delta = TURN[turn]
    return dict(renewal=float(renewal), P0=float(P0), rho=float(rho),
                delta=float(delta))


def fp_of(renewal, P0, turn, seed):
    return fingerprint(seed=seed, S=S, NU=NU, days=DAYS,
                       params={**FIX, **params(renewal, P0, turn)})


def run(renewal, P0, turn, seed):
    """One construct grown to arrest under one corner of the factorial."""
    f = OUT / f"{key(renewal, P0, turn, seed)}.json"
    fp = fp_of(renewal, P0, turn, seed)
    hit = cached(f, fp)
    if hit is not None:
        return hit
    rho, delta = TURN[turn]
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S), 190 * S,
                        days=DAYS, rng=np.random.default_rng(seed), nu=NU,
                        **FIX, **params(renewal, P0, turn))
    c = tr.closure
    r = dict(renewal=renewal, P0=P0, turnover=turn, rho=rho, delta=delta,
             seed=seed, N=len(U), mass=float(m.sum()), spent=float(c["spent"]),
             supply_held=float(c["supply_held"]),
             space_held=float(c["space_held"]), active=float(c["active"]),
             fingerprint=fp)
    json.dump(r, open(f, "w"))
    return r


def _warm(a):
    run(*a)


def jobs(grid):
    """Run whatever is not already cached, in parallel."""
    todo = [g for g in grid
            if cached(OUT / f"{key(*g)}.json", fp_of(*g)) is None]
    if not todo:
        return
    k = max(1, min(int(os.environ.get("OGROW_JOBS", 0)) or (os.cpu_count() or 1),
                   len(todo)))
    print(f"running {len(todo)} simulations on {k} worker(s). Every finished "
          f"run is cached,\nso an interrupted pass resumes where it stopped.",
          flush=True)
    t0, done = time.time(), 0
    step = max(1, len(todo) // 20)
    if k > 1:
        with ProcessPoolExecutor(max_workers=k) as pool:
            for _ in pool.map(_warm, todo):
                done += 1
                if done % step == 0 or done == len(todo):
                    print(f"  {done:>5d}/{len(todo)} after "
                          f"{time.time() - t0:.0f}s", flush=True)
    else:
        for g in todo:
            _warm(g)
            done += 1
            if done % step == 0 or done == len(todo):
                print(f"  {done:>5d}/{len(todo)} after "
                      f"{time.time() - t0:.0f}s", flush=True)


def cell(w, P0, turn):
    v = [run(w, P0, turn, s) for s in SEEDS]
    f = lambda k: float(np.mean([x[k] for x in v]))          # noqa: E731
    rho, delta = TURN[turn]
    return dict(renewal=w, P0=P0, turnover=turn, rho=rho, delta=delta,
                seeds=len(v), N=f("N"),
                N_sd=float(np.std([x["N"] for x in v], ddof=1)),
                spent=f("spent"), supply_held=f("supply_held"),
                space_held=f("space_held"), active=f("active"))


def holds(c):
    """Which bound ends the growth, as the largest share of the pool."""
    k = max(("spent", "supply_held", "space_held", "active"),
            key=lambda n: c[n])
    return {"spent": "spent", "supply_held": "supply", "space_held": "room",
            "active": "still active"}[k]


def elasticity(w, P0, turn):
    """d log N / d log P0 at one corner, paired seed by seed."""
    up = [run(w, P0 * np.exp(H), turn, s) for s in SEEDS]
    dn = [run(w, P0 * np.exp(-H), turn, s) for s in SEEDS]
    e = [(np.log(max(a["N"], 1)) - np.log(max(b["N"], 1))) / (2 * H)
         for a, b in zip(up, dn)]
    return float(np.mean(e)), float(np.std(e, ddof=1) / np.sqrt(len(e)))


if __name__ == "__main__":
    jobs([(w, P0, t, s) for w in RENEWAL for P0 in POOL
          for t, _, _ in TURNOVER for s in SEEDS])

    cells = {}
    print(f"\nWhere the pool stands at arrest, {len(SEEDS)} seeds per cell. The "
          f"last column names the\nlargest of the four shares, which is the "
          f"bound that ended the growth.\n")
    print(f"{'renewal':>8} {'P0':>6} {'turnover':>9} {'rho':>5} {'delta':>6} "
          f"{'N':>7} {'spent':>7} {'supply':>7} {'room':>7} {'active':>7}  "
          f"{'holds':>12}  named")
    for w in RENEWAL:
        for P0 in POOL:
            for t, rho, d in TURNOVER:
                c = cell(w, P0, t)
                cells[(w, P0, t)] = c
                print(f"{w:>8.0f} {P0:>6.0f} {t:>9} {rho:>5.2f} {d:>6.2f} "
                      f"{c['N']:>7.1f} {c['spent']:>7.2f} "
                      f"{c['supply_held']:>7.2f} {c['space_held']:>7.2f} "
                      f"{c['active']:>7.2f}  {holds(c):>12}  "
                      f"{NAMED.get((w, P0, t), '')}", flush=True)
        print()

    # Does renewal change the closure at fixed pool and turnover?
    ncell = len(POOL) * len(TURNOVER)
    print(f"Changing renewal alone, at each of the {ncell} settings of pool and "
          f"turnover:\n")
    print(f"{'P0':>6} {'turnover':>9} {'spent, w=0':>11} {'spent, w=1':>11} "
          f"{'N, w=0':>8} {'N, w=1':>8}  {'holds at w=0':>13} -> "
          f"{'holds at w=1':<13}")
    flips = 0
    for P0 in POOL:
        for t, _, _ in TURNOVER:
            a, b = cells[(0.0, P0, t)], cells[(1.0, P0, t)]
            same = holds(a) == holds(b)
            flips += not same
            print(f"{P0:>6.0f} {t:>9} {a['spent']:>11.2f} {b['spent']:>11.2f} "
                  f"{a['N']:>8.1f} {b['N']:>8.1f}  {holds(a):>13} -> "
                  f"{holds(b):<13}{'' if same else '  (changes)'}", flush=True)
    print(f"\n  renewal changes which bound holds growth in {flips} of "
          f"{ncell} settings of the pool and the turnover.")

    # And does the pool or the turnover change it at fixed renewal?
    for w in RENEWAL:
        seen = {holds(cells[(w, P0, t)]) for P0 in POOL for t, _, _ in TURNOVER}
        print(f"  at renewal {w:.0f}, the {ncell} settings of pool and turnover "
              f"give {len(seen)} distinct bound(s): {', '.join(sorted(seen))}")

    # the elasticity with the pool and the turnover matched
    jobs([(w, P0 * f, t, s) for w in RENEWAL for P0, t in MATCHED
          for f in (np.exp(H), np.exp(-H)) for s in SEEDS])
    print("\nElasticity of the final unit count to the seeded pool, with the "
          "pool and the\nturnover matched between the two sources, paired seed "
          "by seed:\n")
    print(f"{'P0':>6} {'turnover':>9} {'spent, w=0':>18} "
          f"{'renewing, w=1':>18}  gap")
    els = {}
    for P0, t in MATCHED:
        a = elasticity(0.0, P0, t)
        b = elasticity(1.0, P0, t)
        els[f"P{P0:.0f}_{t}"] = dict(P0=P0, turnover=t, rho=TURN[t][0],
                                     delta=TURN[t][1], spent=a, renewing=b)
        gap = a[0] - b[0]
        se = (a[1] ** 2 + b[1] ** 2) ** 0.5
        print(f"{P0:>6.0f} {t:>9} {a[0]:>+12.2f} ± {a[1]:.2f} "
              f"{b[0]:>+12.2f} ± {b[1]:.2f}  "
              f"{gap:+.2f} ({'z = %.1f' % (gap / se) if se > 0 else 'n/a'})",
              flush=True)
    print("\n  The published contrast, +0.71 for a spent source against +0.07 "
          "for a renewing\n  one, is measured at P0 60 and 2600 respectively. "
          "The rows above hold the pool\n  and the turnover fixed, so what is "
          "left is the binary property alone.")

    res = {f"w{w}_P{P0:.0f}_{t}": cells[(w, P0, t)]
           for w in RENEWAL for P0 in POOL for t, _, _ in TURNOVER}
    res["_meta"] = dict(seeds=len(SEEDS), days=DAYS, renewal=list(RENEWAL),
                        pool=list(POOL),
                        turnover=[dict(name=n, rho=r, delta=d)
                                  for n, r, d in TURNOVER],
                        renewal_flips=flips, cells=ncell, H=H)
    res["_elasticity_matched"] = els
    json.dump(res, open(ROOT / "results" / "regimes.json", "w"), indent=1)
