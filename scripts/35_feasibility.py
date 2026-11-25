#!/usr/bin/env python
"""Can a construct be built to a joint specification, and up to what size?

The design rules in this paper say what a supply constraint fixes and what it
leaves free. This asks the question a biofabrication group would actually put
to them: given a target volume and a specification the construct has to meet,
does any setting of the knobs an operator controls meet all of it at once, and
if not, which requirement is the one that cannot be met?

The specification has three requirements, each evaluated over an ensemble
because two of them are ensemble quantities:

    density        the units must occupy at least `fill_min` of the domain's
                   close packing, (R/lam)^3 * 0.74 units
    vascular cost  the tree must use at most `vessel_max` segments
    reproducibility  the unit count must have a coefficient of variation of at
                   most `cv_max` across independent realisations

The knobs are the four things a bioreactor run actually sets: how many
competent cells are seeded (`P0`), the nucleation propensity of the medium
(`nu`), how long the construct is cultured (`days`), and which cell source is
used -- one whose cells are restored when a unit forms, or one whose are
spent. Target volume is not a knob; it is the specification's first line, and
sweeping it is what produces a frontier.

Selection is maximin. Each requirement gives a slack normalised by its own
threshold, a configuration is feasible when all three are non-negative, and
the chosen configuration is the one whose smallest slack is largest. That is
the configuration with the most room to spare on its tightest requirement,
which is what "robustly meets the spec" means when the requirements are not
commensurate.

Search and verification use disjoint seeds: the configuration is chosen on
seeds 500-515 and then re-run on seeds 600-623, which no part of the search
saw. One thing this does NOT do, and the printed output says so again: new
seeds of the same model are new realisations, not new evidence about the
model. This measures whether a configuration chosen on one batch holds up on
another, which is reproducibility. It is not a validation of the model against
a construct, and no claim here should be read as one.

Every run also reports whether it had arrested by the criterion of section
3.3, fewer than 5% more units over 40 further days, read off the daily trace.
A configuration that had not arrested is reported but excluded from selection,
because its numbers are a snapshot of something still growing.

Usage:  python scripts/35_feasibility.py

    OGROW_SEEDS   search seeds per configuration, default 16
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

OUT = ROOT / "results" / "feasibility"
OUT.mkdir(parents=True, exist_ok=True)

LAM, R0, PACK = 7.4, 190.0, 0.74
ARREST_WINDOW, ARREST_TOL = 40, 0.05
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=30000,
           q_set=0.6)

#: the target volumes, as the domain scale; volume goes as S^3
SCALES = (0.4, 0.6, 0.8, 1.0)
#: the knobs an operator sets
POOL = (300.0, 1000.0, 3000.0, 9000.0)
NU = (6.0, 12.0, 24.0)
DAYS = (140, 280)
#: and the cell source, the binary property of section 3.7 with the division
#: and loss rates that go with each published source
SOURCE = {"renewing": dict(renewal=1.0, rho=0.10, delta=0.10),
          "spent": dict(renewal=0.0, rho=0.00, delta=0.05)}

#: the worked specification. Stated here rather than chosen after the fact;
#: the sweep at the end reports how the frontier moves as it changes.
SPEC = dict(fill_min=0.10, vessel_max=2500.0, cv_max=0.15)

NSEARCH = int(os.environ.get("OGROW_SEEDS", 16))
if NSEARCH < 2:
    #: one realisation has no coefficient of variation, and the sample
    #: deviation would come back as nan, which compares false against every
    #: threshold and would report the whole grid infeasible without saying why
    raise SystemExit("OGROW_SEEDS must be at least 2: one of the three "
                     "requirements is a variation across realisations. "
                     "Fewer than 8 makes it too noisy to select on.")
SEARCH = tuple(range(500, 500 + NSEARCH))
HOLDOUT = tuple(range(600, 624))

#: the spec sweep: the largest volume that each pair of requirements admits
FILL_GRID = (0.04, 0.06, 0.08, 0.10, 0.14)
CV_GRID = (0.10, 0.15, 0.20, 0.30, 0.50)


def key(S, P0, nu, days, src, seed):
    return f"S{S}_P{P0:.0f}_nu{nu:g}_d{days}_{src}_s{seed}"


def params(P0, nu, src):
    return dict(nu=float(nu), P0=float(P0), **FIX, **SOURCE[src])


def fp_of(S, P0, nu, days, src, seed):
    return fingerprint(seed=seed, S=S, days=days,
                       params=params(P0, nu, src))


def late_growth(N_day):
    """Fractional gain in unit count over the last `ARREST_WINDOW` days.

    A run shorter than the window has nothing to compare against, so it counts
    as still growing rather than as arrested by default.
    """
    if len(N_day) <= ARREST_WINDOW:
        return float("inf")
    back = N_day[-1 - ARREST_WINDOW]
    return (N_day[-1] - back) / max(back, 1.0)


def run(S, P0, nu, days, src, seed):
    """One construct, grown to `days` at scale `S`, with its readouts."""
    f = OUT / f"{key(S, P0, nu, days, src, seed)}.json"
    fp = fp_of(S, P0, nu, days, src, seed)
    hit = cached(f, fp)
    if hit is not None:
        return hit
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S),
                        R0 * S, days=days, rng=np.random.default_rng(seed),
                        **params(P0, nu, src))
    ceiling = (R0 * S / LAM) ** 3 * PACK
    #: the arrest criterion of section 3.3, off the daily trace
    growth = late_growth(list(tr.N))
    r = dict(S=S, P0=P0, nu=nu, days=days, src=src, seed=seed, N=len(U),
             vessels=len(vt), ceiling=ceiling, fill=len(U) / ceiling,
             mass=float(m.sum()),
             #: clamped so the record stays valid JSON for any reader; a run
             #: shorter than the arrest window is the only way to reach it
             late_growth=float(min(growth, 1e6)),
             arrested=bool(growth < ARREST_TOL), fingerprint=fp)
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


def cell(S, P0, nu, days, src, seeds):
    """The three specified quantities for one configuration over one block."""
    v = [run(S, P0, nu, days, src, s) for s in seeds]
    N = np.array([x["N"] for x in v], float)
    return dict(S=S, P0=P0, nu=nu, days=days, src=src, seeds=len(v),
                N=float(N.mean()),
                fill=float(np.mean([x["fill"] for x in v])),
                vessels=float(np.mean([x["vessels"] for x in v])),
                cv=float(N.std(ddof=1) / max(N.mean(), 1e-9)),
                arrested=bool(all(x["arrested"] for x in v)),
                late_growth=float(np.mean([x["late_growth"] for x in v])))


def slacks(c, spec):
    """Each requirement's slack, normalised by its own threshold."""
    return dict(density=(c["fill"] - spec["fill_min"]) / spec["fill_min"],
                vascular=(spec["vessel_max"] - c["vessels"]) / spec["vessel_max"],
                reproducibility=(spec["cv_max"] - c["cv"]) / spec["cv_max"])


def margin(c, spec):
    """The maximin score: the smallest of the three slacks."""
    return min(slacks(c, spec).values())


def binding(cells, spec):
    """Which requirement no configuration could meet, and the best each got."""
    best = {k: max(slacks(c, spec)[k] for c in cells)
            for k in ("density", "vascular", "reproducibility")}
    return min(best, key=best.get), best


CONFIGS = [(P0, nu, days, src) for P0 in POOL for nu in NU for days in DAYS
           for src in SOURCE]


if __name__ == "__main__":
    jobs([(S, P0, nu, days, src, s) for S in SCALES
          for P0, nu, days, src in CONFIGS for s in SEARCH])

    print(f"\nSearch block: {len(SEARCH)} seeds, {len(CONFIGS)} configurations "
          f"at each of {len(SCALES)} volumes.")
    print(f"Specification: fill >= {SPEC['fill_min']:.2f} of close packing, "
          f"<= {SPEC['vessel_max']:.0f} vessel segments, "
          f"CV <= {SPEC['cv_max']:.2f}.\n")

    search = {S: [cell(S, *c, SEARCH) for c in CONFIGS] for S in SCALES}

    # What any configuration can reach at each volume, before any spec is
    # imposed. The three are reached by different configurations, which is
    # the whole difficulty: they are not jointly attainable.
    print("Best any configuration reaches at each volume, one column at a "
          "time. The three\nbests are different configurations, so no "
          "construct has all three at once.\n")
    print(f"{'S':>5} {'volume':>8} {'best fill':>10} {'fewest vessels':>15} "
          f"{'lowest CV':>10}  {'arrested':>9}")
    reach = {}
    for S in SCALES:
        ok = [c for c in search[S] if c["arrested"]]
        pool_ = ok or search[S]
        bf = max(pool_, key=lambda c: c["fill"])
        bv = min(pool_, key=lambda c: c["vessels"])
        bc = min(pool_, key=lambda c: c["cv"])
        reach[S] = dict(best_fill=bf, fewest_vessels=bv, lowest_cv=bc,
                        arrested=len(ok), configs=len(search[S]))
        print(f"{S:>5.1f} {S ** 3:>8.3f} {100 * bf['fill']:>9.2f}% "
              f"{bv['vessels']:>15.0f} {bc['cv']:>10.3f}  "
              f"{len(ok):>4d}/{len(search[S]):<4d}")

    # The worked specification: choose on the search block, then verify the
    # chosen configuration on seeds the search never saw.
    print("\nThe specification, configuration chosen by maximin on the "
          "search block:\n")
    print(f"{'S':>5} {'feasible':>9} {'P0':>6} {'nu':>5} {'days':>5} "
          f"{'source':>9} {'fill':>7} {'vessels':>8} {'CV':>6} {'margin':>7}")
    chosen, infeasible, held, optimism = {}, {}, {}, []
    for S in SCALES:
        ok = [c for c in search[S] if c["arrested"]]
        feas = [c for c in ok if margin(c, SPEC) >= 0]
        if feas:
            w = max(feas, key=lambda c: margin(c, SPEC))
            chosen[S] = w
            print(f"{S:>5.1f} {'yes':>9} {w['P0']:>6.0f} {w['nu']:>5.0f} "
                  f"{w['days']:>5d} {w['src']:>9} {100 * w['fill']:>6.2f}% "
                  f"{w['vessels']:>8.0f} {w['cv']:>6.3f} "
                  f"{margin(w, SPEC):>+7.2f}")
        else:
            who, best = binding(ok or search[S], SPEC)
            infeasible[S] = dict(binding=who, best=best)
            print(f"{S:>5.1f} {'no':>9} {'-':>6} {'-':>5} {'-':>5} {'-':>9} "
                  f"{'-':>7} {'-':>8} {'-':>6} {'-':>7}   "
                  f"{who} is the requirement none met "
                  f"(best slack {best[who]:+.2f})")

    if chosen:
        jobs([(S, w["P0"], w["nu"], w["days"], w["src"], s)
              for S, w in chosen.items() for s in HOLDOUT])
        print(f"\nThe chosen configurations re-run on {len(HOLDOUT)} held-out "
              f"seeds ({HOLDOUT[0]}-{HOLDOUT[-1]}),\nwhich no part of the "
              f"search saw:\n")
        print(f"{'S':>5} {'fill, search':>13} {'fill, held out':>15} "
              f"{'CV, search':>11} {'CV, held out':>13} {'margin':>15} "
              f"{'still meets spec':>17}")
        for S, w in chosen.items():
            h = cell(S, w["P0"], w["nu"], w["days"], w["src"], HOLDOUT)
            held[S] = h
            ms, mh = margin(w, SPEC), margin(h, SPEC)
            optimism.append(ms - mh)
            print(f"{S:>5.1f} {100 * w['fill']:>12.2f}% "
                  f"{100 * h['fill']:>14.2f}% {w['cv']:>11.3f} "
                  f"{h['cv']:>13.3f} {ms:>+7.2f} -> {mh:>+4.2f} "
                  f"{'yes' if mh >= 0 else 'no':>17}")
        _o = float(np.mean(optimism))
        print(f"\n  The margin measured on the search block "
              f"{'overstates' if _o >= 0 else 'understates'} the held-out "
              f"margin by\n  {abs(_o):.2f} on average: the configuration was "
              f"chosen partly for the realisations that\n  happened to favour "
              f"it, and the coefficient of variation is the noisiest of the\n  "
              f"three quantities to estimate. The held-out column is the "
              f"honest one, and the\n  reason for holding seeds out at all.")
        print("\n  Held-out seeds are new realisations of the same model, not "
              "new evidence\n  about it. This measures whether a configuration "
              "chosen on one batch holds\n  on another. It is not a validation "
              "against a construct.")

    # How the frontier moves with the specification, so that the result is a
    # map and not one hand-picked triple of thresholds.
    print(f"\nLargest volume that stays feasible, as the density and "
          f"reproducibility\nrequirements change, at a vascular budget of "
          f"{SPEC['vessel_max']:.0f} segments. Entries are S;\n'none' means "
          f"not even the smallest volume admits a configuration.\n")
    print(f"{'fill >=':>8} " + " ".join(f"{'CV<=' + format(c, '.2f'):>9}"
                                        for c in CV_GRID))
    frontier = {}
    for fmin in FILL_GRID:
        row = []
        for cmax in CV_GRID:
            sp = dict(fill_min=fmin, vessel_max=SPEC["vessel_max"], cv_max=cmax)
            good = [S for S in SCALES
                    if any(margin(c, sp) >= 0
                           for c in search[S] if c["arrested"])]
            frontier[f"f{fmin}_cv{cmax}"] = max(good) if good else None
            row.append(f"{max(good):>9.1f}" if good else f"{'none':>9}")
        print(f"{fmin:>8.2f} " + " ".join(row))

    # the margins are quoted in the manuscript, so they are stored rather
    # than left to be recomputed by whatever reads this
    for S, c in chosen.items():
        c["margin"] = margin(c, SPEC)
    for S, h in held.items():
        h["margin"] = margin(h, SPEC)

    res = dict(
        spec=SPEC, scales=list(SCALES), search_seeds=list(SEARCH),
        holdout_seeds=list(HOLDOUT), configs=len(CONFIGS),
        knobs=dict(P0=list(POOL), nu=list(NU), days=list(DAYS),
                   source=sorted(SOURCE)),
        reach={str(S): {k: v for k, v in r.items()} for S, r in reach.items()},
        chosen={str(S): c for S, c in chosen.items()},
        held_out={str(S): c for S, c in held.items()},
        infeasible={str(S): v for S, v in infeasible.items()},
        frontier=frontier,
        optimism=float(np.mean(optimism)) if held else None,
        arrest=dict(window=ARREST_WINDOW, tol=ARREST_TOL))
    json.dump(res, open(ROOT / "results" / "feasibility.json", "w"), indent=1)
