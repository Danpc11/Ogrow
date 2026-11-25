#!/usr/bin/env python
"""What changes when a unit that loses its supply can disappear.

`scripts/32_inflow_law.py` reports that a construct held at a fixed inflow
settles in total while its units do not: the residual sum |m - ceiling| over
sum ceiling stays around a third, and most of the units end at the mass floor.
The reason is a mechanism the reference model does not have. A unit whose share
of the delivery collapses has no way to disappear, so it sits at the floor
instead of involuting, and the construct ends as a few large units among many
vestigial ones.

This is the sensitivity that asks what the result looks like when that
mechanism is supplied. `involution` is a rate per day at which a unit whose
supply has fallen below `involution_q` of the setpoint loses mass; below half
the mass a unit is born with it is removed. The reference model is the rate of
zero, so the first row of each table is the published result.

The rate is not calibrated against any measurement. What makes the sensitivity
worth reporting is that the answer barely moves across a tenfold range of it,
so the conclusion rests on the mechanism existing and not on its speed.

Everything is at unit flow, where the published result sits. The demand phi is
measured once on an intact ensemble of a disjoint block of seeds, exactly as in
32, and does not depend on involution because it is read before the injury.

Usage:  python scripts/33_involution.py [REGIME ...]     default: both

    OGROW_SEEDS         replicates, default 24
    OGROW_JOBS          workers
    OGROW_ADOPT_CACHE   see organgrow.runcache
"""
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow, _p
from organgrow.runcache import fingerprint, cached

OUT = ROOT / "results" / "involution"
OUT.mkdir(parents=True, exist_ok=True)
S, NU, CLOSE, END = 0.6, 8.0, 140, 220
RATES = (0.0, 0.02, 0.05, 0.10, 0.20)
# The rate says how fast the mechanism runs. These say what turns it on and
# when a unit is gone, which is a different question: a result that moved with
# them would be a result about the algorithm. REF is the setting the rate sweep
# uses, so each sweep below crosses the rate sweep at one point.
REF = dict(rate=0.05, q=0.1, hold=1, exit=0.5)
TRIGGERS = (("involution_q", "q", (0.03, 0.05, 0.10, 0.20, 0.30)),
            ("involution_hold", "hold", (1, 3, 5, 10)),
            ("involution_exit", "exit", (0.25, 0.50, 0.75, 0.90)))
RECOVER = 0.90          # share of pre-injury mass whose first crossing we time
NSEED = int(os.environ.get("OGROW_SEEDS", 24))
SEEDS = tuple(range(300, 300 + NSEED))          # where the sweep is evaluated
CAL = tuple(range(400, 400 + NSEED))            # where phi is calibrated
REG = {"liver": dict(P0=2600.0, rho=0.10, delta=0.10, q_set=0.6, renewal=1.0),
       "kidney": dict(P0=300.0, rho=0.00, delta=0.05, q_set=0.6, renewal=0.0)}
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=6000)
LOSS, SCALE = 0.3, 1.0
M0 = float(_p("m0", None))


def _kw(reg, rate, phi, q=None, hold=None, exit_=None):
    return dict(injure_day=CLOSE, injury="parenchyma", injury_fraction=LOSS,
                injury_reactivates=1.0, inflow_limited=True,
                inflow_scale=SCALE, phi_demand=float(phi), involution=float(rate),
                involution_q=float(REF["q"] if q is None else q),
                involution_hold=int(REF["hold"] if hold is None else hold),
                involution_exit=float(REF["exit"] if exit_ is None else exit_))


def intact(reg, seed):
    """The intact organ at arrest, which is where phi is measured."""
    f = OUT / f"intact_{reg}_s{seed}.json"
    fp = fingerprint(regime=reg, seed=seed, S=S, NU=NU, days=CLOSE,
                     params={**FIX, **REG[reg]})
    hit = cached(f, fp)
    if hit is not None:
        return hit
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S), 190 * S,
                        days=CLOSE, rng=np.random.default_rng(seed), nu=NU,
                        **FIX, **REG[reg])
    r = dict(regime=reg, seed=seed, mass=float(np.asarray(tr.mass)[-1]),
             Q=float(np.asarray(tr.Q)[-1]), N=len(U), fingerprint=fp)
    r["phi"] = r["Q"] / max(r["mass"], 1e-12)
    json.dump(r, open(f, "w"))
    return r


def run(reg, seed, rate, phi, trig=None):
    """One construct, injured and then followed with involution at `rate`.

    `trig` overrides the conditions that turn the mechanism on: the supply at
    which a unit counts as starving, the consecutive days it must stay there,
    and the mass at which it is removed. Empty means the reference setting.
    """
    trig = dict(trig or {})
    tag = f"{reg}_s{seed}_i{rate}"
    if trig:
        tag += "_" + "_".join(f"{k}{v}" for k, v in sorted(trig.items()))
    f = OUT / f"{tag}.json"
    kw = _kw(reg, rate, phi, **{k + ("_" if k == "exit" else ""): v
                                for k, v in trig.items()})
    fp = fingerprint(regime=reg, seed=seed, S=S, NU=NU, days=END,
                     params={**FIX, **REG[reg], **kw})
    hit = cached(f, fp)
    if hit is not None:
        return hit
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S), 190 * S,
                        days=END, rng=np.random.default_rng(seed), nu=NU,
                        **FIX, **REG[reg], **kw)
    pre = tr.injury["pre"]
    M = np.asarray(tr.mass, float)
    c = tr.final.get("ceiling")
    c = None if c is None else np.asarray(c, float)
    # When the mass first comes back to RECOVER of what it was, counted from
    # the injury. The endpoint can be reached along different trajectories, so
    # a result that is flat in the rate at day 220 need not be flat in time.
    want = RECOVER * float(pre["mass"])
    after = M[CLOSE:]
    hit_day = np.flatnonzero(after >= want)
    # Mass comes back in days, because the surviving units are far below much
    # larger ceilings and the growth clock has just been reset. The trajectory
    # that belongs to this mechanism is the removals: the day by which nine
    # tenths of them have happened, counted from the injury.
    iv = np.asarray(tr.involuted, float)[CLOSE:]
    done = np.flatnonzero(iv >= 0.9 * iv[-1]) if len(iv) and iv[-1] > 0 else []
    r = dict(regime=reg, seed=seed, rate=rate, phi=phi, trig=dict(trig),
             pre_mass=float(pre["mass"]), pre_N=int(pre["N"]),
             pre_Q=float(tr.Q[CLOSE - 1]), mass=float(M[-1]), N=len(U),
             residual=(float("nan") if c is None
                       else float(np.abs(m - c).sum() / max(c.sum(), 1e-12))),
             floor_share=float((m <= 1.05 * M0).mean()) if len(m) else float("nan"),
             drift_per_day=float(abs(M[-1] - M[-11]) / max(M[-1], 1e-12) / 10.0),
             days_to_recover=(int(hit_day[0]) if len(hit_day) else -1),
             days_to_clear=(int(done[0]) if len(done) else -1),
             involuted=int(tr.fate.get("units_involuted", 0)),
             idle_tips=int(tr.idle_tips[-1]) if tr.idle_tips else 0,
             vessels=int(tr.vessels[-1]) if tr.vessels else 0,
             Q_end=float(np.asarray(tr.Q)[-1]), C_end=float(np.asarray(tr.C)[-1]),
             balance=float(tr.fate["balance"]), fingerprint=fp)
    json.dump(r, open(f, "w"))
    return r


def _warm(a):
    (intact(*a[1:]) if a[0] == "intact" else run(*a[1:]))


def jobs(grid):
    if not grid:
        return
    k = max(1, min(int(os.environ.get("OGROW_JOBS", 0)) or (os.cpu_count() or 1),
                   len(grid)))
    print(f"running {len(grid)} simulations on {k} worker(s)", flush=True)
    if k > 1:
        with ProcessPoolExecutor(max_workers=k) as pool:
            for _ in pool.map(_warm, grid):
                pass
    else:
        for g in grid:
            _warm(g)


def summarise(reg, phi, rate, trig=None):
    """Mean behaviour over the seeds at one setting."""
    v = [run(reg, s, rate, phi, trig) for s in SEEDS]
    rel = np.array([r["mass"] / max(r["pre_mass"], 1e-12) for r in v])
    un = np.array([r["N"] / max(r["pre_N"], 1) for r in v])
    day = np.array([r["days_to_recover"] for r in v], float)
    got = day >= 0
    clr = np.array([r.get("days_to_clear", -1) for r in v], float)
    gotc = clr >= 0
    return dict(rate=rate, trig=dict(trig or {}),
                mass_rel=float(rel.mean()), mass_rel_sd=float(rel.std(ddof=1)),
                mass_rel_lo=float(rel.min()), mass_rel_hi=float(rel.max()),
                units_rel=float(un.mean()),
                units_rel_lo=float(un.min()), units_rel_hi=float(un.max()),
                residual=float(np.nanmean([r["residual"] for r in v])),
                floor_share=float(np.nanmean([r["floor_share"] for r in v])),
                drift=float(max(r["drift_per_day"] for r in v)),
                involuted=float(np.mean([r["involuted"] for r in v])),
                idle_frac=float(np.mean([r["idle_tips"] / max(r["vessels"], 1)
                                         for r in v])),
                Q_end=float(np.mean([r["Q_end"] for r in v])),
                recovered=int(got.sum()), seeds=len(v),
                days_to_recover=(float(day[got].mean()) if got.any() else float("nan")),
                days_to_clear=(float(clr[gotc].mean()) if gotc.any() else float("nan")),
                days_to_clear_lo=(float(clr[gotc].min()) if gotc.any() else float("nan")),
                days_to_clear_hi=(float(clr[gotc].max()) if gotc.any() else float("nan")),
                balance=float(max(abs(r["balance"]) for r in v)))


def line(label, d):
    rec = ("-" if d["recovered"] == 0
           else f"{d['days_to_recover']:.0f} d")
    clr = ("-" if d["days_to_clear"] != d["days_to_clear"]
           else f"{d['days_to_clear']:.0f} [{d['days_to_clear_lo']:.0f}, "
                f"{d['days_to_clear_hi']:.0f}]")
    return (f"{label:>14} {d['mass_rel']:>8.3f} [{d['mass_rel_lo']:.2f}, "
            f"{d['mass_rel_hi']:.2f}] {d['units_rel']:>10.2f} "
            f"{d['residual']:>9.3f} {d['floor_share']:>9.2f} "
            f"{d['idle_frac']:>10.2f} {rec:>8} {clr:>16}")


HEAD = (f"{'setting':>14} {'mass / pre-injury':>21} {'units/pre':>10} "
        f"{'residual':>9} {'at floor':>9} {'idle tips':>10} "
        f"{'90 % mass':>8} {'90 % removed':>16}")


if __name__ == "__main__":
    regimes = [a for a in sys.argv[1:] if a in REG] or list(REG)
    res = {}
    for reg in regimes:
        jobs([("intact", reg, s) for s in CAL
              if cached(OUT / f"intact_{reg}_s{s}.json",
                        fingerprint(regime=reg, seed=s, S=S, NU=NU, days=CLOSE,
                                    params={**FIX, **REG[reg]})) is None])
        ref = [intact(reg, s) for s in CAL]
        phi = float(np.mean([r["phi"] for r in ref]))
        print(f"\n{reg}: a 30 % loss at unit flow, demand phi = {phi:.5f} "
              f"calibrated on seeds {CAL[0]} to {CAL[-1]}, evaluated on seeds "
              f"{SEEDS[0]} to {SEEDS[-1]}")

        jobs([("run", reg, s, x, phi) for x in RATES for s in SEEDS
              if cached(OUT / f"{reg}_s{s}_i{x}.json",
                        fingerprint(regime=reg, seed=s, S=S, NU=NU, days=END,
                                    params={**FIX, **REG[reg],
                                            **_kw(reg, x, phi)})) is None])

        print("\n  the rate: how fast the mechanism runs\n")
        print(HEAD)
        rows = [summarise(reg, phi, x) for x in RATES]
        for d in rows:
            print(line(f"rate {d['rate']:.2f}", d), flush=True)

        # The conditions that turn it on, each crossing the rate sweep at the
        # reference setting, so the sweeps agree where they meet.
        trig_rows = {}
        for name, key, vals in TRIGGERS:
            want = [{key: v} for v in vals if v != REF[key]]
            jobs([("run", reg, s_, REF["rate"], phi, t) for t in want for s_ in SEEDS
                  if cached(OUT / (f"{reg}_s{s_}_i{REF['rate']}_"
                                   + "_".join(f"{k}{x}" for k, x in sorted(t.items()))
                                   + ".json"),
                            fingerprint(regime=reg, seed=s_, S=S, NU=NU, days=END,
                                        params={**FIX, **REG[reg],
                                                **_kw(reg, REF["rate"], phi,
                                                      **{k + ("_" if k == "exit" else ""): x
                                                         for k, x in t.items()})})) is None])
            print(f"\n  {name}: what turns it on\n")
            print(HEAD)
            out = []
            for v in vals:
                d = summarise(reg, phi, REF["rate"],
                              None if v == REF[key] else {key: v})
                out.append(d)
                print(line(f"{key} {v}", d), flush=True)
            trig_rows[name] = out

        on = [r for r in rows if r["rate"] > 0]
        res[reg] = dict(phi=phi, seeds=len(SEEDS), rows=rows,
                        reference=dict(REF), triggers=trig_rows,
                        calibration_seeds=[CAL[0], CAL[-1]],
                        evaluation_seeds=[SEEDS[0], SEEDS[-1]],
                        units_rel_on=[min(r["units_rel"] for r in on),
                                      max(r["units_rel"] for r in on)],
                        floor_on=[min(r["floor_share"] for r in on),
                                  max(r["floor_share"] for r in on)],
                        residual_on=[min(r["residual"] for r in on),
                                     max(r["residual"] for r in on)])
        every = on + [d for v in trig_rows.values() for d in v if d["rate"] > 0]
        u = [d["units_rel"] for d in every]
        mm = [d["mass_rel"] for d in every]
        res[reg]["units_rel_any"] = [min(u), max(u)]
        res[reg]["mass_rel_any"] = [min(mm), max(mm)]
        d0 = rows[0]
        print(f"\n  Across every setting of the rate and of the three conditions "
              f"that turn it on,\n  the construct keeps {min(u):.2f} to {max(u):.2f} "
              f"of its units against {d0['units_rel']:.2f} with the\n  mechanism "
              f"absent, and its mass stays between {min(mm):.2f} and {max(mm):.2f} "
              f"of pre-injury\n  against {d0['mass_rel']:.2f}. What the mechanism "
              f"changes is the architecture, not the\n  mass, and that holds "
              f"however it is switched on.", flush=True)

    json.dump(res, open(ROOT / "results" / "involution.json", "w"), indent=1)
