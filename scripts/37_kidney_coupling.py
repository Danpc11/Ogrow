#!/usr/bin/env python
"""An extension: progenitors that make their tip branch, and what ablation does then.

The reference model's kidney-like source does not divide and its tree does
not depend on the pool. Section 3.10 sets it beside two mouse experiments that
each removed about 40 % of nephron progenitors: it agrees with the smaller of
the two deficits, misses the larger by half, and has the mechanism of neither.
This script asks what the model needs for both.

What the experiments say the kidney does, and the reference model does not:

    the progenitors divide       nephron progenitors renew themselves until a
                                 timed stop; here rho > 0 and delta = 0, and
                                 the run ends at a fixed day
    they make the tip branch     fewer cells around a tip, slower branching,
                                 and each branching splits the cap in two
                                 (Cebrian et al. 2014); here a segment branches
                                 only while it holds two caps, `cap_size`
    they may divide faster       when depleted, in one experiment and not in
                                 the other (Muthukrishnan et al. 2018 against
                                 Cebrian et al. 2014); here `cap_response`

The ablation is done as the experiments did it: on one day, each progenitor
dies with probability 0.4 wherever it sits (`ablate_day`, `ablate_fraction`),
and the same seed grown without it is the littermate.

This is an extension and not a recalibration. Nothing in the reference model
uses `cap_size`, the reference run's digest does not move, and no number
reported elsewhere depends on this script. The rules were added after the
comparison of Section 3.10 was seen, so what follows cannot be a prediction of
the endpoint. What it can test is whether one mechanism gives several things
at once that were measured separately: a pool that stays reduced while each
cell divides as before, a cap that recovers on its own, branching that slows,
and a deficit in units; and whether the single difference between the two
experiments, the proliferative response, is enough to separate their outcomes.

The four numbers the extension adds are not measured. The cap size, the
division rate, the window and the day of ablation are fixed below and the two
that bear on whether perfusion limits nucleation, q_set and nu, are swept,
because that turned out to decide the answer.

Usage:  python scripts/37_kidney_coupling.py
        OGROW_JOBS=n   worker count
"""
import json
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow
from organgrow.autotext import write
from organgrow.runcache import fingerprint, cached, save_json, stale

OUT = ROOT / "results" / "kidney_coupling"
DOCS = (ROOT / "manuscript" / "MANUSCRIPT.md", ROOT / "manuscript" / "MANUSCRIPT_A.md")
#: the two settings the manuscript quotes: where perfusion limits nucleation
#: as in the reference model, and where it limits it least
QUOTED = ("reference", "threshold quartered, law slowed")

S, SEEDS = 0.6, tuple(range(300, 324))
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=6000)
#: the coupled source: one cap seeded at the root, dividing, never lost to
#: attrition, spent when it makes a unit
CAP, RHO, WINDOW = 60.0, 0.05, 120
SOURCE = dict(P0=CAP, rho=RHO, delta=0.0, renewal=0.0)
#: the ablation, and when: late enough that the tree has a few sites, early
#: enough that most of the window is still ahead
ABLATE_DAY, ABLATE_FRACTION = 45, 0.4
#: how strongly perfusion limits nucleation. The first is the reference
#: model's own setting; the others lower the threshold and slow the law.
NUCLEATION = {
    "reference": dict(q_set=0.6, nu=8.0),
    "threshold halved": dict(q_set=0.3, nu=8.0),
    "threshold halved, law slowed": dict(q_set=0.3, nu=2.0),
    "threshold quartered, law slowed": dict(q_set=0.15, nu=2.0),
}
#: 0: a depleted cap divides as before. 1: it divides twice as fast, which is
#: what the study that found compensation measured.
RESPONSE = (0.0, 1.0)

#: What was measured, as ablated over control. Typed from the papers.
#: Cebrian et al. 2014, tamoxifen at E12.5: SIX2+ cells per tip 375 against 620
#: at E13.5, 357 against 376 at E14.5, 263 against 253 at E15.5 (Fig. 6C); tips
#: per kidney 55 against 52, 133 against 188, 215 against 364 (Fig. 6F); total
#: cap mesenchyme grew 2.3-fold against 2.2-fold in a day, with no change in
#: mitosis (text); glomeruli 48.3 % at E19.5 and 57.6 % at P50 (Fig. 4B).
#: Muthukrishnan et al. 2018, tamoxifen at E15.5: CITED1+ cells 60 % of control
#: after 24 h and 85 % after 48 h (Fig. 1G); proliferating CITED1+ cells more
#: than doubled (Fig. 1F); glomeruli down 12 % at P3 and 21 % at P10 (text).
MEASURED = {
    "no response (Cebrian 2014)": dict(
        cap=[0.60, 0.95, 1.04], tips=[1.06, 0.71, 0.59], units=[0.483, 0.576]),
    "response (Muthukrishnan 2018)": dict(
        pool=[0.60, 0.85], units=[0.88, 0.79]),
}


def spec(name, response, ablated, seed):
    """The file, the growth arguments and the fingerprint of one run."""
    p = dict(FIX, **SOURCE, **NUCLEATION[name], cap_size=CAP, cap_response=response)
    if ablated:
        p.update(ablate_day=ABLATE_DAY, ablate_fraction=ABLATE_FRACTION)
    tag = name.replace(",", "").replace(" ", "_")
    f = OUT / f"{tag}_r{response:g}_{'abl' if ablated else 'ctl'}_s{seed}.json"
    return f, p, fingerprint(seed=seed, S=S, days=WINDOW, params=p)


def needs(*a):
    f, _, fp = spec(*a)
    return stale(f, fp)


def run(name, response, ablated, seed):
    f, p, fp = spec(name, response, ablated, seed)
    hit = cached(f, fp)
    if hit is not None:
        return hit
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S),
                        190 * S, days=WINDOW, rng=np.random.default_rng(seed), **p)
    return save_json(f, dict(condition=name, response=response, ablated=ablated,
                             seed=seed, N=[int(v) for v in tr.N],
                             P=[float(v) for v in tr.P],
                             caps=[int(v) for v in tr.caps]), fp)


def _warm(a):
    run(*a)
    return a


def over_seeds(ctl, abl, f):
    """f(control runs, ablated runs) with a paired bootstrap over seeds."""
    v = f(ctl, abl)
    rng = np.random.default_rng(2026)
    idx = rng.integers(0, len(ctl), (2000, len(ctl)))
    boot = [f([ctl[i] for i in ix], [abl[i] for i in ix]) for ix in idx]
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return dict(value=float(v), lo=float(lo), hi=float(hi))


def mean_of(runs, key, day):
    return float(np.mean([r[key][day - 1] for r in runs]))


def ratio_at(key, day):
    return lambda c, a: mean_of(a, key, day) / max(mean_of(c, key, day), 1e-9)


def cap_ratio(day):
    """Progenitors per site holding any, ablated over control."""
    def f(c, a):
        ca = mean_of(a, "P", day) / max(mean_of(a, "caps", day), 1e-9)
        cc = mean_of(c, "P", day) / max(mean_of(c, "caps", day), 1e-9)
        return ca / max(cc, 1e-9)
    return f


def units_since(c, a):
    """Units formed after the ablation, ablated over control."""
    d0 = ABLATE_DAY - 1
    made = lambda runs: float(np.mean([r["N"][-1] - r["N"][d0 - 1] for r in runs]))
    return made(a) / max(made(c), 1e-9)


def recovery_day(ctl, abl, level=0.95):
    """Days after the ablation at which the cap is first back to `level`."""
    for d in range(ABLATE_DAY, WINDOW + 1):
        if cap_ratio(d)(ctl, abl) >= level:
            return d - ABLATE_DAY
    return None


def ci(r):
    return f"{r['value']:.2f} [{r['lo']:.2f}, {r['hi']:.2f}]"


def _span(rows, key, word="and"):
    v = [r[key]["value"] for r in rows]
    return f"{min(v):.2f} {word} {max(v):.2f}"


def _inside(r, lo, hi):
    """Whether a measured range and a bootstrap interval overlap."""
    return "inside" if r["lo"] <= hi and r["hi"] >= lo else "outside"


def sentences(res):
    """The strings the manuscript quotes, each built from the results once."""
    c, s = res["conditions"], {}
    off = [c[n]["0"] for n in NUCLEATION]
    on = [c[n]["1"] for n in NUCLEATION]
    ref, low = QUOTED
    s["coupled source, pool after ablation"] = (
        f"between {_span(off, 'pool')} of the control's progenitors")
    s["coupled source, sites after ablation"] = (
        f"between {_span(off, 'sites')} of its branching sites")
    s["coupled source, cap the day after"] = (
        f"cap stands at {_span(off, 'cap_next_day', 'to')} of control")
    s["coupled source, units, reference nucleation"] = (
        f"{ci(c[ref]['0']['units'])} of control where nucleation is set as in "
        f"the reference model")
    s["coupled source, units, perfusion not limiting"] = (
        f"{ci(c[low]['0']['units'])} where the perfusion threshold is a quarter")
    s["coupled source, units since ablation"] = (
        f"{c[ref]['0']['units_since']['value']:.2f} and "
        f"{c[low]['0']['units_since']['value']:.2f} of control")
    s["coupled source with response, pool"] = (
        f"pool ends at {_span(on, 'pool', 'to')} of control")
    s["coupled source with response, units"] = (
        f"units at {ci(c[ref]['1']['units'])} and {ci(c[low]['1']['units'])}")
    return s


def blocks(res):
    """The passage of the manuscript this script owns.

    It carries the numbers and the two statements that depend on them, whether
    each measured range falls inside the model's interval, so that neither can
    be left standing by a run that no longer supports it.
    """
    c, s = res["conditions"], sentences(res)
    ref, low = QUOTED
    m0 = MEASURED["no response (Cebrian 2014)"]
    m1 = MEASURED["response (Muthukrishnan 2018)"]
    lo0, hi0 = min(m0["units"]), max(m0["units"])
    lo1, hi1 = min(m1["units"]), max(m1["units"])
    back = [c[n]["0"]["cap_back_days"] for n in NUCLEATION]
    came = [b for b in back if b is not None]
    body = (
        "\nWith no proliferative response, and in all four settings of the "
        "nucleation law, the ablated construct ends the window with "
        f"{s['coupled source, pool after ablation']} and "
        f"{s['coupled source, sites after ablation']}. The day after the ablation "
        f"its {s['coupled source, cap the day after']}, and it "
        + (f"is back to 0.95 of control within {min(came)} to {max(came)} days in "
           f"{len(came)} of the {len(back)} settings" if came else
           "does not return to 0.95 of control within the window in any setting")
        + ". Its units end at "
        f"{s['coupled source, units, reference nucleation']}, and at "
        f"{s['coupled source, units, perfusion not limiting']} of the reference "
        "and the law four times slower; counting only units formed after the "
        f"ablation, {s['coupled source, units since ablation']}. The measured "
        f"{lo0:.2f} to {hi0:.2f} lies {_inside(c[ref]['0']['units'], lo0, hi0)} the "
        f"first interval and {_inside(c[low]['0']['units'], lo0, hi0)} the second. "
        "With division doubled in a depleted cap, the "
        f"{s['coupled source with response, pool']} and the "
        f"{s['coupled source with response, units']} in the same two settings; "
        f"the measured {lo1:.2f} to {hi1:.2f} lies "
        f"{_inside(c[ref]['1']['units'], lo1, hi1)} the first interval and "
        f"{_inside(c[low]['1']['units'], lo1, hi1)} the second. Ratios of means, "
        f"ablated over control, {len(res['seeds'])} paired seeds, 95 % bootstrap "
        "intervals.\n")
    return {"kidney-coupling": body}


def values(res):
    """Single numbers for MANUSCRIPT_A.md, each quoted between its own markers."""
    c = res["conditions"]
    off = [c[n]["0"] for n in NUCLEATION]
    on = [c[n]["1"] for n in NUCLEATION]
    ref, low = QUOTED
    back = [r["cap_back_days"] for r in off]
    came = [b for b in back if b is not None]
    if not came:
        capback = "in none of the settings within the window"
    elif len(came) == len(back):
        capback = f"within {min(came)} to {max(came)} days in all {len(back)} settings"
    else:
        capback = (f"within {min(came)} to {max(came)} days in {len(came)} of "
                   f"{len(back)} settings")
    return {"k.pool": _span(off, "pool", "to"), "k.sites": _span(off, "sites", "to"),
            "k.capback": capback,
            "k.units.ref": ci(c[ref]["0"]["units"]), "k.units.low": ci(c[low]["0"]["units"]),
            "k.since.low": f"{c[low]['0']['units_since']['value']:.2f}",
            "k.resp.pool": _span(on, "pool", "to"),
            "k.resp.units.low": ci(c[low]["1"]["units"])}


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    grid = [(n, r, a, s) for n in NUCLEATION for r in RESPONSE
            for a in (False, True) for s in SEEDS]
    todo = [g for g in grid if needs(*g)]
    if todo:
        n = max(1, min(int(os.environ.get("OGROW_JOBS", 0)) or (os.cpu_count() or 1),
                       len(todo)))
        print(f"running {len(todo)} simulations on {n} worker(s). Every finished "
              f"run is cached,\nso an interrupted pass resumes where it stopped.",
              flush=True)
        with ProcessPoolExecutor(n) as ex:
            for k, _ in enumerate(ex.map(_warm, todo, chunksize=4), 1):
                if k % max(1, len(todo) // 12) == 0 or k == len(todo):
                    print(f"  {k:5d}/{len(todo)}", flush=True)

    res = dict(seeds=list(SEEDS), cap=CAP, rho=RHO, window=WINDOW,
               ablate_day=ABLATE_DAY, ablate_fraction=ABLATE_FRACTION,
               nucleation=NUCLEATION, measured=MEASURED, conditions={})
    end = WINDOW
    print(f"\nA coupled source: one cap of {CAP:.0f} progenitors at the root, dividing "
          f"at {RHO} a day,\nbranching only from two caps, grown for {WINDOW} days. "
          f"On day {ABLATE_DAY} each progenitor\ndies with probability "
          f"{ABLATE_FRACTION}. Ablated over control, {len(SEEDS)} paired seeds, 95 % "
          f"bootstrap interval.\n")
    head = (f"{'nucleation':>32} {'resp':>5} {'units, control':>15} {'pool':>18} "
            f"{'cap +1 d':>9} {'cap back':>9} {'sites':>18} {'units':>18} "
            f"{'units since':>18}")
    print(head)
    for name in NUCLEATION:
        res["conditions"][name] = {}
        for resp in RESPONSE:
            ctl = [run(name, resp, False, s) for s in SEEDS]
            abl = [run(name, resp, True, s) for s in SEEDS]
            back = recovery_day(ctl, abl)
            row = dict(
                control_units=mean_of(ctl, "N", end),
                control_sites=mean_of(ctl, "caps", end),
                pool=over_seeds(ctl, abl, ratio_at("P", end)),
                cap_next_day=over_seeds(ctl, abl, cap_ratio(ABLATE_DAY + 1)),
                cap_back_days=back,
                sites=over_seeds(ctl, abl, ratio_at("caps", end)),
                units=over_seeds(ctl, abl, ratio_at("N", end)),
                units_since=over_seeds(ctl, abl, units_since),
                days=list(range(1, WINDOW + 1)),
                pool_by_day=[ratio_at("P", d)(ctl, abl) for d in range(1, WINDOW + 1)],
                cap_by_day=[cap_ratio(d)(ctl, abl) for d in range(1, WINDOW + 1)],
                sites_by_day=[ratio_at("caps", d)(ctl, abl) for d in range(1, WINDOW + 1)],
                units_by_day=[ratio_at("N", d)(ctl, abl) for d in range(1, WINDOW + 1)])
            res["conditions"][name][f"{resp:g}"] = row
            print(f"{name:>32} {resp:>5.0f} {row['control_units']:>15.0f} "
                  f"{ci(row['pool']):>18} {row['cap_next_day']['value']:>9.2f} "
                  f"{(str(back) + ' d') if back is not None else 'never':>9} "
                  f"{ci(row['sites']):>18} {ci(row['units']):>18} "
                  f"{ci(row['units_since']):>18}")
        print()

    m0, m1 = MEASURED["no response (Cebrian 2014)"], MEASURED["response (Muthukrishnan 2018)"]
    print("Measured, ablated over control:")
    print(f"  without a proliferative response (Cebrian 2014): cells per tip "
          f"{m0['cap'][0]:.2f}, {m0['cap'][1]:.2f}, {m0['cap'][2]:.2f} one, two and "
          f"three days\n  after tamoxifen; tips {m0['tips'][0]:.2f}, "
          f"{m0['tips'][1]:.2f}, {m0['tips'][2]:.2f}; total progenitors growing as "
          f"fast as control;\n  glomeruli {m0['units'][0]:.2f} at birth and "
          f"{m0['units'][1]:.2f} in the adult.")
    print(f"  with one (Muthukrishnan 2018): progenitors {m1['pool'][0]:.2f} after a "
          f"day and {m1['pool'][1]:.2f} after two;\n  glomeruli {m1['units'][0]:.2f} "
          f"and {m1['units'][1]:.2f} three and ten days after birth.")
    print("\n  Model days are not embryonic days, so the sequence is compared and "
          "the endpoints,\n  not the timing. The extension was written after these "
          "numbers were known.")
    res["strings"] = sentences(res)
    res["values"] = values(res)
    json.dump(res, open(ROOT / "results" / "kidney_coupling.json", "w"), indent=1)
    for _d in [d for d in DOCS if write(d, {**blocks(res), **values(res)})]:
        print(f"\nupdated the passage of {_d.relative_to(ROOT)} that quotes these "
              "numbers.\nCommit that file before `make paper`, or the stamp will "
              "record an uncommitted change.")
