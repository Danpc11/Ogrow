#!/usr/bin/env python
"""Predictions for experiments that have not been done.

Sections 3.10 and 3.11 set the model beside measurements that exist. A model
that is only ever compared with what is known cannot be wrong, so this script
commits it to numbers for experiments nobody has reported, in the same
systems and with the same method as the ones already compared: progenitor
ablation in the mouse embryo, at doses and times that were not used.

    dose      The pancreas-like and liver-like constructs seeded with shares
              of the pool between and beyond the published ones. Published:
              0.92 and 0.39 of pancreatic progenitors, at most 0.35 of hepatic
              ones. Predicted here: 0.75 and 0.60 for the pancreas, 0.75 and
              0.15 for the liver.
    depth     The coupled kidney-like source of scripts/37_kidney_coupling.py
              with 20 % and 60 % of its progenitors killed, where both
              published ablations killed about 40 %.
    timing    The same source with 40 % killed early and late in its window,
              where the two published ablations differ in timing and in
              everything else.

Each prediction is stated with the outcome that would refute it. The coupled
source is read in the setting where perfusion limits nucleation least, which
is the one that came closest to the measured deficit, and with and without a
proliferative response, because which of the two a given experiment will show
is not something the model says.

These numbers are stamped with the commit that produced them. Once one of
the experiments is done, the number that stood here before it is the
prediction, whatever a later version of the model gives.

Usage:  python scripts/38_predictions.py
        OGROW_JOBS=n   worker count
"""
import json
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from organgrow import sdf
from organgrow.autotext import write
from organgrow.config import ROOT
from organgrow.closed import grow
from organgrow.runcache import fingerprint, cached, save_json, stale

OUT = ROOT / "results" / "predictions"
DOCS = (ROOT / "manuscript" / "MANUSCRIPT_A.md",)

S, SEEDS = 0.6, tuple(range(300, 324))
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=6000)

# ---- dose: the reference regimes, as in scripts/36_validation.py -----------
DAYS, NU = 140, 8.0
REG = {"liver": dict(P0=2600.0, rho=0.10, delta=0.10, q_set=0.6, renewal=1.0),
       "pancreas": dict(P0=60.0, rho=0.00, delta=0.00, q_set=0.6, renewal=0.0)}
#: shares of the pool not used in any published ablation
DOSE = {"pancreas": (0.75, 0.60), "liver": (0.75, 0.15)}

# ---- depth and timing: the coupled source, as in scripts/37 ----------------
CAP, RHO, WINDOW = 60.0, 0.05, 120
SOURCE = dict(P0=CAP, rho=RHO, delta=0.0, renewal=0.0)
NUCLEATION = dict(q_set=0.15, nu=2.0)      # perfusion limits nucleation least
RESPONSE = (0.0, 1.0)
#: (day, share killed). The published ablations killed about 0.4; day 45 is
#: the day scripts/37 uses, and 20 and 80 are early and late in the window.
ABLATIONS = ((45, 0.2), (45, 0.4), (45, 0.6), (20, 0.4), (80, 0.4))


def dose_spec(reg, kept, seed):
    p = dict(REG[reg]); p["P0"] = p["P0"] * float(kept)
    f = OUT / f"dose_{reg}_kept{kept:g}_s{seed}.json"
    return f, p, fingerprint(seed=seed, S=S, days=DAYS, params={**FIX, "nu": NU, **p})


def dose_run(reg, kept, seed):
    f, p, fp = dose_spec(reg, kept, seed)
    hit = cached(f, fp)
    if hit is not None:
        return hit
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S),
                        190 * S, days=DAYS, rng=np.random.default_rng(seed),
                        nu=NU, **FIX, **p)
    return save_json(f, dict(N=len(U), mass=float(m.sum())), fp)


def kid_spec(response, day, share, seed):
    p = dict(FIX, **SOURCE, **NUCLEATION, cap_size=CAP, cap_response=response)
    if share > 0:
        p.update(ablate_day=day, ablate_fraction=share)
    f = OUT / f"kidney_r{response:g}_d{day}_f{share:g}_s{seed}.json"
    return f, p, fingerprint(seed=seed, S=S, days=WINDOW, params=p)


def kid_run(response, day, share, seed):
    f, p, fp = kid_spec(response, day, share, seed)
    hit = cached(f, fp)
    if hit is not None:
        return hit
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S),
                        190 * S, days=WINDOW, rng=np.random.default_rng(seed), **p)
    return save_json(f, dict(N=len(U), sites=int(tr.caps[-1])), fp)


def _warm(a):
    (dose_run if a[0] == "dose" else kid_run)(*a[1:])
    return a


def _needs(a):
    f, _, fp = (dose_spec if a[0] == "dose" else kid_spec)(*a[1:])
    return stale(f, fp)


def ratio(ref, cut):
    """Ratio of the means with a paired bootstrap over seeds."""
    ref, cut = np.asarray(ref, float), np.asarray(cut, float)
    rng = np.random.default_rng(2026)
    idx = rng.integers(0, len(ref), (4000, len(ref)))
    lo, hi = np.percentile(cut[idx].mean(1) / ref[idx].mean(1), [2.5, 97.5])
    return dict(value=float(cut.mean() / ref.mean()), lo=float(lo), hi=float(hi))


def ci(r):
    return f"{r['value']:.2f} [{r['lo']:.2f}, {r['hi']:.2f}]"


def sentences(res):
    """The strings the manuscript quotes, each built from the results once."""
    d, k, s = res["dose"], res["kidney"], {}
    s["prediction, pancreas at 0.75 of its progenitors"] = (
        f"a pancreas at {ci(d['pancreas']['0.75'])} of control weight")
    s["prediction, pancreas at 0.60"] = f"and at {ci(d['pancreas']['0.6'])} with 0.60"
    s["prediction, liver at 0.75"] = (
        f"a liver at {ci(d['liver']['0.75'])} of control weight")
    s["prediction, liver at 0.15"] = f"and at {ci(d['liver']['0.15'])} with 0.15"
    for resp, word in (("0", "without"), ("1", "with")):
        r = k[resp]
        s[f"prediction, nephrons after 20 % and 60 % ablation, {word} a response"] = (
            f"{ci(r['45|0.2'])} after 20 % and {ci(r['45|0.6'])} after 60 %")
        s[f"prediction, early against late ablation, {word} a response"] = (
            f"{ci(r['20|0.4'])} early and {ci(r['80|0.4'])} late")
    return s


def blocks(res):
    """The passage of the manuscript this script owns."""
    s, k = sentences(res), res["kidney"]
    deep = k["0"]["45|0.6"]
    prop = 1.0 - 0.6
    body = (
        "\nFor the pancreas, an ablation that leaves 0.75 of the progenitors is "
        f"predicted to give {s['prediction, pancreas at 0.75 of its progenitors']}, "
        f"{s['prediction, pancreas at 0.60']}; a pancreas of normal weight at either "
        "dose refutes it. For the liver, one that leaves 0.75 is predicted to give "
        f"{s['prediction, liver at 0.75']}, {s['prediction, liver at 0.15']}: "
        "compensation that is complete after a shallow ablation and fails after a "
        "deep one. A liver of normal weight after 85 % of its progenitors are removed "
        "refutes it. For the kidney, in the coupled source and where the surviving "
        "progenitors do not divide faster, the units are predicted at "
        f"{s['prediction, nephrons after 20 % and 60 % ablation, without a response']} "
        "of control. "
        + (f"The deficit is not in proportion at depth: a proportional deficit would "
           f"leave {prop:.2f} after 60 %, which is outside the interval. "
           if not (deep["lo"] <= prop <= deep["hi"]) else
           f"A proportional deficit would leave {prop:.2f} after 60 %, which the "
           f"interval includes. ")
        + "Where they do divide faster, the units are predicted at "
        f"{s['prediction, nephrons after 20 % and 60 % ablation, with a response']}. "
        "For the same 40 % killed at two times, the units are predicted at "
        f"{s['prediction, early against late ablation, without a response']} in the "
        "window without a response, and at "
        f"{s['prediction, early against late ablation, with a response']} with one. "
        "Ratios of means to the unablated control, 24 paired seeds, 95 % bootstrap "
        "intervals; model mass for pancreas and liver, units for kidney.\n")
    return {"predictions": body}


def values(res):
    """Single numbers for MANUSCRIPT_A.md, each quoted between its own markers."""
    d, k = res["dose"], res["kidney"]
    return {"p.panc75": ci(d["pancreas"]["0.75"]), "p.panc60": ci(d["pancreas"]["0.6"]),
            "p.liv75": ci(d["liver"]["0.75"]), "p.liv15": ci(d["liver"]["0.15"]),
            "p.kid20": ci(k["0"]["45|0.2"]), "p.kid60": ci(k["0"]["45|0.6"]),
            "p.kid20r": ci(k["1"]["45|0.2"]), "p.kid60r": ci(k["1"]["45|0.6"]),
            "p.early": ci(k["0"]["20|0.4"]), "p.late": ci(k["0"]["80|0.4"])}


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    grid = [("dose", reg, kept, s) for reg in DOSE for kept in (1.0,) + DOSE[reg]
            for s in SEEDS]
    grid += [("kid", r, 0, 0.0, s) for r in RESPONSE for s in SEEDS]
    grid += [("kid", r, d, f, s) for r in RESPONSE for d, f in ABLATIONS for s in SEEDS]
    todo = [g for g in grid if _needs(g)]
    if todo:
        n = max(1, min(int(os.environ.get("OGROW_JOBS", 0)) or (os.cpu_count() or 1),
                       len(todo)))
        print(f"running {len(todo)} simulations on {n} worker(s)", flush=True)
        with ProcessPoolExecutor(n) as ex:
            list(ex.map(_warm, todo, chunksize=4))

    res = dict(seeds=list(SEEDS), dose={}, kidney={}, nucleation=NUCLEATION,
               ablations=[list(a) for a in ABLATIONS])
    print("\nPredictions for ablations that have not been reported. Ratio to the "
          "unablated control,\n24 paired seeds, 95 % bootstrap interval.\n")
    print("Dose: mass of the construct against the share of its pool that is left")
    for reg in ("pancreas", "liver"):
        ref = [dose_run(reg, 1.0, s)["mass"] for s in SEEDS]
        res["dose"][reg] = {}
        for kept in DOSE[reg]:
            r = ratio(ref, [dose_run(reg, kept, s)["mass"] for s in SEEDS])
            res["dose"][reg][f"{kept:g}"] = r
            print(f"  {reg + '-like':>14}  pool left {kept:.2f}   {ci(r)}")
    print("\nDepth and timing: units of the coupled source, perfusion limiting least")
    for resp in RESPONSE:
        ctl = [kid_run(resp, 0, 0.0, s)["N"] for s in SEEDS]
        res["kidney"][f"{resp:g}"] = {}
        for day, share in ABLATIONS:
            r = ratio(ctl, [kid_run(resp, day, share, s)["N"] for s in SEEDS])
            res["kidney"][f"{resp:g}"][f"{day}|{share:g}"] = r
            what = "published dose and day" if (day, share) == (45, 0.4) else "prediction"
            print(f"  response {resp:.0f}  day {day:>3}  killed {share:.0%}   {ci(r)}   {what}")
        print()
    print("  What refutes each: a pancreas of normal weight at either dose; a liver "
          "of normal weight\n  after 85 % of its progenitors are removed; a nephron "
          "deficit in proportion after a\n  60 % ablation without a proliferative "
          "response; the same deficit whenever the\n  ablation is done.")
    res["strings"] = sentences(res)
    res["values"] = values(res)
    json.dump(res, open(ROOT / "results" / "predictions.json", "w"), indent=1)
    changed = [d for d in DOCS if write(d, {**blocks(res), **values(res)})]
    for d in changed:
        print(f"\nupdated the passage of {d.relative_to(ROOT)} that quotes these "
              "numbers.\nCommit that file before `make paper`, or the stamp will "
              "record an uncommitted change.")
