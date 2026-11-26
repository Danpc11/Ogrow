#!/usr/bin/env python
"""Three predictions of the model against measurements that already exist.

Everything else in the paper is the model read against itself. These are the
three places where a number the model gives can be set beside a number someone
measured, without fitting anything.

    ablation  Seed a fraction of the reference pool and read what is left of
              the construct. In mouse, removing pancreatic progenitors leaves
              a pancreas smaller in proportion and removing hepatic ones does
              not limit the liver (Stanger, Tanaka and Melton 2007); removing
              about 40 % of nephron progenitors leaves 10 to 20 % fewer
              nephrons (Cebrian et al. 2014 report a permanent deficit).
    pairs     Two constructs grown from identical parameters differ only by
              their realisation. So do the two kidneys of one person. The
              coefficient of variation across seeds is set beside the
              left-right variation of kidney volume in 48 people (TRUSTED).
    size      The model's spread across seeds grows with the domain. Whether
              the left-right difference grows with kidney size is read off
              the same 48 pairs.

What each comparison can and cannot carry is said where the manuscript uses
it. Two cautions belong here because they are about this script. The three
regimes were named with the first of these experiments in view, so agreeing
with it is consistency and not a prediction; the kidney ablation and the pairs
were not. And the model's reading is a unit count or a mass in model units,
while the measurements are organ weights, nephron counts and CT volumes.

The script also writes the sentences of the manuscript that quote these
numbers, between `<!-- auto:NAME -->` and `<!-- /auto:NAME -->` markers in
manuscript/MANUSCRIPT.md, so that no number there is typed by hand. A second
run writes the same text and leaves the file untouched. Run it once before
committing a change that moves these numbers, so that `make paper` finds the
manuscript already in step and the stamp stays clean.

Usage:  python scripts/36_validation.py
        OGROW_JOBS=n   worker count
"""
import csv
import json
import os
import re
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from organgrow import sdf
from organgrow.config import ROOT
from organgrow.closed import grow
from organgrow.runcache import fingerprint, cached, save_json, stale

OUT = ROOT / "results" / "validation"
DOC = ROOT / "manuscript" / "MANUSCRIPT.md"
PAIRS = ROOT / "source_data" / "Fig6_kidney_invariants.csv"

S, DAYS, NU, SEEDS = 0.6, 140, 8.0, tuple(range(300, 324))
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=6000)
REG = {"liver": dict(P0=2600.0, rho=0.10, delta=0.10, q_set=0.6, renewal=1.0),
       "kidney": dict(P0=300.0, rho=0.00, delta=0.05, q_set=0.6, renewal=0.0),
       "pancreas": dict(P0=60.0, rho=0.00, delta=0.00, q_set=0.6, renewal=0.0)}
#: the share of the reference pool that is seeded. A half and a quarter for
#: the two tissues of the mouse experiment; 0.6 for the kidney because the
#: published ablations removed about 40 %, and a half beside it.
KEPT = {"pancreas": (0.5, 0.25), "liver": (0.5, 0.25), "kidney": (0.6, 0.5)}
NBOOT, BOOT_SEED = 4000, 2026


def spec(reg, kept, seed):
    """The file, the growth parameters and the fingerprint of one run."""
    p = dict(REG[reg]); p["P0"] = p["P0"] * float(kept)
    f = OUT / f"{reg}_kept{kept:g}_s{seed}.json"
    return f, p, fingerprint(seed=seed, S=S, days=DAYS,
                             params={**FIX, "nu": NU, **p})


def needs(*a):
    f, _, fp = spec(*a)
    return stale(f, fp)


def run(reg, kept, seed):
    f, p, fp = spec(reg, kept, seed)
    hit = cached(f, fp)
    if hit is not None:
        return hit
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S),
                        190 * S, days=DAYS, rng=np.random.default_rng(seed),
                        nu=NU, **FIX, **p)
    return save_json(f, dict(regime=reg, kept=kept, seed=seed, N=len(U),
                             mass=float(m.sum())), fp)


def _warm(a):
    run(*a)
    return a


def ratio(ref, cut):
    """Ratio of the means, with a paired bootstrap over seeds.

    The same seeds grow the reference and the reduced pool, so a resample
    draws seeds and takes both runs of each.
    """
    ref, cut = np.asarray(ref, float), np.asarray(cut, float)
    rng = np.random.default_rng(BOOT_SEED)
    idx = rng.integers(0, len(ref), (NBOOT, len(ref)))
    boot = cut[idx].mean(1) / ref[idx].mean(1)
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return dict(value=float(cut.mean() / ref.mean()), lo=float(lo), hi=float(hi))


def cv(a):
    a = np.asarray(a, float)
    return float(a.std(ddof=1) / a.mean())


def paired_variation(left, right):
    """How much two organs of one person differ, as a coefficient of variation.

    The difference of the logs of a pair has twice the variance of one organ's
    own deviation, hence the division by root two. The mean difference between
    sides is anatomy, not chance, and is removed before the spread is taken;
    the value with it left in is returned as well.
    """
    left, right = np.asarray(left, float), np.asarray(right, float)
    d = np.log(left) - np.log(right)
    mean_pair = (left + right) / 2.0
    dev = np.abs(d - d.mean())
    from scipy.stats import spearmanr
    rho, p = spearmanr(mean_pair, dev)
    return dict(pairs=int(len(d)),
                side_offset=float(np.expm1(d.mean())),
                cv=float(d.std(ddof=1) / np.sqrt(2.0)),
                cv_with_offset=float(np.sqrt(np.mean(d ** 2) / 2.0)),
                cv_between=float(((np.log(left) + np.log(right)) / 2).std(ddof=1)),
                size_rho=float(rho), size_p=float(p),
                smallest=float(mean_pair.min()), largest=float(mean_pair.max()))


def read_pairs(path=PAIRS):
    """Left and right volumes per donor, or None when the table is absent."""
    if not path.exists():
        return None
    by = {}
    for r in csv.DictReader(open(path, newline="")):
        by.setdefault(r["donor"], {})[r["side"]] = float(r["volume_cm3"])
    both = [v for v in by.values() if "L" in v and "R" in v]
    return ([v["L"] for v in both], [v["R"] for v in both]) if both else None


def fill(text, name, body):
    """Replace what stands between the two markers of one named block."""
    pat = re.compile(rf"(<!-- auto:{re.escape(name)} -->)(.*?)(<!-- /auto:{re.escape(name)} -->)",
                     re.S)
    if not pat.search(text):
        return text
    return pat.sub(lambda m: m.group(1) + body + m.group(3), text)


def ci(r, nd=2):
    return f"{r['value']:.{nd}f} [{r['lo']:.{nd}f}, {r['hi']:.{nd}f}]"


def sentences(res):
    """The strings the manuscript quotes, each built from the results once.

    `scripts/29_check_manuscript.py` requires every one of them to be present,
    and the two blocks written into the manuscript are assembled from the same
    strings, so the check and the prose cannot disagree.
    """
    a, s = res["ablation"], {}
    s["pancreas-like, half the pool, mass"] = f"{ci(a['pancreas']['0.5']['mass'])} of its mass"
    s["pancreas-like, half the pool, units"] = f"{a['pancreas']['0.5']['N']['value']:.2f} of its units"
    s["liver-like, half the pool, mass"] = f"{ci(a['liver']['0.5']['mass'])} of its mass"
    s["pancreas-like, a quarter of the pool, mass"] = (
        f"keeps {a['pancreas']['0.25']['mass']['value']:.2f} of its mass on a quarter of the pool")
    s["liver-like, a quarter of the pool, mass"] = (
        f"and the liver-like one {a['liver']['0.25']['mass']['value']:.2f}")
    k = a["kidney"]["0.6"]["N"]
    s["kidney-like, 40 % of the pool removed, units"] = f"{ci(k)} of the units"
    s["kidney-like, 40 % of the pool removed, deficit"] = (
        f"a deficit of {100 * (1 - k['value']):.0f} % "
        f"[{100 * (1 - k['hi']):.0f}, {100 * (1 - k['lo']):.0f}]")
    c = res["contingent"]["kidney"]
    s["kidney-like, variation across seeds, mass"] = (
        f"coefficient of variation of {c['mass']:.3f} in mass")
    s["kidney-like, variation across seeds, units"] = f"and of {c['N']:.3f} in unit count"
    p = res.get("pairs")
    if p:
        s["paired kidneys, left-right variation"] = (
            f"coefficient of variation of {p['cv']:.3f} over {p['pairs']} pairs")
        s["paired kidneys, mean difference between sides"] = (
            f"the left smaller by {abs(100 * p['side_offset']):.1f} % on average"
            if p["side_offset"] < 0 else
            f"the left larger by {100 * p['side_offset']:.1f} % on average")
        s["paired kidneys, variation against size"] = (
            f"$\\rho = {p['size_rho']:+.2f}$, $P = {p['size_p']:.2f}$")
    return s


def blocks(res):
    """The two passages of the manuscript this script owns."""
    s, p = sentences(res), res.get("pairs")
    body = (
        "\nSeeded with half its reference pool, the pancreas-like construct ends with "
        f"{s['pancreas-like, half the pool, mass']} and "
        f"{s['pancreas-like, half the pool, units']}, and the liver-like construct with "
        f"{s['liver-like, half the pool, mass']} (ratio of means over 24 paired seeds, "
        "95 % bootstrap interval). The pancreas-like construct "
        f"{s['pancreas-like, a quarter of the pool, mass']}, "
        f"{s['liver-like, a quarter of the pool, mass']}. With 40 % of the "
        "kidney-like pool removed the construct ends with "
        f"{s['kidney-like, 40 % of the pool removed, units']}, "
        f"{s['kidney-like, 40 % of the pool removed, deficit']}. Grown 24 times from "
        "identical parameters, the kidney-like construct has a "
        f"{s['kidney-like, variation across seeds, mass']} "
        f"{s['kidney-like, variation across seeds, units']}.")
    if p:
        body += (
            " In the human cohort the volumes of a person's two kidneys differ with a "
            f"{s['paired kidneys, left-right variation']}, once the mean difference "
            f"between sides is removed ({s['paired kidneys, mean difference between sides']}); "
            "that difference does not rise detectably with kidney size "
            f"({s['paired kidneys, variation against size']}, Spearman, over volumes "
            f"from {p['smallest']:.0f} to {p['largest']:.0f} cm$^3$).")
    body += "\n"
    k = res["ablation"]["kidney"]["0.6"]["N"]
    short = (
        "Halving the seeded pool leaves "
        f"{res['ablation']['pancreas']['0.5']['mass']['value']:.2f} of the mass in the "
        "pancreas-like regime and "
        f"{res['ablation']['liver']['0.5']['mass']['value']:.2f} in the liver-like one, as "
        "progenitor ablation does in mouse; removing 40 % of the kidney-like pool "
        f"costs {100 * (1 - k['value']):.0f} % of the units, against 10 to 20 % of "
        "nephrons measured")
    if p:
        short += (
            "; and two constructs grown alike differ by "
            f"{res['contingent']['kidney']['mass']:.3f} in mass, against "
            f"{p['cv']:.3f} between a person's two kidneys")
    short += "."
    return {"validation": body, "validation-abstract": short}


def write_blocks(res, doc=DOC):
    """Bring the manuscript's two passages in step. True if the file changed."""
    if not doc.exists():
        return False
    old = doc.read_text(encoding="utf-8")
    new = old
    for name, body in blocks(res).items():
        new = fill(new, name, body)
    if new == old:
        return False
    doc.write_text(new, encoding="utf-8")
    return True


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    grid = [(reg, kept, s) for reg in REG for kept in (1.0,) + KEPT[reg] for s in SEEDS]
    todo = [g for g in grid if needs(*g)]
    if todo:
        n = max(1, min(int(os.environ.get("OGROW_JOBS", 0)) or (os.cpu_count() or 1),
                       len(todo)))
        print(f"running {len(todo)} simulations on {n} worker(s)", flush=True)
        with ProcessPoolExecutor(n) as ex:
            list(ex.map(_warm, todo))

    res = dict(seeds=list(SEEDS), kept={k: list(v) for k, v in KEPT.items()},
               ablation={}, contingent={})
    print("\nWhat is left of the construct when a share of the reference pool is "
          "seeded,\nratio of means over 24 paired seeds with a 95 % bootstrap "
          "interval:\n")
    print(f"{'regime':>9} {'pool kept':>10} {'units':>20} {'mass':>20}  elasticity of mass")
    for reg in ("pancreas", "kidney", "liver"):
        ref = [run(reg, 1.0, s) for s in SEEDS]
        res["contingent"][reg] = dict(N=cv([r["N"] for r in ref]),
                                      mass=cv([r["mass"] for r in ref]),
                                      N_mean=float(np.mean([r["N"] for r in ref])))
        res["ablation"][reg] = {}
        for kept in KEPT[reg]:
            cut = [run(reg, kept, s) for s in SEEDS]
            rn = ratio([r["N"] for r in ref], [r["N"] for r in cut])
            rm = ratio([r["mass"] for r in ref], [r["mass"] for r in cut])
            res["ablation"][reg][f"{kept:g}"] = dict(N=rn, mass=rm)
            print(f"{reg:>9} {kept:>10.2f} {ci(rn):>20} {ci(rm):>20}  "
                  f"{np.log(rm['value']) / np.log(kept):+.2f}")
    print("\n  Published: a pancreas smaller in proportion to the progenitors removed "
          "and a liver\n  that is not limited by them (Stanger et al. 2007); 10 to "
          "20 % fewer nephrons after\n  about 40 % of nephron progenitors are "
          "removed.")

    print("\nVariation across 24 seeds grown from identical parameters:\n")
    print(f"{'regime':>9} {'CV of units':>12} {'CV of mass':>11}")
    for reg in ("pancreas", "kidney", "liver"):
        c = res["contingent"][reg]
        print(f"{reg:>9} {c['N']:>12.3f} {c['mass']:>11.3f}")

    lr = read_pairs()
    if lr is None:
        print(f"\nno table of paired kidneys at {PAIRS.relative_to(ROOT)}; "
              "the comparison with them is skipped")
    else:
        p = res["pairs"] = paired_variation(*lr)
        print(f"\nThe two kidneys of one person, {p['pairs']} pairs (TRUSTED, CT "
              f"volume):\n"
              f"  mean difference between sides        {100 * p['side_offset']:+.1f} % (left against right)\n"
              f"  left-right CV, that difference removed {p['cv']:.3f}\n"
              f"  left-right CV, left in                 {p['cv_with_offset']:.3f}\n"
              f"  between people                         {p['cv_between']:.3f}\n"
              f"  left-right difference against size     rho = {p['size_rho']:+.2f}, "
              f"P = {p['size_p']:.2f}, over {p['smallest']:.0f} to {p['largest']:.0f} cm3")
        print("\n  The model's spread across seeds rises with the domain "
              "(PAPER_NUMBERS.md, Organ size).\n  Whether the pairs show the same is "
              "the last line, and it is not the same\n  experiment: people of "
              "different size are not one construct in a larger domain.\n  The CT "
              "volumes also carry segmentation error, which this does not separate.")

    res["strings"] = sentences(res)
    json.dump(res, open(ROOT / "results" / "validation.json", "w"), indent=1)
    if write_blocks(res):
        print(f"\nupdated the passages of {DOC.relative_to(ROOT)} that quote these "
              "numbers.\nCommit that file before `make paper`, or the stamp will "
              "record an uncommitted change.")
