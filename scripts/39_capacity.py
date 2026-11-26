#!/usr/bin/env python
"""Two tests a correlation cannot make.

1. Does nearness to what perfusion supports move the elasticity, or does size?

   Across the factorial of scripts/34_regimes.py the elasticity of the unit
   count to the seeded pool falls with the size a source reaches. Size and
   elasticity come from the same runs, so that is a correlation, and it does
   not say whether the tissue answers less to its pool because it is large or
   because it is near a limit set by perfusion. The two readings part when the
   perfusion threshold is changed inside one source: a higher threshold makes
   the construct smaller and, at the same time, brings it nearer the limit. If
   size were what mattered the elasticity would rise; if nearness to the limit
   is what matters it falls.

   The threshold turns out to move size and limit together, so it separates
   the two readings only weakly. The domain separates them cleanly: the same
   source grown in a smaller or a larger domain has the same pool and a
   different limit. A source that its pool holds far below the limit reaches
   nearly the same size in both, and the question is whether its elasticity
   stays, as size predicts, or moves, as nearness to the limit predicts.

   The limit is measured independently of the source, as the units reached at
   each threshold by a source whose pool cannot run out: a renewing source
   with twenty times the largest pool of the factorial. Whether that is
   enough is checked and not assumed: the same count is taken at half and at
   twice that pool, and a capacity is a capacity only if it has stopped rising
   with the pool. It is also taken at twice the horizon in the smallest and
   the largest domain, because a larger domain takes longer to fill and a
   fixed horizon could favour the small one. It is an operational capacity of
   this domain under these rules, and it moves with the threshold, the tree
   and the domain. It is not a constant of the tissue.

2. Is the coupling of branching to the progenitors needed for the kidney?

   scripts/37_kidney_coupling.py adds three things at once to a kidney-like
   source: division, a fixed window, and branching that waits for its cap.
   The reduced model here keeps the first two and drops the third, so that
   what the coupling alone contributes can be read off: the same ablation, in
   the same two settings of the nucleation law, with the tree free to branch
   whatever the progenitors do.

Usage:  python scripts/39_capacity.py
        OGROW_JOBS=n   worker count
        OGROW_SEEDS=n  fewer seeds, for a quick look; the paper uses 24
"""
import json
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.stats import spearmanr

from organgrow import sdf
from organgrow.autotext import write
from organgrow.config import ROOT
from organgrow.closed import grow
from organgrow.runcache import fingerprint, cached, save_json, stale

OUT = ROOT / "results" / "capacity"
DOC = ROOT / "manuscript" / "MANUSCRIPT_A.md"

S, DAYS, NU, H = 0.6, 140, 8.0, 0.3
SEEDS = tuple(range(300, 300 + int(os.environ.get("OGROW_SEEDS", 24))))
FIX = dict(R_unit=650.0, units_per_branch=650 / 600.0, max_vessels=6000)
REG = {"liver": dict(P0=2600.0, rho=0.10, delta=0.10, renewal=1.0),
       "kidney": dict(P0=300.0, rho=0.00, delta=0.05, renewal=0.0),
       "pancreas": dict(P0=60.0, rho=0.00, delta=0.00, renewal=0.0)}
#: a pool that cannot run out: the renewing source at twenty times its pool.
#: At ten times the count was still rising with the pool; see POOL_CHECK.
SATURATING = dict(P0=52000.0, rho=0.10, delta=0.10, renewal=1.0)
#: the same source at half and at twice that pool, at the reference threshold
POOL_CHECK = (0.5, 1.0, 2.0)
#: and at twice the horizon, in the smallest and the largest domain
LONG_DAYS = 280
THRESHOLDS = (0.3, 0.4, 0.5, 0.6, 0.7)
#: the domain, as a linear scale; 0.6 is the reference. The liver-like source
#: is left out: it sits near its capacity at every scale and its elasticity is
#: within noise of zero, so it cannot show a shift either way.
SCALES = (0.4, 0.6, 0.8)
DOMAIN_SOURCES = ("kidney", "pancreas")

# ---- the kidney, with and without the coupling ------------------------------
CAP, RHO, WINDOW, ABLATE_DAY, ABLATE_FRACTION = 60.0, 0.05, 120, 45, 0.4
SOURCE = dict(P0=CAP, rho=RHO, delta=0.0, renewal=0.0)
NUCLEATION = {"reference": dict(q_set=0.6, nu=8.0),
              "threshold quartered, law slowed": dict(q_set=0.15, nu=2.0)}


def spec(kind, *a):
    if kind == "size":
        reg, q, factor, seed = a[:4]
        scale = a[4] if len(a) > 4 else S
        days = a[5] if len(a) > 5 else DAYS
        p = dict(SATURATING if reg == "saturating" else REG[reg], nu=NU, q_set=q)
        p["P0"] = p["P0"] * float(factor)
        tag = ("" if scale == S else f"_S{scale:g}") + ("" if days == DAYS else f"_d{days}")
        f = OUT / f"size_{reg}_q{q:g}_f{factor:.4f}{tag}_s{seed}.json"
        return f, dict(FIX, **p), days, fingerprint(seed=seed, S=scale, days=days, params=dict(FIX, **p))
    name, coupled, ablated, seed = a
    p = dict(FIX, **SOURCE, **NUCLEATION[name])
    if coupled:
        p.update(cap_size=CAP)
    if ablated:
        p.update(ablate_day=ABLATE_DAY, ablate_fraction=ABLATE_FRACTION)
    tag = name.replace(",", "").replace(" ", "_")
    f = OUT / f"kidney_{tag}_{'coupled' if coupled else 'free'}_{'abl' if ablated else 'ctl'}_s{seed}.json"
    return f, p, WINDOW, fingerprint(seed=seed, S=S, days=WINDOW, params=p)


def run(kind, *a):
    f, p, days, fp = spec(kind, *a)
    hit = cached(f, fp)
    if hit is not None:
        return hit
    sc = a[4] if (kind == "size" and len(a) > 4) else S
    seed = a[3]
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / sc, y / sc, z / sc),
                        190 * sc, days=days, rng=np.random.default_rng(seed), **p)
    n_day = list(tr.N)
    late = float((n_day[-1] - n_day[-41]) / max(n_day[-41], 1)) if len(n_day) > 40 else float("nan")
    return save_json(f, dict(N=len(U), P=float(tr.P[-1]), tips=int(len(vt.tips())),
                             caps=int(tr.caps[-1]), vessels=int(tr.vessels[-1]),
                             late=late), fp)


def _warm(a):
    run(*a)
    return a


def _needs(a):
    f, _, _, fp = spec(*a)
    return stale(f, fp)


def boot(f, n, B=4000):
    """A statistic of seed indices, with a 95 % bootstrap interval."""
    rng = np.random.default_rng(2026)
    v = f(np.arange(n))
    b = [f(ix) for ix in rng.integers(0, n, (B, n))]
    lo, hi = np.percentile(b, [2.5, 97.5])
    return dict(value=float(v), lo=float(lo), hi=float(hi))


def ci(r, nd=2):
    return f"{r['value']:.{nd}f} [{r['lo']:.{nd}f}, {r['hi']:.{nd}f}]"


def sci(r):
    return f"{r['value']:+.2f} [{r['lo']:+.2f}, {r['hi']:+.2f}]"


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    up, dn = float(np.exp(H)), float(np.exp(-H))
    grid = [("size", "saturating", q, 1.0, s) for q in THRESHOLDS for s in SEEDS]
    grid += [("size", reg, q, f, s) for reg in REG for q in THRESHOLDS
             for f in (1.0, up, dn) for s in SEEDS]
    grid += [("size", "saturating", 0.6, 1.0, s, sc) for sc in SCALES if sc != S for s in SEEDS]
    grid += [("size", "saturating", 0.6, pf, s) for pf in POOL_CHECK if pf != 1.0 for s in SEEDS]
    grid += [("size", "saturating", 0.6, 1.0, s, sc, LONG_DAYS) for sc in (SCALES[0], SCALES[-1])
             for s in SEEDS]
    grid += [("size", reg, 0.6, f, s, sc) for reg in DOMAIN_SOURCES for sc in SCALES if sc != S
             for f in (1.0, up, dn) for s in SEEDS]
    grid += [("kidney", n, c, a, s) for n in NUCLEATION for c in (True, False)
             for a in (False, True) for s in SEEDS]
    todo = [g for g in grid if _needs(g)]
    if todo:
        n = max(1, min(int(os.environ.get("OGROW_JOBS", 0)) or (os.cpu_count() or 1),
                       len(todo)))
        print(f"running {len(todo)} simulations on {n} worker(s). Every finished run "
              f"is cached,\nso an interrupted pass resumes where it stopped.", flush=True)
        with ProcessPoolExecutor(n) as ex:
            for k, _ in enumerate(ex.map(_warm, todo, chunksize=4), 1):
                if k % max(1, len(todo) // 10) == 0 or k == len(todo):
                    print(f"  {k:5d}/{len(todo)}", flush=True)

    def N(reg, q, f, sc=S):
        return np.array([run("size", reg, q, f, s, sc)["N"] if sc != S
                         else run("size", reg, q, f, s)["N"] for s in SEEDS], float)

    res = dict(seeds=list(SEEDS), thresholds=list(THRESHOLDS), H=H, saturating=SATURATING,
               points=[], capacity={})
    print(f"\n1. The perfusion threshold changed inside each source, {len(SEEDS)} paired "
          f"seeds.\n   Capacity is the units reached at that threshold by a source whose "
          f"pool cannot run out.\n")
    print(f"{'source':>9} {'q_set':>6} {'units':>7} {'capacity':>9} {'units/capacity':>15} "
          f"{'elasticity to the pool':>24}")
    for q in THRESHOLDS:
        res["capacity"][f"{q:g}"] = float(N("saturating", q, 1.0).mean())

    # ---- is it a capacity? the pool doubled and halved, and the horizon doubled
    chk = {f"{pf:g}": float(N("saturating", 0.6, pf).mean()) for pf in POOL_CHECK}
    late = float(np.mean([run("size", "saturating", 0.6, 1.0, s)["late"] for s in SEEDS]))
    hor = {}
    for sc in (SCALES[0], SCALES[-1]):
        short = np.array([run("size", "saturating", 0.6, 1.0, s, sc)["N"] if sc != S
                          else run("size", "saturating", 0.6, 1.0, s)["N"] for s in SEEDS], float)
        long_ = np.array([run("size", "saturating", 0.6, 1.0, s, sc, LONG_DAYS)["N"] for s in SEEDS], float)
        hor[f"{sc:g}"] = dict(short=float(short.mean()), long=float(long_.mean()),
                              ratio=float(long_.mean() / short.mean()))
    res["checks"] = dict(pool=chk, pools=[SATURATING["P0"] * pf for pf in POOL_CHECK],
                         late_gain=late, horizon=hor, long_days=LONG_DAYS)
    print("   is it a capacity? units of the saturating source at the reference threshold, "
          "with its pool\n   at " + ", ".join(f"{SATURATING['P0'] * pf:,.0f}" for pf in POOL_CHECK)
          + ": " + ", ".join(f"{chk[f'{pf:g}']:.0f}" for pf in POOL_CHECK)
          + f". Gain over its last 40 days: {late:.1%}.")
    print(f"   at {LONG_DAYS} days instead of {DAYS}: "
          + "; ".join(f"scale {k_}, {v_['short']:.0f} to {v_['long']:.0f} "
                      f"({v_['ratio'] - 1:+.0%})" for k_, v_ in hor.items()) + ".\n")
    for reg in REG:
        for q in THRESHOLDS:
            n0, nu_, nd_ = N(reg, q, 1.0), N(reg, q, up), N(reg, q, dn)
            cap = N("saturating", q, 1.0)
            e = (np.log(np.maximum(nu_, 1)) - np.log(np.maximum(nd_, 1))) / (2 * H)
            pt = dict(source=reg, q_set=q, units=float(n0.mean()), capacity=float(cap.mean()),
                      fill=float(n0.mean() / cap.mean()), eps=float(e.mean()),
                      eps_se=float(e.std(ddof=1) / np.sqrt(len(e))))
            res["points"].append(pt)
            print(f"{reg:>9} {q:>6.1f} {pt['units']:>7.1f} {pt['capacity']:>9.1f} "
                  f"{pt['fill']:>15.2f} {pt['eps']:>+16.2f} ± {pt['eps_se']:.2f}")
        print()
    P = res["points"]
    eps = [p["eps"] for p in P]
    rf, pf = spearmanr(eps, [p["fill"] for p in P])
    rn, pn = spearmanr(eps, [p["units"] for p in P])
    res["across"] = dict(n=len(P), rho_fill=float(rf), p_fill=float(pf),
                         rho_units=float(rn), p_units=float(pn))
    print(f"   across all {len(P)} points: elasticity against units/capacity, Spearman "
          f"{rf:+.2f} (P = {pf:.2g});\n   against units alone, {rn:+.2f} (P = {pn:.2g}).")

    # inside one source: raise the threshold from its lowest to its highest
    # value. Units fall. If size were what mattered the elasticity would rise.
    res["within"] = {}
    lo_q, hi_q = THRESHOLDS[0], THRESHOLDS[-1]
    print(f"\n   inside one source, threshold raised from {lo_q} to {hi_q}:")
    for reg in REG:
        def change(ix, reg=reg):
            def el(q):
                a, b = N(reg, q, up)[ix], N(reg, q, dn)[ix]
                return float(np.mean((np.log(np.maximum(a, 1)) - np.log(np.maximum(b, 1))) / (2 * H)))
            return el(hi_q) - el(lo_q)
        d = boot(change, len(SEEDS), B=2000)
        lo_pt = next(p for p in P if p["source"] == reg and p["q_set"] == lo_q)
        hi_pt = next(p for p in P if p["source"] == reg and p["q_set"] == hi_q)
        res["within"][reg] = dict(change=d, units=[lo_pt["units"], hi_pt["units"]],
                                  fill=[lo_pt["fill"], hi_pt["fill"]],
                                  eps=[lo_pt["eps"], hi_pt["eps"]])
        print(f"   {reg + '-like':>14}: units {lo_pt['units']:.0f} to {hi_pt['units']:.0f}, "
              f"units/capacity {lo_pt['fill']:.2f} to {hi_pt['fill']:.2f}, elasticity "
              f"{lo_pt['eps']:+.2f} to {hi_pt['eps']:+.2f}, change {sci(d)}")
    print("   A source made smaller by a higher threshold is nearer its capacity. Size "
          "alone\n   predicts a rise of the elasticity there, nearness to capacity a fall.")

    # ---- the domain changed, the pool and the threshold kept
    res["domain"] = {}
    print("\n   the same source in a smaller and a larger domain (threshold 0.6, same pool):\n")
    print(f"{'source':>9} {'scale':>6} {'units':>7} {'capacity':>9} {'units/capacity':>15} "
          f"{'elasticity to the pool':>24}")
    for reg in DOMAIN_SOURCES:
        res["domain"][reg] = {}
        for sc in SCALES:
            n0, cap = N(reg, 0.6, 1.0, sc), N("saturating", 0.6, 1.0, sc)
            e = (np.log(np.maximum(N(reg, 0.6, up, sc), 1)) - np.log(np.maximum(N(reg, 0.6, dn, sc), 1))) / (2 * H)
            pt = dict(source=reg, scale=sc, units=float(n0.mean()), capacity=float(cap.mean()),
                      fill=float(n0.mean() / cap.mean()), eps=float(e.mean()),
                      eps_se=float(e.std(ddof=1) / np.sqrt(len(e))))
            res["domain"][reg][f"{sc:g}"] = pt
            print(f"{reg:>9} {sc:>6.1f} {pt['units']:>7.1f} {pt['capacity']:>9.1f} "
                  f"{pt['fill']:>15.2f} {pt['eps']:>+16.2f} ± {pt['eps_se']:.2f}")

        def dchange(ix, reg=reg):
            def el(sc):
                a_, b_ = N(reg, 0.6, up, sc)[ix], N(reg, 0.6, dn, sc)[ix]
                return float(np.mean((np.log(np.maximum(a_, 1)) - np.log(np.maximum(b_, 1))) / (2 * H)))
            return el(SCALES[-1]) - el(SCALES[0])
        res["domain"][reg]["change"] = boot(dchange, len(SEEDS), B=2000)
        print(f"{'':>9} elasticity, largest domain minus smallest: "
              f"{sci(res['domain'][reg]['change'])}\n")
    allp = P + [res["domain"][r][f"{sc:g}"] for r in DOMAIN_SOURCES for sc in SCALES if sc != S]
    ra, pa = spearmanr([p_["eps"] for p_ in allp], [p_["fill"] for p_ in allp])
    rb, pb = spearmanr([p_["eps"] for p_ in allp], [p_["units"] for p_ in allp])
    res["all"] = dict(n=len(allp), rho_fill=float(ra), p_fill=float(pa),
                      rho_units=float(rb), p_units=float(pb))
    print(f"   over all {len(allp)} points, thresholds and domains: elasticity against "
          f"units/capacity,\n   Spearman {ra:+.2f} (P = {pa:.2g}); against units alone, "
          f"{rb:+.2f} (P = {pb:.2g}).")

    # ---- the kidney with and without the coupling
    print(f"\n2. The kidney-like dividing source with {ABLATE_FRACTION:.0%} of its "
          f"progenitors killed on day {ABLATE_DAY},\n   with branching coupled to the cap "
          f"and with the tree free. Ablated over control.\n")
    print(f"{'nucleation':>32} {'branching':>10} {'pool':>18} {'tips of the tree':>18} "
          f"{'sites with a cap':>18} {'units':>18}")
    res["kidney"] = {}
    for name in NUCLEATION:
        res["kidney"][name] = {}
        for coupled in (True, False):
            c = [run("kidney", name, coupled, False, s) for s in SEEDS]
            a = [run("kidney", name, coupled, True, s) for s in SEEDS]
            row = {k: boot(lambda ix, k=k: float(np.mean([a[i][k] for i in ix])
                                                 / max(np.mean([c[i][k] for i in ix]), 1e-9)),
                           len(SEEDS), B=2000) for k in ("P", "tips", "caps", "N")}
            res["kidney"][name]["coupled" if coupled else "free"] = row
            print(f"{name:>32} {'coupled' if coupled else 'free':>10} {ci(row['P']):>18} "
                  f"{ci(row['tips']):>18} {ci(row['caps']):>18} {ci(row['N']):>18}")
        # the same seeds grow both trees, so the difference is paired
        R4 = {(cpl, abl): [run("kidney", name, cpl, abl, s) for s in SEEDS]
              for cpl in (True, False) for abl in (False, True)}

        def delta(ix, k):
            r = lambda cpl: (np.mean([R4[(cpl, True)][i][k] for i in ix])        # noqa: E731
                             / max(np.mean([R4[(cpl, False)][i][k] for i in ix]), 1e-9))
            return float(r(True) - r(False))
        res["kidney"][name]["difference"] = {k: boot(lambda ix, k=k: delta(ix, k), len(SEEDS), B=2000)
                                             for k in ("tips", "N")}
        dd = res["kidney"][name]["difference"]
        print(f"{'':>32} coupled minus free: tips {sci(dd['tips'])}, units {sci(dd['N'])}\n")
    print("   Measured without a proliferative response (Cebrian et al. 2014): progenitors "
          "stay reduced,\n   tips 0.59, glomeruli 0.48 to 0.58.")

    low = "threshold quartered, law slowed"
    k, w, ac, dm = res["kidney"], res["within"], res["all"], res["domain"]
    lo_s, hi_s = f"{SCALES[0]:g}", f"{SCALES[-1]:g}"
    res["values"] = {
        "c.n": f"{ac['n']}",
        "c.rho.fill": f"{ac['rho_fill']:+.2f}",
        "c.rho.units": f"{ac['rho_units']:+.2f}",
        "c.cap.lo": f"{res['capacity'][f'{THRESHOLDS[-1]:g}']:.0f}",
        "c.cap.hi": f"{res['capacity'][f'{THRESHOLDS[0]:g}']:.0f}",
        **{f"c.{r}.units": f"{w[r]['units'][0]:.0f} to {w[r]['units'][1]:.0f}" for r in REG},
        **{f"c.{r}.fill": f"{w[r]['fill'][0]:.2f} to {w[r]['fill'][1]:.2f}" for r in REG},
        **{f"c.{r}.eps": f"{w[r]['eps'][0]:+.2f} to {w[r]['eps'][1]:+.2f}" for r in REG},
        **{f"c.{r}.change": sci(w[r]["change"]) for r in REG},
        **{f"c.dom.{r}.units": f"{dm[r][lo_s]['units']:.0f} and {dm[r][hi_s]['units']:.0f}" for r in DOMAIN_SOURCES},
        **{f"c.dom.{r}.fill": f"{dm[r][lo_s]['fill']:.2f} and {dm[r][hi_s]['fill']:.2f}" for r in DOMAIN_SOURCES},
        **{f"c.dom.{r}.eps": f"{dm[r][lo_s]['eps']:+.2f} and {dm[r][hi_s]['eps']:+.2f}" for r in DOMAIN_SOURCES},
        **{f"c.dom.{r}.change": sci(dm[r]["change"]) for r in DOMAIN_SOURCES},
        "c.kid.units.delta": sci(k[low]["difference"]["N"]),
        "c.kid.tips.delta": sci(k[low]["difference"]["tips"]),
        "c.sat.pools": ", ".join(f"{x:,.0f}" for x in res["checks"]["pools"]),
        "c.sat.units": ", ".join(f"{res['checks']['pool'][f'{pf:g}']:.0f}" for pf in POOL_CHECK),
        "c.hor.small": f"{res['checks']['horizon'][lo_s]['ratio'] - 1:+.0%}".replace("%", " %"),
        "c.hor.large": f"{res['checks']['horizon'][hi_s]['ratio'] - 1:+.0%}".replace("%", " %"),
        "c.kid.tips.coupled": ci(k[low]["coupled"]["tips"]),
        "c.kid.sites.coupled": ci(k[low]["coupled"]["caps"]),
        "c.kid.tips.free": ci(k[low]["free"]["tips"]),
        "c.kid.units.coupled": ci(k[low]["coupled"]["N"]),
        "c.kid.units.free": ci(k[low]["free"]["N"]),
        "c.kid.pool.free": ci(k[low]["free"]["P"]),
    }
    json.dump(res, open(ROOT / "results" / "capacity.json", "w"), indent=1)
    if write(DOC, res["values"]):
        print(f"\nupdated the values of {DOC.relative_to(ROOT)} that this script writes."
              "\nCommit that file before `make paper`, or the stamp will record an "
              "uncommitted change.")
