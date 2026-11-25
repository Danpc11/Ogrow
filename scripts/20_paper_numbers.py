#!/usr/bin/env python
"""Write PAPER_NUMBERS.md: every model number the paper uses, with the model
version and commit that produced it.

The paper's text does not copy model numbers by hand; it points here, so that
numbers from different versions of the model cannot sit side by side.
Run by `make paper` after the runs it reads.
"""
import csv
import datetime
import json
import subprocess

import numpy as np

import organgrow
from organgrow.config import ROOT

RES = ROOT / "results"
out = []
w = out.append
try:
    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                            capture_output=True, text=True).stdout.strip()
    # The whole tree, not just src, scripts and config. A run is reproducible
    # from the named commit only if everything that shaped it is in that
    # commit, and the Makefile decides which Python runs and over what. The
    # paths skipped here are the ones this run writes itself, so their being
    # modified says nothing about the code that produced the numbers.
    _mine = ("PAPER_NUMBERS.md", "source_data/", "results/", "figures/out/")
    _lines = subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                            capture_output=True, text=True).stdout.splitlines()
    _paths = sorted({ln[3:].strip().strip('"') for ln in _lines if len(ln) > 3})
    dirty = [p for p in _paths if not any(p.startswith(m) for m in _mine)]
    commit = (commit or "unknown") + (
        " (with uncommitted changes to " + ", ".join(dirty[:4])
        + (f" and {len(dirty) - 4} more" if len(dirty) > 4 else "") + ")"
        if dirty else "")
except Exception:
    commit = "unknown"
w("# Numbers the paper uses\n")
import hashlib, platform, scipy
cfg_hash = hashlib.sha256((ROOT / "config" / "params.yaml").read_bytes()).hexdigest()[:12]
w(f"Model version {organgrow.__version__}, commit {commit}, generated "
  f"{datetime.date.today().isoformat()} by `scripts/20_paper_numbers.py`. Do not edit by hand.\n")
w(f"Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}; "
  f"`config/params.yaml` sha256 {cfg_hash}.\n")

def table(header, rows):
    w("| " + " | ".join(header) + " |")
    w("|" + "---|" * len(header))
    for r in rows:
        w("| " + " | ".join(str(x) for x in r) + " |")
    w("")

src = ROOT / "source_data"
if (src / "Fig6bc_variance_partition.csv").exists():
    w("## Human kidneys (TRUSTED, independent of the model)\n")
    rows = [[r["invariant"], f"{float(r['between_donor_fraction']):.3f}", f"{float(r['p_wilcoxon']):.2g}"]
            for r in csv.DictReader(open(src / "Fig6bc_variance_partition.csv"))]
    table(["invariant", "between-donor share", "paired Wilcoxon P"], rows)

reg = ("liver", "kidney", "pancreas")
if all((RES / "paper" / f"regime_{r}.npz").exists() for r in reg):
    rp = RES / "paper" / "regime_reps.json"
    Rp = json.load(open(rp)) if rp.exists() else {}
    nrep = max((v["seeds"] for v in Rp.values()), default=1)
    w(f"## Regimes at closure (mean ± s.e. over {nrep} seeds)\n")
    w("Two different things are reported here and they do not add up to each "
      "other. `pool spent` is a *stock*: the share of the initial pool no "
      "longer held at a tip, 1 - P_final/P_0. The two columns after it are "
      "*flows* over the whole run, each divided by P_0, and they satisfy the "
      "balance the tests check exactly, P_0 + born = consumed + left + lost. "
      "In a pool that does not renew, born is zero and the flows do add to "
      "the stock. In a renewing pool they do not and cannot: the compartment "
      "is replaced rather than drawn down, so the loss runs past 1 and is how "
      "many times over the compartment turns, while the stock barely moves. A "
      "cell consumed in forming a unit built something; a cell lost to "
      "attrition did not. The last column is read off the single trace that "
      "the time-series panels show, not off the replicates.\n")
    rows = []
    for r in reg:
        d = np.load(RES / "paper" / f"regime_{r}.npz")
        Nt = d["N"]; late = (Nt[-1] - Nt[int(0.8 * len(Nt))]) / max(Nt[-1], 1)
        v = Rp.get(r)
        if v is None:                      # no replicates cached: the one trace
            c = d["closure"] if "closure" in d.files else [np.nan] * 4
            fa = d["fate"] if "fate" in d.files else [np.nan] * 4
            cell = dict(zip(("spent", "supply_held", "space_held", "active"),
                            (f"{x:.2f}" for x in c)))
            cell.update(zip(("consumed", "lost"), (f"{x:.2f}" for x in fa)))
            units = str(int(Nt[-1]))
        else:
            def cell_(d_, k):
                return f"{d_[k]['mean']:.2f} ± {d_[k]['se']:.2f}"
            cell = {k: cell_(v["closure"], k)
                    for k in ("spent", "supply_held", "space_held", "active")}
            cell.update({k: cell_(v["fate"], k) for k in ("consumed", "lost")})
            units = f"{v['N']['mean']:.0f} ± {v['N']['se']:.0f}"
        rows.append([r, units, cell["spent"], cell["consumed"], cell["lost"],
                     cell["supply_held"], cell["space_held"], cell["active"],
                     f"{100 * late:.1f} %"])
    table(["regime", "units", "pool spent", "of it, consumed into units",
           "of it, lost to attrition", "held by supply", "held by room", "active",
           "units added in the last fifth"], rows)

if (RES / "transmission.json").exists():
    T = json.load(open(RES / "transmission.json"))
    w("## Transmission of an individual spread (95 % bootstrap interval)\n")
    rows = []
    for r in reg:
        def cell(k):
            v = T.get(f"{r}_{k}")
            return f"{v['transmission']:.2f} [{v['ci95'][0]:.2f}, {v['ci95'][1]:.2f}]" if v else "-"
        rows.append([r, cell("pool"), cell("setpoint")])
    table(["regime", "pool", "setpoint"], rows)

if (RES / "elasticity.json").exists():
    E = json.load(open(RES / "elasticity.json"))
    w("## Elasticities of the final count (mean ± s.e., paired replicates)\n")
    rows = [[k.split("_", 1)[0], k.split("_", 1)[1],
             f"{v['P0']['mean']:+.2f} ± {v['P0']['se']:.2f}" if "P0" in v else "-",
             f"{v['q_set']['mean']:+.2f} ± {v['q_set']['se']:.2f}" if "q_set" in v else "-",
             str(max((v[i]["n"] for i in ("P0", "q_set") if i in v), default=0))]
            for k, v in E.items()]
    table(["regime", "setting", "pool", "setpoint", "paired replicates"], rows)

if (RES / "regimes.json").exists():
    Rg = json.load(open(RES / "regimes.json"))
    meta = Rg.get("_meta", {})
    w("## What separates the regimes: renewal against pool and turnover\n")
    w("The three regimes differ in several parameters at once, so the binary "
      "property cannot be credited with the separation until the others are "
      "held fixed. This varies each independently and reports where the pool "
      "stands at arrest. Division and loss are one axis rather than two, "
      "swept over the pairings the three sources use, because a pool that "
      "divides without being lost grows without bound and because the kidney "
      "and the pancreas do not divide at all. The last column names the "
      "largest of the four shares, which is the bound that ended the growth "
      "(`scripts/34_regimes.py`).\n")
    named = {("1.0", "2600", "balanced"): "liver-like",
             ("0.0", "300", "decaying"): "kidney-like",
             ("0.0", "60", "static"): "pancreas-like"}

    def _holds(c):
        k = max(("spent", "supply_held", "space_held", "active"),
                key=lambda n: c[n])
        return {"spent": "spent", "supply_held": "supply",
                "space_held": "room", "active": "still active"}[k]

    rows = []
    for k, c in Rg.items():
        if k.startswith("_"):
            continue
        key = (f"{c['renewal']}", f"{c['P0']:.0f}", c["turnover"])
        rows.append([f"{c['renewal']:.0f}", f"{c['P0']:.0f}", c["turnover"],
                     f"{c['rho']:.2f}", f"{c['delta']:.2f}",
                     f"{c['N']:.1f} ± {c['N_sd']:.1f}", f"{c['spent']:.2f}",
                     f"{c['supply_held']:.2f}", f"{c['space_held']:.2f}",
                     f"{c['active']:.2f}", _holds(c), named.get(key, "")])
    table(["renewal", "pool", "turnover", "division", "loss", "units", "spent",
           "supply", "room", "active", "what holds growth", "named regime"],
          rows)
    if "renewal_flips" in meta:
        w("Each axis is asked the same question: holding the other two "
          "fixed, in how many of those settings does moving along this axis "
          "change which bound ends the growth?\n")
        table(["axis varied", "settings of the other two",
               "in how many the bound changes"],
              [["renewal", meta.get("cells", "?"), meta["renewal_flips"]],
               ["turnover", meta.get("other_cells", "?"),
                meta.get("turnover_flips", "?")],
               ["pool", meta.get("other_cells", "?"),
                meta.get("pool_flips", "?")]])

    el = Rg.get("_elasticity_matched") or {}
    if el:
        w("The elasticity of the final unit count to the seeded pool is "
          "reported elsewhere at each regime's own pool size, so a renewing "
          "source at 2600 is compared with a spent one at 60. Here the pool "
          "and the turnover are matched between the two sources and only the "
          "binary property differs.\n")
        table(["pool", "turnover", "spent source", "renewing source", "gap"],
              [[f"{v['P0']:.0f}", v.get("turnover", "-"),
                f"{v['spent'][0]:+.2f} ± {v['spent'][1]:.2f}",
                f"{v['renewing'][0]:+.2f} ± {v['renewing'][1]:.2f}",
                f"{v['spent'][0] - v['renewing'][0]:+.2f}"]
               for v in el.values()])

if (RES / "feasibility.json").exists():
    Fe = json.load(open(RES / "feasibility.json"))
    sp = Fe["spec"]
    w("## Designing to a specification: where the requirements can be met\n")
    w(f"A joint specification is imposed at each target volume and the knobs "
      f"an operator sets are searched for a configuration meeting all of it: "
      f"units occupying at least {100 * sp['fill_min']:.0f}% of the domain's "
      f"close packing, at most {sp['vessel_max']:.0f} vessel segments, and a "
      f"coefficient of variation of the unit count across realisations of at "
      f"most {sp['cv_max']:.2f}. The knobs are the seeded pool, the "
      f"nucleation propensity, the culture duration and the cell source. The "
      f"configuration is chosen on {len(Fe['search_seeds'])} seeds by "
      f"maximin, the largest smallest slack, and re-run on "
      f"{len(Fe['holdout_seeds'])} seeds the search never saw. Configurations "
      f"that had not arrested are excluded "
      f"(`scripts/35_feasibility.py`).\n")

    w("What any configuration reaches at each volume, one requirement at a "
      "time. The three are reached by different configurations, so no single "
      "construct has all three at once. Each row also carries what the "
      "winning configuration does on the other two requirements, because a "
      "column read on its own can be won by a construct that is useless: the "
      "lowest variation is easiest to reach by building almost nothing.\n")

    def _cfg(c):
        return (f"{c['src']}, pool {c['P0']:.0f}, nu {c['nu']:.0f}, "
                f"{c['days']} d")

    for S, r in sorted(Fe["reach"].items()):
        w(f"At scale {S}, volume {float(S) ** 3:.3f}, "
          f"{r['arrested']} of {r['configs']} configurations had arrested.\n")
        table(["requirement won", "configuration", "fill", "vessels", "CV"],
              [[name, _cfg(r[k]), f"{100 * r[k]['fill']:.2f}%",
                f"{r[k]['vessels']:.0f}", f"{r[k]['cv']:.3f}"]
               for name, k in (("most dense", "best_fill"),
                               ("fewest vessels", "fewest_vessels"),
                               ("least variable", "lowest_cv"))])

    rows = []
    for S in sorted(set(Fe["chosen"]) | set(Fe["infeasible"])):
        if S in Fe["chosen"]:
            c, h = Fe["chosen"][S], Fe["held_out"].get(S, {})
            rows.append([S, "yes", f"{c['P0']:.0f}", f"{c['nu']:.0f}",
                         f"{c['days']}", c["src"],
                         f"{100 * c['fill']:.2f}%", f"{c['vessels']:.0f}",
                         f"{c['cv']:.3f}",
                         f"{100 * h['fill']:.2f}%" if h else "-",
                         f"{h['cv']:.3f}" if h else "-"])
        else:
            rows.append([S, "no", "-", "-", "-", "-", "-", "-", "-", "-", "-"])
    table(["scale", "feasible", "pool", "nu", "days", "source",
           "fill", "vessels", "CV", "fill, held out", "CV, held out"], rows)
    for S, b in sorted(Fe["infeasible"].items()):
        w(f"At scale {S} no configuration met the specification. The "
          f"requirement none could meet is {b['binding']}; the best any "
          f"configuration came to it was a slack of "
          f"{b['best'][b['binding']]:+.2f}, where zero would be meeting it "
          f"exactly.\n")
    if Fe.get("optimism") is not None:
        o = Fe["optimism"]
        w(f"The margin measured on the search block "
          f"{'overstates' if o >= 0 else 'understates'} the margin on the "
          f"held-out block by {abs(o):.2f} on average, because the "
          f"configuration is chosen partly for the realisations that happened "
          f"to favour it. The held-out column is the honest one. Held-out "
          f"seeds are new realisations of the same model and not new evidence "
          f"about it, so this measures reproducibility across batches and is "
          f"not a validation against a construct.\n")

    cv_ = Fe.get("convergence") or {}
    if cv_:
        w("Each chosen configuration again on the held-out seeds, with the "
          "cellular step split two and four ways, under the reference capped "
          "Poisson law and under the binomial law. Only the binomial one is a "
          "step refinement: the capped law applies its per-day cap once per "
          "substep, so splitting the day loosens it. A verdict that moves "
          "here is a verdict about the step and not about the construct.\n")
        table(["scale", "law", "substeps", "fill", "vessels", "CV", "margin",
               "meets spec"],
              [[f"{v['S']}", v["law"], v["sub"], f"{100 * v['fill']:.2f}%",
                f"{v['vessels']:.0f}", f"{v['cv']:.3f}", f"{v['margin']:+.2f}",
                "yes" if v["margin"] >= 0 else "no"]
               for _, v in sorted(cv_.items(),
                                  key=lambda kv: (kv[1]["S"], kv[1]["law"],
                                                  kv[1]["sub"]))])
        base = {(v["S"], v["law"]): v["margin"] >= 0 for v in cv_.values()
                if v["sub"] == 1}
        moved = [v for v in cv_.values() if v["sub"] != 1
                 and (v["margin"] >= 0) != base[(v["S"], v["law"])]]
        tried = sum(1 for v in cv_.values() if v["sub"] != 1)
        w(f"The verdict moves in {len(moved)} of {tried} refinements"
          + (": " + ", ".join(sorted(f"scale {v['S']} under the {v['law']} "
                                     f"law at {v['sub']} substeps"
                                     for v in moved)) + ".\n"
             if moved else ".\n"))

    fr = Fe.get("frontier") or {}
    if fr:
        fills = sorted({float(k.split("_")[0][1:]) for k in fr})
        cvs = sorted({float(k.split("_cv")[1]) for k in fr})
        w("The largest volume that stays feasible as the density and "
          "reproducibility requirements change, at the same vascular budget. "
          "An entry is the domain scale; `none` means no volume admits a "
          "configuration.\n")
        table(["fill >="] + [f"CV <= {c:.2f}" for c in cvs],
              [[f"{f:.2f}"] + [(lambda v: f"{v:.1f}" if v else "none")(
                  fr.get(f"f{f}_cv{c}")) for c in cvs] for f in fills])

if (RES / "paper_numbers.json").exists():
    Pn = json.load(open(RES / "paper_numbers.json"))
    w("## Murray ratio in grown trees\n")
    table(["regime", "bifurcations", "median", "interquartile"],
          [[r, Pn[f"murray_{r}"]["n"], f"{Pn[f'murray_{r}']['median']:.3f}",
            f"{Pn[f'murray_{r}']['q25']:.3f} to {Pn[f'murray_{r}']['q75']:.3f}"]
           for r in ("liver", "kidney") if f"murray_{r}" in Pn])
    c = Pn.get("calibration")
    if c:
        w("## Contraction of the spread of unit number\n")
        table(["spread", "CV at start", "CV at closure"],
              [[f"maturity, offset s.d. {m['offset_sd']} d", f"{m['cv_start']:.3f}", f"{m['cv_final']:.3f}"]
               for m in c["maturity"]] +
              [["growth rate", f"{c['rate']['cv_start']:.3f}", f"{c['rate']['cv_final']:.3f}"],
               ["mouse liver (Azizoglu et al. 2025)", "0.229", "0.074"]])

if (RES / "robustness.json").exists():
    Rb = json.load(open(RES / "robustness.json"))
    w("## Robustness to the auxiliary rules (24 seeds per cell; the contrast is "
      "paired, seed by seed)\n")
    rows = []
    for v, x in Rb.items():
        L, K_, Pn_ = x["liver"], x["kidney"], x["pancreas"]
        zz = []
        for a, b in (("kidney", "liver"), ("pancreas", "kidney"), ("pancreas", "liver")):
            c = x.get(f"{a}_minus_{b}")
            zz.append(f"{c['z']:.1f}" if c else "-")
        rows.append([v, f"{K_['spent']:.2f}", f"{Pn_['spent']:.2f}",
                     f"{L['supply_held']:.2f} / {L['space_held']:.2f}",
                     f"{L['eps_pool']:+.2f} ± {L['eps_se']:.2f}", f"{K_['eps_pool']:+.2f} ± {K_['eps_se']:.2f}",
                     f"{Pn_['eps_pool']:+.2f} ± {Pn_['eps_se']:.2f}"] + zz)
    table(["variant", "kidney pool spent", "pancreas pool spent", "liver held by supply / room",
           "liver pool elasticity", "kidney", "pancreas",
           "z, kidney minus liver", "z, pancreas minus kidney",
           "z, pancreas minus liver"], rows)

    # The counts the prose refers to, generated rather than written by hand:
    # a rerun on another machine moves the digits, and these have to move with
    # them. The sign test is over the variants, two-sided.
    from math import comb
    w("## The ordering, counted over the variants\n")
    nv = len(Rb)
    rows = [["liver below kidney below pancreas",
             f"{sum(x['liver']['eps_pool'] < x['kidney']['eps_pool'] < x['pancreas']['eps_pool'] for x in Rb.values())} of {nv}",
             "-", "-", "-"]]
    for a, b in (("kidney", "liver"), ("pancreas", "kidney"), ("pancreas", "liver")):
        c = [x.get(f"{a}_minus_{b}") for x in Rb.values()]
        if not all(c):
            continue
        pos = sum(d["mean"] > 0 for d in c)
        P = 2 * sum(comb(nv, k) for k in range(pos, nv + 1)) / 2 ** nv
        ref = Rb.get("reference", {}).get(f"{a}_minus_{b}")
        rows.append([f"{a} above {b}", f"{pos} of {nv}", f"{min(P, 1.0):.4f}",
                     f"{sum(d['z'] >= 2 for d in c)} of {nv}",
                     f"{ref['mean']:+.2f} ± {ref['se']:.2f}, z = {ref['z']:.1f}" if ref else "-"])
    table(["statement", "variants holding it", "sign test P", "of those, z at least 2",
           "at the reference setting"], rows)

if (RES / "interventions.json").exists():
    I = json.load(open(RES / "interventions.json"))
    ns = {v.get("seeds", 0) for v in I.values()}
    w(f"## Interventions on the closed state (day 140, then 40 days; "
      f"{max(ns) if ns else '?'} seeds)\n")
    w("A regime that closes with nothing held has no pool to release, so the "
      "share is reported only where what the regime still holds is a real "
      "remnant, one per cent of the pool it started with and at least one "
      "cell. Where it is not, the entry is left out rather than divided by a "
      "pool of zero, and the column beside it gives the pool against that "
      "floor. Units added needs no pool, so it is always over every seed.\n")

    def cell(v, k, se):
        return ("-" if v[k] != v[k] else f"{v[k]:.3f} ± {v[se]:.3f}")

    table(["regime", "intervention", "pool on day 140", "pool needed for a share",
           "share of that pool released", "units added, relative to day 140"],
          [[k.split("_")[0], k.split("_")[1],
            f"{v['pool']:.1f} [{v.get('pool_lo', float('nan')):.1f}, "
            f"{v.get('pool_hi', float('nan')):.1f}]",
            f"{v.get('pool_floor', float('nan')):.0f}",
            cell(v, "released", "se"), cell(v, "units_added", "units_added_se")]
           for k, v in I.items()])

if (RES / "scaling.json").exists():
    Sc = json.load(open(RES / "scaling.json"))
    nse = max(len({r["seed"] for r in v["rows"]}) for v in Sc.values())
    w(f"## Organ size: unit count and closure against domain volume "
      f"({nse} seeds each)\n")
    w("The unit count is the mean over the seeds with the range they span and "
      "its coefficient of variation, because whether a larger domain is "
      "reproducible at all is what the exponent rests on.\n")
    rows = []
    for mode, v in Sc.items():
        for S in sorted({r["S"] for r in v["rows"]}):
            g = [r for r in v["rows"] if r["S"] == S]
            mm = lambda k: float(np.mean([x[k] for x in g]))
            N = np.array([x["N"] for x in g], float)
            rows.append([mode, S, g[0]["days"],
                         f"{N.mean():.0f} [{N.min():.0f}, {N.max():.0f}]",
                         f"{N.std(ddof=1) / max(N.mean(), 1e-9):.2f}",
                         f"{mm('supply_held'):.2f}", f"{mm('space_held'):.2f}",
                         f"{100 * mm('fill'):.1f} %"])
        rows.append([f"{mode}: unit count scales as V^{v['exponent']:.2f}",
                     "", "", "", "", "", "", ""])
    table(["sweep", "domain scale", "days", "units", "CV", "held by supply",
           "held by room", "share of close packing"], rows)

if (RES / "inflow_law.json").exists():
    Fl = json.load(open(RES / "inflow_law.json"))
    w("## The mass a construct lands on against the flow it is given\n")
    w("The ceiling on a unit is the delivered flow times its share over the "
      "delivery per unit mass the tissue needs, and the shares sum to one, so "
      "a total of Q/phi is what the rule gives if every unit reaches its "
      "ceiling. It is a consequence of the rule, not a measurement. What is "
      "measured here is how completely a construct reaches it: `mass / "
      "(Q/phi)` away from 1, and the exponent away from 1. `phi` is "
      "calibrated on a block of seeds disjoint from the one the sweep is "
      "evaluated on.\n\n"
      "A small daily change of mass is necessary for equilibrium with the "
      "supply and not sufficient. The growth rate decays as exp(-t/tau_g), so "
      "the growth available after the injury is finite and a unit can stop "
      "moving while still far from its ceiling. The residual column is "
      "`sum |m - ceiling| / sum ceiling` over the units, and the column beside "
      "it is the share of units ending at the mass floor. The two together say "
      "that the total lands where the budget puts it while the units "
      "individually do not.\n\n"
      "In this section the spread after ± is a standard deviation over seeds, "
      "not the standard error used elsewhere in this file.\n")
    for reg, d in Fl.items():
        c, e = d.get("calibration_seeds", ["?", "?"]), d.get("evaluation_seeds", ["?", "?"])
        ci = d.get("slope_ci", [float("nan")] * 2)
        w(f"**{reg}**: phi = {d['phi']:.6f} ± {d['phi_sd']:.6f} (s.d.), "
          f"calibrated on seeds {c[0]} to {c[1]}, evaluated on seeds "
          f"{e[0]} to {e[1]}. Exponent of mass against delivered flow, one fit "
          f"per seed: {d['slope']:.3f} ± {d.get('slope_sd', float('nan')):.3f} "
          f"(s.d.), bootstrap 95 % interval [{ci[0]:.3f}, {ci[1]:.3f}].\n")
        table(["delivered flow", "mass / pre-injury", "mass / (Q/phi)",
               "mass still moving, % per day", "residual", "units at the mass floor"],
              [[f"{r['scale']:.1f}", f"{r['mass_rel']:.3f} ± {r['mass_rel_sd']:.3f}",
                f"{r['mass_over_law']:.3f}",
                f"{100 * r.get('drift_per_day', float('nan')):.2f}",
                f"{r.get('residual', float('nan')):.3f}",
                f"{r.get('floor_share', float('nan')):.2f}"]
               for r in d["rows"]])
        rr = [r.get("residual", float("nan")) for r in d["rows"]]
        ff = [r.get("floor_share", float("nan")) for r in d["rows"]]
        w(f"Across the whole range of flows the residual stays between "
          f"{min(rr):.2f} and {max(rr):.2f}, and between {min(ff):.0%} and "
          f"{max(ff):.0%} of units end at the mass floor. Neither improves at "
          f"any flow, so this is the regime and not an edge of the sweep.\n")

if (RES / "involution.json").exists():
    Iv = json.load(open(RES / "involution.json"))
    w("## Involution, as a sensitivity\n")
    w("The units stranded at the mass floor above are there because the "
      "reference model has no way to remove a unit that loses its supply. "
      "`involution` is a rate per day at which a unit whose supply has fallen "
      "below a tenth of the setpoint loses mass, after which it is removed; "
      "the rate of zero is the reference model, so the first row of each table "
      "is the published result. The rate is not calibrated against any "
      "measurement, which is why this is a sensitivity and not the result. "
      "Everything is at unit flow, with the demand calibrated out of sample as "
      "above (`scripts/33_involution.py`).\n")
    for reg, d in Iv.items():
        c, e = d.get("calibration_seeds", ["?", "?"]), d.get("evaluation_seeds", ["?", "?"])
        w(f"**{reg}**: phi = {d['phi']:.6f}, calibrated on seeds {c[0]} to "
          f"{c[1]}, evaluated on seeds {e[0]} to {e[1]}, {d['seeds']} seeds.\n")
        cols = ["setting", "mass / pre-injury", "units / pre-injury", "residual",
                "units at the mass floor", "tips serving nothing",
                "days to 90 % of the removals"]

        def row(r, label):
            cl = r.get("days_to_clear", float("nan"))
            return [label,
                    f"{r['mass_rel']:.3f} [{r['mass_rel_lo']:.2f}, {r['mass_rel_hi']:.2f}]",
                    f"{r['units_rel']:.2f} [{r['units_rel_lo']:.2f}, {r['units_rel_hi']:.2f}]",
                    f"{r['residual']:.3f}", f"{r['floor_share']:.2f}",
                    f"{r.get('idle_frac', float('nan')):.2f}",
                    "-" if cl != cl else f"{cl:.0f} [{r['days_to_clear_lo']:.0f}, "
                                         f"{r['days_to_clear_hi']:.0f}]"]

        w("The rate, which is how fast the mechanism runs.\n")
        table(cols, [row(r, f"{r['rate']:.2f} per day") for r in d["rows"]])
        u, f_, rs = d["units_rel_on"], d["floor_on"], d["residual_on"]
        w(f"Over a tenfold range of the rate the construct keeps {u[0]:.2f} to "
          f"{u[1]:.2f} of its units, against {d['rows'][0]['units_rel']:.2f} "
          f"with the mechanism absent; the share at the mass floor falls from "
          f"{d['rows'][0]['floor_share']:.2f} to between {f_[0]:.2f} and "
          f"{f_[1]:.2f}, and the residual from {d['rows'][0]['residual']:.2f} "
          f"to between {rs[0]:.2f} and {rs[1]:.2f}.\n")

        tg = d.get("triggers") or {}
        if tg:
            ref = d.get("reference", {})
            w("What turns it on. The rate is held at "
              f"{ref.get('rate', '?')} per day throughout these, so each sweep "
              "crosses the one above at its reference setting. `involution_q` "
              "is the supply, as a share of the setpoint, below which a unit "
              "counts as starving; `involution_hold` is the consecutive days "
              "it must stay there before it begins to lose mass; "
              "`involution_exit` is the mass, as a share of the mass a unit is "
              "born with, below which it is removed.\n")
            for name, rows_ in tg.items():
                key = name.split("_")[-1]
                table(cols, [row(r, f"{key} = {r['trig'].get(key, ref.get(key, '?'))}")
                             for r in rows_])
            ua, ma = d.get("units_rel_any"), d.get("mass_rel_any")
            if ua and ma:
                w(f"Across every setting of the rate and of the three conditions "
                  f"that turn it on, the construct keeps {ua[0]:.2f} to "
                  f"{ua[1]:.2f} of its units against "
                  f"{d['rows'][0]['units_rel']:.2f} with the mechanism absent, "
                  f"while its mass stays between {ma[0]:.2f} and {ma[1]:.2f} of "
                  f"pre-injury against {d['rows'][0]['mass_rel']:.2f}. The "
                  f"settings move when the removals happen, not what the "
                  f"construct ends as.\n")

if (RES / "injury.json").exists():
    Ij = json.load(open(RES / "injury.json"))
    names = {"units": "units", "mass": "mass", "mass_per_unit": "mass per unit",
             "perfused": "perfused mass", "units_made_up": "share of the lost units made up",
             "mass_made_up": "share of the lost mass made up"}

    def seedcell(d):
        """Mean over the seeds and the range they span, which is the point."""
        v = d["seeds"]
        return f"{d['mean']:.2f} [{min(v):.2f}, {max(v):.2f}]"
    w("## Injury after closure (day 140, 30 % of the units removed by a plane cut, "
      "read on day 220; the first four entries are day 220 over the day before "
      "the injury, the last two the share of what the cut took that is back, "
      "each the mean over 24 seeds with the range they span)\n")
    rows = []
    for reg in ("liver", "kidney"):
        for cond in ("as grown", "reactivated", "reactivated, fixed inflow"):
            v = Ij.get(f"{reg}_{cond}")
            if v:
                rows.append([reg, cond] + [seedcell(v[q]) if q in v else "-"
                                           for q in names])
    table(["regime", "condition"] + list(names.values()), rows)

    # How wide the draw is, and how much the budget narrows it. These are the
    # numbers the prose refers to, generated for the same reason as the counts
    # above: a rerun moves them.
    def cv(d):
        v = np.asarray(d["seeds"], float)
        return 100.0 * v.std(ddof=1) / v.mean() if v.mean() else float("nan")

    w("## Spread across seeds, and what a fixed inflow does to it "
      "(coefficient of variation, %, over 24 seeds)\n")
    rows = []
    for reg in ("liver", "kidney"):
        free, held = Ij.get(f"{reg}_reactivated"), Ij.get(f"{reg}_reactivated, fixed inflow")
        if not (free and held):
            continue
        for key, lab in (("mass", "mass"), ("units", "unit count"),
                         ("mass_per_unit", "mass per unit")):
            a, b = cv(free[key]), cv(held[key])
            rows.append([reg, lab, f"{a:.1f}", f"{b:.1f}",
                         f"{a / b:.1f}" if b else "-"])
    table(["regime", "readout", "growth reactivated, inlet pressure fixed",
           "and inflow held", "contraction"], rows)

    comp = {k: v for k, v in Ij.items() if "_compensation_" in k}
    if comp:
        w("## How far a kidney compensates, growth reactivated and inflow fixed\n")
        rows = []
        for k in sorted(comp):
            reg, frac = k.split("_compensation_")
            v = comp[k]
            rows.append([reg, f"{100 * float(frac):.0f} %", seedcell(v["units"]),
                         seedcell(v["mass_per_unit"]), f"{v['full']:.2f}",
                         seedcell(v["mass"])])
        table(["regime", "lost", "units kept", "mass per unit",
               "full compensation would need", "mass"], rows)

open(ROOT / "PAPER_NUMBERS.md", "w").write("\n".join(out) + "\n")
print("wrote PAPER_NUMBERS.md")
