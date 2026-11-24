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
    dirty = subprocess.run(["git", "status", "--porcelain", "src", "scripts", "config"],
                           cwd=ROOT, capture_output=True, text=True).stdout.strip()
    commit = (commit or "unknown") + (" (with uncommitted changes)" if dirty else "")
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
    w("`pool spent` is the share of the initial pool no longer held, which is "
      "the two columns after it added together. A cell consumed in forming a "
      "unit built something; a cell lost to attrition did not. In a renewing "
      "pool the loss runs past 1, because the compartment is replaced rather "
      "than drawn down: the liver's figure is how many times over it turns. "
      "The last column is read off the single trace that the time-series "
      "panels show, not off the replicates.\n")
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
    w("A regime that closes with nothing held has no pool to release, so its "
      "share is left out rather than divided by a pool of zero; the last "
      "column says how many seeds carried one.\n")

    def cell(v, k, se):
        return ("-" if v[k] != v[k] else f"{v[k]:.3f} ± {v[se]:.3f}")

    table(["regime", "intervention", "pool on day 140", "share of that pool released",
           "units added, relative to day 140", "seeds with a pool"],
          [[k.split("_")[0], k.split("_")[1], f"{v['pool']:.0f}",
            cell(v, "released", "se"), cell(v, "units_added", "units_added_se"),
            f"{v.get('seeds_used', 0)} of {v.get('seeds', 0)}"] for k, v in I.items()])

if (RES / "scaling.json").exists():
    Sc = json.load(open(RES / "scaling.json"))
    nse = max(len({r["seed"] for r in v["rows"]}) for v in Sc.values())
    w(f"## Organ size: unit count and closure against domain volume "
      f"({nse} seeds each)\n")
    rows = []
    for mode, v in Sc.items():
        for S in sorted({r["S"] for r in v["rows"]}):
            g = [r for r in v["rows"] if r["S"] == S]
            mm = lambda k: float(np.mean([x[k] for x in g]))
            rows.append([mode, S, g[0]["days"], f"{mm('N'):.0f}", f"{mm('supply_held'):.2f}",
                         f"{mm('space_held'):.2f}", f"{100 * mm('fill'):.1f} %"])
        rows.append([f"{mode}: unit count scales as V^{v['exponent']:.2f}", "", "", "", "", "", ""])
    table(["sweep", "domain scale", "days", "units", "held by supply", "held by room",
           "share of close packing"], rows)

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
