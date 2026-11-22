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
    w("## Regimes at closure (one organ each)\n")
    rows = []
    for r in reg:
        d = np.load(RES / "paper" / f"regime_{r}.npz")
        c = d["closure"] if "closure" in d.files else [np.nan] * 4
        Nt = d["N"]; late = (Nt[-1] - Nt[int(0.8 * len(Nt))]) / max(Nt[-1], 1)
        rows.append([r, int(Nt[-1]), f"{c[0]:.2f}", f"{c[1]:.2f}", f"{c[2]:.2f}", f"{c[3]:.2f}",
                     f"{100 * late:.1f} %"])
    table(["regime", "units", "pool spent", "held by supply", "held by room", "active",
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
             f"{v['q_set']['mean']:+.2f} ± {v['q_set']['se']:.2f}" if "q_set" in v else "-"]
            for k, v in E.items()]
    table(["regime", "setting", "pool", "setpoint"], rows)

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
    w("## Robustness to the auxiliary rules (six seeds per cell)\n")
    rows = []
    for v, x in Rb.items():
        L, K_, Pn_ = x["liver"], x["kidney"], x["pancreas"]
        d = Pn_["eps_pool"] - L["eps_pool"]; s = float(np.hypot(Pn_["eps_se"], L["eps_se"]))
        rows.append([v, f"{K_['spent']:.2f}", f"{Pn_['spent']:.2f}",
                     f"{L['supply_held']:.2f} / {L['space_held']:.2f}",
                     f"{L['eps_pool']:+.2f} ± {L['eps_se']:.2f}", f"{K_['eps_pool']:+.2f} ± {K_['eps_se']:.2f}",
                     f"{Pn_['eps_pool']:+.2f} ± {Pn_['eps_se']:.2f}", f"{d / s:.1f}"])
    table(["variant", "kidney pool spent", "pancreas pool spent", "liver held by supply / room",
           "liver pool elasticity", "kidney", "pancreas", "z, pancreas minus liver"], rows)

if (RES / "interventions.json").exists():
    I = json.load(open(RES / "interventions.json"))
    w("## Interventions on the closed state (day 140, then 40 days; six seeds)\n")
    table(["regime", "intervention", "pool on day 140", "share of that pool released",
           "units added, relative to day 140"],
          [[k.split("_")[0], k.split("_")[1], f"{v['pool']:.0f}", f"{v['released']:.3f} ± {v['se']:.3f}",
            f"{v['units_added']:.3f} ± {v['units_added_se']:.3f}"] for k, v in I.items()])

if (RES / "scaling.json").exists():
    Sc = json.load(open(RES / "scaling.json"))
    w("## Organ size: unit count and closure against domain volume (two seeds each)\n")
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
    w("## Injury after closure (day 140, 30 % of the organ, read on day 220; three seeds)\n")
    names = {"N": "units", "mass": "mass", "perfused": "perfused mass", "Q": "root flow"}
    rows = []
    for k, v in Ij.items():
        reg, kind = k.rsplit("_", 1)
        if kind == "control":
            continue
        for q, nm in names.items():
            c = v[q]
            rec = "not lowered" if c["recovered"] != c["recovered"] else \
                  f"{c['recovered']:.2f} ± {c['se']:.2f}"
            rows.append([reg, kind, nm, rec, f"{c['of_control']:.2f}"])
    table(["regime", "injury", "readout", "share of the loss recovered",
           "day 220, relative to uninjured"], rows)

open(ROOT / "PAPER_NUMBERS.md", "w").write("\n".join(out) + "\n")
print("wrote PAPER_NUMBERS.md")
