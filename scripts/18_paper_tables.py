#!/usr/bin/env python
"""Source data for the closed-model figures, one CSV per table."""
import csv
import json

import numpy as np

from organgrow.config import ROOT

RES, OUT = ROOT / "results", ROOT / "source_data"
P = RES / "paper"


def write(name, header, rows):
    with open(OUT / name, "w", newline="") as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)
    print(f"  {name}: {len(rows)} rows")


T = json.load(open(RES / "transmission.json"))
write("Fig3closed_transmission.csv", ["regime", "input_spread", "replicate", "final_units"],
      [[k.split("_")[0], k.split("_")[1], i, n] for k, v in T.items() for i, n in enumerate(v["N"])])
rows = []
for reg in ("liver", "kidney", "pancreas"):
    d = np.load(P / f"regime_{reg}.npz")
    for t in range(len(d["N"])):
        rows.append([reg, t + 1, d["N"][t], d["q"][t], d["P"][t], d["vessels"][t]])
write("Fig2closed_regimes.csv", ["regime", "day", "units", "supply_per_unit", "pool", "vessels"], rows)
rows = []
for grp in ("cohort", "rate"):
    for k in range(8):
        d = np.load(P / f"{grp}_k{k}.npz")
        rows += [[grp, k, float(d["nu"]), t + 1, d["N"][t]] for t in range(len(d["N"]))]
write("Fig4closed_trajectories.csv", ["group", "individual", "nu", "day", "units"], rows)
rows = []
for reg in ("liver", "kidney"):
    rows += [[reg, float(v)] for v in np.load(P / f"regime_{reg}.npz")["murray"]]
write("Fig5closed_murray.csv", ["regime", "ratio_rp_x_over_sum_rc_x"], rows)
