#!/usr/bin/env python
"""Export the source data behind every panel, one CSV per table.

These are the numbers plotted, not the raw simulation state: a reader can
redraw any panel from them without running the engine.
"""
import csv
import json
from itertools import combinations
import numpy as np
from organgrow.config import CFG, ROOT

RES = ROOT / "results"
OUT = ROOT / "source_data"
OUT.mkdir(exist_ok=True)


def write(name, header, rows):
    with open(OUT / name, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print(f"  {name:<34} {len(rows):>5} rows")


# ---------------------------------------------------------------- Figure 1
cost = json.load(open(RES / "cost.json"))
write("Fig1e_flow_solve_cost.csv",
      ["edges", "tree_traversal_ms", "sparse_solve_ms"],
      [[c["E"], c["tree"] * 1e3,
        "" if c["sparse"] != c["sparse"] else c["sparse"] * 1e3] for c in cost])

# ---------------------------------------------------------------- Figure 2
E = np.load(RES / "ens_main.npz")
organs = CFG["organs"]["growth_series"]
rows = []
for k in organs:
    N, H = E[f"{k}_N"], E[f"{k}_hypo"]
    for d in range(N.shape[1]):
        rows.append([k, d + 1, N[:, d].mean(), np.percentile(N[:, d], 10),
                     np.percentile(N[:, d], 90), np.median(H[:, d])])
write("Fig2efgh_growth_trajectories.csv",
      ["organ", "day", "units_mean", "units_p10", "units_p90",
       "under_perfused_fraction_median"], rows)

# ---------------------------------------------------------------- Figure 3
M = np.load(RES / "f3_main.npz")
rows = []
for k in M.files:
    A = M[k]
    for d in range(A.shape[1]):
        rows.append([k, d + 1, A[:, d].mean(), A[:, d].std(ddof=1)])
write("Fig3bc_unit_number_dispersion.csv",
      ["organ", "day", "units_mean", "units_sd_across_replicates"], rows)

rows = []
for k in ("liver", "kidney"):
    S = np.load(RES / f"f3_sweep_{k}.npz")
    for i, sp in enumerate(S["sig"]):
        A = S["N"][i]
        rows.append([k, float(sp), A[:, 0].std(ddof=1), A[:, -1].std(ddof=1),
                     A[:, -1].std(ddof=1) / A[:, 0].std(ddof=1),
                     A[:, -1].std(ddof=1) / A[:, -1].mean()])
write("Fig3ef_parameter_disorder_sweep.csv",
      ["closure_mode_organ", "sigma_param", "sd_day1", "sd_day70",
       "sd_ratio", "terminal_cv"], rows)

# ---------------------------------------------------------------- Figure 4
rows = []
for tag in ("nodisorder", "disorder"):
    A = json.load(open(RES / f"scaling_{tag}.json"))
    for r in A:
        rows.append([tag, r["days"], r["N"], r["cv_N"], r["cv_mean_unit"],
                     r["cv_mass"], r["cv_dissipation"],
                     r["vessel_identity"], r["density_discord"]])
write("Fig4bd_scaling_with_organ_size.csv",
      ["parameter_disorder", "days", "units", "cv_unit_number",
       "cv_mean_unit_mass", "cv_total_mass", "cv_dissipation",
       "vessel_identity_disagreement", "density_pattern_discord"], rows)

A0 = json.load(open(RES / "scaling_nodisorder.json"))
N0 = np.array([r["N"] for r in A0])
rows = []
for key, nm in [("cv_N", "unit number"), ("cv_mean_unit", "mean unit mass"),
                ("cv_mass", "total mass"),
                ("vessel_identity", "vessel identity"),
                ("density_discord", "local density pattern")]:
    y = np.array([r[key] for r in A0])
    g, cov = np.polyfit(np.log(N0), np.log(y), 1, cov=True)
    rows.append([nm, -g[0], float(np.sqrt(cov[0, 0]))])
write("Fig4c_fitted_exponents.csv", ["observable", "exponent_a", "std_error"], rows)

REF = np.load(RES / "ref_cohort.npy")
inv = ["unit number", "mean unit mass", "unit-mass spread", "mass per volume",
       "unit spacing", "spacing spread", "supply distance", "segment length"]
C = np.corrcoef(REF.T)
write("Fig4e_invariant_correlation.csv", [""] + inv,
      [[inv[i]] + list(np.round(C[i], 4)) for i in range(len(inv))])

# ---------------------------------------------------------------- Figure 5
I = json.load(open(RES / "ident.json"))
write("Fig5c_condition_numbers.csv",
      ["observable_set", "rank", "n_parameters", "condition_number"],
      [[k, v["rank"], len(I["params"]), v["cond"]] for k, v in I["sets"].items()])
rows = []
for k, v in I["sets"].items():
    for p, se, n in zip(I["params"], v["se"], v["nrep"]):
        rows.append([k, p, se, n])
write("Fig5de_parameter_precision.csv",
      ["observable_set", "parameter", "standard_error_log", "replicates_for_20pct"],
      rows)
S = np.load(RES / "sens.npz")
days, keys, pars = S["days"], [str(x) for x in S["keys"]], [str(x) for x in S["params"]]
cols = [f"{k}_day{d}" for k in keys for d in days]
write("Fig5b_sensitivity_matrix.csv", ["parameter"] + cols,
      [[pars[i]] + list(np.round(S["J"][i], 5)) for i in range(len(pars))])

# ---------------------------------------------------------------- Figure 6
D = json.load(open(RES / "trusted_split.json"))
ks = sorted(D)
write("Fig6_kidney_invariants.csv",
      ["kidney", "donor", "side", "volume_cm3", "surface_cm2", "sphericity",
       "elongation", "flatness", "surface_per_volume"],
      [[k, D[k]["pid"], D[k]["side"], D[k]["volume_cm3"], D[k]["surface_cm2"],
        D[k]["sphericity"], D[k]["elong"], D[k]["flat"], D[k]["surf_vol"]]
       for k in ks])

P = json.load(open(RES / "partition.json"))
write("Fig6bc_variance_partition.csv",
      ["invariant", "between_donor_fraction", "within_donor_fraction",
       "left_minus_right_pct", "n_left_greater", "n_pairs", "p_sign",
       "p_wilcoxon"],
      [[nm, v["between_frac"], 1 - v["between_frac"], v["lr_rel_diff"],
        v["n_left_greater"], v["n_pairs"], v["p_sign"], v["p_wilcoxon"]]
       for nm, v in P["per_invariant"].items()])
write("Fig6a_within_between_distances.csv",
      ["quantity", "value"],
      [["kidneys", P["kidneys"]], ["donors", P["donors"]],
       ["within_donor_mean", P["intra"]], ["between_donor_mean", P["inter"]],
       ["ratio", P["ratio"]], ["permutation_null_mean", P["null_mean"]],
       ["permutation_null_sd", P["null_sd"]], ["p_value", P["p_permutation"]]])

print(f"\nsource data in {OUT.relative_to(ROOT)}/")
