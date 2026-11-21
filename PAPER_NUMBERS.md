# Numbers the paper uses

Model version 4.9.0, commit bfded3f, generated 2026-11-21 by `scripts/20_paper_numbers.py`. Do not edit by hand.

Python 3.9.6, numpy 1.26.4, scipy 1.13.1; `config/params.yaml` sha256 4d4c28e4b024.

## Human kidneys (TRUSTED, independent of the model)

| invariant | between-donor share | paired Wilcoxon P |
|---|---|---|
| volume_cm3 | 0.887 | 0.00067 |
| surface_cm2 | 0.898 | 0.02 |
| sphericity | 0.196 | 0.023 |
| elong | 0.423 | 0.029 |
| flat | 0.242 | 5.4e-06 |
| surf_vol | 0.760 | 4.2e-05 |

## Regimes at closure (one organ each)

| regime | units | pool spent | held by supply | held by room | active | units added in the last fifth |
|---|---|---|---|---|---|---|
| liver | 231 | 0.00 | 0.91 | 0.24 | 0.00 | 2.6 % |
| kidney | 77 | 1.00 | 0.00 | 0.00 | 0.00 | 0.0 % |
| pancreas | 41 | 0.68 | 0.32 | 0.00 | 0.00 | 0.0 % |

## Transmission of an individual spread (95 % bootstrap interval)

| regime | pool | setpoint |
|---|---|---|
| liver | 0.18 [0.00, 0.34] | 1.87 [1.15, 2.48] |
| kidney | 0.38 [0.20, 0.49] | 1.31 [0.75, 1.89] |
| pancreas | 0.90 [0.61, 1.11] | 1.24 [0.62, 1.84] |

## Elasticities of the final count (mean ± s.e., paired replicates)

| regime | setting | pool | setpoint |
|---|---|---|---|
| liver | ref | +0.03 ± 0.17 | -1.45 ± 0.57 |
| liver | nu_lo | -0.09 ± 0.16 | -1.83 ± 0.40 |
| liver | nu_hi | +0.15 ± 0.08 | -1.63 ± 0.60 |
| liver | q_lo | +0.11 ± 0.07 | -1.18 ± 0.09 |
| liver | q_hi | +0.23 ± 0.20 | -2.76 ± 0.07 |
| kidney | ref | +0.36 ± 0.07 | -1.49 ± 0.12 |
| kidney | nu_lo | +0.38 ± 0.11 | -1.45 ± 0.09 |
| kidney | nu_hi | +0.34 ± 0.07 | -1.52 ± 0.12 |
| kidney | q_lo | +0.47 ± 0.09 | -0.56 ± 0.06 |
| kidney | q_hi | +0.42 ± 0.08 | -2.37 ± 0.09 |
| pancreas | ref | +0.65 ± 0.04 | -0.81 ± 0.12 |
| pancreas | nu_lo | +0.65 ± 0.10 | -1.06 ± 0.19 |
| pancreas | nu_hi | +0.64 ± 0.09 | -0.80 ± 0.22 |
| pancreas | q_lo | +0.66 ± 0.07 | -0.31 ± 0.05 |
| pancreas | q_hi | +0.64 ± 0.09 | -2.59 ± 0.08 |

## Murray ratio in grown trees

| regime | bifurcations | median | interquartile |
|---|---|---|---|
| liver | 200 | 0.875 | 0.496 to 0.932 |
| kidney | 51 | 0.919 | 0.883 to 0.940 |

## Contraction of the spread of unit number

| spread | CV at start | CV at closure |
|---|---|---|
| maturity, offset s.d. 2 d | 0.371 | 0.077 |
| maturity, offset s.d. 4 d | 0.399 | 0.077 |
| maturity, offset s.d. 6 d | 0.415 | 0.077 |
| maturity, offset s.d. 8 d | 0.416 | 0.077 |
| maturity, offset s.d. 10 d | 0.429 | 0.077 |
| growth rate | 0.175 | 0.269 |
| mouse liver (Azizoglu et al. 2025) | 0.229 | 0.074 |

## Robustness to the auxiliary rules (six seeds per cell)

| variant | kidney pool spent | pancreas pool spent | liver held by supply / room | liver pool elasticity | kidney | pancreas | z, pancreas minus liver |
|---|---|---|---|---|---|---|---|
| reference | 1.00 | 0.80 | 0.74 / 0.18 | +0.03 ± 0.17 | +0.36 ± 0.07 | +0.65 ± 0.04 | 3.6 |
| root_surface | 1.00 | 0.77 | 0.72 / 0.32 | +0.07 ± 0.07 | +0.28 ± 0.10 | +0.64 ± 0.12 | 4.2 |
| root_random | 1.00 | 0.56 | 0.64 / 0.31 | -0.06 ± 0.16 | +0.15 ± 0.13 | +0.64 ± 0.24 | 2.5 |
| no_deficit | 1.00 | 0.76 | 0.79 / 0.23 | +0.08 ± 0.13 | +0.36 ± 0.07 | +0.64 ± 0.06 | 3.9 |
| no_load | 1.00 | 0.58 | 1.02 / 0.03 | -0.27 ± 1.12 | +0.43 ± 0.24 | +0.23 ± 0.32 | 0.4 |
| no_crowding | 1.00 | 0.77 | 0.74 / 0.30 | +0.12 ± 0.09 | +0.36 ± 0.06 | +0.61 ± 0.11 | 3.3 |
| no_mesenchyme | 1.00 | 0.85 | 0.74 / 0.23 | +0.19 ± 0.07 | +0.37 ± 0.15 | +0.84 ± 0.06 | 7.1 |
| territory | 1.00 | 0.78 | 0.73 / 0.27 | +0.09 ± 0.09 | +0.32 ± 0.08 | +0.65 ± 0.05 | 5.2 |
| binomial | 1.00 | 0.66 | 0.69 / 0.32 | +0.15 ± 0.07 | +0.40 ± 0.04 | +0.80 ± 0.08 | 6.0 |
| binomial_half_day | 1.00 | 0.80 | 0.67 / 0.26 | +0.26 ± 0.07 | +0.26 ± 0.06 | +0.75 ± 0.11 | 3.7 |
| binomial_quarter_day | 1.00 | 0.78 | 0.76 / 0.21 | +0.13 ± 0.07 | +0.40 ± 0.09 | +0.89 ± 0.07 | 7.9 |
| competing_fates | 1.00 | 0.80 | 0.75 / 0.20 | +0.25 ± 0.05 | +0.39 ± 0.05 | +0.65 ± 0.04 | 6.4 |
| no_self_renewal | 1.00 | 0.80 | 0.71 / 0.18 | +0.28 ± 0.08 | +0.36 ± 0.07 | +0.65 ± 0.04 | 4.2 |

## Interventions on the closed state (day 140, then 40 days; six seeds)

| regime | intervention | pool on day 140 | share of that pool released | units added, relative to day 140 |
|---|---|---|---|---|
| liver | none | 2398 | 0.002 ± 0.001 | 0.030 ± 0.009 |
| liver | supply | 2398 | 0.013 ± 0.002 | 0.154 ± 0.026 |
| liver | room | 2398 | 0.027 ± 0.002 | 0.323 ± 0.029 |
| kidney | none | 0 | 0.000 ± 0.000 | 0.000 ± 0.000 |
| kidney | supply | 0 | 0.000 ± 0.000 | 0.000 ± 0.000 |
| kidney | room | 0 | 0.000 ± 0.000 | 0.000 ± 0.000 |
| pancreas | none | 12 | 0.019 ± 0.019 | 0.003 ± 0.003 |
| pancreas | supply | 12 | 0.123 ± 0.063 | 0.032 ± 0.022 |
| pancreas | room | 12 | 0.037 ± 0.037 | 0.007 ± 0.007 |

## Organ size: unit count and closure against domain volume (two seeds each)

| sweep | domain scale | days | units | held by supply | held by room | share of close packing |
|---|---|---|---|---|---|---|
| geometry | 0.4 | 140 | 116 | 0.97 | 0.02 | 14.5 % |
| geometry | 0.6 | 140 | 206 | 0.80 | 0.20 | 7.6 % |
| geometry | 0.8 | 140 | 246 | 0.78 | 0.31 | 3.8 % |
| geometry | 1.0 | 140 | 442 | 0.70 | 0.42 | 3.5 % |
| geometry: unit count scales as V^0.44 |  |  |  |  |  |  |
| allometric | 0.4 | 103 | 98 | 0.96 | 0.14 | 12.3 % |
| allometric | 0.6 | 140 | 206 | 0.80 | 0.20 | 7.6 % |
| allometric | 0.8 | 174 | 458 | 0.82 | 0.14 | 7.1 % |
| allometric | 1.0 | 205 | 430 | 0.60 | 0.43 | 3.4 % |
| allometric: unit count scales as V^0.57 |  |  |  |  |  |  |

