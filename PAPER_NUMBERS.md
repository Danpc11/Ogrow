# Numbers the paper uses

Model version 4.5.0, commit unknown, generated 2026-11-20 by `scripts/20_paper_numbers.py`. Do not edit by hand.

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

| regime | units | pool spent | held by supply | held by room | active |
|---|---|---|---|---|---|
| liver | 202 | 0.02 | 0.79 | 0.19 | 0.00 |
| kidney | 79 | 1.00 | 0.00 | 0.00 | 0.00 |
| pancreas | 41 | 0.68 | 0.32 | 0.00 | 0.00 |

## Transmission of an individual spread (95 % bootstrap interval)

| regime | pool | setpoint |
|---|---|---|
| liver | 0.12 [0.00, 0.28] | 1.61 [1.08, 2.14] |
| kidney | 0.39 [0.23, 0.50] | 1.29 [0.73, 1.85] |
| pancreas | 0.91 [0.63, 1.12] | 1.25 [0.63, 1.85] |

## Elasticities of the final count (mean ± s.e., paired replicates)

| regime | setting | pool | setpoint |
|---|---|---|---|
| liver | ref | +0.22 ± 0.05 | -1.89 ± 0.31 |
| liver | nu_lo | +0.14 ± 0.08 | -2.01 ± 0.14 |
| liver | nu_hi | +0.22 ± 0.03 | -1.82 ± 0.33 |
| liver | q_lo | +0.26 ± 0.08 | -0.88 ± 0.16 |
| liver | q_hi | +0.13 ± 0.20 | -2.79 ± 0.04 |
| kidney | ref | +0.35 ± 0.06 | -1.51 ± 0.12 |
| kidney | nu_lo | +0.36 ± 0.09 | -1.57 ± 0.06 |
| kidney | nu_hi | +0.37 ± 0.04 | -1.51 ± 0.12 |
| kidney | q_lo | +0.40 ± 0.07 | -0.50 ± 0.07 |
| kidney | q_hi | +0.38 ± 0.11 | -2.39 ± 0.08 |
| pancreas | ref | +0.61 ± 0.04 | -0.87 ± 0.11 |
| pancreas | nu_lo | +0.63 ± 0.11 | -1.00 ± 0.19 |
| pancreas | nu_hi | +0.61 ± 0.09 | -0.83 ± 0.22 |
| pancreas | q_lo | +0.70 ± 0.07 | -0.27 ± 0.06 |
| pancreas | q_hi | +0.63 ± 0.11 | -2.56 ± 0.10 |

## Murray ratio in grown trees

| regime | bifurcations | median | interquartile |
|---|---|---|---|
| liver | 148 | 0.900 | 0.808 to 0.936 |
| kidney | 53 | 0.923 | 0.898 to 0.942 |

## Contraction of the spread of unit number

| spread | CV at start | CV at closure |
|---|---|---|
| maturity, offset s.d. 2 d | 0.398 | 0.050 |
| maturity, offset s.d. 4 d | 0.429 | 0.050 |
| maturity, offset s.d. 6 d | 0.449 | 0.050 |
| maturity, offset s.d. 8 d | 0.448 | 0.050 |
| maturity, offset s.d. 10 d | 0.461 | 0.050 |
| growth rate | 0.183 | 0.167 |
| mouse liver (Azizoglu et al. 2025) | 0.229 | 0.074 |

## Robustness to the auxiliary rules (six seeds per cell)

| variant | kidney pool spent | pancreas pool spent | liver held by supply / room | liver pool elasticity | kidney | pancreas | z, pancreas minus liver |
|---|---|---|---|---|---|---|---|
| reference | 1.00 | 0.79 | 0.75 / 0.15 | +0.22 ± 0.05 | +0.35 ± 0.06 | +0.61 ± 0.04 | 6.0 |
| root_surface | 1.00 | 0.77 | 0.65 / 0.25 | +0.19 ± 0.07 | +0.32 ± 0.05 | +0.60 ± 0.12 | 2.9 |
| root_random | 1.00 | 0.55 | 0.63 / 0.32 | +0.28 ± 0.13 | +0.09 ± 0.10 | +0.64 ± 0.23 | 1.4 |
| no_deficit | 1.00 | 0.76 | 0.74 / 0.17 | +0.28 ± 0.05 | +0.28 ± 0.07 | +0.65 ± 0.06 | 5.1 |
| no_load | 1.00 | 0.58 | 0.61 / 0.25 | -0.55 ± 0.55 | +0.43 ± 0.24 | +0.23 ± 0.32 | 1.2 |
| no_crowding | 1.00 | 0.76 | 0.73 / 0.21 | +0.33 ± 0.05 | +0.37 ± 0.07 | +0.59 ± 0.10 | 2.3 |
| no_mesenchyme | 1.00 | 0.88 | 0.73 / 0.18 | +0.11 ± 0.06 | +0.40 ± 0.08 | +0.81 ± 0.06 | 8.3 |
| territory | 1.00 | 0.79 | 0.66 / 0.23 | +0.06 ± 0.08 | +0.35 ± 0.08 | +0.67 ± 0.05 | 6.6 |

