# Numbers the paper uses

Model version 4.4.0, commit unknown, generated 2026-11-20 by `scripts/20_paper_numbers.py`. Do not edit by hand.

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
| liver | 242 | 0.05 | 0.56 | 0.39 | 0.00 |
| kidney | 82 | 1.00 | 0.00 | 0.00 | 0.00 |
| pancreas | 41 | 0.68 | 0.32 | 0.00 | 0.00 |

## Transmission of an individual spread (95 % bootstrap interval)

| regime | pool | setpoint |
|---|---|---|
| liver | 0.24 [0.00, 0.34] | 1.60 [1.05, 2.16] |
| kidney | 0.37 [0.16, 0.51] | 1.32 [0.76, 1.89] |
| pancreas | 0.94 [0.65, 1.13] | 1.23 [0.56, 1.86] |

## Elasticities of the final count (mean ± s.e., paired replicates)

| regime | setting | pool | setpoint |
|---|---|---|---|
| liver | ref | +0.33 ± 0.10 | -1.93 ± 0.09 |
| liver | nu_lo | +0.26 ± 0.11 | -2.04 ± 0.09 |
| liver | nu_hi | +0.22 ± 0.08 | -1.94 ± 0.15 |
| liver | q_lo | +0.25 ± 0.09 | -0.98 ± 0.06 |
| liver | q_hi | +0.03 ± 0.16 | -2.80 ± 0.05 |
| kidney | ref | +0.36 ± 0.07 | -1.46 ± 0.12 |
| kidney | nu_lo | +0.39 ± 0.06 | -1.59 ± 0.04 |
| kidney | nu_hi | +0.34 ± 0.09 | -1.51 ± 0.15 |
| kidney | q_lo | +0.44 ± 0.06 | -0.53 ± 0.07 |
| kidney | q_hi | +0.36 ± 0.11 | -2.40 ± 0.09 |
| pancreas | ref | +0.64 ± 0.04 | -0.79 ± 0.11 |
| pancreas | nu_lo | +0.67 ± 0.10 | -0.92 ± 0.17 |
| pancreas | nu_hi | +0.62 ± 0.10 | -0.82 ± 0.23 |
| pancreas | q_lo | +0.71 ± 0.08 | -0.32 ± 0.06 |
| pancreas | q_hi | +0.53 ± 0.07 | -2.58 ± 0.10 |

## Murray ratio in grown trees

| regime | bifurcations | median | interquartile |
|---|---|---|---|
| liver | 187 | 0.896 | 0.526 to 0.934 |
| kidney | 55 | 0.925 | 0.903 to 0.941 |

## Contraction of the spread of unit number

| spread | CV at start | CV at closure |
|---|---|---|
| maturity, offset s.d. 2 d | 0.399 | 0.072 |
| maturity, offset s.d. 4 d | 0.428 | 0.072 |
| maturity, offset s.d. 6 d | 0.447 | 0.072 |
| maturity, offset s.d. 8 d | 0.445 | 0.072 |
| maturity, offset s.d. 10 d | 0.457 | 0.072 |
| growth rate | 0.183 | 0.145 |
| mouse liver (Azizoglu et al. 2025) | 0.229 | 0.074 |

## Robustness to the auxiliary rules (six seeds per cell)

| variant | kidney pool spent | pancreas pool spent | liver held by supply / room | liver pool elasticity | kidney | pancreas | z, pancreas minus liver |
|---|---|---|---|---|---|---|---|
| reference | 1.00 | 0.81 | 0.54 / 0.37 | +0.33 ± 0.10 | +0.36 ± 0.07 | +0.64 ± 0.04 | 2.9 |
| root_surface | 1.00 | 0.76 | 0.42 / 0.45 | +0.36 ± 0.12 | +0.35 ± 0.04 | +0.61 ± 0.14 | 1.4 |
| root_random | 1.00 | 0.64 | 0.46 / 0.47 | -0.68 ± 0.53 | +0.45 ± 0.08 | +0.82 ± 0.07 | 2.8 |
| no_deficit | 1.00 | 0.76 | 0.53 / 0.35 | +0.24 ± 0.09 | +0.31 ± 0.07 | +0.66 ± 0.05 | 4.3 |
| no_load | 1.00 | 0.61 | 0.60 / 0.24 | -0.26 ± 0.68 | +0.30 ± 0.28 | +0.45 ± 0.46 | 0.9 |
| no_crowding | 1.00 | 0.77 | 0.58 / 0.35 | +0.28 ± 0.04 | +0.33 ± 0.07 | +0.60 ± 0.14 | 2.3 |
| no_mesenchyme | 1.00 | 0.91 | 0.43 / 0.45 | +0.19 ± 0.08 | +0.41 ± 0.06 | +0.90 ± 0.03 | 8.5 |
| territory | 1.00 | 0.77 | 0.62 / 0.33 | +0.19 ± 0.05 | +0.35 ± 0.08 | +0.60 ± 0.04 | 6.7 |

