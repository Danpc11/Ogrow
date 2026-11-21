# Numbers the paper uses

Model version 4.7.0, commit 2536cd6, generated 2026-11-21 by `scripts/20_paper_numbers.py`. Do not edit by hand.

Python 3.9.6, numpy 1.26.4, scipy 1.13.1; `config/params.yaml` sha256 04415cee25e4.

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
| liver | 214 | 0.16 | 0.70 | 0.14 | 0.00 | 6.5 % |
| kidney | 77 | 1.00 | 0.00 | 0.00 | 0.00 | 0.0 % |
| pancreas | 41 | 0.68 | 0.32 | 0.00 | 0.00 | 0.0 % |

## Transmission of an individual spread (95 % bootstrap interval)

| regime | pool | setpoint |
|---|---|---|
| liver | 0.14 [0.00, 0.30] | 1.56 [1.01, 2.10] |
| kidney | 0.38 [0.20, 0.49] | 1.31 [0.75, 1.89] |
| pancreas | 0.90 [0.61, 1.11] | 1.24 [0.62, 1.84] |

## Elasticities of the final count (mean ± s.e., paired replicates)

| regime | setting | pool | setpoint |
|---|---|---|---|
| liver | ref | +0.28 ± 0.08 | -2.08 ± 0.14 |
| liver | nu_lo | +0.17 ± 0.03 | -2.12 ± 0.13 |
| liver | nu_hi | +0.25 ± 0.06 | -2.18 ± 0.09 |
| liver | q_lo | +0.24 ± 0.05 | -0.94 ± 0.06 |
| liver | q_hi | +0.24 ± 0.15 | -2.73 ± 0.05 |
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
| liver | 165 | 0.890 | 0.523 to 0.938 |
| kidney | 51 | 0.919 | 0.883 to 0.940 |

## Contraction of the spread of unit number

| spread | CV at start | CV at closure |
|---|---|---|
| maturity, offset s.d. 2 d | 0.397 | 0.087 |
| maturity, offset s.d. 4 d | 0.425 | 0.087 |
| maturity, offset s.d. 6 d | 0.444 | 0.087 |
| maturity, offset s.d. 8 d | 0.442 | 0.087 |
| maturity, offset s.d. 10 d | 0.454 | 0.087 |
| growth rate | 0.155 | 0.147 |
| mouse liver (Azizoglu et al. 2025) | 0.229 | 0.074 |

## Robustness to the auxiliary rules (six seeds per cell)

| variant | kidney pool spent | pancreas pool spent | liver held by supply / room | liver pool elasticity | kidney | pancreas | z, pancreas minus liver |
|---|---|---|---|---|---|---|---|
| reference | 1.00 | 0.80 | 0.71 / 0.18 | +0.28 ± 0.08 | +0.36 ± 0.07 | +0.65 ± 0.04 | 4.2 |
| root_surface | 1.00 | 0.77 | 0.66 / 0.29 | +0.23 ± 0.09 | +0.28 ± 0.10 | +0.64 ± 0.12 | 2.8 |
| root_random | 1.00 | 0.56 | 0.52 / 0.35 | +0.01 ± 0.16 | +0.15 ± 0.13 | +0.64 ± 0.24 | 2.2 |
| no_deficit | 1.00 | 0.76 | 0.68 / 0.27 | +0.22 ± 0.08 | +0.36 ± 0.07 | +0.64 ± 0.06 | 4.2 |
| no_load | 1.00 | 0.58 | 0.61 / 0.25 | -0.55 ± 0.55 | +0.43 ± 0.24 | +0.23 ± 0.32 | 1.2 |
| no_crowding | 1.00 | 0.77 | 0.55 / 0.34 | +0.25 ± 0.07 | +0.36 ± 0.06 | +0.61 ± 0.11 | 2.7 |
| no_mesenchyme | 1.00 | 0.85 | 0.71 / 0.30 | +0.15 ± 0.05 | +0.37 ± 0.15 | +0.84 ± 0.06 | 8.7 |
| territory | 1.00 | 0.78 | 0.64 / 0.24 | +0.14 ± 0.09 | +0.32 ± 0.08 | +0.65 ± 0.05 | 4.9 |
| binomial | 1.00 | 0.66 | 0.63 / 0.24 | +0.14 ± 0.05 | +0.40 ± 0.04 | +0.80 ± 0.08 | 6.8 |
| binomial_half_day | 1.00 | 0.80 | 0.57 / 0.31 | +0.14 ± 0.05 | +0.26 ± 0.06 | +0.75 ± 0.11 | 5.0 |
| binomial_quarter_day | 1.00 | 0.78 | 0.69 / 0.23 | +0.07 ± 0.08 | +0.40 ± 0.09 | +0.89 ± 0.07 | 7.7 |
| competing_fates | 1.00 | 0.80 | 0.75 / 0.16 | +0.19 ± 0.04 | +0.39 ± 0.05 | +0.65 ± 0.04 | 8.8 |

## Interventions on the closed state (day 140, then 40 days; six seeds)

| regime | intervention | pool on day 140 | share of that pool released | units added, relative to day 140 |
|---|---|---|---|---|
| liver | none | 2316 | 0.003 ± 0.000 | 0.037 ± 0.006 |
| liver | supply | 2316 | 0.012 ± 0.003 | 0.130 ± 0.024 |
| liver | room | 2316 | 0.026 ± 0.002 | 0.278 ± 0.023 |
| kidney | none | 0 | 0.000 ± 0.000 | 0.000 ± 0.000 |
| kidney | supply | 0 | 0.000 ± 0.000 | 0.000 ± 0.000 |
| kidney | room | 0 | 0.000 ± 0.000 | 0.000 ± 0.000 |
| pancreas | none | 12 | 0.019 ± 0.019 | 0.003 ± 0.003 |
| pancreas | supply | 12 | 0.123 ± 0.063 | 0.032 ± 0.022 |
| pancreas | room | 12 | 0.037 ± 0.037 | 0.007 ± 0.007 |

