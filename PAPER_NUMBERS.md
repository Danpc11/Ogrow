# Numbers the paper uses

Model version 4.46.0, commit 567bfb7, generated 2026-11-24 by `scripts/20_paper_numbers.py`. Do not edit by hand.

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

## Regimes at closure (mean ± s.e. over 24 seeds)

Two different things are reported here and they do not add up to each other. `pool spent` is a *stock*: the share of the initial pool no longer held at a tip, 1 - P_final/P_0. The two columns after it are *flows* over the whole run, each divided by P_0, and they satisfy the balance the tests check exactly, P_0 + born = consumed + left + lost. In a pool that does not renew, born is zero and the flows do add to the stock. In a renewing pool they do not and cannot: the compartment is replaced rather than drawn down, so the loss runs past 1 and is how many times over the compartment turns, while the stock barely moves. A cell consumed in forming a unit built something; a cell lost to attrition did not. The last column is read off the single trace that the time-series panels show, not off the replicates.

| regime | units | pool spent | of it, consumed into units | of it, lost to attrition | held by supply | held by room | active | units added in the last fifth |
|---|---|---|---|---|---|---|---|---|
| liver | 208 ± 4 | 0.06 ± 0.02 | 0.00 ± 0.00 | 13.66 ± 0.15 | 0.78 ± 0.03 | 0.19 ± 0.02 | 0.00 ± 0.00 | 2.6 % |
| kidney | 89 ± 2 | 1.00 ± 0.00 | 0.30 ± 0.01 | 0.70 ± 0.01 | 0.00 ± 0.00 | 0.00 ± 0.00 | 0.00 ± 0.00 | 0.0 % |
| pancreas | 45 ± 1 | 0.75 ± 0.02 | 0.75 ± 0.02 | 0.00 ± 0.00 | 0.25 ± 0.02 | 0.00 ± 0.00 | 0.00 ± 0.00 | 0.0 % |

## Transmission of an individual spread (95 % bootstrap interval)

| regime | pool | setpoint |
|---|---|---|
| liver | 0.18 [0.00, 0.34] | 1.87 [1.15, 2.48] |
| kidney | 0.38 [0.20, 0.49] | 1.31 [0.75, 1.89] |
| pancreas | 0.90 [0.61, 1.11] | 1.24 [0.62, 1.84] |

## Elasticities of the final count (mean ± s.e., paired replicates)

| regime | setting | pool | setpoint | paired replicates |
|---|---|---|---|---|
| liver | ref | +0.07 ± 0.06 | -1.83 ± 0.22 | 24 |
| liver | nu_lo | +0.08 ± 0.06 | -1.54 ± 0.20 | 24 |
| liver | nu_hi | +0.10 ± 0.05 | -1.91 ± 0.23 | 24 |
| liver | q_lo | +0.09 ± 0.04 | -1.16 ± 0.04 | 24 |
| liver | q_hi | +0.37 ± 0.14 | -2.90 ± 0.03 | 24 |
| kidney | ref | +0.37 ± 0.03 | -1.65 ± 0.06 | 24 |
| kidney | nu_lo | +0.33 ± 0.04 | -1.54 ± 0.06 | 24 |
| kidney | nu_hi | +0.36 ± 0.03 | -1.61 ± 0.06 | 24 |
| kidney | q_lo | +0.40 ± 0.03 | -0.64 ± 0.03 | 24 |
| kidney | q_hi | +0.38 ± 0.04 | -2.44 ± 0.04 | 24 |
| pancreas | ref | +0.71 ± 0.05 | -0.77 ± 0.07 | 24 |
| pancreas | nu_lo | +0.76 ± 0.05 | -1.11 ± 0.09 | 24 |
| pancreas | nu_hi | +0.72 ± 0.04 | -0.66 ± 0.08 | 24 |
| pancreas | q_lo | +0.80 ± 0.04 | -0.35 ± 0.04 | 24 |
| pancreas | q_hi | +0.70 ± 0.06 | -2.52 ± 0.06 | 24 |

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

## Robustness to the auxiliary rules (24 seeds per cell; the contrast is paired, seed by seed)

| variant | kidney pool spent | pancreas pool spent | liver held by supply / room | liver pool elasticity | kidney | pancreas | z, kidney minus liver | z, pancreas minus kidney | z, pancreas minus liver |
|---|---|---|---|---|---|---|---|---|---|
| reference | 1.00 | 0.75 | 0.78 / 0.19 | +0.07 ± 0.06 | +0.37 ± 0.03 | +0.71 ± 0.05 | 4.1 | 5.7 | 7.8 |
| root_surface | 1.00 | 0.73 | 0.77 / 0.27 | +0.10 ± 0.04 | +0.34 ± 0.04 | +0.79 ± 0.06 | 4.8 | 6.0 | 10.0 |
| root_random | 1.00 | 0.66 | 0.64 / 0.33 | +0.05 ± 0.06 | +0.27 ± 0.05 | +0.71 ± 0.08 | 3.0 | 4.9 | 6.4 |
| no_deficit | 1.00 | 0.72 | 0.79 / 0.21 | +0.12 ± 0.05 | +0.37 ± 0.03 | +0.73 ± 0.05 | 4.9 | 6.3 | 9.2 |
| no_load | 1.00 | 0.72 | 0.93 / 0.08 | -0.49 ± 0.45 | +0.26 ± 0.12 | +0.36 ± 0.16 | 1.5 | 0.6 | 1.8 |
| no_crowding | 1.00 | 0.76 | 0.76 / 0.26 | +0.12 ± 0.04 | +0.38 ± 0.03 | +0.77 ± 0.05 | 5.4 | 6.2 | 8.7 |
| no_mesenchyme | 1.00 | 0.86 | 0.70 / 0.30 | +0.14 ± 0.04 | +0.34 ± 0.05 | +0.83 ± 0.03 | 3.3 | 8.0 | 13.3 |
| territory | 1.00 | 0.72 | 0.74 / 0.28 | +0.08 ± 0.06 | +0.27 ± 0.04 | +0.72 ± 0.04 | 2.5 | 8.2 | 8.9 |
| binomial | 1.00 | 0.71 | 0.73 / 0.25 | +0.15 ± 0.04 | +0.42 ± 0.03 | +0.80 ± 0.05 | 5.7 | 7.9 | 9.8 |
| binomial_half_day | 1.00 | 0.78 | 0.75 / 0.27 | +0.11 ± 0.04 | +0.37 ± 0.04 | +0.75 ± 0.05 | 4.8 | 6.0 | 10.2 |
| binomial_quarter_day | 1.00 | 0.80 | 0.78 / 0.22 | +0.14 ± 0.05 | +0.35 ± 0.04 | +0.91 ± 0.04 | 3.4 | 9.0 | 12.7 |
| competing_fates | 1.00 | 0.75 | 0.78 / 0.21 | +0.07 ± 0.04 | +0.36 ± 0.03 | +0.71 ± 0.05 | 5.7 | 5.9 | 10.5 |
| no_self_renewal | 1.00 | 0.75 | 0.70 / 0.19 | +0.24 ± 0.04 | +0.37 ± 0.03 | +0.71 ± 0.05 | 2.9 | 5.7 | 7.6 |

## The ordering, counted over the variants

| statement | variants holding it | sign test P | of those, z at least 2 | at the reference setting |
|---|---|---|---|---|
| liver below kidney below pancreas | 13 of 13 | - | - | - |
| kidney above liver | 13 of 13 | 0.0002 | 12 of 13 | +0.31 ± 0.08, z = 4.1 |
| pancreas above kidney | 13 of 13 | 0.0002 | 12 of 13 | +0.34 ± 0.06, z = 5.7 |
| pancreas above liver | 13 of 13 | 0.0002 | 12 of 13 | +0.65 ± 0.08, z = 7.8 |

## Interventions on the closed state (day 140, then 40 days; 24 seeds)

A regime that closes with nothing held has no pool to release, so the share is reported only where what the regime still holds is a real remnant, one per cent of the pool it started with and at least one cell. Where it is not, the entry is left out rather than divided by a pool of zero, and the column beside it gives the pool against that floor. Units added needs no pool, so it is always over every seed.

| regime | intervention | pool on day 140 | pool needed for a share | share of that pool released | units added, relative to day 140 |
|---|---|---|---|---|---|
| liver | none | 2534.8 [2014.0, 3306.0] | 26 | 0.002 ± 0.000 | 0.028 ± 0.003 |
| liver | supply | 2534.8 [2014.0, 3306.0] | 26 | 0.018 ± 0.002 | 0.226 ± 0.027 |
| liver | room | 2534.8 [2014.0, 3306.0] | 26 | 0.027 ± 0.002 | 0.323 ± 0.017 |
| kidney | none | 0.1 [0.0, 1.0] | 3 | - | 0.000 ± 0.000 |
| kidney | supply | 0.1 [0.0, 1.0] | 3 | - | 0.000 ± 0.000 |
| kidney | room | 0.1 [0.0, 1.0] | 3 | - | 0.000 ± 0.000 |
| pancreas | none | 15.0 [4.0, 23.0] | 1 | 0.008 ± 0.006 | 0.002 ± 0.001 |
| pancreas | supply | 15.0 [4.0, 23.0] | 1 | 0.170 ± 0.030 | 0.062 ± 0.014 |
| pancreas | room | 15.0 [4.0, 23.0] | 1 | 0.015 ± 0.010 | 0.003 ± 0.002 |

## Organ size: unit count and closure against domain volume (24 seeds each)

The unit count is the mean over the seeds with the range they span and its coefficient of variation, because whether a larger domain is reproducible at all is what the exponent rests on.

| sweep | domain scale | days | units | CV | held by supply | held by room | share of close packing |
|---|---|---|---|---|---|---|---|
| geometry | 0.4 | 140 | 130 [108, 158] | 0.09 | 0.89 | 0.09 | 16.2 % |
| geometry | 0.6 | 140 | 208 [174, 241] | 0.09 | 0.78 | 0.19 | 7.7 % |
| geometry | 0.8 | 140 | 289 [219, 525] | 0.23 | 0.71 | 0.30 | 4.5 % |
| geometry | 1.0 | 140 | 520 [274, 1074] | 0.47 | 0.63 | 0.37 | 4.2 % |
| geometry: unit count scales as V^0.45 |  |  |  |  |  |  |  |
| allometric | 0.4 | 103 | 107 [81, 125] | 0.09 | 0.86 | 0.15 | 13.3 % |
| allometric | 0.6 | 140 | 208 [174, 241] | 0.09 | 0.78 | 0.19 | 7.7 % |
| allometric | 0.8 | 174 | 338 [259, 651] | 0.25 | 0.72 | 0.29 | 5.3 % |
| allometric | 1.0 | 205 | 572 [308, 1052] | 0.40 | 0.64 | 0.37 | 4.6 % |
| allometric: unit count scales as V^0.58 |  |  |  |  |  |  |  |

## The mass a construct lands on against the flow it is given

The ceiling on a unit is the delivered flow times its share over the delivery per unit mass the tissue needs, and the shares sum to one, so a total of Q/phi is what the rule gives if every unit reaches its ceiling. It is a consequence of the rule, not a measurement. What is measured here is how completely a construct reaches it: `mass / (Q/phi)` away from 1, and the exponent away from 1. `phi` is calibrated on a block of seeds disjoint from the one the sweep is evaluated on.

A small daily change of mass is necessary for equilibrium with the supply and not sufficient. The growth rate decays as exp(-t/tau_g), so the growth available after the injury is finite and a unit can stop moving while still far from its ceiling. The residual column is `sum |m - ceiling| / sum ceiling` over the units, and the column beside it is the share of units ending at the mass floor. The two together say that the total lands where the budget puts it while the units individually do not.

In this section the spread after ± is a standard deviation over seeds, not the standard error used elsewhere in this file.

**liver**: phi = 0.000752 ± 0.000034 (s.d.), calibrated on seeds 400 to 423, evaluated on seeds 300 to 323. Exponent of mass against delivered flow, one fit per seed: 0.899 ± 0.024 (s.d.), bootstrap 95 % interval [0.891, 0.909].

| delivered flow | mass / pre-injury | mass / (Q/phi) | mass still moving, % per day | residual | units at the mass floor |
|---|---|---|---|---|---|
| 0.6 | 0.645 ± 0.021 | 1.760 | 0.02 | 0.377 | 0.74 |
| 0.8 | 0.833 ± 0.030 | 1.278 | 0.01 | 0.349 | 0.73 |
| 1.0 | 1.020 ± 0.041 | 1.001 | 0.02 | 0.335 | 0.72 |
| 1.2 | 1.203 ± 0.051 | 0.819 | 0.03 | 0.334 | 0.73 |
| 1.4 | 1.382 ± 0.056 | 0.692 | 0.03 | 0.328 | 0.73 |

Across the whole range of flows the residual stays between 0.33 and 0.38, and between 72% and 74% of units end at the mass floor. Neither improves at any flow, so this is the regime and not an edge of the sweep.

**kidney**: phi = 0.000905 ± 0.000047 (s.d.), calibrated on seeds 400 to 423, evaluated on seeds 300 to 323. Exponent of mass against delivered flow, one fit per seed: 0.932 ± 0.023 (s.d.), bootstrap 95 % interval [0.923, 0.941].

| delivered flow | mass / pre-injury | mass / (Q/phi) | mass still moving, % per day | residual | units at the mass floor |
|---|---|---|---|---|---|
| 0.6 | 0.603 ± 0.031 | 1.674 | 0.01 | 0.243 | 0.49 |
| 0.8 | 0.787 ± 0.043 | 1.228 | 0.01 | 0.218 | 0.43 |
| 1.0 | 0.970 ± 0.055 | 0.968 | 0.02 | 0.206 | 0.40 |
| 1.2 | 1.150 ± 0.066 | 0.798 | 0.03 | 0.199 | 0.36 |
| 1.4 | 1.329 ± 0.077 | 0.677 | 0.04 | 0.195 | 0.35 |

Across the whole range of flows the residual stays between 0.19 and 0.24, and between 35% and 49% of units end at the mass floor. Neither improves at any flow, so this is the regime and not an edge of the sweep.

## Involution, as a sensitivity

The units stranded at the mass floor above are there because the reference model has no way to remove a unit that loses its supply. `involution` is a rate per day at which a unit whose supply has fallen below a tenth of the setpoint loses mass, after which it is removed; the rate of zero is the reference model, so the first row of each table is the published result. The rate is not calibrated against any measurement, which is why this is a sensitivity and not the result. Everything is at unit flow, with the demand calibrated out of sample as above (`scripts/33_involution.py`).

**liver**: phi = 0.000752, calibrated on seeds 400 to 423, evaluated on seeds 300 to 323, 24 seeds.

| involution, per day | mass / pre-injury | units / pre-injury | residual | units at the mass floor |
|---|---|---|---|---|
| 0.00 | 1.020 [0.91, 1.08] | 0.91 [0.72, 1.04] | 0.335 | 0.72 |
| 0.02 | 0.954 [0.84, 1.01] | 0.35 [0.29, 0.45] | 0.231 | 0.14 |
| 0.05 | 0.966 [0.89, 1.07] | 0.37 [0.28, 0.50] | 0.233 | 0.13 |
| 0.10 | 0.965 [0.83, 1.03] | 0.40 [0.25, 0.77] | 0.234 | 0.15 |
| 0.20 | 0.977 [0.89, 1.05] | 0.39 [0.28, 0.46] | 0.227 | 0.14 |

Over a tenfold range of the rate the construct keeps 0.35 to 0.40 of its units, against 0.91 with the mechanism absent; the share at the mass floor falls from 0.72 to between 0.13 and 0.15, and the residual from 0.33 to between 0.23 and 0.23. The rate sets how fast, not where it lands.

**kidney**: phi = 0.000905, calibrated on seeds 400 to 423, evaluated on seeds 300 to 323, 24 seeds.

| involution, per day | mass / pre-injury | units / pre-injury | residual | units at the mass floor |
|---|---|---|---|---|
| 0.00 | 0.970 [0.87, 1.05] | 0.70 [0.69, 0.70] | 0.206 | 0.40 |
| 0.02 | 0.966 [0.86, 1.05] | 0.51 [0.30, 0.63] | 0.177 | 0.04 |
| 0.05 | 0.969 [0.86, 1.05] | 0.52 [0.33, 0.63] | 0.178 | 0.04 |
| 0.10 | 0.971 [0.86, 1.06] | 0.52 [0.30, 0.63] | 0.176 | 0.03 |
| 0.20 | 0.973 [0.86, 1.06] | 0.53 [0.38, 0.64] | 0.177 | 0.03 |

Over a tenfold range of the rate the construct keeps 0.51 to 0.53 of its units, against 0.70 with the mechanism absent; the share at the mass floor falls from 0.40 to between 0.03 and 0.04, and the residual from 0.21 to between 0.18 and 0.18. The rate sets how fast, not where it lands.

## Injury after closure (day 140, 30 % of the units removed by a plane cut, read on day 220; the first four entries are day 220 over the day before the injury, the last two the share of what the cut took that is back, each the mean over 24 seeds with the range they span)

| regime | condition | units | mass | mass per unit | perfused mass | share of the lost units made up | share of the lost mass made up |
|---|---|---|---|---|---|---|---|
| liver | as grown | 0.95 [0.77, 1.01] | 0.79 [0.67, 0.86] | 0.84 [0.69, 0.92] | 0.81 [0.67, 0.90] | 0.83 [0.24, 1.05] | 0.18 [0.04, 0.30] |
| liver | reactivated | 0.95 [0.80, 1.02] | 1.82 [1.54, 2.07] | 1.93 [1.73, 2.19] | 1.97 [1.76, 2.18] | 0.82 [0.34, 1.06] | 4.37 [2.79, 6.32] |
| liver | reactivated, fixed inflow | 0.92 [0.72, 1.05] | 1.00 [0.95, 1.03] | 1.10 [0.90, 1.32] | 1.09 [1.02, 1.13] | 0.72 [0.06, 1.18] | 1.01 [0.83, 1.12] |
| kidney | as grown | 0.70 [0.69, 0.70] | 0.71 [0.60, 0.79] | 1.02 [0.85, 1.13] | 0.73 [0.56, 0.82] | 0.00 [0.00, 0.00] | -0.00 [-0.00, 0.00] |
| kidney | reactivated | 0.70 [0.69, 0.70] | 1.85 [1.53, 2.21] | 2.64 [2.18, 3.18] | 2.02 [1.66, 2.41] | 0.00 [0.00, 0.00] | 4.09 [2.31, 6.62] |
| kidney | reactivated, fixed inflow | 0.70 [0.69, 0.70] | 0.97 [0.93, 0.99] | 1.39 [1.34, 1.43] | 1.04 [1.00, 1.08] | 0.00 [0.00, 0.00] | 0.88 [0.68, 0.96] |

## Spread across seeds, and what a fixed inflow does to it (coefficient of variation, %, over 24 seeds)

| regime | readout | growth reactivated, inlet pressure fixed | and inflow held | contraction |
|---|---|---|---|---|
| liver | mass | 6.3 | 2.3 | 2.7 |
| liver | unit count | 5.9 | 6.9 | 0.9 |
| liver | mass per unit | 6.1 | 6.7 | 0.9 |
| kidney | mass | 8.2 | 1.4 | 5.8 |
| kidney | unit count | 0.5 | 0.5 | 1.0 |
| kidney | mass per unit | 8.3 | 1.5 | 5.6 |

## How far a kidney compensates, growth reactivated and inflow fixed

| regime | lost | units kept | mass per unit | full compensation would need | mass |
|---|---|---|---|---|---|
| kidney | 30 % | 0.70 [0.69, 0.70] | 1.39 [1.34, 1.43] | 1.43 | 0.97 [0.93, 0.99] |
| kidney | 50 % | 0.50 [0.50, 0.52] | 1.88 [1.78, 1.94] | 2.00 | 0.95 [0.92, 0.98] |

