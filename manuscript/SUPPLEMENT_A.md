# Supplemental information

**Coupled local constraints modulate organ size and growth arrest**

Hugo Tovar, Genaro Vázquez-Victorio, Miguel Ángel Peña-Rico, Angélica S. Martínez-Ramírez, Itzel Escobedo-Avila, Enrique Hernández-Lemus, Daniel Pérez-Calixto

This file contains Figure S1 (provided as a separate file), Notes S1 and S2, and Tables S1 to S8. Every table is a section of `PAPER_NUMBERS.md` as stamped, and is rebuilt by `make paper`.

Model version 4.68.0, commit d67f5bd, generated 2026-11-26 by `scripts/20_paper_numbers.py`. Do not edit by hand.

Python 3.9.6, numpy 1.26.4, scipy 1.13.1; `config/params.yaml` sha256 4d4c28e4b024.

Reference run digest `05e86dcb755f30be` (136 units, 372 vessels) on `Darwin arm64`. Reference ensemble: 133.25 ± 8.75 units over 24 seeds (mean ± s.d.).

---

## Figure S1

**Figure S1. Two organs grown alike, related to Figure 3**

Deviation of one organ from another grown like it: log mass of the kidney-like construct across 24 seeds, and log volume of each kidney relative to its pair in 48 people [S1], with the mean difference between sides removed. The standard deviation of each group is given above it.

## Note S1. Variation between two organs grown alike

Two constructs grown from identical parameters differ only by chance, as do the two kidneys of one person. Across seeds, the standard deviation of the log mass of the kidney-like construct is 0.104. In CT volumes of 48 pairs of human kidneys [S1], the within-pair standard deviation of log volume, after removal of the mean difference between sides, is 0.098 (Figure S1). Nothing in the model was set by the cohort. The agreement is nonetheless limited: part of the measured variation is segmentation error; a mass in model units is compared with a volume; and the variation in the model grows with the size of its domain, whereas within-pair variation in the cohort does not rise detectably with kidney volume (Spearman ρ = +0.25, P = 0.09). We report the comparison for completeness and draw no conclusion from it: the two kidneys of one person share a body, and no hypothesis here connects the variation of the model across seeds to bilateral variation in people.

## Note S2. How to read the tables

All values are from 24 independent seeds of the stochastic model unless stated. Counts and shares are means with their standard errors; contrasts between a perturbed construct and its reference use the same seeds and are ratios of means with 95% bootstrap intervals. The three cell sources are named liver-like, kidney-like, and pancreas-like for the tissues whose cell source they resemble; no parameter was optimized against the measurements with which the model is compared.

**How to read the shares of the pool.** The columns *supply*, *room* and *active* divide the progenitors a construct still holds at arrest among tips whose perfusion has fallen to the threshold, tips above it with no free site, and tips that can still form a unit. They are mutually exclusive, and each is given as a share of the pool that was seeded. *Spent* is the part of the seeded pool no longer held, one minus the pool held over the pool seeded, and is set to zero when the pool has grown. In a source that is spent the four add to one. In a source that divides, the pool held can exceed the pool seeded, so the three held shares can add to more than one while *spent* is zero: 0.21 + 1.03 + 0.18 in one row is a pool that has grown to about 1.2 times what was seeded, not an error of accounting. The shares are averaged over seeds after *spent* is set to zero where the pool has grown, so a row can show a spent share above zero beside held shares that add to more than one. The bound named in the last column is the largest of the four.

## Table S1. The factorial of renewal, seeded pool, and turnover, with the elasticity of each source, related to Figure 2

The three regimes differ in several parameters at once, so the binary property cannot be credited with the separation until the others are held fixed. This varies each independently and reports where the pool stands at arrest. Division and loss are one axis rather than two, swept over the pairings the three sources use, because a pool that divides without being lost grows without bound and because the kidney and the pancreas do not divide at all. The last column names the largest of the four shares, which is the bound that ended the growth (`scripts/34_regimes.py`).

| renewal | pool | turnover | division | loss | units | spent | supply | room | active | what holds growth | elasticity of units to the pool | named regime |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 60 | static | 0.00 | 0.00 | 45.0 ± 4.6 | 0.75 | 0.25 | 0.00 | 0.00 | spent | +0.71 ± 0.05 | pancreas-like |
| 0 | 60 | decaying | 0.00 | 0.05 | 39.5 ± 3.0 | 1.00 | 0.00 | 0.00 | 0.00 | spent | +0.62 ± 0.03 |  |
| 0 | 60 | balanced | 0.10 | 0.10 | 43.3 ± 6.1 | 0.74 | 0.28 | 0.00 | 0.00 | spent | +0.73 ± 0.08 |  |
| 0 | 300 | static | 0.00 | 0.00 | 121.9 ± 8.7 | 0.41 | 0.49 | 0.10 | 0.00 | supply | +0.42 ± 0.05 |  |
| 0 | 300 | decaying | 0.00 | 0.05 | 89.0 ± 7.8 | 1.00 | 0.00 | 0.00 | 0.00 | spent | +0.37 ± 0.03 | kidney-like |
| 0 | 300 | balanced | 0.10 | 0.10 | 116.0 ± 11.9 | 0.41 | 0.51 | 0.08 | 0.00 | supply | +0.48 ± 0.06 |  |
| 0 | 2600 | static | 0.00 | 0.00 | 212.6 ± 15.1 | 0.08 | 0.74 | 0.18 | 0.00 | supply | +0.14 ± 0.04 |  |
| 0 | 2600 | decaying | 0.00 | 0.05 | 154.8 ± 13.1 | 1.00 | 0.00 | 0.00 | 0.00 | spent | +0.24 ± 0.03 |  |
| 0 | 2600 | balanced | 0.10 | 0.10 | 211.0 ± 16.8 | 0.12 | 0.70 | 0.19 | 0.00 | supply | +0.24 ± 0.04 |  |
| 1 | 60 | static | 0.00 | 0.00 | 87.8 ± 18.8 | 0.00 | 0.75 | 0.25 | 0.01 | supply | +0.44 ± 0.12 |  |
| 1 | 60 | decaying | 0.00 | 0.05 | 51.5 ± 7.7 | 1.00 | 0.00 | 0.00 | 0.00 | spent | +0.51 ± 0.07 |  |
| 1 | 60 | balanced | 0.10 | 0.10 | 77.4 ± 21.0 | 0.21 | 1.03 | 0.18 | 0.00 | supply | +0.62 ± 0.11 |  |
| 1 | 300 | static | 0.00 | 0.00 | 143.7 ± 21.0 | 0.00 | 0.78 | 0.22 | 0.00 | supply | +0.18 ± 0.06 |  |
| 1 | 300 | decaying | 0.00 | 0.05 | 98.4 ± 8.4 | 1.00 | 0.00 | 0.00 | 0.00 | spent | +0.34 ± 0.05 |  |
| 1 | 300 | balanced | 0.10 | 0.10 | 140.3 ± 12.9 | 0.06 | 0.78 | 0.32 | 0.01 | supply | +0.23 ± 0.07 |  |
| 1 | 2600 | static | 0.00 | 0.00 | 223.7 ± 22.9 | 0.00 | 0.80 | 0.20 | 0.00 | supply | +0.13 ± 0.04 |  |
| 1 | 2600 | decaying | 0.00 | 0.05 | 152.0 ± 14.8 | 1.00 | 0.00 | 0.00 | 0.00 | spent | +0.22 ± 0.04 |  |
| 1 | 2600 | balanced | 0.10 | 0.10 | 208.1 ± 17.9 | 0.06 | 0.78 | 0.19 | 0.00 | supply | +0.07 ± 0.06 | liver-like |

Does the bound that holds the pool predict the elasticity? Across the 18 sources it sorts them coarsely: spent, 8 sources, +0.22 to +0.73; supply, 10 sources, +0.07 to +0.62. The Spearman correlation of the elasticity with the share of the pool spent is +0.52 (P = 0.028), and with the units reached -0.92 (P = 8.4e-08). The size a source reaches orders the elasticity; the label of its bound does not.

Each axis is asked the same question: holding the other two fixed, in how many of those settings does moving along this axis change which bound ends the growth?

| axis varied | settings of the other two | in how many the bound changes |
|---|---|---|
| renewal | 9 | 2 |
| turnover | 6 | 5 |
| pool | 6 | 2 |

The elasticity of the final unit count to the seeded pool is reported elsewhere at each regime's own pool size, so a renewing source at 2600 is compared with a spent one at 60. Here the pool and the turnover are matched between the two sources and only the binary property differs.

| pool | turnover | spent source | renewing source | gap |
|---|---|---|---|---|
| 2600 | balanced | +0.24 ± 0.04 | +0.07 ± 0.06 | +0.17 |
| 300 | decaying | +0.37 ± 0.03 | +0.34 ± 0.05 | +0.04 |
| 60 | static | +0.71 ± 0.05 | +0.44 ± 0.12 | +0.28 |

## Table S3. Elasticities of the final unit count under alternative settings of the nucleation law, related to Figure 2

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

## Table S4. Robustness of the three sources to the auxiliary rules of the model, related to Figures 1 and 2

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

**The ordering, counted over the variants**

| statement | variants holding it | sign test P | of those, z at least 2 | at the reference setting |
|---|---|---|---|---|
| liver below kidney below pancreas | 13 of 13 | - | - | - |
| kidney above liver | 13 of 13 | 0.0002 | 12 of 13 | +0.31 ± 0.08, z = 4.1 |
| pancreas above kidney | 13 of 13 | 0.0002 | 12 of 13 | +0.34 ± 0.06, z = 5.7 |
| pancreas above liver | 13 of 13 | 0.0002 | 12 of 13 | +0.65 ± 0.08, z = 7.8 |

## Table S5. The reference model beside published progenitor ablations and paired human kidneys, related to Figure 3

The only numbers here that are set beside a measurement. A share of the reference pool is seeded and what is left of the construct is read, as a ratio of means over 24 paired seeds with a 95 % bootstrap interval. The elasticity is log(mass ratio) / log(pool kept).

| regime | pool kept | units | mass | elasticity of mass |
|---|---|---|---|---|
| pancreas | 0.92 | 0.90 [0.86, 0.93] | 0.81 [0.72, 0.90] | +2.57 |
| pancreas | 0.5 | 0.58 [0.56, 0.61] | 0.53 [0.42, 0.64] | +0.91 |
| pancreas | 0.39 | 0.46 [0.43, 0.48] | 0.44 [0.37, 0.52] | +0.87 |
| pancreas | 0.25 | 0.32 [0.30, 0.33] | 0.17 [0.12, 0.23] | +1.27 |
| kidney | 0.6 | 0.80 [0.77, 0.83] | 0.84 [0.79, 0.88] | +0.35 |
| kidney | 0.5 | 0.74 [0.71, 0.77] | 0.76 [0.73, 0.80] | +0.39 |
| liver | 0.5 | 0.92 [0.88, 0.97] | 0.92 [0.90, 0.95] | +0.11 |
| liver | 0.35 | 0.87 [0.82, 0.93] | 0.86 [0.84, 0.89] | +0.14 |
| liver | 0.25 | 0.80 [0.76, 0.84] | 0.77 [0.75, 0.79] | +0.19 |

What the experiments measured, as a share of the littermate control, with the model at the same share of the pool. The model column is mass for pancreas and liver and units for the kidney, to match what was measured. The numbers are typed from the papers; `scripts/36_validation.py` says which figure each is read from.

| regime | pool left | measured | what, and by whom | model |
|---|---|---|---|---|
| pancreas | 0.92 | 0.77 [0.53, 1.02] | weight at birth, Stanger 2007 | 0.81 [0.72, 0.90] |
| pancreas | 0.39 | 0.36 [0.19, 0.53] | weight at birth, Stanger 2007 | 0.44 [0.37, 0.52] |
| kidney | 0.60 | 0.88 | glomeruli at P3, Muthukrishnan 2018 | 0.80 [0.77, 0.83] |
| kidney | 0.60 | 0.79 | glomeruli at P10, Muthukrishnan 2018 | 0.80 [0.77, 0.83] |
| kidney | 0.60 | 0.48 | glomeruli at E19.5, Cebrian 2014 | 0.80 [0.77, 0.83] |
| kidney | 0.60 | 0.58 | glomeruli at P50, Cebrian 2014 | 0.80 [0.77, 0.83] |
| liver | 0.35 | 0.96 | weight four days after the ablation ended, Stanger 2007 | 0.86 [0.84, 0.89] |

The regime parameters were fixed before this comparison was computed and were not adjusted to it. The two kidney studies removed the same share of progenitors and differ twofold in what was left; the model agrees with one and has the mechanism of neither.

Variation across 24 seeds grown from identical parameters, which is what two organs of one individual share.

| regime | CV of units | CV of mass | s.d. of log mass |
|---|---|---|---|
| pancreas | 0.103 | 0.168 | 0.172 |
| kidney | 0.088 | 0.096 | 0.104 |
| liver | 0.086 | 0.057 | 0.058 |

The two kidneys of one person, CT volume, 48 pairs (TRUSTED): left-right standard deviation of log volume, per kidney, 0.098 with the mean difference between sides (-7.4 %) removed and 0.111 with it left in; between people 0.322. The left-right difference against kidney size gives Spearman rho = +0.25, P = 0.09, over 109 to 553 cm3. The volumes carry segmentation error, which is not separated here.

## Table S6. The extension in which progenitors drive branching, in every setting examined, related to Figure 4 and Table 2

Not the reference model, and no number above depends on it. A source of one cap of 60 progenitors that divides at 0.05 a day, is spent in making a unit, and lets a segment branch only while it holds two caps, grown for 120 days. On day 45 each progenitor dies with probability 0.4. Ablated over control, 24 paired seeds, 95 % bootstrap intervals. The response is the fractional rise of the division rate in a cap below its size. The rules were written after the measurements below were known.

| nucleation | q_set | nu | response | units, control | pool | cap, next day | cap back to 0.95 | sites | units | units since the ablation |
|---|---|---|---|---|---|---|---|---|---|---|
| reference | 0.6 | 8.0 | 0 | 429 | 0.54 [0.49, 0.58] | 0.65 | 50 d | 0.55 [0.50, 0.60] | 0.81 [0.72, 0.94] | 0.77 [0.68, 0.92] |
| reference | 0.6 | 8.0 | 1 | 462 | 0.68 [0.66, 0.70] | 0.67 | 9 d | 0.71 [0.68, 0.75] | 0.85 [0.75, 0.95] | 0.84 [0.72, 0.94] |
| threshold halved | 0.3 | 8.0 | 0 | 384 | 0.53 [0.48, 0.59] | 0.61 | 10 d | 0.52 [0.44, 0.60] | 0.73 [0.64, 0.86] | 0.65 [0.52, 0.82] |
| threshold halved | 0.3 | 8.0 | 1 | 764 | 0.65 [0.62, 0.68] | 0.68 | 11 d | 0.67 [0.65, 0.70] | 0.88 [0.85, 0.90] | 0.86 [0.84, 0.89] |
| threshold halved, law slowed | 0.3 | 2.0 | 0 | 421 | 0.54 [0.49, 0.60] | 0.69 | 3 d | 0.60 [0.54, 0.66] | 0.79 [0.71, 0.91] | 0.74 [0.63, 0.88] |
| threshold halved, law slowed | 0.3 | 2.0 | 1 | 786 | 0.67 [0.65, 0.70] | 0.68 | 25 d | 0.72 [0.69, 0.75] | 0.88 [0.84, 0.91] | 0.87 [0.83, 0.91] |
| threshold quartered, law slowed | 0.15 | 2.0 | 0 | 384 | 0.57 [0.53, 0.62] | 0.90 | 7 d | 0.58 [0.51, 0.66] | 0.64 [0.57, 0.70] | 0.52 [0.41, 0.61] |
| threshold quartered, law slowed | 0.15 | 2.0 | 1 | 886 | 0.66 [0.63, 0.68] | 0.67 | 14 d | 0.71 [0.68, 0.73] | 0.86 [0.83, 0.89] | 0.85 [0.81, 0.88] |

Measured, ablated over control. Without a proliferative response (Cebrian et al. 2014): cells per tip 0.60, 0.95 and 1.04 one, two and three days after tamoxifen; tips 1.06, 0.71 and 0.59; total progenitors growing as fast as in control; glomeruli 0.48 at birth and 0.58 in the adult. With one (Muthukrishnan et al. 2018): progenitors 0.60 after a day and 0.85 after two; glomeruli 0.88 and 0.79 three and ten days after birth. Model days are not embryonic days.

## Table S8. Predictions for ablations that have not been performed, related to Figure 5

Numbers for experiments nobody has reported, stamped with the commit above. Once one is done, the number that stood here before it is the prediction, whatever a later version gives. Ratio to the unablated control, 24 paired seeds, 95 % bootstrap intervals.

Dose, reference model: mass of the construct against the share of its pool left.

| regime | pool left | mass | what refutes it |
|---|---|---|---|
| pancreas | 0.75 | 0.70 [0.58, 0.81] | a pancreas of normal weight |
| pancreas | 0.6 | 0.61 [0.51, 0.70] | a pancreas of normal weight |
| liver | 0.75 | 0.97 [0.94, 1.00] | a deficit after a shallow ablation, or none after a deep one |
| liver | 0.15 | 0.70 [0.68, 0.72] | a deficit after a shallow ablation, or none after a deep one |

Depth and timing, the extension, perfusion limiting nucleation least: units of the coupled source.

| response | day of the ablation | share killed | units | status |
|---|---|---|---|---|
| 0 | 45 | 0.2 | 0.82 [0.75, 0.89] | prediction |
| 0 | 45 | 0.4 | 0.64 [0.57, 0.70] | compared with a measurement |
| 0 | 45 | 0.6 | 0.56 [0.50, 0.63] | prediction |
| 0 | 20 | 0.4 | 0.67 [0.56, 0.79] | prediction |
| 0 | 80 | 0.4 | 0.75 [0.71, 0.79] | prediction |
| 1 | 45 | 0.2 | 0.93 [0.90, 0.97] | prediction |
| 1 | 45 | 0.4 | 0.86 [0.83, 0.89] | compared with a measurement |
| 1 | 45 | 0.6 | 0.79 [0.77, 0.82] | prediction |
| 1 | 20 | 0.4 | 0.93 [0.86, 1.05] | prediction |
| 1 | 80 | 0.4 | 0.87 [0.84, 0.90] | prediction |

What refutes the kidney predictions: a deficit in proportion after a 60 % ablation without a proliferative response, and the same deficit whenever the ablation is done.

## Supplemental references

[S1] Ndzimbong, W., Fourniol, C., Themyr, L., Thome, N., Keeza, Y., Sauer, B., Piéchaud, P.-T., Méjean, A., Marescaux, J., George, D., et al. (2025). TRUSTED: the paired 3D transabdominal ultrasound and CT human data for kidney segmentation and registration research. Sci. Data *12*, 615. https://doi.org/10.1038/s41597-025-04467-1.
