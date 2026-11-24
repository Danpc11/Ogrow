# Figures, with legends

The numbering is that of `MANUSCRIPT.md`. File names are historical and do not
follow it, so this table is the mapping. The legends below carry the
manuscript number in their heading.

| Manuscript | File | Script | What it carries |
| --- | --- | --- | --- |
| Figure 1 | `Fig1_closed.pdf` | `figures/fig_framework.py` | the framework, and the three ways growth stops |
| Figure 2 | `Fig5_pool.pdf` | `figures/fig_framework.py` | where each pool is held, and what seeding buys |
| Figure 3 | `Fig3_interventions.pdf` | `figures/fig_design.py` | what each relief does to a construct that has stopped |
| Figure 4 | `Fig8_injury.pdf` | `figures/fig_injury.py` | injury after arrest, and what the boundary condition does |
| Figure 5 | `Fig9_transcriptome.pdf` | `figures/fig_transcriptome.py` | when three developing tissues lose the pool they are named for |
| Figure 6 | `Fig6_scaling.pdf` | `figures/fig_design.py` | what domain size costs in density and in reproducibility |

The repository also builds these, which the current manuscript does not use.
They belong to the earlier framing kept in `PAPER.md` and are available as
supplementary material.

| File | Script | What it carries |
| --- | --- | --- |
| `Fig2_closed.pdf` | `figures/fig_closed_paper.py` | the three regimes as they grow |
| `Fig3_closed.pdf` | `figures/fig3_closed.py` | an individual's spread, and how much reaches size |
| `Fig6.pdf` | `figures/fig6.py` | 96 human kidneys: size between donors, shape within |
| `Fig5_closed.pdf` | `figures/fig_closed_paper.py` | the vascular rules |
| `Fig4_closed.pdf` | `figures/fig_closed_paper.py` | contraction of the spread |

**Figure 1. A growing tree whose tips hold supply and progenitors.**
(A) A liver grown in the model at closure, projected on one plane: vessels
(width scaled by radius), units, and tips that hold progenitors with supply
below the threshold. (B) The rule at each tip. (C) Supply per unit relative
to the threshold against the fraction of the progenitor pool left, at
closure, for 24 replicates of each regime.

**Figure 2. What holds each organ's progenitors, and how much of their
number reaches its size.** (A) Share of the initial progenitor pool spent,
held at tips short of supply, held for lack of room (estimated by trial
placements), and still active at closure; mean over the reference seeds of
`scripts/21_robustness.py`. (B) Elasticity of the final unit count to the
progenitor pool at five settings per regime; bars, s.e. over the paired
replicates (`PAPER_NUMBERS.md`, Elasticities). (C) Difference in
that elasticity between pancreas and liver under thirteen variants of the
auxiliary rules, of the progenitor process and of self-renewal, with z;
bars, s.e.

**Supplementary. Three regimes from one law.** (A) Units, (B) supply per unit
relative to the threshold and (C) progenitor pool relative to its initial
value, against day, for one liver, kidney and pancreas.

**Supplementary. Which individual spread reaches the final unit count.**
(A) Final unit count per replicate, relative to the mean without added
spread, for each regime with no added spread, an individual progenitor pool
or an individual supply threshold (24 replicates each); bars mark means.
(B) Transmission of each spread to the final count.

**Supplementary. Size and shape partition in 96 human kidneys from 48 donors.**
(A) Design: distances between donors and between the left and right kidney of
one donor; within-to-between distance ratio, permutation P and the
permutation null (mean ± s.d.). (B) Share of the variance of six invariants
between donors and within a donor. (C) Left minus right difference, as a
percentage of the mean, for each invariant, with the number of donors whose
left kidney is larger and the paired Wilcoxon P; asterisks mark invariants
that survive Bonferroni correction (α = 0.0083). (D) Flatness of the left
against the right kidney, one point per donor, with identity (dashed) and
trend line, slope and r. (E) The three sources of variation.

**Figure 4. Injury after arrest, and the two things recovery needs.** An
organ loses 30 % of its parenchyma on day 140, after closure, and is read on
day 220; mean over the seeds of `scripts/24_injury.py`. (A) Liver mass against day, relative to the day
before the injury, as grown (the growth term has decayed), with growth
reactivated at a fixed inlet pressure, and with growth reactivated and the
root flow held at its pre-injury value; the dashed line marks the injury.
(B) Liver units, mass, mass per unit and perfused mass on day 220, relative
to the day before the injury, under the three conditions. (C) Kidney after
30 % and 50 % loss, with growth reactivated and the inflow held: nephrons
kept, mass per nephron, and the mass per nephron that full compensation would
require (hatched).

**Supplementary. Vascular rules.** (A) Ratio $r_p^{2.9}/\sum r_c^{2.9}$ at each
bifurcation of a grown liver and kidney; bars mark medians. (B) Wall shear
relative to its peak against radius relative to the critical radius, for a
segment in series with its tissue. (C) Conductance of two daughters relative
to their parent against the wall cost exponent $b$; the vertical line marks
$b = 0.9$.

**Supplementary. Contraction of the spread of unit number.** (A) Final unit count
against nucleation rate for eight livers with individual rates, with trend
line, Spearman correlation and P. (B) CV of unit number at closure against
its CV at the start of the measured interval, for a spread in maturity
(offsets of 2 to 10 days), a spread in rate, and the mouse liver.

**Figure 5. When three developing tissues lose the pool they are named for.** (A) Share of the
nephron lineage in each of three states, cap mesenchyme still uncommitted, a
pretubular aggregate or renal vesicle being induced, and a unit already
formed, at five gestational weeks of human fetal kidney (GSE114530); the
number above each bar is the cells of the lineage scored at that week.
(B) Progenitor share of the lineage and number induced per progenitor still
held, against gestational week, on a log scale, one point per sample; the
open square is the progenitor share of an independent week-14 kidney
(GSE139280) read the same way, and the title gives Spearman's rho against
stage and its P for each ratio, over the five samples of the first series. (C) Where the initial progenitor pool goes in the model, for
the two regimes whose pool is a budget, as the mean over the replicates of
`scripts/17_paper_runs.py reps`; the liver is not shown because it renews, so
a share of what it started with is not a denominator for it.
(D) Progenitor share of the lineage against embryonic day in mouse pancreas,
multipotent cells by the authors' annotation (GSE132188), and mouse liver,
cells still running the hepatoblast programme among DLK1-sorted cells
(GSE90047); bars are 95 % intervals over cells, and points at the floor of
the axis are zero. (E) The same three organs on one axis: progenitor share
relative to the first stage measured, against the fraction of the organ's
unit-forming window that has passed, with the window taken from the
literature as weeks 5 to 35 for human nephrogenesis, E9.5 to E15.5 for the
mouse pancreatic multipotent window, and E10.5 to P20 for mouse hepatocyte
number; the dashed line marks the end of the window.

**Figure 3. What each relief does to a construct that has stopped.** Growth is
run to arrest on day 140, one bound is relieved by 20 %, and growth continues
to day 180; 24 seeds per bar, error bars s.e. (A) Units added by day 180
relative to day 140, under no relief, a 20 % lower perfusion setpoint, and a
20 % shorter minimum spacing. (B) Share of the pool held on day 140 that
becomes units by day 180; left blank for a regime that holds nothing.
(C) Progenitors held on day 140, log scale, with the range the seeds span; the
dashed line is the pool below which a released share is not reported.

**Figure 6. What domain size costs.** The renewing regime run at four domain
scales over a 15-fold range of volume, 24 seeds per scale. (A) Units at arrest
against domain volume, one point per seed, with the fitted power of each sweep
and the $V^{1}$ that a conserved unit density would require. (B) Share of the
pool held by perfusion and held by room against domain scale, domain scaled
alone. (C) Coefficient of variation of the unit count across seeds against
domain scale, for both sweeps.
