# Figures, with legends

The numbering is that of `MANUSCRIPT.md`. File names are historical and do not
follow it, so this table is the mapping. The legends below carry the
manuscript number in their heading.

| Manuscript | File | Script | What it carries |
| --- | --- | --- | --- |
| Figure 1 | `Fig1_model.pdf` | `figures/fig_paper.py` | the law, what can run out, the one property, a grown construct, where each pool ends |
| Figure 2 | `Fig2_bounds.pdf` | `figures/fig_paper.py` | what seeding buys, what the threshold buys, the ordering, what each relief does |
| Figure 3 | `Fig8_injury.pdf` | `figures/fig_injury.py` | injury after arrest, and what the boundary condition does |
| Figure 4 | `Fig9_transcriptome.pdf` | `figures/fig_transcriptome.py` | when three developing tissues lose the pool they are named for |
| Figure 5 | `Fig6_scaling.pdf` | `figures/fig_design.py` | what domain size costs in density and in reproducibility |
| Figure 6 | `Fig10_validation.pdf` | `figures/fig_validation.py` | three numbers of the model beside published measurements |

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

**Figure 1. The model, and a construct grown under it.** Panels A to C and F
are schematic and carry no data; proportions in them are not to scale. Colour
in those panels marks the three resources rather than the three regimes:
supply purple, competent cells teal, units and the room between them ink.
(A) A vascular tip, the tissue it serves, the supply per unit it delivers, the
competent cells it holds and the least spacing between units already built,
with the nucleation law. (B) The three quantities that can run out at a tip,
each with the state it ends in. (C) The one property varied between regimes:
whether the cell that forms a unit is spent in forming it. (D) One renewing
construct at arrest, projected on a plane: vessels with calibre as line width,
units, and the tips that still hold competent cells below the supply
threshold. (E) Disposition of the initial pool at arrest for each regime, mean
over 24 seeds; the annotation splits what the spent-with-loss regime spent
into the share that built units and the share lost to attrition. (F) The bound
that holds each regime at arrest, and whether growth can restart.

**Figure 2. What the binding resource is worth.** 24 seeds throughout, error
bars s.e. (A) Elasticity of the final unit count to the seeded pool, and
(B) to the perfusion setpoint, at five parameter settings per regime; points
are the mean over paired replicates. (C) Paired difference in pool elasticity
between the spent-no-loss and the renewing regime under thirteen variants of
the auxiliary rules, with z. (D) Units added by day 180 relative to day 140
after relieving one bound on day 140, (E) share of the day-140 pool released
by the same relief, left blank for a regime that holds nothing, and
(F) progenitors held on day 140 with the range the seeds span, log scale; the
dashed line is the pool below which a released share is not reported.

**Figure 3. Injury after arrest, and the two things recovery needs.** An
organ loses 30 % of its parenchyma on day 140, after closure, and is read on
day 220; mean over the seeds of `scripts/24_injury.py`. (A) Liver mass against day, relative to the day
before the injury, as grown (the growth term has decayed), with growth
reactivated at a fixed inlet pressure, and with growth reactivated and the
root flow held at its pre-injury value; the dashed line marks the injury.
(B) Liver units, mass, mass per unit and perfused mass on day 220, relative
to the day before the injury, under the three conditions. (C) Kidney after
30 % and 50 % loss, with growth reactivated and the inflow held: nephrons
kept, mass per nephron, and the mass per nephron that full compensation would
require (hatched). (D) Day-220 mass against the flow the
circulation delivers, both relative to the day before the injury, with the
demand per unit mass measured once on the intact ensemble and held fixed; one
point per seed, mean with s.d., the fitted exponent, and the dashed line the
mass the supply budget would give if every unit reached its ceiling.

**Figure 4. When three developing tissues lose the pool they are named for.** (A) Share of the
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

**Figure 5. What domain size costs.** The renewing regime run at four domain
scales over a 15-fold range of volume, 24 seeds per scale. (A) Units at arrest
against domain volume, one point per seed, with the fitted power of each sweep
and the $V^{1}$ that a conserved unit density would require. (B) Share of the
pool held by perfusion and held by room against domain scale, domain scaled
alone. (C) Coefficient of variation of the unit count across seeds against
domain scale, for both sweeps.

**Figure 6. Three numbers of the model beside measurements.** Filled symbols
are the model, open symbols are measurements, and nothing is fitted. (A) Share
of the organ that results against the share of the progenitor pool left.
Model: mass of the pancreas-like and liver-like constructs at day 140 as a
ratio to the reference, 24 paired seeds, bars the 95 % bootstrap interval.
Measured: weight of the mouse pancreas at birth after 8 % and 61 % of its
progenitors were removed, with 95 % intervals, and of the liver four days
after at least 65 % of its progenitors were removed [2]. The dashed line is an
organ in proportion to its pool. (B) The same for the kidney-like construct,
in units, against glomerular number in two mouse studies that each removed
about 40 % of nephron progenitors: through a CITED1 driver, counted three and
ten days after birth [15], and through a Gdnf driver, counted at birth and in
the adult [14]. The two are offset sideways to be told apart. (C) How far one
organ sits from the other grown like it: the kidney-like construct across 24
seeds, as the deviation of log mass from the mean, and the two kidneys of each
of 48 people, as the deviation of log volume from the pair with the mean
difference between sides removed [16]. Bars span one standard deviation
either side of zero.

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
