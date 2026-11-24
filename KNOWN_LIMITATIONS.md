# What this work does not establish

Stated once, here, rather than hedged throughout.

**Closure by capacity and the measured contraction are not compatible.** The
same feedback that closes the window is the one that contracts the variance,
so they cannot be asked for separately. Parameter choices that let the liver
reach its capacity ceiling within the run contract the CV of unit number by 8
to 52-fold; the choice that reproduces the 3.09-fold contraction measured in
mouse leaves the organ at about a fifth of its ceiling (N = 351 at day 70
against N_inf = C / c_set of roughly 1800). If the mouse liver both closes its
window in a week and contracts only 3.09-fold, capacity-limited nucleation
cannot be the only closure mechanism for that organ.

This is also where `THEORY.md` is internally inconsistent: the section on the
window eliminates sigma(t) for the liver by replacing it with capacity, while
the table of closure modes keeps signal decay (venous sprouting, about one
postnatal week). The engine shows the two are not equivalent. Resolving it is
a decision about the theory, not about calibration, and it is open.

**The engine is calibrated against one measurement.** nu and sigma_param are
set so that the liver reproduces the postnatal contraction measured in mouse
(Azizoglu et al. 2025): unit number grows 1.82-fold while the CV of lobule
number contracts 3.09-fold. With nu = 400, an earlier value chosen to make the
mass curve fit, the engine contracted 7.6-fold, more than twice what the animal
does. One anchor is not a validation; the other three organs are not
independently calibrated.

**The engine is three times more reproducible than the animal.** At the
calibration point the terminal CV of unit number is 0.025 against 0.074
measured. Two candidate explanations were tested and rejected: varying the
organ domain between individuals changes nothing, because supply feedback
absorbs it, and the vascular floor is too small at this unit count. The
observed 0.074 sits just below the Poisson limit of counting 146 lobules
(0.083), so part of it may be measurement, but that cannot be shown from the
published data. The gap is open.

**The contraction law has a regime of validity.** The bound
CV_i/CV_f = (N_f/N_i)^2 describes dilution and fits the mouse. It is the weak
feedback limit of a more general form,

    sd^2(t) = sd_0^2 exp(-2 lambda t) + sd_inf^2 (1 - exp(-2 lambda t))

which fits the engine with R^2 between 0.92 and 0.99. lambda is set by the
closure mode and sd_inf by the vasculature. Where feedback is strong the
contraction exceeds the bound, so the bound is not an upper limit on
contraction but the floor of it.

**Flow is solved on a tree, not a network.** Anastomoses, pulsatility and
non-Newtonian rheology are absent, and pressure is implicit. The shear setpoint
rises as p^(b/2) with b = 0.9, shared by both engines and derived from the cost of mural mass
(Taber 1998; Pries, Reglin & Secomb), not fitted. An earlier version used an
exponent of 2.0, chosen so that supply per unit would fall; that produced a
sign reversal in the central result which was an artefact of the choice.

**What in the closed model is set rather than measured.** The paper's
results come from the closed model (`src/organgrow/closed.py`). Its nucleation
rate (nu = 8), supply threshold (0.6), the three regimes' pools and turnover,
and the density, reach and kill distance of the mesenchymal signal are chosen,
not fitted to data; the junction volume is geometric ($4\pi r^3/3$). The liver
regime uses a pool much larger than the units it makes (2,600 progenitors for
about 200 units), so that its pool does not bind. The one comparison with data
is the mouse liver contraction, which the model exceeds. Results come from one
domain scale (0.6); the transmission uses 24 replicates per condition with
bootstrap intervals, the elasticities 24 paired replicates at five settings
per regime, and the robustness checks 24 seeds in each of thirteen variants.

**Three errors in earlier versions, now fixed, that invalidated results.**
Nucleation truncated its rate to an integer each day, which turned the law
into a staircase of four units per day and hid saturation entirely; the
fractional part is now carried over. The pressure field saturated at zero for
81 % of segments, so the p^(b/2) setpoint operated on a safety clamp rather
than its own form; pressure is now rescaled by the root-to-tip span. And the
per-unit supply factor used the flow its branch carried, which is the mass
hanging from it, making a unit's ceiling proportional to its own mass; that
term is removed and the tree now acts only through the organ conductance.
Figures produced before these fixes are not comparable.

**Spatial exclusion is enforced against the previous step, not within it.**
Candidates in one colonisation iteration are tested against the units that
already exist, not against each other, so a small fraction end up closer than
lambda.

**The simulated and measured invariant sets are not the same.** The engine
produces unit number, unit mass, spacing and supply distance. The human
cohorts, at 100 um to 1 mm per voxel, yield shape and intensity descriptors
only, so a simulated organ cannot yet be scored against them. Closing this
needs the 1 to 6 um zoom volumes.

**The exponent a = 1 for unit number under exhaustion is dilution.** With no
parameter disorder, exhaustion adds a deterministic number of units each day,
so the absolute spread sd(N) is fixed by the seeded start and the CV falls as
1/N exactly. That is the ambiguity THEORY §7 warns of, not a concentration of
the distribution. Source data now carry sd(N) beside CV(N) (Fig4bd, and the
exponent of sd in Fig4c) so the two can be told apart: self-averaging in the
sense of §8 is a falling sd, and a CV exponent equal to the sd exponent plus
one is pure dilution.

**The vascular dispersion at the fastest calibre rate is seed-sensitive.**
Near a relaxation rate of 1 the calibre loop is at the edge of stability, and
CV(C) there changes by an order of magnitude between replicate sets.
`results/vascular_disorder.json` reports a bootstrap interval for every rate;
`THEORY.md` quotes the range, not a point.

**Self-averaging depends on the closure mode, not only on the observable.**
Under exhaustion, unit number concentrates with a = 1.00 exactly, which is
dilution of a constant spread rather than averaging; under saturation the
spread itself falls. Reporting a single exponent per observable hides that
difference.

**Three of six invariants survive correction.** The left-right asymmetry in
human kidneys is tested by paired Wilcoxon with Bonferroni over six invariants
(alpha = 0.0083). Volume, flatness and surface-per-volume survive; surface,
sphericity and elongation do not.

**The human cohorts carry disease.** Kidney donors are elderly and
hypertensive; the lung cohort is COVID-19. Reported variability is an upper
bound on developmental variability.

**The pairwise distance in invariant space is self-normalising.** Dividing by
the cohort's own spread makes the mean squared distance 2 for any
exchangeable sample (and the mean distance at most sqrt(2)), so that value is
a property of the metric, not a measurement.

**Absolute timings in Figure 1E are machine-dependent.**

**Errors found by the external review of 20 November 2026.** The critical
radius was written and implemented as $(3\ell/R_t)^{1/4}$ instead of
$(\ell/3R_t)^{1/4}$, and the flow split at a node ignored tissue hanging from
that node. Both are fixed and covered by tests in `tests/test_smoke.py`; the
second never applied in practice. The time step is one day; the progenitor
processes converge with one, two and four steps a day with the tree held fixed
within the day (THEORY.md §13). The coupled vascular and tissue dynamics were
not run with a smaller step. Sensitivity to where the vascular root
enters and to each sprouting rule is reported in `THEORY.md` §13, Robustness.

**Units are assigned to the nearest tip in space.** Each day a unit is fed by
the terminal branch whose end is closest to it, not by the branch it is
closest to along the tree or by hydraulic resistance. Where branches pass
close to each other a unit can change territory. Assignment by vascular
territory (`assign_mode="segment"`, exact projection on segments) moves 19 to
25 % of the liver's units and
leaves the closure modes and the pancreas-liver contrast unchanged.

**What self-renewal does and does not settle.** Version 4.8 makes the liver's
competent compartment self-renewing (renewal = 1), because mononuclear
hepatocytes are their own progenitors, and keeps the kidney's and the
pancreas's consumable (renewal = 0). This replaces a pool thirteen times
larger than the units it makes with a mechanism, and it reconciles the liver
regime with fetal transcriptomes, in which the hepatoblast programme is fully
converted before birth while lobules keep forming. It does not change the
results: the liver's elasticity to the initial pool is 0.13 with it and 0.22
without, both far below the pancreas at 0.74. The difference between organs therefore rests on the size of the pool
relative to what supply and room allow, not on self-renewal, and the initial
pool cannot be made small, because a tip makes at most as many units per step
as it holds progenitors.

**The spatial constraint is local, not the organ's geometry.** Units are
excluded within one spacing of each other and are placed near a vascular tip,
so what the model calls lack of room is crowding at the vascular front. Across
a 15-fold range of domain volume the units occupy 3 to 17 % of close packing,
and that fraction falls as the domain grows, so the packing limit of the organ
never binds. Unit number scales as $V^{0.37}$ with the domain alone and
$V^{0.74}$ with the allometric scalings of pool, window and rates, against the
$V^{1}$ that a conserved unit size requires. Replicate spread grows with scale.
All regimes in this work are reported at one domain scale.

**How far compensation goes is not right.** Given an injury that restores
the ability of units to enlarge and an inflow held at its pre-injury value,
the model comes close to making up the whole loss: a kidney that loses 30 %
of its units ends with about 1.4 times the mass per unit, against the 1.43 of
full compensation, and about 1.9 against 2.00 after losing half. A real kidney that
loses half grows the remnant by about a half. Nothing in the model limits how
large a unit may become. The vascular and combined injuries were examined only
under the original condition, in which no organ recovers, and have not been
rerun with reactivation and a fixed inflow.

**A budget cannot be spent to the last cell.** This is structural, not a
choice of parameters. A progenitor is tied to the tip it sits at, so where
supply there falls short it strands: with no attrition, no initial pool
between 80 and 300 leaves less than a third of the pool held, and a larger
pool strands more. Two repairs fail for reasons worth stating. An attrition
term empties the pool without building anything, which is what the kidney
regime does now, a quarter of its pool consumed in nephrons and three quarters
lost. A residual drive below the threshold empties a budget but wrecks a
renewing organ, taking the liver from 201 units to 950 and its perfused share
from 0.70 to 0.36; making that drive organ-specific would add a second
property differing exactly as `renewal` does, which is the claim under test.
The variant is kept as `final_induction`, and `scripts/27_pool_fate.py`
regenerates every figure quoted here. What
the model lacks is a per-progenitor induction rate that rises as the pool
falls, which is what the human fetal kidney shows.

**The digits are not reproducible to two decimals across machines.** Every
run is seeded, so one machine repeats a run exactly. The events are discrete
and sit on floating-point comparisons, though, so the same seed on a different
platform takes a different path once two quantities that differ in the last
bits are compared. Run on two machines with three seeds, the kidney's mass per
nephron after a 30 % loss came out 1.36 and 1.39. Both lie inside the range
that 24 seeds span on one machine, 1.35 to 1.42, so what looked like a
platform artefact was the spread of the draw. The tables report the range
beside the mean and the per-seed values are kept in the source data. Summary
statistics that a rerun can move, the coefficients of variation, the paired
contrasts and the counts of variants that clear a threshold, are generated
into `PAPER_NUMBERS.md` rather than written into the prose, so that the
manuscript and the machine that regenerated it cannot disagree.

**What the injury experiment is and is not.** Growing an organ to closure,
destroying a part of it and running on measures recovery rather than further
growth, which the earlier interventions on an intact organ did not. It is not
a validation: it has not been compared with a measured regeneration, and the
comparison that would be natural, a partial hepatectomy, restores mass largely
by hypertrophy and division of the hepatocytes that remain rather than by
building lobules, so the measured observable is mass and not unit count. The
cut is a plane through the organ, the vascular injury is occlusion rather than
excision, and the fraction destroyed was not varied.
