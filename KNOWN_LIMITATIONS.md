# What this work does not establish

Stated once, here, rather than hedged throughout.

**Closure by capacity and the measured contraction are not compatible.** The
same feedback that closes the window is the one that contracts the variance,
so they cannot be asked for separately. Parameter choices that let the liver
reach its capacity ceiling within the run contract the CV of unit number by 8
to 52-fold; the choice that reproduces the 3.09-fold contraction measured in
mouse leaves the organ at half its ceiling. If the mouse liver both closes its
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
rises as p^(b/2) with b = 1, which is derived from the cost of mural mass
(Taber 1998; Pries, Reglin & Secomb), not fitted. An earlier version used an
exponent of 2.0, chosen so that supply per unit would fall; that produced a
sign reversal in the central result which was an artefact of the choice.

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
the cohort's own spread makes the expected distance sqrt(2) for any
exchangeable sample, so that value is a property of the metric, not a
measurement.

**Absolute timings in Figure 1E are machine-dependent.**
