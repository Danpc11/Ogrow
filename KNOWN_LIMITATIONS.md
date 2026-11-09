# What this work does not establish

Stated once, here, rather than hedged throughout.

**Flow is solved on a tree, not a network.** Anastomoses, pulsatility and
non-Newtonian rheology are absent, and pressure is implicit rather than carried
as a variable. Murray's law emerges rather than being imposed, which is a check
on the implementation, not a validation against measured calibres.

**The central result changed when hemodynamics were added.** With a lumped
conductance proportional to mass, capacity saturation destroyed variance. With
calibres that relax toward a shear setpoint, it amplifies it, because the tree
disperses supply between replicates and saturation is the only closure mode
that reads supply. Results from the earlier engine are not comparable.

**Spatial exclusion is enforced against the previous step, not within it.**
Candidates in one colonisation iteration are tested against the units that
already exist, not against each other, so a small fraction end up closer than
lambda: about 1.5 % in liver and 4.5 % in kidney at day 70. The stated rule is
"only if", so this is a known gap between the equation and the code.

**The simulated and measured invariant sets are not the same.** The engine
produces unit number, unit mass, spacing and supply distance. The human
cohorts, at 100 um to 1 mm per voxel, yield shape and intensity descriptors
only. A simulated organ therefore cannot yet be scored against the human
cohort. Closing this needs the 1 to 6 um zoom volumes.

**Self-averaging depends on the closure mode, not only on the observable.**
Under exhaustion, unit number concentrates with a = 1.00; under capacity
saturation the same quantity disperses, a = -2.45 +/- 0.77. The difference is
not an artefact of the fit: with ten size points over a factor of three in N,
the saturating case gives a = -3.64 +/- 0.37 at R^2 = 0.92. Figure 4 therefore
reports two closure modes rather than one, and the claim that additive
quantities self-average holds only where the organ does not read a supply that
the vasculature disperses.

**Finite differences run on a stochastic process with integer truncation.**
Nucleation truncates to whole units, so the Jacobian, and with it the rank and
condition number, carry truncation noise. The condition number for
N + per-unit supply moves between 65 and 119 across seeds;
`scripts/04_identifiability.py` reports that spread.

**The pairwise distance in invariant space is self-normalising.** Dividing by
the cohort's own spread makes the expected distance sqrt(2) for any
exchangeable sample, so that value is a property of the metric, not a
measurement. `scripts/06_exchangeability.py` reports it against a random
reference to make this explicit. What does carry information are quantities in
a fixed metric: raw dispersions, the within-to-between ratio, and the annotator
disagreement.

**Segmentation is not validated against ground truth.** The boundary-gradient
method is gated by `qc.py` and an outlier rule, and every rejection is recorded
in `results/hoa_rejected.csv`, but it has not been compared with manual
annotation.

**The human cohorts carry disease.** Kidney donors are elderly and
hypertensive; the lung cohort is COVID-19. Reported variability is an upper
bound on developmental variability.

**Three of six invariants survive correction.** The left-right asymmetry is
tested by paired Wilcoxon with Bonferroni over six invariants
(alpha = 0.0083). Volume, flatness and surface-per-volume survive; surface,
sphericity and elongation do not.

**Absolute timings in Figure 1E are machine-dependent.** The ordering and the
linear slope reproduce; the milliseconds do not.
