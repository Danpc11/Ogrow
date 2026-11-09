# What this work does not establish

Stated once, here, rather than hedged throughout.

**The engine is not hemodynamic.** It carries no segment radii, no Poiseuille
conductance and no shear-driven calibre relaxation. Supply is a single lumped
conductance attenuated by distance to the nearest vessel (`THEORY.md` Eq. 5).
The benchmark in Figure 1E shows that a tree solve would be affordable at this
scale; it is not part of the growth loop. Adding radii and a shear setpoint is
the natural next version and would change which architectures are reachable.

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

**Self-averaging exponents rest on five size points.** The fitted exponent for
unit number is 0.96 +/- 0.19; the ordering is robust, the values are not
precise.

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
