# The closed model (work in progress)

`src/organgrow/closed.py` is a rewrite of the engine, kept beside the one the
figures use rather than replacing it. The figures and the audit still run on
`engine.py`; nothing here affects them.

## Why a rewrite

The old engine carried three loops that were open in ways that decided the
answers before any measurement:

- **Flow was demand-driven.** A branch carried the mass hanging from it, so
  the tree's conductance limited nothing, Murray's relation followed from the
  rule rather than from the physics, and the tree could not be a source of
  disorder because nothing fed back into it.
- **Supply was imposed.** `C_root(t)` was the body's allometry entered as an
  input, so much of the growth curve was the datum rather than the result.
- **The three closure modes did not compete.** Saturation had feedback;
  exhaustion and pool were decay curves with none. Asking which mode destroys
  variance had the answer built in: an open loop conserves the spread exactly
  because it reads nothing.
- **Nucleation had no noise**, so all terminal variance entered through
  sigma_param, which is why the calibration against mouse never fitted: in the
  animal the floor is demographic.

## What is in place and verified

**Pressure-driven flow.** Conductance is reduced leaves-to-root, the root
draws Q = dP * C, and that flow splits among daughters in proportion to their
conductances. Mass is conserved to 7e-15. A branch now carries what its
calibre allows.

**One nucleation law.** At each tip, with supply per unit q_i,

    lambda_i = nu * P_i * (q_i - q_set)_+      units made
    mu_i     = sprout * (q_set - q_i)_+        vessels sprouted
    k_i      ~ Poisson(lambda_i)               demographic noise
    dP_i/dt  = (rho - delta) * P_i - lambda_i

One signal with two responses of opposite sign: excess supply makes tissue, a
shortfall makes vessels. The three closure modes are two numbers,
(rho - delta, P_0), not a branch in the code, so which bound binds first is
an outcome rather than an assumption.

**Scales are derived, not fitted.** Pressure dP = 1, length lam = 1, radius
r_root = 1; terminal conductance per unit of mass is 1 / N_scale, so supply
per unit is of order 1 at the design size. Every run starts at q = 1.00
exactly, which is the check that the scales are right.

**A critical radius falls out of the physics.** Shear on a segment in series
with its tissue is not monotone in the radius: tau = Q/r^3 rises while the
segment's own resistance dominates and falls once the tissue's does, peaking
at

    r_c = (3 l / R_t)^(1/4)

Below r_c the feedback is positive and the vessel closes; above it the vessel
holds. This is a pruning criterion derived rather than imposed, and it
predicts that a calibre distribution has a lower cut-off rather than a tail
running to zero. It belongs in the theory.

**A fixed point exists.** With q_set at or above 1 the system settles: supply
per unit holds at 0.933 and radii at 0.461, with no clamping.

## Three loops closed, one window left

**The calibre loop closes.** An isolated vessel with no demand regresses from
r = 1.0 to the floor, so `relax` was never at fault. The engine was forcing
radii above r_crit immediately after relaxing, which cancelled exactly the
regression the rule implies; vessels that fall below r_crit are now pruned
rather than propped up.

**The terminal-conductance loop closes.** A unit's conductance rose with the
mass it had reached, while its mass was capped by the flow that conductance
brought: more mass bought more flow bought more mass. Demand per unit mass is
fixed by the tissue, so conductance per unit is fixed too. Root flow fell from
78 to 0.03, which is the design scale.

**The setpoint is now derived from the coupling point.** tau_set is the shear
a segment carries when its own conductance matches the tissue's below it,

    C_t = units_per_tip / R0,  r_match = (l C_t)^(1/4),
    tau_set = (dP C_t / 2) / r_match^3

Deriving it from the peak of the shear curve instead left it orders of
magnitude low, so vessels widened until they reached it. Radii now settle
between 0.66 and 0.95 rather than at the clamp.

**What is left** is a narrow window in q_set, around 0.75 and 0.9, where the
system still runs away. It converges on either side, including at 0.8, so it
is not a monotone failure: it is a narrow resonance between vessel sprouting
and nucleation, and it appears exactly where the tree branches to three
vessels. Everything else in the swept range converges with supply per unit
between 0.76 and 1.10.

Next: find what sets that window, then calibrate the closed engine against the
mouse contraction the way `engine.py` is calibrated, and only then move the
scripts over.

## What this changes in the manuscript if it works

The central sentence stops being "identity is a closure parameter" and becomes
"identity is which bound binds first, and that choice fixes what fraction of
the final variance is specifiable". The phase diagram in (nu, rho - delta)
replaces Figure 3, and the three noise sources of the model (parametric,
domain, demographic) map one to one onto the systemic, positional and
contingent partition measured in human kidneys, so Figure 6 becomes a
prediction of the engine rather than a separate result.
