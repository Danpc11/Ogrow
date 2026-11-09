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

## What is still broken

When nucleation stops, the vessel loop does not close. A trajectory shows the
unit count frozen at 27 while tree conductance climbs from 0.0015 to 754: with
no new tissue to feed, each vessel's shear falls, which widens it, which
lowers the shear again. The shrinking tendency does not win even at eight
times its nominal value, which suggests the term is not being applied as
intended rather than being too weak.

First thing to do next: isolate `relax` on a single vessel with no demand and
check that it actually regresses.

## What this changes in the manuscript if it works

The central sentence stops being "identity is a closure parameter" and becomes
"identity is which bound binds first, and that choice fixes what fraction of
the final variance is specifiable". The phase diagram in (nu, rho - delta)
replaces Figure 3, and the three noise sources of the model (parametric,
domain, demographic) map one to one onto the systemic, positional and
contingent partition measured in human kidneys, so Figure 6 becomes a
prediction of the engine rather than a separate result.
