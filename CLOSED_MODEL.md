# The closed model (runs; not yet calibrated)

`src/organgrow/closed.py` is a rewrite of the engine, kept beside the one the
figures use rather than replacing it. The figures and the audit still run on
`engine.py`; nothing here affects them. `make test` exercises it.

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
conductances. Mass is conserved at every branch. A tip with no tissue below it
is a dead end with zero conductance, not an open drain (see below). A branch
carries what its calibre allows.

**Supply is read from pressure, not from the tissue's own draw.** The supply
a unit can take at segment i is what a mature unit would draw at the outlet
pressure the tree delivers there:

    q_i = p_out,i / dP

The flow a branch carries is the tissue's demand, but the pressure it
delivers is the tree's doing, so a unit's ceiling does not depend on the mass
already hanging beside it. The earlier draft used flow per unit, which is
mass per unit, and froze every organ at a few dozen units.

**One nucleation law.** At each segment i that carries tissue or progenitors,

    lambda_i = nu * P_i * (q_i - q_set)_+      units made
    k_i      ~ Poisson(lambda_i)               demographic noise
    dP_i/dt  = (rho - delta) * P_i - k_i       progenitors spent on units made
    mu_i     = sprout * (q_set - q_i)_+        vessels sprouted toward tissue short of supply
             + sprout_room * (placements refused for lack of space)

One signal with two responses of opposite sign: excess supply makes tissue, a
shortfall makes vessels. The second sprouting term is geometric: when a
segment has supply and progenitors but no room within the spacing rule, the
tree reaches into empty space. The three closure modes are two numbers,
(rho - delta, P_0), not a branch in the code.

**Tissue hangs along a vessel, not only at its end.** A segment that has
branched keeps the units and the pool it carried, and a sprout toward a
unit takes with it the units now closer to its end, with the pool split in
the same proportion. Without this, progenitors were stranded on segments the
nucleation loop no longer visited.

**Scales are derived, not fitted.** Pressure dP = 1, length lam = 1, radius
r_root = 1; a mature unit has conductance 1 / N_scale, so an organ of N_scale
mature units matches the tree. Every run starts at q = 1.00 exactly.

**A critical radius falls out of the physics.** Shear on a segment in series
with the conductance below it is not monotone in the radius: tau = Q/r^3
rises while the segment's own resistance dominates and falls once the
tissue's does, peaking at

    r_c = (3 l / R_below)^(1/4)

Below r_c the feedback is positive and the vessel closes; above it the vessel
holds. It is computed for every segment at every step from the real
conductance below it (`VascularTree.r_crit`) and is the floor the segment may
not cross. It predicts a lower cut-off in measured calibre distributions.

**Every parameter other than identity is in `config/params.yaml`,** section
`closed`. A keyword passed to `grow` overrides the file.

## What the earlier draft got wrong, and the fix

The draft reported a calibre runaway when nucleation stopped (conductance
climbing from 0.0015 to 754 with no new tissue) and attributed it to the
shrinking tendency being too weak. The cause was elsewhere: `perfuse` gave a
tip with nothing below it its own full conductance, as if it discharged to
zero pressure. Every sprouted vessel without a unit yet was an open drain;
its shear rose with its radius, so the loop ran away. With dead ends at zero
conductance the tree's conductance is bounded by the tissue's, which the
smoke test now asserts.

## What is still open

- **Calibration.** nu, q_set, rho, delta and P_0 are not fitted to anything.
  The first anchor is the mouse liver (N grows 1.82-fold, CV of lobule number
  contracts 3.09-fold, terminal CV 0.074 against a Poisson floor of 0.083 for
  146 lobules). With Poisson nucleation the floor is demographic, so the
  contraction should come out near the measured one without tension with
  closure; if it does not, the model is wrong, not the calibration.
- **Growth rate.** With the default scales the organ is supply-limited early:
  the pressure drop per segment at the shear setpoint is tau_set * l / r, and
  with q_set near 1 a path of a few segments already falls short. q_set and
  tau_set_fraction have to be set together; the smoke test uses q_set = 0.5.
- **The old engine's analyses are not ported.** Closure comparison,
  self-averaging and identifiability still run on `engine.py`.

## What this changes in the manuscript if it works

The central sentence stops being "identity is a closure parameter" and becomes
"identity is which bound binds first, and that choice fixes what fraction of
the final variance is specifiable". The phase diagram in (nu, rho - delta)
replaces Figure 3, and the three noise sources of the model (parametric,
domain, demographic) map one to one onto the systemic, positional and
contingent partition measured in human kidneys, so Figure 6 becomes a
prediction of the engine rather than a separate result.
