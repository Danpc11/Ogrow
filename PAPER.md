# The paper: framework, claims and evidence

Target: *Frontiers in Cell and Developmental Biology*, Hypothesis and Theory.
Reference model: version 4.14.0 (`src/organgrow/closed.py`). Every number from
the model that the paper uses is in `PAPER_NUMBERS.md`, written by
`make paper` with the version and commit that produced it; this file does not
copy model numbers by hand. Numbers given here are either data (the human
cohort) or analytic. Model numbers depend on the random stream, which numpy
does not keep fixed across versions, so they are produced with the versions
pinned in `requirements.txt` and recorded in the stamp.

## The question

Some organs rebuild themselves after damage and others never do. The usual
account of that difference is about stem cells: an organ regenerates if it
keeps a reserve. We argue that the difference is set earlier, by whatever
stopped the organ growing in the first place, and that the mechanism of
arrest and the capacity to regenerate are the same mechanism seen twice.

## The hypothesis

An organ stops making functional units when a resource runs out at the tips of
its growing vascular tree. Three resources are local to a tip: perfusion per
unit, cells competent to form a unit, and room to place one. Which of them
binds first is what distinguishes a liver from a kidney, it decides which
individual variation reaches final size, and it decides whether growth can
ever start again.

## One property separates the organs

Organs are not given different rules here. They are given one difference:
whether the cell that forms a unit is spent in forming it.

In the liver it is not. Mononuclear hepatocytes are their own progenitors, so
the compartment is restored by the very cells it produces and forming a
lobule does not draw it down. In the kidney the cap mesenchyme is a separate
population, consumed and not replaced, which is why no nephrons form after
birth; in the pancreas the multipotent progenitors are likewise spent. One
parameter carries this, and everything below follows from it together with the
pool each organ starts with.

## Three ways to stop

With that single difference the three organs close in three ways. The kidney
spends its pool to the last cell. The pancreas spends most of it and keeps a
remainder at tips short of supply. The liver keeps nearly all of it, held
partly where perfusion has fallen to threshold and partly where existing units
leave no room. These are read from the distribution of the pool over tips at
closure, with the room at each tip tested by trial placements independent of
the dynamics, not from an average. They survive thirteen variants of the rules
that are not the nucleation law. `PAPER_NUMBERS.md`, Regimes and Robustness;
Figures 1 and 2A, and Figure 3 for the trajectories.

## In development: which resource binds decides what reaches size

If the pool binds, the number of progenitors an individual starts with sets
the organ it ends with. If supply binds, it does not. The elasticity of final
unit number to the initial pool is accordingly high in the pancreas and low in
the liver, in every variant where it is measurable.

This is the pattern Stanger, Tanaka and Melton found in mice (*Nature* 2007):
halving the pancreatic progenitors halved the pancreas, which never
compensated, while the liver compensated fully. The model was not fitted to
that experiment. It does not separate the kidney from the liver robustly, so
it orders two organs and not three. `PAPER_NUMBERS.md`, Transmission,
Elasticities and Robustness; Figures 2B, 2C and 4.

## In people: size is systemic, shape is positional

In 96 kidneys from 48 donors (TRUSTED, clinical CT), 89 % of the variance in
volume lies between donors and 80 % of the variance in sphericity within a
donor; the within-to-between distance ratio is 1.54 (permutation
$P$ = 0.0005). Left and right differ systematically: by paired Wilcoxon with
Bonferroni over six invariants, volume ($P$ = 0.0007), flatness and
surface-per-volume ($P$ < 0.0001) survive. What an individual carries sets how
large their kidneys are; where each kidney sits sets its shape. Figure 5.

## In regeneration: the way an organ stopped decides whether it can restart

Closure is stable in all three regimes: forty more days add little. It differs
in whether it can be undone. Relieving perfusion, by lowering the threshold,
or room, by shortening the least spacing, restarts the liver's growth on the
day it closed. Nothing restarts the kidney's, because its pool is spent and no
relief returns it. The liver's arrest is a standoff between resources it still
has; the kidney's is an exhaustion. That the adult liver regrows and the adult
kidney makes no new nephrons is the correspondence this predicts, and it is a
correspondence rather than a test. `PAPER_NUMBERS.md`, Interventions.

## But restarting is not restoring, and recovery needs two things

An organ that can restart is not an organ that recovers. Grown to closure,
then stripped of 30 % of its units by a plane cut and left for eighty more
days, the liver makes back 85 % of the units it lost but only a sixth of the
mass, and its perfused function tracks the mass. The reason is in the model's
own rules:
the rate at which a unit grows decays through development, so by closure a
unit can no longer enlarge, and an organ that rebuilds after closure fills
with units that never mature.

Recovery therefore needs the first of two things that growth does not supply.
**Injury must restore a unit's ability to enlarge**, which is what the growth
factors of regeneration do. Given that, both organs recover their mass.

It also needs the second. **The circulation must cap how far they go.** Under
the model's own boundary condition, a fixed pressure at the inlet, the tree
grows its conductance as the tissue demands more and the supply per unit never
falls, so the organ overshoots to nearly twice the mass it had. Holding the
inflow instead at what the intact organ was receiving, which is what a gut
does to a resected liver, each unit's share is its own conductance times the
pressure its tip delivers, so delivery per unit mass falls as the organ grows.
The organ then restores its mass and stops. The rule stays local, a unit grows while its own
delivery per mass exceeds a threshold; the ceiling on the whole organ follows
from it. `PAPER_NUMBERS.md`, Injury; Figure 6.

## The organ that cannot rebuild units enlarges the ones it kept

How an organ recovers follows from what it spent in development. The kidney
closed with no competent cells left, so it cannot make nephrons and can only
enlarge what it has: after losing 30 % of its units it ends with 0.70 of the
count and about 1.4 times the mass per unit, its function restored. The
constraint also decides how repeatable that landing is, and the way it does
says what the constraint is. Across 24 seeds, holding the inflow contracts the
spread of final mass about 3.7 fold in both organs, from 8.5 % of the mean to
2.3 % in the liver and from 6.1 % to 1.7 % in the kidney, with bootstrap 95 %
intervals of 2.5 to 5.6 and 2.6 to 5.2. It leaves the spread of unit count
untouched: 7.2 % against 7.9 % in the liver, 0.4 % against 0.4 % in the
kidney. A fixed budget fixes the mass and says nothing about how an organ
divides it. The kidney's count is already pinned by its own biology, so its
size per nephron is forced as well, 6.0 % down to 1.6 %; the liver's count is
free to absorb the variation, so its size per unit stays loose, 9.6 % against
7.2 %. The property that separates the two organs in growth separates how
they recover. That is
compensatory renal hypertrophy as it is described after a nephrectomy, where
no new nephrons form, the remaining ones enlarge and function is maintained.
The liver, whose competent compartment renews itself, has both routes open.

The model compensates further than a kidney does. It comes close to making up
the whole loss, reaching a mass per nephron of about 1.4 after a 30 % loss and
about 1.9 after 50 %, against the 1.43 and 2.00 that full compensation
requires, whereas a kidney that loses half its mass grows the remnant by about
a half rather than doubling it.
Nothing in the model limits how large a unit may become; in a kidney that
limit exists and passing it is hyperfiltration injury.

## The machinery underneath

**One wall-cost exponent sets the vascular rules.** Minimising
$\mathcal{C} = Q^2/r^4 + k\,(p\,r^2)^b$ per unit length, with wall
cross-section $\propto p\,r^2$ by Laplace's law, gives branching
$r_p^{b+2} = \sum r_c^{b+2}$ and wall shear $\propto p^{b/2} r^{b-1}$. Two
daughters conduct $2\cdot 2^{-4/(b+2)}$ of their parent, below one for every
$b < 2$, so a branch buys reach and not conductance. Shear on a segment in
series with its tissue peaks at $r_c = (\ell/3R_t)^{1/4}$, the radius below
which a vessel regresses. `PAPER_NUMBERS.md`, Murray; Figure 7.

**Capacity closure contracts the spread of unit number.** If every individual
gains the same increment of $N^2$, the ratio of initial to final CV is
$(N_f/N_i)^2$, which gives 3.31 for the mouse liver against 3.09 measured
(Azizoglu et al. 2025). `PAPER_NUMBERS.md`, Contraction; Figure 8.

## What this is for

Two consequences follow for work that does not use this model.

**Counting units is the wrong readout for regeneration.** A tissue that
rebuilds after its growth window has closed recovers its unit number long
before, or without ever, recovering its perfused mass. A study that measures
recovery by counting would call such a tissue healed.

**Perfusion comes before cells when building an organ.** In an organ whose
competent compartment renews itself, adding progenitors barely changes final
size, while the perfusion threshold passes almost entirely into it. The
requirements that follow, in physical units for human liver, are in the
Discussion.

## What would test it

The results above are predictions, not measurements. Each of these would
decide one of them.

- **Count is the wrong readout for regeneration.** A tissue that rebuilds
  units after its growth window has closed, without restoring their ability to
  enlarge, recovers its unit number long before, or without ever, recovering
  its perfused mass. Measuring recovery by counting units would call such a
  tissue healed.
- **What stops compensatory growth is outside the organ.** Holding the inflow
  to a regenerating organ at its pre-injury value should let it restore its
  mass and stop there; raising the inflow with the demand should let it
  overshoot. The brake is the delivery the circulation provides, not a size
  the tissue remembers.
- **The vasculature of lost tissue is worth closing.** After parenchymal
  loss, removing or occluding the vessels that served the lost region should
  improve the function the remaining tissue recovers, because they compete for
  the same inflow.
- **A budget organ tracks its progenitors.** Pancreas size should scale with
  progenitor number, as observed, while the liver's should not, because the
  liver's competent compartment is restored by the cells it makes; a kidney whose
  progenitors are also lost during induction should track them weakly.
- **Paired kidneys.** With the model's elasticities, the share of size
  variance between individuals is
  $(\varepsilon_q^2\sigma_{q,s}^2 + \varepsilon_P^2\sigma_{P,s}^2)$ over that
  plus the local and contingent terms; paired nephron counts and split renal
  function in healthy donors would turn this into a test.

## Not claimed

- The figures of the open engine (`make all`), in which exhaustion and
  limiting pool conserved spread by construction.
- That the pool waiting just below the threshold sets the response to the
  threshold: that formula is the direct, frozen-state term, and feedback
  through the tree dominates once the critical radius is correct.
- That capacity closure treats a spread in maturity and a spread in growth
  rate differently. The difference first came from a defect that let a tip
  make more units than it had progenitors; it disappeared when that was fixed
  and reappears in the current version, where the final count tracks the
  nucleation rate. A result that changes with versions of the model is not
  claimed in either direction.
- A branching exponent for human liver vessels, whose radius reads the minor
  axis of flattened ex vivo vessels.
- Glomerular density and spacing in human kidneys, which need a validated
  detector.

## Limitations

The model's parameters are chosen, not fitted; the one calibration is the
mouse contraction. Results come from one domain scale, so the model's
contingent spread cannot stand in for a human organ's. The spatial constraint
is the local exclusion between neighbouring units, not the packing limit of
the organ, which never binds: units occupy 3 to 17 % of close packing across a
15-fold range of domain volume, and unit count scales as domain volume to the
0.74 under the allometric scalings of pool, window and rates, against the
exponent 1 that a conserved unit size requires. The progenitor processes converge with one, two
and four steps a day with the tree held fixed within the day; the coupled
vascular and tissue dynamics run once a day, and their convergence with a
smaller step was not tested. The two interventions are not matched in
strength, so their effects show that both constraints hold the liver, not
which holds it more. The critical radius is the floor of calibre
adaptation, so simulated trees cannot show calibres below it. The pool elasticity of the liver is not
measurable without load-driven sprouting (s.e. 1.12), so one of the thirteen
robustness variants neither supports nor contradicts the contrast. The room left at a tip is
estimated from trial placements in the final state, which measures residual
capacity, not whether growth could resume later. The
human data give shape and size but no unit counts, so model and cohort meet
only through the partition formula, which awaits independent input spreads.

## Before submission

1. Commit, then `make test` and `make paper` on that commit, so that
   `PAPER_NUMBERS.md` names the current version and a clean commit.
2. Write the manuscript from this file, `PAPER_NUMBERS.md` and `FIGURES.md`.
