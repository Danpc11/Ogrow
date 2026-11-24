# Perfusion-limited growth sets the size a vascularised tissue reaches, whether it can be restarted, and how poorly it scales

Daniel Pérez Calixto

Instituto Nacional de Medicina Genómica (INMEGEN), Mexico City, Mexico
Facultad de Ciencias, Universidad Nacional Autónoma de México, Mexico City, Mexico

Target journal: *Biofabrication* (IOP). Structure follows IOP's Introduction,
Methods, Results, Discussion. Abstract is under 300 words.

Reference model version 4.29.0 (`src/organgrow/closed.py`). Every model number
quoted below is taken from `PAPER_NUMBERS.md`, which `make paper` regenerates
with the version and commit that produced it, and is checked against the
result files by `scripts/29_check_manuscript.py`. No model number in this file
is written by hand.

---

## Abstract

Engineered tissues stop growing before they reach a useful size, and making a
construct larger yields functional units more slowly than it yields volume.
Which resource runs out, and whether anything can restart growth once it has
stopped, are not addressed by current design tools, which generate a vascular
geometry and then treat growth as a separate problem. We study a single
nucleation law on a perfused vascular tree: a functional unit forms at a tip
at a rate proportional to the competent cells held there and to the excess of
local perfusion over a setpoint, bounded locally by perfusion per unit, by
competent cells, and by room. Three tissue-like regimes follow from one binary
property of the cell source, whether the cell that forms a unit is spent in
forming it. The regimes differ in what seeding buys: doubling the seeded pool
changes the final unit count by an elasticity of +0.07 ± 0.06 when the source
renews and +0.71 ± 0.05 when it does not, while lowering the perfusion
setpoint is the stronger lever in every regime. They differ in whether arrest
can be undone: relieving perfusion or room restarts a renewing construct and
nothing restarts one whose pool is spent. After a 30 % loss of units, how the
perfusion circuit is driven decides the outcome; at a fixed inlet pressure the
construct overshoots to 1.82 times its mass, and holding the root flow at its
pre-injury value returns it to 1.01 and stops it. Scaling the domain
15-fold raises unit count only as V^0.45, and past an intermediate size the
unit count also stops being reproducible, its coefficient of variation rising
from 0.09 to 0.47. The model is not fitted to any construct and its
predictions are untested; we state the measurements that would decide each.

**Keywords:** tissue engineering, vascularisation, perfusion, computational
model, scale-up, organ growth, design rules

---

## 1. Introduction

Two problems limit the size of an engineered tissue. The first is that
constructs stop growing. The second is that making them larger does not make
them proportionally more functional: unit density falls with size, and the
reason is usually attributed to transport without being resolved into which
resource actually binds.

Model-guided design has made real progress on the geometry of the supply
network. Algorithms now generate organ-spanning vascular trees for arbitrary
shapes in minutes, and such trees have been printed and perfused [1]. That
work designs a static network for a given target shape. It does not ask what
happens while the tissue grows together with its vasculature, which is the
regime in which a construct either reaches a useful size or stops short.

The question of what stops growth has an established answer in developmental
biology, framed around stem cell reserves: a tissue regenerates if it keeps
one. We take a different starting point. We ask which locally available
resource runs out at the tips of a growing vascular tree, and we show that the
answer also determines whether growth can be restarted, what seeding more
cells buys, and how the construct behaves when scaled up.

The model contains one nucleation law and three local bounds. It contains no
term that measures the size of the whole construct. Size is an outcome, not a
setpoint, and the question "how does the tissue know when it is big enough"
does not arise in it: the tissue stops when a local resource is exhausted, and
which resource that is follows from one binary property of the cell source.

That property is whether the cell that forms a functional unit is spent in
forming it. In a liver-like source the cell is not spent, because mononuclear
hepatocytes are their own progenitors and the compartment is restored by the
cells it produces. In a kidney-like or pancreas-like source it is spent, and
the compartment is drawn down. We carry this as a single parameter and show
that three qualitatively different growth regimes follow from it, together
with the size of the pool each starts with.

We report five results and their reading for biofabrication. First, the three
regimes and what distinguishes them at arrest. Second, that what seeding more
cells buys depends on which bound is active, and that the perfusion setpoint
is the stronger lever in all three. Third, that arrest is stable, and whether
it can be undone is set by the same property. Fourth, that after injury the
boundary condition on the perfusion circuit, not anything the tissue
remembers, decides the final size. Fifth, that scaling the domain up loses
unit density and, past an intermediate size, loses reproducibility as well.

None of these predictions has been tested against a construct. We state in the
Discussion the measurements that would decide each of them.

## 2. Methods

### 2.1 The nucleation law

A functional unit is nucleated at a vascular tip $i$ at rate

$$\lambda_i = \nu\, P_i\, (q_i - q_{\mathrm{set}})_+$$

where $P_i$ is the number of competent cells held at that tip, $q_i$ is the
local perfusion per unit of tissue served by the tip, $q_{\mathrm{set}}$ is a
threshold below which no unit forms, $\nu$ is a rate constant, and
$(\cdot)_+$ denotes the positive part. The law is local: no term in it refers
to the state of the construct as a whole.

Nucleation is bounded at each tip by three quantities, each of which can run
out independently: the perfusion available per unit, the competent cells
present, and the room to place a unit at least a minimum spacing
$\lambda_{\min}$ from its neighbours.

### 2.2 The cell source

One binary parameter, `renewal` ($w$), states whether the cell that forms a
unit is spent. With $w = 1$ the compartment is restored by the cells it
produces and forming a unit does not draw it down. With $w = 0$ the cell is
consumed. Three regimes are defined by this parameter together with the
initial pool $P_0$ and a loss rate $\delta$:

| regime | $w$ | $P_0$ | $\delta$ |
|---|---|---|---|
| liver-like, renewing source | 1 | 2600 | 0.10 |
| kidney-like, spent source with loss | 0 | 300 | 0.05 |
| pancreas-like, spent source without loss | 0 | 60 | 0.00 |

No other parameter differs between the three. All runs use
$\nu = 8$, $q_{\mathrm{set}} = 0.6$, and the same vascular and geometric
settings.

### 2.3 The vascular tree and its boundary condition

The tree grows by sprouting at tips and adapts its calibres under a wall-cost
objective $\mathcal{C} = Q^2/r^4 + k\,(p\,r^2)^b$ minimised per unit length,
with wall cross-section proportional to $p\,r^2$ by Laplace's law. This gives
the branching relation $r_p^{\,b+2} = \sum r_c^{\,b+2}$ and wall shear
proportional to $p^{b/2} r^{\,b-1}$. Shear on a segment in series with its
tissue peaks at $r_c = (\ell / 3R_t)^{1/4}$, the radius below which a vessel
regresses.

The network is linear in the driving pressure, so two boundary conditions are
available and we use both. Under a **fixed inlet pressure**, the standard
condition, the tree grows its conductance as the tissue demands more, and
perfusion per unit need not fall as the construct grows. Under a **fixed root
flow**, the whole solve is rescaled so that the flow entering the root equals a
chosen value; each unit's share is then its own conductance times the pressure
its tip delivers,

$$w_j \propto m_j^{\,\xi}\, q_{\mathrm{tip}(j)},$$

normalised over units, so delivery per unit mass falls as the construct grows.
The second condition is what a circulation imposes on a tissue that is part of
a body, and it is what a perfusion bioreactor imposes when it is run at
constant flow rather than constant pressure.

### 2.4 Readouts, replicates and robustness

Growth runs for 140 days, by which point all three regimes have stopped. At
arrest we report the final unit count and the disposition of the initial pool,
separated into the share consumed into units, the share lost to attrition, the
share held at tips short of supply, the share held for lack of room, and the
share still active. Room at a tip is measured by trial placements in the final
state, independently of the dynamics.

Every reported quantity is the mean over 24 independent seeds, with the
standard error or the range the seeds span. Contrasts between regimes are
paired seed by seed.

Auxiliary rules were varied in 13 ways that do not touch the nucleation law:
root placement at the surface or at random, removal of deficit-driven
sprouting, of load-driven sprouting, of crowding, and of the mesenchymal term,
territory-based rather than nearest-tip assignment, binomial nucleation at
one, two and four substeps per day, competing fates, and removal of
self-renewal.

### 2.5 Developmental data

Four public single-cell series were used to ask whether the progenitor pool
each tissue is named for is present or absent at a given stage, relative to
the window over which that tissue is still forming units: human fetal kidney
(GSE114530, weeks 9, 11, 13, 16, 18; and GSE139280, week 14, as an independent
preparation), mouse pancreas (GSE132188, E12.5 to E15.5) and mouse liver
(GSE90047, E10.5 to E17.5, DLK1-sorted). Cell states were assigned by marker
score, z-scored within each sample, with a two-component Gaussian mixture used
to separate states. Only presence and absence are read across studies, for
reasons given in Section 3.6.

### 2.6 Code and data availability

All code, parameters and generated tables are at
`https://github.com/Danpc11/Ogrow`. `make paper` regenerates every number in
`PAPER_NUMBERS.md` from scratch and records the model version and the commit
that produced it. `make test` runs the conservation and consistency tests,
including exact progenitor balance in every realisation. The developmental
series are public under the accessions above.

## 3. Results

### 3.1 One property of the cell source separates three growth regimes

With the single difference of Section 2.2, the three sources arrest in three
different states (`PAPER_NUMBERS.md`, Regimes at closure; Figure 1).

The renewing source reaches 208 ± 4 units and keeps essentially its whole
pool, held partly where perfusion has fallen to threshold (0.78 ± 0.03 of the
pool) and partly where existing units leave no room (0.19 ± 0.02). It is
limited by its vasculature and its geometry, not by cells.

The spent source without loss reaches 45 ± 1 units, consumes 0.75 ± 0.02 of
its pool into units, and strands the remaining 0.25 ± 0.02 at tips short of
supply. It is limited by cells, with a remainder it cannot use.

The spent source with loss reaches 89 ± 2 units and ends holding nothing at
all. Its emptiness is not what the word "spent" suggests: of the pool it
starts with, 0.30 ± 0.01 is consumed into units and 0.70 ± 0.01 is lost to
attrition. Both leave the compartment and only the first built anything. More
than twice as much of the pool is lost as is used.

That last figure is a property of the model and not a claim about a kidney. It
is also, read as engineering, a statement of a failure mode: seeded cells that
end up at tips the vasculature does not reach never build anything and are
simply lost. We return to it in Section 4.3.

The ordering of the three regimes survives all 13 variants of the auxiliary
rules (`PAPER_NUMBERS.md`, Robustness).

### 3.2 What seeding more cells buys depends on which bound is active

The elasticity of the final unit count to the initial pool, measured by paired
runs at $P_0 e^{\pm 0.3}$ over 24 seeds, separates the regimes sharply
(`PAPER_NUMBERS.md`, Elasticities; Figure 2):

| source | elasticity to the seeded pool | elasticity to the perfusion setpoint |
|---|---|---|
| renewing | +0.07 ± 0.06 | −1.83 ± 0.22 |
| spent, with loss | +0.37 ± 0.03 | −1.65 ± 0.06 |
| spent, without loss | +0.71 ± 0.05 | −0.77 ± 0.07 |

Two readings follow. First, doubling the number of cells seeded into a
renewing construct changes its final unit count by a few per cent, while
doubling it into a spent-source construct changes it by most of a factor of
two. Which of those a given construct does is set by one property of its cell
source, and it can be determined before any cells are seeded.

Second, and in every regime, the perfusion setpoint is the stronger lever. Its
elasticity is negative and of order one to two, against a pool elasticity
between 0.07 and 0.71. Lowering the local perfusion a unit requires buys more
than seeding more cells, in all three regimes, by a margin of between two-fold
and twenty-fold.

The ordering of pool elasticities, renewing below spent-with-loss below
spent-without-loss, holds in 13 of 13 variants of the auxiliary rules, each a
sign test at $P = 0.0002$, and the separations reach $z \ge 2$ in 12 of the 13
(`PAPER_NUMBERS.md`, The ordering counted over the variants). The one variant
in which it is not resolved is the one without load-driven sprouting, where
the renewing source's elasticity has a standard error of the order of its own
range.

This ordering is also the pattern reported in mice by Stanger, Tanaka and
Melton [2], who found that halving the pancreatic progenitor pool halved the
pancreas, which never compensated, while the liver compensated fully. The
model was not fitted to that experiment.

### 3.3 Arrest is stable, and whether it can be undone follows from the same property

At day 140 growth has stopped. Letting it run 40 days further adds 2.8 ± 0.3 %
more units to the renewing construct and nothing to either spent source, so
arrest is a state and not a slowdown (`PAPER_NUMBERS.md`, Interventions;
Figure 3).

Arrest can then be probed by relieving one bound at a time. Lowering the
perfusion setpoint by 20 %, which raises every tip's availability relative to
threshold without touching the pool, adds 22.6 ± 2.7 % to the renewing
construct. Shortening the minimum spacing by 20 %, which frees room without
touching perfusion or the pool, adds 32.3 ± 1.7 %. Both bounds hold it at
once, and either relief restarts it.

Neither relief does anything to the spent-with-loss construct, which closed
holding nothing: it has no pool for a relieved bound to release, and in 24
seeds its pool on day 140 is 0.1 cells out of 300 seeded. The spent-without-loss
construct, which closed holding a quarter of its budget, gains 6.2 ± 1.4 % from
the perfusion relief and 0.3 ± 0.2 % from the room relief, which is what a
construct stranded for want of supply rather than space should do.

The engineering statement is that a stalled construct is diagnosable. If
relieving perfusion restarts it, perfusion was binding; if relieving spacing
restarts it, room was; if neither does, the cell source is exhausted and no
change to the circuit or the scaffold will help.

### 3.4 After injury, how the circuit is driven decides the final size

A construct grown to arrest, then stripped of 30 % of its units by a plane
cut and followed for 80 more days, does not recover on its own
(`PAPER_NUMBERS.md`, Injury; Figure 4). Left as grown, the renewing construct
returns to 0.95 [0.77, 1.01] of its unit count but only 0.79 [0.67, 0.86] of
its mass, and its perfused mass tracks the mass at 0.81 [0.67, 0.90]. The
reason is internal to the model: the rate at which a unit grows decays through
development, so by arrest a unit can no longer enlarge, and a construct that
rebuilds after arrest fills with units that never mature.

Recovery therefore needs the ability of a unit to enlarge to be restored, which
is what the growth factors of regeneration do. With growth reactivated, both
constructs recover their mass. They also overshoot it badly: at the model's
default boundary condition, a fixed inlet pressure, the renewing construct
reaches 1.82 [1.54, 2.07] times its pre-injury mass and the spent-source
construct 1.85 [1.53, 2.21], because the tree grows its conductance as the
tissue demands more and perfusion per unit never falls.

Holding the root flow instead at its pre-injury value changes the outcome
qualitatively. Each unit's share is then its own conductance times the pressure
its tip delivers, so delivery per unit mass falls as the construct grows. The
renewing construct restores its mass to 1.01 [0.97, 1.04] and stops. The
spent-source construct reaches 0.97 [0.94, 0.99]. The rule stays local
throughout: a unit grows while its own delivery per mass exceeds a threshold,
and the ceiling on the whole construct follows from that, with nothing
measuring the whole.

This is the sharpest prediction in the paper and the cheapest to test. A
perfusion bioreactor run at constant pressure and one run at constant flow
should give different final sizes from identical constructs, overshoot in the
first case and self-limitation in the second.

Holding the inflow also makes the outcome more repeatable. Over 24 seeds the
coefficient of variation of final mass falls from 6.3 % to 2.1 % in the
renewing construct and from 8.2 % to 1.3 % in the spent-source construct, while
the spread of unit count is left essentially untouched, 5.9 % to 7.0 % and
0.5 % to 0.5 % (`PAPER_NUMBERS.md`, Spread across seeds). A fixed supply budget
fixes the mass and says nothing about how a construct divides it.

### 3.5 A construct that cannot rebuild units enlarges the ones it kept

The spent-source construct arrested with no competent cells, so it cannot make
new units and can only enlarge those it has. After losing 30 % of its units it
ends with 0.70 [0.69, 0.70] of the count and 1.39 [1.35, 1.43] times the mass
per unit, against the 1.43 that full compensation would require, and its mass
returns to 0.97 [0.94, 0.99]. After a 50 % loss it reaches 1.88 [1.79, 1.95]
against the 2.00 required. The model therefore compensates almost completely by
enlargement alone.

Nothing in the model limits how large a unit may become. In a real tissue that
limit exists, and this is where the model and the biology part company
(Section 5).

### 3.6 The one property, checked against three developing tissues

The model's central claim is about the cell source, so the one check available
without new experiments is whether the progenitor pool each tissue is named for
is present or absent while that tissue is still forming units
(`PAPER_NUMBERS.md`; Figure 5).

All three tissues deplete the progenitor they are named for, which on its own
settles nothing, since differentiation empties a named compartment whether or
not that compartment is what stops the tissue. What separates them is when,
measured against the window over which units are still being formed. In mouse
pancreas the multipotent pool is at 0.0001 of the lineage at E15.5, the end of
the window in which that organ specifies its units: 100 % of the window
elapsed. In human fetal kidney the progenitors are still at 0.1689 of the
lineage at week 18, 43 % of the way through a window that runs to weeks 34 to
36. In mouse liver the hepatoblast programme is at 0.0000 by E17.5, 24 % of the
way through the window over which hepatocyte number is still rising, and the
liver goes on to gain more than fivefold in mass by P20 by dividing cells it
has already made.

Two tissues are limited by the pool they deplete and one is not. That is what
the single property predicts, and it is what counting progenitors alone would
miss.

Only presence and absence are read across studies, and the reason is
measurable. An independent week-14 kidney has 0.662 of its nephron lineage in
progenitors, against 0.159 to 0.332 at weeks 13 to 16 in the first series. That
spread at the same stage is not closed by batch correction: centring each gene
within its sample, which is the location adjustment a batch correction makes,
leaves the two studies apart, because an expression-level correction aligns
shared states and does not change which cells a preparation contains. A share
is a property of the preparation as much as of the tissue. An absent
compartment is not, since no protocol recovers one that is not there.

The law's other prediction cannot be decided with these data. While the supply
term stays open, induction should be proportional to the pool, so the number
induced per progenitor still held should hold steady as the pool falls. Across
the five kidney samples that ratio does not order with stage (Spearman
$\rho = +0.30$, $P = 0.62$), and neither does the progenitor share
($\rho = -0.50$, $P = 0.39$). With one sample per stage, a change over
development and a difference between samples cannot be told apart. Deciding it
needs a series with replicate donors per stage.

### 3.7 Scaling up loses unit density, and past an intermediate size loses reproducibility

Scaling the domain over a 15-fold range of volume, with the pool, the window
and the rates held at the reference, raises the unit count only as
$V^{0.45}$. Carrying the standard allometric scalings of pool, window and rates
as well raises it as $V^{0.58}$ (`PAPER_NUMBERS.md`, Organ size; Figure 6).
Both are far below the $V^{1}$ that a conserved unit size would require, so
unit density falls as the construct grows. Doubling the volume buys about 1.4
times the units.

Which bound is active also changes with size. The share of the pool held for
lack of room rises from 0.09 at the smallest domain to 0.37 at the largest,
while the share held by perfusion falls from 0.89 to 0.63. That reads as
geometry starting to bind, and it is not: the units occupy 16.2 % of close
packing at the smallest domain and 4.2 % at the largest, so the construct is
emptier at size, not fuller. What rises is local crowding at the vascular
front, where a unit must be placed within a fixed distance of a tip. Scaling up
is limited by how fast the tree spreads, not by how much room the container
has.

The more consequential result is that the unit count stops being reproducible.
Its coefficient of variation over 24 seeds is 0.09 at the two smallest domains,
then rises to 0.23 and 0.47; at the largest domain the seeds span 274 to 1074
units under identical rules and parameters. The onset lies between the two
middle scales, and it coincides with the crossover at which room overtakes
perfusion as the binding constraint.

We state the consequence plainly rather than working around it. The exponents
above are fitted across points whose spread at the large end is nearly
fourfold, so they describe a trend and are not precise constants. And if the
effect is not an artefact of the model's discretisation, it says something a
construct designer needs to know: under a purely local placement rule, two
identical large builds need not give the same number of functional units. Below
the crossover, they do.

The exponent is measured off the sweep rather than derived, and it is stated at
one nucleation rate. Whether a rate matched to the speed of the vascular front
would carry it towards 1 is the question `scripts/25_front.py` asks; we do not
claim its answer here.

## 4. Discussion

### 4.1 Four design rules, stated as the model gives them

**Drive the perfusion circuit at constant flow, not constant pressure, when
the construct should stop at a target size.** At a fixed inlet pressure the
model's construct overshoots to 1.82 times its mass after injury; at fixed root
flow it returns to 1.01 and stops. Nothing in the tissue measures its own size
in either case. The brake is the supply the circuit provides. Constant flow also
tightens the spread of the final mass by a factor of three to six.

**Decide whether your cell source renews before deciding how many cells to
seed.** If it renews, seeding more barely changes the final unit count
(elasticity +0.07 ± 0.06) and the vasculature sets the size. If it is spent in
forming units, the seeded pool sets the size directly (+0.71 ± 0.05) and there
is no second chance once it is gone.

**Lower the perfusion a unit requires before seeding more cells.** The
setpoint elasticity is between −0.77 and −1.83 across the three regimes,
against a pool elasticity between +0.07 and +0.71. In every regime the
threshold is the stronger lever, and in the renewing regime it is the only one
that matters.

**Match the seeding rate to the speed of the vascular front, rather than
maximising it.** This one is a hypothesis the model raises and does not yet
answer. Its basis is that what strands a progenitor is sitting at a tip the
tree has not reached, that 0.70 of the pool is lost that way in the
spent-with-loss regime, and that scaling up is limited by the spread of the
front rather than by room in the container. If that reading is right, the
scaling exponent is a ratio that can be tuned, not a constant.

### 4.2 Relation to model-guided vascular design

Existing model-guided design generates a vascular geometry for a target shape
and validates it by printing and perfusing it [1]. The two approaches are
complementary rather than competing. A generated tree answers where the vessels
should go; the present model answers what happens to the tissue while it grows
with them, and says that even a well-designed tree does not by itself fix the
number of units that get built. The rate at which units are seeded relative to
the advance of the front, and the boundary condition on the circuit, are
separate controls with separate effects.

### 4.3 What the model cannot do

The model cannot spend a budget to the last cell, and this is structural
rather than a matter of parameter choice. A progenitor is tied to the tip it
sits at, so where supply there falls short it strands. With no attrition, no
initial pool in the range tested leaves less than a third of the pool held, and
a larger pool strands more. Two repairs were measured and both fail for stated
reasons. An attrition term empties the pool without building anything, which is
what the spent-with-loss regime does now. A residual drive below the threshold
does empty a budget, but it takes the renewing construct to nearly five times
its unit count and halves the share of its mass the circulation reaches;
making the drive source-specific would add a second property differing exactly
as `renewal` does, which is the claim under test. The variant is kept as
`final_induction`.

For a construct this is not only a defect. It names a failure mode: cells
placed where perfusion will not arrive are lost, and the question for a seeding
protocol is not how many cells to use but how to keep them where the front will
reach them.

### 4.4 What would test this

Each of the following would decide one result. None has been done.

- **Constant flow against constant pressure.** Two identical perfused
  constructs, one driven at constant inlet pressure and one at constant flow,
  followed past the point where growth would arrest. The first should overshoot
  and the second should stop at a reproducible mass with a tighter
  construct-to-construct spread.
- **Count against function.** Unit number and perfused mass measured in the
  same construct as it rebuilds after arrest. The count should return first and
  the perfused mass lag or never arrive. A construct whose count and function
  recover together would contradict the model.
- **Seeding against source.** Final unit count against seeded pool, in one
  source that renews and one that is spent. The elasticities above are the
  prediction, and they differ by an order of magnitude.
- **Diagnosing a stalled construct.** Relieving perfusion and relieving
  spacing, separately, on a construct that has stopped. Which relief restarts
  it identifies the binding constraint; neither restarting identifies an
  exhausted source.
- **Reproducibility at scale.** The same build repeated at a small and at a
  large domain. The model predicts that the spread of unit count is small at
  small size and large at large size, with a crossover.
- **Seeding rate against front speed.** Unit density at fixed size, across a
  range of seeding or induction rates. The model predicts a rate that maximises
  it rather than monotone improvement.

## 5. Limitations

The parameters are chosen, not fitted; the single calibration is the
contraction of the spread of unit number in mouse liver. No prediction in this
paper has been tested against a construct or a tissue.

The spatial constraint is local exclusion between neighbouring units, not the
packing limit of the container, which never binds.

The model compensates after injury further than a kidney does, reaching 1.39
times the mass per unit after a 30 % loss against the 1.43 of full
compensation, whereas a kidney that loses half its mass grows the remnant by
about a half rather than doubling it. Nothing in the model limits how large a
unit may become; in a tissue that limit exists and passing it is injury.

The two interventions of Section 3.3 are not matched in strength, so they show
that both bounds hold the renewing construct, not which holds it more. Room at
a tip is estimated from trial placements in the final state, which measures
residual capacity rather than whether growth could resume.

The progenitor processes converge with one, two and four steps per day with the
tree held fixed within the day; the coupled vascular and tissue dynamics run
once per day and their convergence at a smaller step was not tested. The
critical radius is the floor of calibre adaptation, so simulated trees cannot
show calibres below it.

The developmental series carry one sample per stage, so they support presence
and absence and not the magnitude of any ratio.

## 6. Conclusion

The size a perfused tissue reaches, whether its arrest can be undone, and what
seeding more cells buys all follow from one binary property of the cell source
and from the boundary condition on the perfusion circuit. Neither is a
property of the scaffold, and neither requires anything in the tissue to
measure its own size. Scaling the construct up costs unit density and, past an
intermediate size, costs reproducibility as well. These are predictions from a
model that has not been fitted to a construct, and each has a measurement that
would decide it.

## Acknowledgements

To be completed: funding (SECIHTI), computing (INMEGEN cluster), and author
contributions.

## References

[1] Sexton Z I et al. 2023 Rapid model-guided design of organ-scale synthetic
vasculature for biomanufacturing *Science* **382** eadj6152

[2] Stanger B Z, Tanaka A J and Melton D A 2007 Organ size is limited by the
number of embryonic progenitor cells in the pancreas but not the liver
*Nature* **445** 886–91

[3] Murray C D 1926 The physiological principle of minimum work. I. The
vascular system and the cost of blood volume *Proc. Natl Acad. Sci. USA*
**12** 207–14

<!--
REFERENCES STILL NEEDED. These claims are in the text and currently rest on
the repository's notes rather than on a citation I have verified. Each must be
given a source before submission, or the claim removed:

  a. Azizoglu et al. 2025, the mouse liver CV at start and at closure
     (0.229 and 0.074). Used as the single calibration.
  b. The human nephrogenesis window, weeks 5 to 35 or 36, ending with the
     progenitors exhausted.
  c. The mouse pancreatic multipotent window, E9.5 to E15.5.
  d. The mouse liver window over which hepatocyte number is still rising, to
     P20, with the more than fivefold gain in mass.
  e. Mononuclear hepatocytes acting as their own progenitors.
  f. Compensatory renal hypertrophy after nephrectomy: no new nephrons, the
     remaining ones enlarge, function maintained, and the figure that a kidney
     losing half its mass grows the remnant by about a half.
  g. The four GEO accessions' original publications.
-->

## Figures

Legends are in `FIGURES.md` and are regenerated with the figures by
`make paper`. File names in the repository are historical; the mapping is:

| here | file | script |
|---|---|---|
| Figure 1 | `Fig1_closed` | `figures/fig_framework.py` |
| Figure 2 | `Fig5_pool` | `figures/fig_framework.py` |
| Figure 3 | `Fig3_interventions` | `figures/fig_design.py` |
| Figure 4 | `Fig8_injury` | `figures/fig_injury.py` |
| Figure 5 | `Fig9_transcriptome` | `figures/fig_transcriptome.py` |
| Figure 6 | `Fig6_scaling` | `figures/fig_design.py` |

`FIGURES.md` also lists the figures the repository builds that this manuscript
does not use. They belong to the earlier framing kept in `PAPER.md` and are
available as supplementary material.
