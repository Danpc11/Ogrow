# Perfusion-limited growth sets the size a vascularised tissue reaches, whether it can be restarted, and how poorly it scales

Hugo Tovar^1^, Genaro Vázquez-Victorio^2^, Miguel Ángel Peña-Rico^3^,
Angélica S. Martínez-Ramírez^3^, Itzel Escobedo-Avila^4^,
Enrique Hernández-Lemus^1^, Daniel Pérez-Calixto^1,2^\*

^1^ Instituto Nacional de Medicina Genómica (INMEGEN), Mexico City, Mexico

^2^ Facultad de Ciencias, Universidad Nacional Autónoma de México, Mexico City,
Mexico

^3^ Centro de Investigaciones Científicas, Instituto de Biotecnología,
Universidad del Papaloapan, Oaxaca, Mexico

^4^ Instituto de Fisiología Celular, Universidad Nacional Autónoma de México,
Ciudad Universitaria, Mexico City, Mexico

\* Corresponding author: Daniel Pérez-Calixto, dperez@inmegen.gob.mx

Target journal: *Biofabrication* (IOP). Structure follows IOP's Introduction,
Methods, Results, Discussion. Abstract is under 300 words.

Reference model version 4.55.0 (`src/organgrow/closed.py`). Every model number
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
competent cells, and by room. Three tissue-like regimes follow from how the cell
source is maintained: whether the cell that forms a unit is spent in forming
it, and, where it is spent, how fast what remains is lost. The regimes differ in what seeding buys: doubling the seeded pool
changes the final unit count by an elasticity of +0.07 ± 0.06 when the source
renews and +0.71 ± 0.05 when it does not, while lowering the perfusion
setpoint is the stronger lever in every regime. They differ in whether arrest
can be undone: relieving perfusion or room restarts a renewing construct and
nothing restarts one whose pool is spent. After a 30 % loss of units, how the
perfusion circuit is driven decides the outcome; at a fixed inlet pressure the
construct overshoots to 1.82 times its mass without settling, while one held
at a fixed root flow tracks a supply budget in the aggregate, its mass rising
as the delivered flow to an exponent of 0.90 while most of its units stop
short of their own share of the supply. That gap is the paper's sharpest
result: a supply constraint fixes how much tissue a construct carries and
leaves what the tissue is made of undetermined. Giving the model the
involution it lacks returns the same mass on roughly a third of the units, so
mass, unit number and function come apart, and a construct optimised to
recover mass need not recover the architecture it was designed around.
Scaling the domain
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
which resource that is follows from how the cell source is maintained.

That property is whether the cell that forms a functional unit is spent in
forming it. In a liver-like source the cell is not spent, because mononuclear
hepatocytes are their own progenitors [4] and the compartment is restored by
the cells it produces. In a kidney-like or pancreas-like source it is spent, and
the compartment is drawn down. We carry this as a single parameter and show
that three qualitatively different growth regimes follow from it together with
the size of the pool each starts with and the rates at which that pool divides
and is lost. Section 3.7 separates the three contributions and finds the
turnover the largest of them.

We report six results and their reading for biofabrication. First, the three
regimes and what distinguishes them at arrest. Second, that what seeding more
cells buys depends on which bound is active, and that the perfusion setpoint
is the stronger lever in all three. Third, that arrest is stable, and whether
it can be undone is set by what the pool has left, which the turnover decides
more often than the binary property does. Fourth, that after injury the
boundary condition on the perfusion circuit, not anything the tissue
remembers, decides the final size. Fifth, that scaling the domain up loses
unit density and, past an intermediate size, loses reproducibility as well.
Sixth, that a joint specification on density, vascular cost and reproducibility
is met up to an intermediate volume and not beyond it, and that what fails is
density, while reproducibility at size is recovered by the choice of cell
source.

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
$\lambda_{\min}$ from its neighbours. The law, the three bounds and the one
property of Section 2.2 are drawn in Figure 1a to 1c.

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
in the sense of Murray [3], with wall cross-section proportional to $p\,r^2$
by Laplace's law. This gives the branching relation
$r_p^{\,b+2} = \sum r_c^{\,b+2}$ and wall shear
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

Under it a unit's ceiling is $m_j^{*} = (Q/\phi)\,w_j/\sum_k w_k$, where $Q$
is the delivered flow and $\phi$ the delivery per unit mass the tissue
requires. The shares sum to one, so if every unit sits at its ceiling then

$$\sum_j m_j^{*} = Q/\phi .$$

A ceiling on the whole construct therefore follows from a rule that refers
only to a tip. This is a consequence of the rule and not a measurement, and we
report it as one: setting $\phi$ to a construct's own pre-injury delivery per
unit mass forces that construct back to its pre-injury mass by arithmetic, for
every seed, which is why Section 3.4 sets $\phi$ once on an intact ensemble
and then varies $Q$. With $\xi = 1$ a unit's ceiling moves with its own mass,
so the identity above is approached rather than reached, and how closely is
measured in Section 3.4.

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
(GSE114530 [6], weeks 9, 11, 13, 16, 18; and GSE139280 [7], week 14, as an
independent preparation), mouse pancreas (GSE132188 [8], E12.5 to E15.5) and
mouse liver (GSE90047 [9], E10.5 to E17.5, DLK1-sorted). Cell states were assigned by marker
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

### 3.1 The three cell sources arrest in three different states

The three sources arrest in three different states (`PAPER_NUMBERS.md`,
Regimes at closure; Figure 1). They differ in three parameters at once, and
Section 3.7 separates them: the binary property is one of the three and the
weakest, while how the pool turns over is what most often decides which bound
ends the growth.

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
two. Which of those a given construct does is set by its cell source and by the
pool it starts with, both of which are fixed before any cells are seeded;
Section 3.7 reports how much of the difference each contributes.

Second, and in every regime, the perfusion setpoint is the stronger lever. Its
elasticity is negative and of order one to two, against a pool elasticity
between 0.07 and 0.71. Lowering the local perfusion a unit requires buys more
than seeding more cells, in all three regimes, by a margin of between two-fold
and twenty-fold.

The ordering of pool elasticities, renewing below spent-with-loss below
spent-without-loss, holds in all 13 variants of the auxiliary rules, and the
separation reaches $z \ge 2$ in 12 of them (`PAPER_NUMBERS.md`, The ordering
counted over the variants). The contrast within each variant is paired seed by
seed and its $z$ is a test; the count across variants is not, because the
variants share a model, a parameter set and the same 24 seeds. It says that
the ordering does not depend on which auxiliary rule is in force, not that it
has been confirmed 13 times over.

The variant in which it is not resolved is the one without load-driven
sprouting, where the three contrasts fall to $z$ of 1.5, 0.6 and 1.8. That
rule is therefore not auxiliary: sprouting driven by load carries part of the
separation between regimes, and removing it costs both the size of the
contrast and its precision.

This ordering is also the pattern reported in mice by Stanger, Tanaka and
Melton [2], who found that halving the pancreatic progenitor pool halved the
pancreas, which never compensated, while the liver compensated fully. The
model was not fitted to that experiment.

### 3.3 Arrest is stable, and whether it can be undone follows from what the pool has left

We call a construct arrested when letting it run a further 40 days, nearly a
third again of its growth, adds less than 5 % to its unit count. That is an
operational criterion and not a claim that no cell could ever divide again:
what it rules out is a construct still climbing slowly enough to be mistaken
for a stopped one. By it, all three regimes are arrested at day 140. The
renewing construct adds 2.8 ± 0.3 % more units over those 40 days and either
spent source adds nothing, so for the spent sources the arrest is also
terminal, the pool having no competent cells left to spend
(`PAPER_NUMBERS.md`, Interventions; Figure 2). The distinction matters for what
follows: a bound that still holds a pool can be relieved, and an exhausted
pool cannot.

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
(`PAPER_NUMBERS.md`, Injury; Figure 3). Left as grown, the renewing construct
returns to 0.95 [0.77, 1.01] of its unit count but only 0.79 [0.67, 0.86] of
its mass, and its perfusion-weighted tissue mass tracks the mass at
0.81 [0.67, 0.90]. The
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

Holding the root flow fixed instead changes the outcome qualitatively, and
the reason is the budget of Section 2.3. Each unit's share is its own
conductance times the pressure its tip delivers, the shares sum to one, so if
every unit sits at its ceiling the whole construct sits at $Q/\phi$. A ceiling
on the construct follows from a rule that refers only to a tip, with nothing
measuring the construct. That is the claim; the value of the ceiling is a
consequence of the rule and not an independent result, and we report it as
such.

Two things follow that are not consequences of the rule. The first is what the
construct actually reaches, since with $\xi = 1$ a unit's ceiling moves with
its own mass and the approach is partial. We calibrate $\phi$ once on an
intact ensemble and evaluate on a disjoint block of seeds, so the demand is
not fitted on the constructs it is applied to, and then vary the delivered
flow between 0.6 and 1.4 times its pre-injury value
(`scripts/32_inflow_law.py`; `PAPER_NUMBERS.md`, The mass a construct lands
on). Mass rises as the delivered flow to an exponent of 0.90, fitted once per
seed. Against the budget itself the construct overshoots when asked to shrink,
reaching about 1.8 times $Q/\phi$ at the lowest flow, and falls short when
asked to grow, reaching about 0.7 of it at the highest.

By day 220 the mass is changing by well under a tenth of a per cent per day.
That is necessary for equilibrium with the supply and it is not sufficient,
and in this model it is not what has happened. The rate at which a unit grows
decays as $e^{-t/\tau_g}$, so the growth available after the injury is finite,
$g_0\tau_g$, and a unit can stop moving while still far from its ceiling.
Summing that distance over the units, $\sum_j |m_j - m_j^*| / \sum_j m_j^*$
reaches as much as 0.38 for the renewing construct, which ends with up to 74 %
of its units at the mass floor, and as much as 0.24 for the spent source, with
up to 49 % at the floor (`PAPER_NUMBERS.md`, The mass a construct lands on).
Neither improves at any flow in the sweep, so this is the regime and not an
edge of it. The total lands where the budget puts it, because the shares sum
to one; the units individually do not reach their ceilings.

What separates the two is a mechanism the model does not have. A unit whose
share of the delivery collapses has no way to disappear, so it stays at the
mass floor rather than involuting, and the construct ends as a few large units
among many vestigial ones.

Supplying that mechanism, as a rate at which a unit whose supply has fallen
below a tenth of the setpoint loses mass and is then removed, changes what the
construct is made of and not how much of it there is
(`scripts/33_involution.py`; `PAPER_NUMBERS.md`, Involution, as a
sensitivity). At a rate of 0.05 per day the renewing construct still returns to
0.97 of its pre-injury mass, against 1.02 without the mechanism, but it does so
with 0.37 of its pre-injury unit count rather than 0.91: the same mass carried
by fewer and larger units, which is the behaviour Section 3.5 reports for the
construct that cannot make new ones. The spent source goes from 0.70 of its
units to 0.52 at an unchanged 0.97 of its mass. The share left at the mass
floor falls from 0.72 to 0.13 in the renewing construct and from 0.40 to 0.04
in the spent source. Over a tenfold range of the rate, 0.02 to 0.20 per day,
the renewing construct keeps between 0.35 and 0.40 of its units and the spent
source between 0.51 and 0.53. What that supports is narrow and worth stating
narrowly: the final unit fraction is insensitive to the rate across that range
at this horizon, which is 80 days after the injury. The same endpoint can be
reached along different trajectories, and we have not measured them.

What the unit selected for removal is starving of also needs saying, because
the criterion and the growth rule could have been reading different
quantities. They are not. A unit's delivery per unit mass under a fixed inflow
is $Q s_j / m_j$, and with $s_j \propto m_j^{\xi} q_j$ at $\xi = 1$ the mass
cancels exactly: $Q s_j / m_j = Q q_j / \sum_k m_k q_k$. The delivery per unit
mass is therefore the tip supply $q_j$ times a quantity common to every unit
in the construct, so a threshold on either selects the same units in the same
order and the two differ only in where the cut falls. The criterion is the
supply the circulation actually provides, not a proxy for it.

Where the cut falls is a separate question, because $Q / \sum_k m_k q_k$ is a
property of the construct and could drift as it rebuilds or as the delivered
flow is changed, so that a fixed threshold on $q_j$ would not be a fixed
requirement in delivery per unit mass. Measured across three delivered flows
and three seeds at day 220, it spans a factor of 1.04. The threshold sweep
covers a factor of ten and moves the result very little, so a four per cent
drift in the conversion between the two decides nothing here. The argument
does not carry to $\xi \neq 1$, where the mass no longer cancels and the two
criteria would have to be compared again.

The residual falls with it, to 0.23 in the renewing construct and to 0.18 in
the spent source, and does not approach zero. Clearing the stranded units
does not give a surviving unit
back the growth it would need to reach its ceiling, so these are two separate
defects and involution addresses one of them. We keep it out of the result
above because its rate is calibrated against nothing, and what this section
reports should not rest on a number we chose.

All of this is structure, and none of it is function. A construct that carries
0.97 of its pre-injury mass on 0.37 of its units has recovered its mass and
has not been shown to recover anything a liver does. The model has no measure
of function at all (Section 5), so what is reported here is structural
recovery and perfusion-weighted mass, and the distance between those and
function is the point rather than an aside: mass, unit number and function are
three quantities, and a supply constraint fixes only the first.

A construct fed at constant flow therefore tracks its supply in the aggregate,
and does so asymmetrically: it grows into a larger budget far more readily
than it gives up mass to a smaller one. The spent-source construct behaves the same way, at
its own $\phi$ and its own exponent, so the asymmetry is a property of the
rule and not of either cell source. Calibrating the demand out of sample barely moves where the construct lands
at unit flow and roughly doubles the range the seeds span. That is the
clearest statement of what the earlier figure was: not the mean, which the
rule gives anyway, but the tightness, which came from handing every construct
its own pre-injury state as its target.

The second is the contrast with constant pressure, which carries no reference
to the pre-injury state at all. There the ceiling is the supply at the unit's
tip with no budget, the tree grows its conductance as the tissue demands more,
and the construct overshoots to 1.82 times its mass without settling.

For completeness, with $\phi$ instead taken from each construct's own
pre-injury delivery per unit mass, which is what `scripts/24_injury.py` does,
the renewing construct returns to 1.00 [0.95, 1.03] and the spent-source
construct to 0.97 [0.93, 0.99] (`PAPER_NUMBERS.md`, Injury). Those numbers,
and the contraction of the spread of final mass that comes with them from
6.3 % to 2.3 % and from 8.2 % to 1.4 % (`PAPER_NUMBERS.md`, Spread across
seeds), are that choice of $\phi$ written out. They are not evidence that the
construct regulates its own size.

The experiment this asks for is a bioreactor run at several constant flows.
The model predicts that final mass rises with the delivered flow, close to
proportionally, and that a construct held at constant pressure does not settle
at all. Both are cheap to measure and either would falsify the rule.

### 3.5 A construct that cannot rebuild units enlarges the ones it kept

The spent-source construct arrested with no competent cells, so it cannot make
new units and can only enlarge those it has. After losing 30 % of its units it
ends with 0.70 [0.69, 0.70] of the count and 1.39 [1.34, 1.43] times the mass
per unit, against the 1.43 that full compensation would require, and its mass
returns to 0.97 [0.93, 0.99]. After a 50 % loss it reaches 1.88 [1.78, 1.94]
against the 2.00 required. The model therefore compensates almost completely by
enlargement alone.

Nothing in the model limits how large a unit may become. In a real tissue that
limit exists, and this is where the model and the biology part company
(Section 5).

### 3.6 The cell source, checked against three developing tissues

The model's central claim is about the cell source, so the one check available
without new experiments is whether the progenitor pool each tissue is named for
is present or absent while that tissue is still forming units
(`PAPER_NUMBERS.md`; Figure 4).

All three tissues deplete the progenitor they are named for, which on its own
settles nothing, since differentiation empties a named compartment whether or
not that compartment is what stops the tissue. What separates them is when,
measured against the window over which units are still being formed. In mouse
pancreas the multipotent pool is at 0.0001 of the lineage at E15.5, the end of
the window in which that organ specifies its units: 100 % of the window
elapsed. In human fetal kidney the progenitors are still at 0.1689 of the
lineage at week 18, 43 % of the way through a window that runs to weeks 34 to
36 and closes when those progenitors are exhausted [5]. When nephrogenesis
stops varies between individuals, and in a series spanning weeks 20 to 41 one
infant was still forming nephrons at week 37 [11], so the end of the window is
a typical value and not a constant. In mouse liver the hepatoblast programme
is at 0.0000 by E17.5, 24 % of the way through the window over which hepatocyte
number is still rising. That fraction is an upper bound. We take the window to
end at P20, which is conservative: the one quantification of hepatocyte number
against postnatal age in mouse puts the hyperplastic phase between one and two
months and liver growth ending after two months [12], which would make the
fraction a tenth rather than a quarter. Either way the liver goes on building
after its named programme is gone, by dividing cells it has already made.

Two tissues are limited by the pool they deplete and one is not. That is what
the model predicts from the cell source together with its turnover
(Section 3.7), and it is what counting progenitors alone would miss.

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

### 3.7 Which of the three parameters separates the regimes

The three regimes differ in the binary property $w$, in the size of the
seeded pool $P_0$, and in how the pool turns over, all at once, so none of the
three can be credited with the separation from the comparison alone. Turnover
is a division rate and a loss rate together, and the two are swept as one axis
over the pairings the three sources use: static ($\rho = \delta = 0$),
decaying ($\rho = 0$, $\delta = 0.05$) and balanced ($\rho = \delta = 0.10$).
Separating them would name cells after the kidney and the pancreas that are
not those regimes, since neither divides, and would pair a dividing pool with
no loss, which grows without bound because a refused nucleation does not spend
a progenitor. We vary each axis independently over a grid of two values of
$w$, three of $P_0$ and three of the turnover, and read where the pool stands
at arrest in every cell: spent, held where perfusion has fallen to threshold,
held where units leave no room, or still able to nucleate. Those four shares
are what a regime is here, so a regime that survives a change of $P_0$ and of
the turnover at fixed $w$ is one the binary property owns
(`PAPER_NUMBERS.md`, What separates the regimes; `scripts/34_regimes.py`).

The grid does not credit the binary property. Changing it alone, with the pool
and the turnover held fixed, changes which bound ends the growth in 2 of the 9
settings of the other two, and both are at the smallest pool. Changing the
turnover alone changes it in 5 of the 6 settings of renewal and the pool, and
changing the pool alone in 2 of 6. The binary property is one of three things
that decide which bound ends growth, and on this grid the weakest of them.

The pattern behind those counts reads off the table. A decaying turnover, loss
without division, ends with the pool gone in all six of its cells, whatever
the binary property and whatever the pool size. A renewing source is not
spared: a cell returned when a unit forms is still lost at rate $\delta$ on
every day it forms none. Under a static or a balanced turnover the pool
survives and perfusion binds instead, except at the smallest pool with a spent
source, which is the pancreas-like corner. The binary property shifts where
the boundary between those two outcomes falls; the turnover decides which
outcomes are available at all.

The same grid settles a confound in the elasticity of Section 3.2. That
elasticity is measured at each regime's own pool size, so a renewing source at
$P_0 = 2600$ is compared with a spent one at 60, and a difference between them
could be a difference in what renewal does or in where on the curve each
regime sits. We measure it again with $P_0$ and the turnover matched between
the two sources, leaving the binary property as the only difference. Matched,
the gap is +0.17 at the liver-like corner and +0.28 at the pancreas-like one,
each about two standard errors from zero, and +0.04 at the kidney-like corner,
which is not. The unmatched comparison reports +0.71 against +0.07, a gap of
0.64. So under half of the published contrast is the binary property and the
rest is the difference between a pool of 60 and a pool of 2600.

Three sources that differ in three things at once cannot be read as a
demonstration about one of them. The binary property does change the final
unit count, the elasticity and, at a small pool, which bound ends growth, and
those effects are real and measured. They are also smaller than the comparison
between the three regimes makes them look, and the turnover is what most often
decides where a construct arrests.

### 3.8 Scaling up loses unit density, and past an intermediate size loses reproducibility

Scaling the domain over a 15-fold range of volume, with the pool, the window
and the rates held at the reference, raises the unit count only as
$V^{0.45}$. Carrying the standard allometric scalings of pool, window and rates
as well raises it as $V^{0.58}$ (`PAPER_NUMBERS.md`, Organ size; Figure 5).
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

### 3.9 A joint specification is met up to an intermediate volume, and it is density that fails

The design rules above say what a supply constraint fixes. Put to them as a
specification, they can also be asked to fail. We require of a construct, at a
target volume, that its units occupy at least 10 % of the domain's close
packing, that the tree use at most 2500 segments, and that the unit count vary
by at most 0.15 across independent realisations. We then search the four
things an operator sets, the seeded pool, the nucleation propensity of the
medium, the culture duration and the cell source, for a setting that meets all
three at once. Selection is maximin on the slack each requirement has in its
own units, so the chosen configuration is the one with the most room to spare
on its tightest requirement. It is chosen on 16 seeds and re-run on 24 that
the search never saw, and configurations that had not arrested by the
criterion of Section 3.3 are excluded (`PAPER_NUMBERS.md`, Designing to a
specification; `scripts/35_feasibility.py`).

The three requirements are met by three different configurations, which is the
whole of the difficulty. At the largest domain the best density any
configuration reaches is 5.10 % of close packing, the fewest segments any
reaches is 1468, and the lowest variation any reaches is 0.057; no single
construct has all three.

The specification is met at the two smaller volumes and at neither of the two
larger. What fails is density: at the third volume the best any configuration
comes to the requirement is a slack of −0.34, and at the largest a slack of
−0.49. Both solutions put the knobs at their stops, the largest pool and the
longest culture. The first has a margin of +0.48 on the search seeds and +0.54
held out. The second sits on its requirement rather than above it, with a
margin of +0.01 on the search seeds and +0.01 held out: it meets the
specification, and nothing about it is comfortable. On these runs the search
block understated the held-out margin by 0.03 on average rather than
overstating it, so the configurations did not owe their selection to the
realisations they were chosen on.

Sweeping the specification rather than quoting one triple, no combination of
five density requirements and five reproducibility requirements admits the
largest volume. Loosening reproducibility all the way to 0.50 buys a larger
volume at one density requirement only; at every other, density alone decides.

That is the result worth carrying into Section 4.1, and it sharpens Section
3.8 rather than repeating it. The reproducibility that collapses with size
there, a coefficient of variation of 0.47 at the largest domain, is a property
of the reference configuration and not of the volume: searching the knobs
recovers 0.057 at that same volume, with a spent source. Reproducibility at
size is recoverable by choosing the cell source. Density at size is not
recoverable by any setting we searched.

Held-out seeds are new realisations of the same model and not new evidence
about it. What this measures is whether a configuration chosen on one batch
holds on another, which is reproducibility across batches. It is not a
validation of the model against a construct, and the boundary it draws is a
boundary of the model.

## 4. Discussion

### 4.1 Five design rules, stated as the model gives them

**A constant flow gives a tunable brake on total mass that a constant pressure
does not.** At a fixed inlet pressure the model's construct overshoots to 1.82
times its mass after injury and does not settle; at fixed root flow its mass
rises as the delivered flow to an exponent of 0.90 over a range of 0.6 to 1.4.
Nothing in the tissue measures its own size in either case, so the brake is the
supply the circuit provides rather than a property of the tissue. Two
qualifications travel with it. The brake sets the total mass and not the number
of units that carry it: in this model the construct reaches its budget in
aggregate while most of its units stop short of their own ceilings, for the
reason given in Section 3.4. And the comparison is between two boundary
conditions in a model, not between two bioreactors: maximum admissible
pressure, shear damage, and the capacity of the perfusion circuit are
constraints the model does not carry, and any of them could decide the question
in practice.

**What seeding more cells buys depends on the source and on the pool you
are already at, and more on the second than on the first.** Across the three
regimes as published the elasticity of the final unit count to the seeded pool
runs from +0.07 ± 0.06 to +0.71 ± 0.05, and it is tempting to read the range as
the difference between a renewing source and a spent one. Matched for pool and
turnover, under half of that range is the source (Section 3.7). The operational
rule is the narrower one: a small pool repays seeding and a large one does not,
whichever source it is, and a source that is spent in forming units offers no
second chance once the pool is gone.

**Lower the perfusion a unit requires before seeding more cells.** The
setpoint elasticity is between −0.77 and −1.83 across the three regimes,
against a pool elasticity between +0.07 and +0.71. In every regime the
threshold is the stronger lever, and in the renewing regime it is the only one
that matters.

**Check that your specifications are jointly reachable at the volume you
want, before optimising any of them.** Density, vascular cost and
reproducibility are reached by different settings of the same knobs, so a
construct that is best at one is not best at the others. In this model a
requirement of 10 % of close packing, 2500 segments and a coefficient of
variation of 0.15 is met up to an intermediate volume and at no larger one,
and what fails is density: at the largest domain no setting we searched
exceeds 5.10 % of close packing. Reproducibility is the requirement that looks
hopeless at size and is not, since choosing a spent source recovers a
coefficient of variation of 0.057 where the reference configuration gives 0.47.
The boundary is a boundary of this model and not a measured limit of any
construct, and it is the kind of statement a design study should produce
whether or not it falls where one hoped.

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

- **Constant flow against constant pressure, and final mass against flow.**
  Identical perfused constructs driven at constant inlet pressure and at
  several constant flows, followed past the point where growth would arrest.
  The first should not settle; the others should settle at masses that rise
  close to proportionally with the delivered flow. A construct whose final
  mass does not depend on the flow it is given would falsify the rule.
- **Count against perfused tissue.** Unit number and perfusion-weighted
  tissue mass measured in the same construct as it rebuilds after arrest. The
  count should return first and the perfusion-weighted mass lag or never
  arrive. A construct whose count and perfused tissue recover together would
  contradict the model.
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

No parameter in this paper is fitted to data. The values are chosen, the same
values are used for all three regimes except the four that define a regime,
and what we report are ratios, orderings and contrasts rather than absolute
sizes, because those are what a chosen parameter set can support. During the
development of the model an earlier version was tuned so that its renewing
organ reproduced a postnatal contraction of the spread of lobule number
measured in mouse [13]; that tuning is not in the pipeline this paper runs,
no number reported here depends on it, and we mention it only because the
repository keeps it. It would in any case be thin evidence: that study counted
lobules at two ages, in three animals at one week and six at eight weeks, and
reported the change as roughly twofold without giving a spread in the text.

No prediction in this paper has been tested against a construct or a tissue,
and the model carries no measure of biological function, so every recovery
reported here is structural.

The spatial constraint is local exclusion between neighbouring units, not the
packing limit of the container, which never binds.

**The reference model excludes involution; an extended model introduces it as
an uncalibrated sensitivity.** In the reference model a unit that loses its
supply stops growing and stays at the mass floor, where a tissue would atrophy
and resorb it, and that is why a construct held at a fixed inflow reaches its
budget in total while most of its units sit far below their own ceilings. It
is the clearest place where the model is missing a process rather than a
parameter. The extended model supplies the process, and Section 3.4 reports
what it does; none of the four numbers it runs on is measured, which is why
every result above is the reference model's. The two also differ in what the
mass floor is: the reference model holds every unit at or above the mass a
unit is born with, and the extended one lets a shrinking unit pass below it on
its way out.

The model compensates after injury further than a kidney does. After removal of
one of two kidneys the adult mammalian remnant forms no new nephrons and gains
about 30 to 50 % in mass, recovering roughly 70 to 80 % of the filtration rate
rather than all of it [10]. The model, given the same 50 % loss, reaches 1.88
times the mass per unit against the 2.00 of full compensation: it reproduces
the mechanism, enlargement of the units that remain with no new ones formed,
but overshoots its extent. Nothing in the model limits how large a unit may
become; in a tissue that limit exists and passing it is injury.

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

The perfusion-weighted tissue mass is a hydraulic quantity, the mass of each
unit weighted by its supply against the threshold. It is not a measurement of
any biological function, and nothing here relates it quantitatively to one.

Recovery after injury requires the ability of a unit to enlarge to be restored
by hand, through a parameter that resets the growth clock. The model does not
derive that restoration from the injury, and the consequences of restoring it
only in part are not swept.

The involution of Section 3.4 carries four choices and all four are swept: the
rate, the supply at which a unit counts as starving, the consecutive days it
must stay there before it begins to lose mass, and the mass at which it is
removed rather than left to shrink (`PAPER_NUMBERS.md`, Involution, as a
sensitivity). What they move is when the removals happen, not what the
construct ends as. That is a sweep and not a calibration: none of the four is
measured, and a mechanism whose every setting gives the same answer is still a
mechanism nobody has yet observed in a construct.

The vascular tree does not lose a branch when the tissue under it goes.
Calibres adapt and a vessel below its critical radius regresses, so a branch
that loses its tissue loses its flow and narrows over the following days, but
that is slower than the loss it answers and the model does not prune. The
removal does leave the tree coherent, with flow conserved and every surviving
unit still hanging from a segment that exists, and the share of tips serving
no unit is reported alongside the rest.

## 6. Conclusion

The size a perfused tissue reaches, whether its arrest can be undone, and what
seeding more cells buys all follow from how the cell source is maintained and
from the boundary condition on the perfusion circuit. Under a fixed inlet
pressure there is no ceiling; under a fixed flow the ceiling is the delivered
flow over what the tissue requires per unit mass, which is a consequence of
the local rule rather than a separate assumption. Neither is a
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

[4] Wang B et al. 2015 Self-renewing diploid Axin2+ cells fuel homeostatic
renewal of the liver *Nature* **524** 180–5

[5] Perl A J, Schuh M P and Kopan R 2022 Regulation of nephron progenitor cell
lifespan and nephron endowment *Nat. Rev. Nephrol.* **18** 683–95

[6] Hochane M et al. 2019 Single-cell transcriptomics reveals gene expression
dynamics of human fetal kidney development *PLoS Biol.* **17** e3000152

[7] Lindström N O et al. 2021 Spatial transcriptional mapping of the human
nephrogenic program *Dev. Cell* **56** 2381–98

[8] Bastidas-Ponce A et al. 2019 Comprehensive single cell mRNA profiling
reveals a detailed roadmap for pancreatic endocrinogenesis *Development*
**146** dev173716

[9] Yang L et al. 2017 A single-cell transcriptomic analysis reveals precise
pathways and regulatory mechanisms underlying hepatoblast differentiation
*Hepatology* **66** 1387–401

[10] Qu X, Shi Y, Li L and Liu T 2026 Compensatory nephron hypertrophy:
mechanisms, pathophysiology dynamics, and clinical implications in acute kidney
injury and chronic kidney disease *Kidney Res. Clin. Pract.* advance online
publication, doi:10.23876/j.krcp.25.405

[11] Ryan D et al. 2018 Development of the human fetal kidney from mid to late
gestation in male and female infants *EBioMedicine* **27**,
doi:10.1016/j.ebiom.2017.12.016

[12] Delone G V, Uryvaeva I V, Koretskiĭ V F and Brodskiĭ V Ya 1987 Analysis of
the postnatal growth of the mouse liver by estimating the hepatocyte count,
their weight and ploidy *Ontogenez* **18** 304–8 (in Russian, English abstract)

[13] Azizoglu D B, Perez K, Zheng S L, Rahman S, Rim E Y, Anbarchian T, Fish M,
Loh K M, Red-Horse K and Nusse R 2025 Vascularization of neonatal liver lobules
presages adult liver size *Nat. Commun.* **16** 9989

## Notes for the authors, to be removed before submission

**Still to confirm on the title page.** The spelling of Enrique
Hernández-Lemus, and whether he should also carry a Centro de Ciencias de la
Complejidad, UNAM, affiliation. Whether Miguel Ángel Peña-Rico carries the
Universidad del Papaloapan affiliation alone or shares another, since it was
given for him and Angélica S. Martínez-Ramírez together. Where the three
authors added last should sit in the order: they are placed here after Genaro
Vázquez-Victorio and before the senior author, which keeps the corresponding
author last, but the group decides. ORCID identifiers are recommended by IOP
and are not yet here.

**Every claim that needed a citation now has one.** References [4] to [13]
were located and verified against the publishers' own pages. The author lists
are given as *et al.* and must be expanded to the form Biofabrication asks for,
and DOIs added, at submission. Two fields are still missing: the volume and
pages of [11], of which only the volume was confirmed, and whether [10] has
been assigned a volume since it appeared online.

The three that were open are closed as follows, and two of them changed what
the manuscript says rather than only adding a reference.

1. **The calibration was not a calibration of this paper.** The tuning against
   mouse lobule number lives in `scripts/12_calibrate_closed.py`, which sits
   under the `closed` target, is labelled there as exploratory and with
   parameters from earlier versions, and is not run by `make paper`. Nothing in
   `PAPER_NUMBERS.md` depends on it: nu = 8 is a chosen value used throughout,
   and `sigma_param` belongs to the open engine, which this manuscript does not
   use. Section 5 claimed a single calibration and now says what is true, that
   no parameter is fitted and that the results are ratios and orderings. The
   mouse tuning is mentioned once, with [13], and with the reason it would be
   thin evidence anyway: that study counted lobules at two ages only, three
   animals at one week and six at eight weeks, and gives no spread in its text.
   The numbers the repository tuned against, 0.229, 0.074 and 1.82-fold, are
   read from Fig. 1g and are not in the text; they stay in the repository as
   the record of what was done and no longer appear as a claim.
2. **The weeks of the human nephrogenesis window** are now [11], which also
   supplies something better than a constant: across weeks 20 to 41 the stop
   varies between individuals and one infant was still forming nephrons at
   week 37. The manuscript says so.
3. **The mouse liver window** is now stated as an upper bound rather than a
   figure. P20 is kept in `WINDOW` as the short end of what the literature
   supports, and [12] is cited for the hyperplastic phase running between one
   and two months, which would put the fraction nearer a tenth. Keeping the
   conservative window means no checked number moves and the claim is the
   weaker of the two, which is the one worth making.

The mouse pancreatic window E9.5 to E15.5 is covered by [2] for the claim the
manuscript actually makes about it. If a source for the window itself is
wanted, Zhou Q et al. 2007 *Dev. Cell* **13** 103–14 is the candidate and has
not been verified here.

**Closed since the fourth review.** The involution sensitivity now sweeps all
four of its choices, the rate and the three conditions that turn it on, and
reports the day by which nine tenths of the removals have happened alongside
the endpoint (`PAPER_NUMBERS.md`, Involution, as a sensitivity). The tree is
checked for coherence after a removal in `tests/test_smoke.py`.

**From the fifth review.** One of its findings did not survive checking: it
read the mass floor as a per-timestep decision, so that a surviving unit's
mass would depend on whether some other unit happened to be removed that day.
It is not. The `else` that re-applies the floor belongs to `if involution >
0.0`, not to the removal branch, so the floor is a property of which model is
running and not of the day. The reference model applies it every day and the
extended model never does. Two checks in `tests/test_smoke.py` now pin both
halves; one of them ends a run on a day that removed nothing and finds units
below the floor, which could not happen if the reading were right.

Three of its findings are in the text above: arrest now carries a stated
criterion (Section 3.3), the limitations no longer say the model has no
involution while Section 3.4 reports it, and the drift of the factor between
the two starvation criteria is measured rather than assumed (Section 3.4).

Its strongest finding is answered in Section 3.7, and the answer went against
the framing this manuscript was written in. The three regimes differ in the
binary property, the seeded pool and the turnover at once, so the separation
could not be credited to the first until the others were held fixed. Held
fixed, the binary property changes which bound ends growth in 2 of 9 settings
against the turnover's 5 of 6, and under half of the elasticity contrast
survives matching. Sections 3.1, 3.3 and 3.6 were headed "one property" and
are no longer: what they report is unchanged, and what it is evidence for is
narrower than those headings claimed.

**Four results from the earlier framing are not in this manuscript** and are
listed in `FIGURES.md` as supplementary: the 96 human kidneys from 48 donors,
the Murray ratio in grown trees, the contraction of the spread of unit number
with its calibration against mouse liver, and the transmission of an
individual spread. The calibration has the strongest case for returning to the
main text, since it is the only place the model is anchored to a measurement.


## Figures

Legends are in `FIGURES.md` and are regenerated with the figures by
`make paper`. File names in the repository are historical; the mapping is:

| here | file | script |
|---|---|---|
| Figure 1 | `Fig1_model` | `figures/fig_paper.py` |
| Figure 2 | `Fig2_bounds` | `figures/fig_paper.py` |
| Figure 3 | `Fig8_injury` | `figures/fig_injury.py` |
| Figure 4 | `Fig9_transcriptome` | `figures/fig_transcriptome.py` |
| Figure 5 | `Fig6_scaling` | `figures/fig_design.py` |

Figures 1 and 2 each carry a schematic and the measurement it explains, so the
argument can be followed from the figures alone.

`FIGURES.md` also lists the figures the repository builds that this manuscript
does not use. They belong to the earlier framing kept in `PAPER.md` and are
available as supplementary material.
