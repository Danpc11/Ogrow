# Which resource runs out decides whether an organ's size depends on its progenitors

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

Target journal: *Cell Systems*. This is the second of two manuscripts drawn
from the same model and results. It asks what limits the size of a growing
tissue and how the tissue then answers a loss of progenitors; `MANUSCRIPT.md`
carries the scaling and the design specification and is the base of the
engineering paper. The order here is Summary, Introduction, Results,
Discussion, Methods; highlights and the remaining Cell Press items are in the
notes at the end.

Reference model version 4.66.0 (`src/organgrow/closed.py`). Every model number
quoted below is taken from `PAPER_NUMBERS.md`, which `make paper` regenerates
with the version and commit that produced it, and is checked against the
result files by `scripts/29_check_manuscript.py`, which holds this manuscript
and `MANUSCRIPT.md` to the same results. The passages between `auto:` markers
are written by the scripts that compute their numbers.

---

## Summary

Why does removing progenitors shrink a pancreas in proportion and leave a
liver its size? We study one growth law on a perfused branching tree: units
form at tips at a rate set by the progenitors there and by perfusion above a
threshold, until perfusion, progenitors or room runs out. Which runs out
decides how the tissue answers a loss. Seeded with the progenitors that mouse
ablations left,
<!-- auto:validation-brief -->a spent source keeps 0.44 of its mass, against 0.36 measured in pancreas, and a renewing one 0.86, against about 0.96 in liver<!-- /auto:validation-brief -->.
A kidney-like source fails until its progenitors are made to drive branching;
then pool, cap and tips behave as measured, and the nephron deficit
approaches proportion only where perfusion is not limiting. The model names
the measurement that identifies the limiting resource and gives the outcome of
ablations not yet done.

**Highlights**

- One local nucleation law gives three ways a growing tissue stops
- Which resource runs out sets how organ size answers progenitor loss
- Pancreas and liver sizes after progenitor ablation come out without fitting
- Kidney needs progenitor-driven branching; untested ablations are predicted

**Keywords:** organ size, progenitor pool, regeneration, perfusion, branching
morphogenesis, nephron endowment, computational model

## 1. Introduction

In the mouse embryo, killing a share of the progenitors of the pancreas leaves
an animal with a pancreas smaller in the same proportion for life, and killing
a larger share of the progenitors of the liver leaves a liver of normal size
four days later [2]. Killing two fifths of the progenitors of the kidney
leaves a fifth fewer nephrons in one study and half as many in another [14,
15]. These experiments ask one question of three organs, how the size of an
organ depends on the cells it is built from, and get three answers. The usual
reading sorts organs into those whose size is set intrinsically, by a number
of progenitors each good for a fixed amount of tissue, and those whose size is
regulated from outside [2]. That names the difference without saying where it
comes from, and it has no place for a kidney that answers both ways.

Kidney, liver and pancreas are among the organs for which replacement tissue
is most needed. Kidney and liver are the two most transplanted organs, 110 467
and 42 497 of the 173 727 solid organ transplants reported worldwide in 2024,
an activity estimated to cover at most a tenth of the need [17]; the pancreas
is transplanted far less often, 2 066 in the same year, and what is sought
there is mostly its endocrine tissue. The three answer a loss in opposite
ways. A liver regains its mass after most of it is removed. A kidney forms no
nephron once nephrogenesis has ended, and a pancreas recovers little. In
liver and pancreas no dedicated stem-cell niche does the work of repair, which
falls to differentiated resident cells [18], so the difference between them is
not one of having or lacking a stem cell. Why they differ, and what would
have to be supplied to a tissue of each kind for it to grow or to grow back,
are the questions this paper puts to a model.

We take a different starting point from a count of progenitors or a reserve of
stem cells. We ask which locally available resource runs out at the tips of a
growing, perfused tree, and show that the answer decides how the tissue
responds when progenitors are removed, whether its growth can be restarted,
and what sets its size after an injury.

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
and is lost. Section 2.2 separates the three contributions and finds the
turnover the largest of them.

The paper has three parts. The first describes the three regimes, what holds
each at arrest, and which parameter of the cell source separates them
(Sections 2.1 to 2.3). The second sets the model beside measurements that
exist: the size of pancreas and liver after progenitor ablation, the nephron
deficit after two ablations of the kidney, and the difference between the two
kidneys of one person (Section 2.4). The kidney comparison fails, and an
extension in which progenitors drive branching shows what the model was
missing and what the failure says about that kidney (Section 2.5). The third
follows the tissue past arrest, through restart and injury (Sections 2.6 to
2.8), and commits the model to the outcome of ablations that have not been
done (Section 2.9). No parameter is fitted, and no prediction has yet been
tested.

## 2. Results

### 2.1 The three cell sources arrest in three different states

The three sources arrest in three different states (`PAPER_NUMBERS.md`,
Regimes at closure; Figure 1). They differ in three parameters at once, and
Section 2.2 separates them: the binary property is one of the three and the
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
simply lost. We return to it in Section 3.3.

The ordering of the three regimes survives all 13 variants of the auxiliary
rules (`PAPER_NUMBERS.md`, Robustness).

### 2.2 Which of the three parameters separates the regimes

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

The same grid settles a confound in the elasticity of Section 2.3. That
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

### 2.3 How size answers the pool depends on which bound is active

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
Section 2.2 reports how much of the difference each contributes.

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
model was not fitted to that experiment; Section 2.4 puts the same
perturbation to the model and reads the sizes.

### 2.4 Three numbers set beside measurements

Everything above is the model read against itself. Three of its numbers can be
set beside measurements that already exist, with nothing fitted (Figure 3;
`PAPER_NUMBERS.md`, Against published perturbations and paired human
kidneys).

<!-- auto:validation -->
Seeded with 39 % of its reference pool, the share of pancreatic progenitors the mouse ablation left, the pancreas-like construct ends with 0.44 [0.37, 0.52] of its mass; the newborn pancreas weighed 0.36 of control (95 % interval 0.19 to 0.53). Seeded with 92 %, it ends with 0.81 [0.72, 0.90], against 0.77 (0.53 to 1.02) measured. Seeded with 35 %, the liver-like construct ends with 0.86 [0.84, 0.89] of its mass; the liver, at about a third of control weight when its ablation ended, stood at about 0.96 four days later. Model values are ratios of means over 24 paired seeds with a 95 % bootstrap interval. The pancreas-like construct keeps 0.17 of its mass on a quarter of the pool, and the liver-like one 0.77. With 40 % of the kidney-like pool removed the construct ends with 0.80 [0.77, 0.83] of the units, a deficit of 20 % [17, 23]; the two published ablations of that size left 12 % and 21 % fewer glomeruli in one study and 42 % and 52 % fewer in the other. Grown 24 times from identical parameters, the kidney-like construct has a coefficient of variation of 0.096 in mass and of 0.088 in unit count; on the scale the human pairs are read on, the standard deviation of its log mass is 0.104. In the human cohort the volumes of a person's two kidneys differ with a standard deviation of log volume of 0.098 per kidney, over 48 pairs, once the mean difference between sides is removed (the left smaller by 7.4 % on average); that difference does not rise detectably with kidney size ($\rho = +0.25$, $P = 0.09$, Spearman, over volumes from 109 to 553 cm$^3$).
<!-- /auto:validation -->

*Removing progenitors from pancreas and liver.* Stanger and colleagues killed
pancreatic progenitors in the mouse embryo over a controlled window and
weighed the organ at birth [2]. Removing 61 % of the pool left a pancreas at
0.36 of control weight; removing 8 % left one at 0.77, not distinguishable
from control; and the small pancreas stayed proportionately small to eleven
weeks of age. Seeded with the same shares of its pool, the pancreas-like
construct ends close to both (Figure 3A). The same method removed at least
65 % of hepatic progenitors: the liver weighed about a third of control when
the ablation ended and was back to control weight four days later [2]. The
liver-like construct recovers most of the way and not all of it, and its
shortfall grows as the pool is cut further. Three qualifications. The regime
parameters were fixed before this comparison was computed and were not
adjusted to it, but the experiment was known when the regimes were described
(Section 2.3), so the agreement is a check passed and not a blind prediction.
The embryonic liver was already built when its progenitors were removed, which
makes that experiment closer to an injury than to a smaller seeding, so only
the outcome is compared. And the liver values are read from a fitted curve in
the published figure.

*Removing nephron progenitors.* Two groups removed about 40 % of nephron
progenitors in the mouse embryo with the same toxin and found different
things. Ablating the CITED1-expressing subset at E15.5 left about 12 % fewer
glomeruli three days after birth and 21 % fewer at ten days; the surviving
progenitors more than doubled their proliferation and refilled the niche
[15]. Ablating cap mesenchyme through a Gdnf driver at E12.5 left 52 % fewer
glomeruli at birth and 42 % fewer in the adult, a deficit in proportion to the
cells removed; there the survivors did not divide faster, the ureteric tree
branched more slowly, and nephrogenesis ended on schedule [14]. The
kidney-like construct with 40 % of its pool removed loses a share of its units
that sits at the edge of the first range and is less than half of the second
(Figure 3B). The model has the mechanism of neither experiment. Its
kidney-like source does not divide, so it cannot refill a niche, and its tree
does not depend on the pool, so branching cannot slow. We therefore do not
claim that the kidney-like regime predicts a nephron deficit. What the second
experiment shows is that where progenitors also drive the branching that
makes room for them, the deficit is proportional, and that coupling is absent
from the model.

*Two organs grown alike.* Two constructs grown from identical parameters
differ only in their realisation, which is the situation of the two kidneys of
one person. The spread of the kidney-like construct across seeds and the
left-right spread of kidney volume are close (Figure 3C), and nothing in the
model was set by the cohort. Four things limit what that carries. A mass in
model units is compared with a CT volume. The two raters of the cohort agreed
to a Dice score above 93 % [16], so part of the left-right spread is
segmentation and the biological part is smaller than the figure quoted. The
agreement is at the one domain scale all regimes are reported at. And the
model's spread rises with the domain (companion manuscript), while the pairs show no
clear rise with kidney size; people of different size are not one construct in
a larger domain, so this does not refute the model, but it does not support
that part of it.

### 2.5 An extension: progenitors that make their tip branch

The kidney comparison of Section 2.4 fails for a reason that can be named. In
the embryonic kidney the progenitors around a tip are what makes the tip
branch [14], and the reference model has no path from the pool to the tree.
We asked what the model needs for both ablation experiments by adding that
path and nothing else (Section 4.5). The rules were written after the
measurements were known, so the endpoint is not a prediction. What can be
asked is whether one path gives several things at once that were measured
separately, and whether the single difference between the two experiments,
whether the surviving progenitors divide faster, is enough to separate their
outcomes (Figure 4; `PAPER_NUMBERS.md`, An extension).

<!-- auto:kidney-coupling -->
With no proliferative response, and in all four settings of the nucleation law, the ablated construct ends the window with between 0.53 and 0.57 of the control's progenitors and between 0.52 and 0.60 of its branching sites. The day after the ablation its cap stands at 0.61 to 0.90 of control, and it is back to 0.95 of control within 3 to 50 days in 4 of the 4 settings. Its units end at 0.81 [0.72, 0.94] of control where nucleation is set as in the reference model, and at 0.64 [0.57, 0.70] where the perfusion threshold is a quarter of the reference and the law four times slower; counting only units formed after the ablation, 0.77 and 0.52 of control. The measured 0.48 to 0.58 lies outside the first interval and inside the second. With division doubled in a depleted cap, the pool ends at 0.65 to 0.68 of control and the units at 0.85 [0.75, 0.95] and 0.86 [0.83, 0.89] in the same two settings; the measured 0.79 to 0.88 lies inside the first interval and inside the second. Ratios of means, ablated over control, 24 paired seeds, 95 % bootstrap intervals.
<!-- /auto:kidney-coupling -->

*The branching side comes out of the coupling alone.* After ablation without a
proliferative response, the progenitors of the mouse kidney went on dividing
as fast as in control, so the deficit in cells persisted; the cells around
each tip stood at 0.60 of control a day after tamoxifen and at 0.95 and 1.04
on the next two days; and the tips, unchanged at first, fell to 0.71 and then
0.59 of control [14]. The model gives the same three things in every setting
of the nucleation law: the pool stays reduced, the cap returns on its own, and
the sites settle near the share of the pool that was left. The reason is the
one the experiment proposed. A tip that has lost part of its cap waits to
branch until the cap has regrown, so the cap is restored at the cost of the
tree. The model's cap recovers less cleanly than the measured one,
overshooting before it settles, on a tree of only a few sites.

*Fewer sites become proportionally fewer units only where perfusion does not
limit nucleation.* With the nucleation law set as in the reference model, the
units of the ablated construct recover much of the way: each remaining site is
better supplied and makes more units. That is the central mechanism of this
paper acting as a compensation, and the embryonic kidney does not show it.
Read through the model, a nephron deficit in proportion to the progenitors
removed says that nephron formation in that kidney was limited by its
progenitors and its branching, and not by perfusion.

*One difference between the experiments is enough to separate them.* Doubling
the division rate in a depleted cap, which is what the second study measured
[15], brings the units of the ablated construct most of the way to control. It
does so with a pool that recovers less than the measured one, which stood at
0.85 of control two days after the ablation [15], so the route is right in
kind and not in detail.

The extension carries numbers that are not measured: the size of a cap, the
division rate, the length of the window and the day of the ablation, which
were fixed, and the two that set how far perfusion limits nucleation, which
were swept and decided the second result. It states what a kidney-like tissue
needs that the reference model lacks. It is not a calibrated model of
nephrogenesis, and every other result in this paper is the reference model's.

### 2.6 Arrest is stable, and whether it can be undone follows from what the pool has left

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

### 2.7 After injury, how the circuit is driven decides the final size

A construct grown to arrest, then stripped of 30 % of its units by a plane
cut and followed for 80 more days, does not recover on its own
(`PAPER_NUMBERS.md`, Injury; Figure 5). Left as grown, the renewing construct
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
the reason is the budget of Section 4.3. Each unit's share is its own
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
by fewer and larger units, which is the behaviour Section 2.8 reports for the
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
of function at all (Section 3.4), so what is reported here is structural
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

### 2.8 A construct that cannot rebuild units enlarges the ones it kept

The spent-source construct arrested with no competent cells, so it cannot make
new units and can only enlarge those it has. After losing 30 % of its units it
ends with 0.70 [0.69, 0.70] of the count and 1.39 [1.34, 1.43] times the mass
per unit, against the 1.43 that full compensation would require, and its mass
returns to 0.97 [0.93, 0.99]. After a 50 % loss it reaches 1.88 [1.78, 1.94]
against the 2.00 required. The model therefore compensates almost completely by
enlargement alone.

Nothing in the model limits how large a unit may become. In a real tissue that
limit exists, and this is where the model and the biology part company
(Section 3.4).

### 2.9 Predictions for ablations that have not been done

Sections 2.4 and 2.5 compare the model with what is known. Here it is
committed to numbers for experiments nobody has reported, in the same systems
and by the same method: progenitor ablation in the mouse embryo, at doses and
at times that were not used (Figure 6; `PAPER_NUMBERS.md`, Predictions). The
dose predictions are the reference model's. Those for the kidney are the
extension's, in the setting where perfusion limits nucleation least, and are
given with and without a proliferative response because the model does not
say which a given experiment will show.

<!-- auto:predictions -->
For the pancreas, an ablation that leaves 0.75 of the progenitors is predicted to give a pancreas at 0.70 [0.58, 0.81] of control weight, and at 0.61 [0.51, 0.70] with 0.60; a pancreas of normal weight at either dose refutes it. For the liver, one that leaves 0.75 is predicted to give a liver at 0.97 [0.94, 1.00] of control weight, and at 0.70 [0.68, 0.72] with 0.15: compensation that is complete after a shallow ablation and fails after a deep one. A liver of normal weight after 85 % of its progenitors are removed refutes it. For the kidney, in the coupled source and where the surviving progenitors do not divide faster, the units are predicted at 0.82 [0.75, 0.89] after 20 % and 0.56 [0.50, 0.63] after 60 % of control. The deficit is not in proportion at depth: a proportional deficit would leave 0.40 after 60 %, which is outside the interval. Where they do divide faster, the units are predicted at 0.93 [0.90, 0.97] after 20 % and 0.79 [0.77, 0.82] after 60 %. For the same 40 % killed at two times, the units are predicted at 0.67 [0.56, 0.79] early and 0.75 [0.71, 0.79] late in the window without a response, and at 0.93 [0.86, 1.05] early and 0.87 [0.84, 0.90] late with one. Ratios of means to the unablated control, 24 paired seeds, 95 % bootstrap intervals; model mass for pancreas and liver, units for kidney.
<!-- /auto:predictions -->

Two of these go against what the experiments so far would lead one to expect.
The published liver ablation removed about two thirds of the progenitors and
the organ recovered [2]; the model recovers most of the way at that depth and
says the recovery fails beyond it, so that a liver is not an organ whose size
is free of its progenitors but one with a wide margin. And the kidney
ablation that gave a deficit in proportion removed two fifths of the
progenitors [14]; the extension says the proportion does not hold at depth,
because the sites that remain are each better supplied, which is the
compensation of Section 2.5 returning as the tree thins. Either is a single
experiment with tools that exist.

## 3. Discussion

### 3.1 Five criteria, as the model gives them

**Find the bound that ended growth before intervening.** Growth stops when one
of three local resources runs out, and each leaves its own signature at
arrest: a pool still held at tips short of supply in the liver-like regime, a
pool that is gone in the kidney-like regime, a budget part spent and part
stranded in the pancreas-like one (Section 2.1). What helps is different in
each case, and Section 3.2 gives the two reliefs that tell them apart in a
construct.

**A tissue can resume only if it stopped with competent cells left.**
Relieving perfusion or room restarts the liver-like construct, and nothing
restarts one whose pool is spent (Section 2.6). Regenerative capacity is not a
separate property in this model. It is what the pool has left at arrest, which
the turnover of the source decides more often than whether its cells are spent
(Section 2.2). For kidney-like and pancreas-like tissue the criterion
therefore applies when the tissue is built: competent cells have to be kept or
supplied, because nothing done to perfusion afterwards makes a new unit. The
extension of Section 2.5 adds a requirement for kidney-like tissue that the
reference model cannot state: where progenitors also make the tree branch,
the number of sites follows the number of progenitors, so they have to be
supplied in proportion to the units wanted, and better perfusion makes up a
shortfall only where perfusion was the limit.


**Hold the delivered flow, not the pressure, to set the size regained after
injury.** A constant flow gives a tunable brake on total mass that a constant
pressure does not. At a fixed inlet pressure the model's construct overshoots to 1.82
times its mass after injury and does not settle; at fixed root flow its mass
rises as the delivered flow to an exponent of 0.90 over a range of 0.6 to 1.4.
Nothing in the tissue measures its own size in either case, so the brake is the
supply the circuit provides rather than a property of the tissue. Two
qualifications travel with it. The brake sets the total mass and not the number
of units that carry it: in this model the construct reaches its budget in
aggregate while most of its units stop short of their own ceilings, for the
reason given in Section 2.7. And the comparison is between two boundary
conditions in a model, not between two bioreactors: maximum admissible
pressure, shear damage, and the capacity of the perfusion circuit are
constraints the model does not carry, and any of them could decide the question
in practice.

**Judge recovery by units, not only by mass.** A construct held at a fixed
flow regains its mass on roughly a third of its units once involution is
allowed (Section 2.7), and a kidney-like construct regains mass by enlarging
the units it kept (Section 2.8). Mass, unit number and perfused tissue recover
apart, so a criterion of success written in mass alone can be met by a tissue
that has lost the architecture it was built for.

**Lower the perfusion a unit requires before seeding more cells.** What
seeding more cells buys depends on the source and on the pool you are already
at, and more on the second than on the first. Across the three
regimes as published the elasticity of the final unit count to the seeded pool
runs from +0.07 ± 0.06 to +0.71 ± 0.05, and it is tempting to read the range as
the difference between a renewing source and a spent one. Matched for pool and
turnover, under half of that range is the source (Section 2.2). The operational
rule is the narrower one: a small pool repays seeding and a large one does not,
whichever source it is, and a source that is spent in forming units offers no
second chance once the pool is gone.

The threshold is the stronger lever. The
setpoint elasticity is between −0.77 and −1.83 across the three regimes,
against a pool elasticity between +0.07 and +0.71. In every regime the
threshold is the stronger lever, and in the renewing regime it is the only one
that matters.

### 3.2 What would test this

The predictions of Section 2.9 are stated with what refutes each. The
following would each decide one of the other results, and none has been done.

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

### 3.3 What the model cannot do

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

### 3.4 Limitations of the study

No parameter in this paper is fitted to data. The values are chosen, the same
values are used for all three regimes except the four listed in Section 4.2,
and what we report are ratios, orderings and contrasts rather than absolute
sizes, because those are what a chosen parameter set can support. During the
development of the model an earlier version was tuned so that its renewing
organ reproduced a postnatal contraction of the spread of lobule number
measured in mouse [13]; that tuning is not in the pipeline this paper runs,
no number reported here depends on it, and we mention it only because the
repository keeps it. It would in any case be thin evidence: that study counted
lobules at two ages, in three animals at one week and six at eight weeks, and
reported the change as roughly twofold without giving a spread in the text.

No criterion in this paper has been tested in a construct. Three numbers are
set beside existing measurements in Section 2.4, and that comparison has the
limits stated there: the pancreas and liver experiment was known when the
regimes were described, the liver-like construct recovers less completely
than a liver does, the kidney-like deficit matches the smaller of two
published deficits and is less than half the larger, the model has the
mechanism of neither kidney experiment, and the paired kidneys are CT volumes
that carry segmentation error. The model carries no measure of
biological function, so every recovery reported here is structural.

The extension of Section 2.5 was written after the measurements it is set
against were known, carries six numbers of which none is measured, and was
swept in two. It is run in a domain of fixed size, where an embryonic kidney
makes its own, and its time is not calibrated to embryonic days, so sequences
and endpoints are compared and timings are not.

All three regimes are grown in one domain shape and at one set of geometric
settings. They are liver-like, kidney-like and pancreas-like in their cell
source and in nothing else, and no statement here is about the anatomy of an
organ.

The spatial constraint is local exclusion between neighbouring units, not the
packing limit of the container, which never binds.

**The reference model excludes involution; an extended model introduces it as
an uncalibrated sensitivity.** In the reference model a unit that loses its
supply stops growing and stays at the mass floor, where a tissue would atrophy
and resorb it, and that is why a construct held at a fixed inflow reaches its
budget in total while most of its units sit far below their own ceilings. It
is the clearest place where the model is missing a process rather than a
parameter. The extended model supplies the process, and Section 2.7 reports
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

The two interventions of Section 2.6 are not matched in strength, so they show
that both bounds hold the renewing construct, not which holds it more. Room at
a tip is estimated from trial placements in the final state, which measures
residual capacity rather than whether growth could resume.

The progenitor processes converge with one, two and four steps per day with the
tree held fixed within the day; the coupled vascular and tissue dynamics run
once per day and their convergence at a smaller step was not tested. The
critical radius is the floor of calibre adaptation, so simulated trees cannot
show calibres below it.

The perfusion-weighted tissue mass is a hydraulic quantity, the mass of each
unit weighted by its supply against the threshold. It is not a measurement of
any biological function, and nothing here relates it quantitatively to one.

Recovery after injury requires the ability of a unit to enlarge to be restored
by hand, through a parameter that resets the growth clock. The model does not
derive that restoration from the injury, and the consequences of restoring it
only in part are not swept.

The involution of Section 2.7 carries four choices and all four are swept: the
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

## 4. Methods

### 4.1 The nucleation law

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
property of Section 4.2 are drawn in Figure 1a to 1c.

### 4.2 The cell source

One binary parameter, `renewal` ($w$), states whether the cell that forms a
unit is spent. With $w = 1$ the compartment is restored by the cells it
produces and forming a unit does not draw it down. With $w = 0$ the cell is
consumed. Three regimes are defined by this parameter together with the
initial pool $P_0$ and the rates at which that pool divides and is lost,
$\rho$ and $\delta$:

| regime | $w$ | $P_0$ | $\rho$ | $\delta$ |
|---|---|---|---|---|
| liver-like, renewing source | 1 | 2600 | 0.10 | 0.10 |
| kidney-like, spent source with loss | 0 | 300 | 0.00 | 0.05 |
| pancreas-like, spent source without loss | 0 | 60 | 0.00 | 0.00 |

Four parameters differ and no others. The division rate $\rho$ belongs with
the loss rate: the liver-like source divides as fast as it is lost and holds
its pool, while neither of the other two divides at all. Section 2.2 treats
the pair as one axis for that reason, and reports that it is the axis that
most often decides which bound ends the growth. All runs use
$\nu = 8$, $q_{\mathrm{set}} = 0.6$, and the same vascular and geometric
settings.

### 4.3 The vascular tree and its boundary condition

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
every seed, which is why Section 2.7 sets $\phi$ once on an intact ensemble
and then varies $Q$. With $\xi = 1$ a unit's ceiling moves with its own mass,
so the identity above is approached rather than reached, and how closely is
measured in Section 2.7.

### 4.4 Readouts, replicates and robustness

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

### 4.5 Published measurements, the extension and the predictions

Three numbers of the model are set beside published measurements (Section
2.4). The seeded pool of each regime was reduced to a half and a quarter of
its reference value, and to the shares that three mouse experiments left after
removing progenitors with diphtheria toxin: 0.92 and 0.39 of pancreatic
progenitors and at most 0.35 of hepatic ones [2], and 0.6 of nephron
progenitors [14, 15]. Each reduced pool was grown over the same 24 seeds as
the reference, and mass and unit count at day 140 are reported as ratios of
means with a paired bootstrap interval. They are compared with organ weight
for pancreas and liver and with glomerular number for kidney, each as a share
of littermate controls, read from the figures and text of those papers. The
variation across seeds grown from identical parameters is compared with the
left-right variation of kidney volume in 48 people, 96 kidneys segmented from
clinical CT [16]: the standard deviation of the log ratio of the two sides,
divided by $\sqrt{2}$, after the mean difference between sides is removed.

The reference model has no path from the progenitor pool to the tree, and
Section 2.5 reports an extension that adds one. In it a source divides and is
not lost to attrition, a segment of the tree branches only while it holds two
caps' worth of progenitors and hands one to the new branch, a cap below its
size may divide faster by a set factor, and the run ends on a fixed day. On
one day each progenitor dies with probability 0.4, and the same seed grown
without that is the control. The extension is switched off in every other
result of this paper, and its settings are in `PAPER_NUMBERS.md`, An
extension.

The predictions of Section 2.9 use the same two procedures at doses and times
that were not published. The pancreas-like and liver-like constructs were
seeded with 0.75 and 0.60, and with 0.75 and 0.15, of their reference pools.
The extended source was run in the setting of the nucleation law in which
perfusion limits least, with 20 %, 40 % and 60 % of its progenitors killed on
day 45 of a window of 120, and with 40 % killed on days 20 and 80, each with
and without a doubled division rate in a depleted cap. They are stamped with
the commit that produced them (`PAPER_NUMBERS.md`, Predictions).

### 4.6 Code and data availability

All code, parameters and generated tables are at
`https://github.com/Danpc11/Ogrow`. `make paper` regenerates every number in
`PAPER_NUMBERS.md` from scratch and records the model version and the commit
that produced it. `make test` runs the conservation and consistency tests,
including exact progenitor balance in every realisation. The developmental
series are public under the accessions above.

The numbers in `PAPER_NUMBERS.md` are exact on the platform recorded there
(macOS, arm64). With the same package versions Linux gives different single
realisations and the same means within sampling error, which `make test`
checks against a reference ensemble. Every claim in this paper is about a mean
or a contrast over seeds, not about a single realisation.

## References

[1] Sexton Z A et al. 2025 Rapid model-guided design of organ-scale synthetic
vasculature for biomanufacturing *Science* **388** 1198–204

[2] Stanger B Z, Tanaka A J and Melton D A 2007 Organ size is limited by the
number of embryonic progenitor cells in the pancreas but not the liver
*Nature* **445** 886–91

[3] Murray C D 1926 The physiological principle of minimum work. I. The
vascular system and the cost of blood volume *Proc. Natl Acad. Sci. USA*
**12** 207–14

[4] Wang B, Zhao L, Fish M, Logan C Y and Nusse R 2015 Self-renewing diploid
Axin2+ cells fuel homeostatic renewal of the liver *Nature* **524** 180–5

[5] Perl A J, Schuh M P and Kopan R 2022 Regulation of nephron progenitor cell
lifespan and nephron endowment *Nat. Rev. Nephrol.* **18** 683–95

[6] Hochane M, van den Berg P R, Fan X, Bérenger-Currias N, Adegeest E,
Bialecka M, Nieveen M, Menschaart M, Chuva de Sousa Lopes S M and Semrau S
2019 Single-cell transcriptomics reveals gene expression dynamics of human
fetal kidney development *PLoS Biol.* **17** e3000152

[7] Lindström N O et al. 2021 Spatial transcriptional mapping of the human
nephrogenic program *Dev. Cell* **56** 2381–98

[8] Bastidas-Ponce A et al. 2019 Comprehensive single cell mRNA profiling
reveals a detailed roadmap for pancreatic endocrinogenesis *Development*
**146** dev173849

[9] Yang L et al. 2017 A single-cell transcriptomic analysis reveals precise
pathways and regulatory mechanisms underlying hepatoblast differentiation
*Hepatology* **66** 1387–401

[10] Qu X, Shi Y, Li L and Liu T 2026 Compensatory nephron hypertrophy:
mechanisms, pathophysiology dynamics, and clinical implications in acute kidney
injury and chronic kidney disease *Kidney Res. Clin. Pract.* advance online
publication, doi:10.23876/j.krcp.25.405

[11] Ryan D et al. 2018 Development of the human fetal kidney from mid to late
gestation in male and female infants *EBioMedicine* **27** 275–83

[12] Delone G V, Uryvaeva I V, Koretskiĭ V F and Brodskiĭ V Ya 1987 Analysis of
the postnatal growth of the mouse liver by estimating the hepatocyte count,
their weight and ploidy *Ontogenez* **18** 304–8 (in Russian, English abstract)

[13] Azizoglu D B, Perez K, Zheng S L, Rahman S, Rim E Y, Anbarchian T, Fish M,
Loh K M, Red-Horse K and Nusse R 2025 Vascularization of neonatal liver lobules
presages adult liver size *Nat. Commun.* **16** 9989

[14] Cebrian C, Asai N, D'Agati V and Costantini F 2014 The number of fetal
nephron progenitor cells limits ureteric branching and adult nephron endowment
*Cell Rep.* **7** 127–37

[15] Muthukrishnan S D, Ryzhov S, Karolak M and Oxburgh L 2018 Nephron
progenitor cell death elicits a limited compensatory response associated with
interstitial expansion in the neonatal kidney *Dis. Model. Mech.* **11**
dmm030544

[16] Ndzimbong W, Fourniol C, Themyr L, Thome N, Keeza Y, Sauer B, Piéchaud
P-T, Méjean A, Marescaux J, George D, Mutter D, Hostettler A and Collins T
2025 TRUSTED: the paired 3D transabdominal ultrasound and CT human data for
kidney segmentation and registration research *Sci. Data* **12** 615; data at
https://doi.org/10.6084/m9.figshare.27981050

[17] Martin F, Carmona M, Mahillo B, Alvarez M, Luengo A, Chatzixiros E,
López-Fraga M, Domínguez-Gil B and Tullius S G 2026 Organ donation and
transplantation worldwide: the Global Observatory on Donation and
Transplantation 2024 report *Transplantation* **110** e655–69

[18] Sljukic A, Jenkinson J G, Niksic A, Prior N and Huch M 2026 Advances in
liver and pancreas organoids: how far we have come and where we go next *Nat.
Rev. Gastroenterol. Hepatol.* **23** 44–64

## Notes for the authors, to be removed before submission

**Two manuscripts.** This one is for *Cell Systems* and `MANUSCRIPT.md` is the
base of the engineering paper for *Biofabrication*. They share the model, the
three regimes and Methods 4.1 to 4.4, and nothing else should appear in both:
the engineering paper keeps the scaling, the joint specification, the design
rules and the developmental series, and should cite this one for the model
and drop Sections 3.10 and 3.11 once this is public. Tell each journal about
the other manuscript at submission. `scripts/29_check_manuscript.py` holds
both to the same results.

**What is still the full manuscript's text.** Sections 2.1 to 2.3 and 2.6 to
2.8, the criteria and Methods are carried over with only their cross
references renumbered. Section 2.7 is long for this paper: the involution
sweep belongs in a supplement, with the thirteen robustness variants. The
reference list is the full manuscript's, unpruned and in its numbering, so
that every citation resolves; renumber by order of appearance at formatting.

**The passages a script writes.** The clause of the Summary and the paragraphs
of Sections 2.4, 2.5 and 2.9 between `auto:` markers are written by
`scripts/36_validation.py`, `37_kidney_coupling.py` and `38_predictions.py`.
Run them and commit this file before `make paper`. Read what they wrote before
trusting the paragraphs after them: in Section 2.5 the interval where
perfusion limits least only touches the measured range at its lower end, and
the text says "approaches" for that reason.

**Before submission.** Send a presubmission inquiry with the Summary. Confirm
the two liver values of [2] against its Fig. 6b and supplement. Say in Section
2.4 whether the regime parameters were chosen before anyone in the group had
read [2]. Cell Press items not done: the eTOC paragraph, the key resources
table, and Methods in STAR form. The Summary is at the 150-word limit once its
clause is filled.

## Figures

Figure 1, the model and its three bounds (`Fig1_model.pdf`). Figure 2, the
regimes at arrest (`Fig2_bounds.pdf`). Figure 3, three numbers beside
measurements (`Fig10_validation.pdf`). Figure 4, the extension
(`Fig11_coupling.pdf`). Figure 5, restart and injury (`Fig8_injury.pdf`).
Figure 6, predictions (`Fig12_predictions.pdf`). Legends are in
`FIGURES_A.md`.
