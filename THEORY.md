# Theory

Every symbol used in the code is defined here, in the order the model builds
them up. Equations are numbered so the scripts can refer to them.

---

## 1. Scope

The framework covers **exchange organs**: organs whose function is to pass
something across an epithelial–endothelial barrier, and whose functional unit
is an element of that barrier. Four qualify: liver, kidney, pancreas and lung.

The heart does not. It is mesenchyme and fibre performing mechanical work; its
unit is not a barrier element and its limiting variable is tension per fibre,
not surface per flow. Including it as a degenerate case would describe it
wrongly, so it is excluded by construction rather than by lack of data.

---

## 2. State variables

| Symbol | Meaning | Units |
| --- | --- | --- |
| $N(t)$ | number of functional units | count |
| $m_j(t)$ | mass of unit $j$ | arbitrary mass |
| $M(t)$ | organ mass | arbitrary mass |
| $P_i(t)$ | progenitor pool at tip $i$ | count |
| $\mathbb{C}(t)$ | effective supply conductance of the organ | as $w_k$ |
| $c_u(t)$ | supply per unit, $\mathbb{C}/N$ | as $w_k$ |
| $\phi(\mathbf{x})$ | signed distance to the tissue boundary | length |

**Equation 1 — two fields.**

$$M(t) = N(t)\,\overline{m}(t), \qquad \overline{m} = \frac{1}{N}\sum_j m_j$$

The separation into a count and a size is forced by data, not chosen for
convenience: a single density field fits the mass curve and misses the unit
count by a factor of six. The two fields also close at different times. $N$
freezes first, $\overline m$ keeps rising for weeks, and $M$ later still.

---

## 3. The functional unit

The unit is **not** a diffusion-limited object. The oxygen diffusion length
sets the capillary spacing (~10 µm), not the lobule (~250 µm) or the nephron
(1–2 mm).

> The functional unit is the territory of one terminal branch of the **primary
> tree**, and which tree is primary depends on the organ.

| Organ | Primary tree | Unit |
| --- | --- | --- |
| Liver | vascular | lobule |
| Kidney | epithelial (ureteric bud) | nephron |
| Pancreas | epithelial (ductal) | acinus, islet |
| Lung | epithelial (bronchial) | acinus |

In the three epithelial-tree organs the vasculature **modulates** rather than
organises: epithelium branches without endothelial cells, vascular ablation
changes stereotypy but not branching rate, and heterotypic recombination shows
the mesenchyme reprogramming the pattern. The hierarchy is epithelium →
mesenchyme → vasculature.

The liver is the one organ where the organising tree **is** the supply tree,
which is why the supply bound bites there directly.

The lung is the one organ where the organising tree carries the **exchange
medium** (air) rather than the supply, so a second tree runs alongside it. In
the code this is the parameter `vessels`: `"demand"` for the other three,
`"bundle"` for the lung.

---

## 4. Supply: flow, shear and calibre on the tree

Sections 4 to 8 describe the open engine, `src/organgrow/engine.py`, which
produces the figures. Where the closed model behaves differently it is said in
place and developed in §13.

A tree has analytic resistance, so no linear system is needed. Demand is
accumulated from the leaves to the root in one pass and conductance combined
in another, both at cost $\mathcal{O}(E)$ in the number of edges
(`hemodynamics.py`).

**Equation 2 — Poiseuille conductance.**

$$w_k = \frac{r_k^{4}}{\ell_k}$$

**Equation 3 — series-parallel reduction.** Walking leaves to root, each
segment's downstream conductance is its own in series with the parallel sum of
its daughters',

$$\mathbb{C}_k = \left(\frac{1}{w_k} + \frac{1}{\sum_{c\,\in\,\text{daughters}(k)} \mathbb{C}_c}\right)^{-1}$$

and $\mathbb{C}_0$ is the organ's distal conductance. The root is in series
with it,

$$\mathbb{C}(t) = \left(\frac{1}{C_{\text{root}}(t)} + \frac{1}{\beta\,\mathbb{C}_0(t)}\right)^{-1}$$

so the root caps the total however wide the distal tree becomes.

**Equation 4 — the root grows with the body.**

$$C_{\text{root}}(t) = C_0\left(1 + \frac{t}{\tau_b}\right)^{\alpha}, \qquad \alpha = 1.155$$

The exponent is the measured developmental allometry of organ mass against
body mass, not the 3/4 of metabolic scaling.

**Equation 5 — flow and wall shear.** Demand drawn at each tip is the mass of
the units nearest to it; flows accumulate to the root. For Poiseuille flow,

$$\tau_k \;\propto\; \frac{f_k}{r_k^{3}}$$

**Equation 6 — calibre relaxation toward a shear setpoint.** The radius that
would put a segment at its setpoint is $(f_k/\tau_{\text{set}})^{1/3}$, and
each calibre moves a fixed fraction $\eta$ of the way there per day,

$$r_k \leftarrow r_k + \eta\left[\left(\frac{f_k}{\tau_{\text{set}}}\right)^{1/3} - r_k\right]$$

Murray's law is **not imposed**. It emerges: at the fixed point every segment
carries the same shear, so $r^3 \propto f$, and since flow is conserved at a
bifurcation, $r_p^{3} = \sum_c r_c^{3}$. Measured over 52 bifurcations of a synthetic tree relaxed to the fixed point
under uniform demand, the ratio $r_p^{3}/\sum_c r_c^{3}$ is 1.00. That is a
check that the implementation reaches its own fixed point, not a validation
against measured calibres: at a uniform setpoint the relation is an algebraic
consequence of the rule.

Murray's exponent 3 is the special case b = 1 of a more general law (§13).
With a wall cost exponent b, the branching exponent is b + 2, so b between 0.7
and 1.0 gives 2.7 to 3.0, the range measured in real trees.

**Equation 7 — supply per unit.**

$$c_u(t) = \frac{\mathbb{C}(t)}{N(t)}$$

**Equation 8 — what reaches one unit.** Supply falls with distance $d_j$ to
the nearest vessel and rises with that vessel's calibre,

$$c_{u,j} = c_u \, e^{-d_j/\ell_a} \sqrt{r_{v(j)}/r_{\text{root}}}$$

This holds for the open engine. In the closed model supply per unit is read
from the outlet pressure of the terminal branch that feeds the unit,
$q = p_{\text{out}}/\Delta P$, with no proximity term.

**Equation 9 — homogenised shear.** For a population of calibres,

$$r_{\text{eff}} = \frac{\langle r^{5}\rangle}{\langle r^{4}\rangle}$$

exact, with no closure assumption, but dependent on the whole distribution and
not its mean. For a log-normal of width $\sigma$, $r_{\text{eff}} =
r_m e^{9\sigma^{2}/2}$: at $\sigma = 1$ the factor is 90, so the width
dominates any exponent.

### A critical radius, below which a vessel cannot hold

Shear on a segment in series with the tissue it feeds is not monotone in the
radius. With $Q = \Delta P / (\ell/r^4 + R_t)$ and $\tau = Q/r^3$, shear rises
while the segment's own resistance dominates and falls once the tissue's does,
so that $\tau \propto r/(\ell + R_t r^4)$, whose derivative vanishes where
$\ell = 3R_t r^4$, that is at

$$r_c = \left(\frac{\ell}{3R_t}\right)^{1/4}$$

(An earlier version wrote $(3\ell/R_t)^{1/4}$, $\sqrt{3}$ times too large, and
the closed model used it as the floor of every calibre; all closed-model
results were recomputed after the correction.)

Below $r_c$ the shear feedback is positive and the vessel closes; above it the
feedback is negative and the vessel holds. The criterion is derived from the
adaptation rule. The closed model, however, uses $r_c$ as the floor of the
adaptation, so the absence of calibres below $r_c$ in simulated trees is
imposed, not predicted. For real trees the derivation predicts a lower
cut-off in the calibre distribution rather than a tail running to zero, which
is a test on measured calibres, not a result of the model.

### How supply scales with unit count

Write $\mathbb{C} \sim N^{\gamma}$. Since $c_u = \mathbb{C}/N \sim
N^{\gamma-1}$, the sign of $\gamma - 1$ decides whether adding units raises or
lowers the supply each one gets, and so whether a capacity-limited window can
close. Measured along grown trajectories:
$\gamma$ = 0.79 (liver).
$\gamma$ = 0.42 (kidney).
$\gamma$ = 0.60 (pancreas).
$\gamma$ = 0.50 (lung).
All are below 1, so the feedback is stabilising and no runaway occurs. They are
also well below the 4/3 that uniform-shear adaptation alone would give, because
the root is in series and caps the total. That cap is a property of the open
engine: in the closed model the root takes under 3 % of the pressure drop and
its radius follows the shear law, growing as $Q^{1/3}$.

### The tree adds little disorder of its own

With calibres frozen, replicates that reach the same unit count differ in
conductance by a CV of order $10^{-4}$: the topology of the tree carries
almost no disorder. Letting the calibres relax raises the dispersion, but it
stays below about 0.01 for relaxation rates up to 0.4 per day. At the
fastest rate tested the calibre loop is at the edge of stability and the
dispersion becomes seed-sensitive, varying by an order of magnitude between
replicate sets; `results/vascular_disorder.json` carries a bootstrap interval
for every rate and `scripts/11_vascular_disorder.py` prints it. An earlier
version of this engine reported a resonant peak near CV = 0.25; that was an
artefact of three errors since fixed (integer truncation in nucleation, a
pressure field clamped at zero, and a circular per-unit supply term), and it
does not survive their correction.

The vasculature is therefore a conduit here, not a source of quenched
disorder, at least at the scale this engine resolves and with flow driven by
demand. In the closed model the opposite held: the topology of the tree was
the main source of spread between replicates, through where its root entered
and how it branched (§13).

## 5. The fabrication window and how it closes

**Organ identity is one parameter: the closure mode.** (This is the open
engine's formulation. In the closed model, §13, identity is which of two local
bounds, supply per unit or the progenitor pool, binds first, and the three
modes below are its limits.) Everything else in the
engine is shared.

> **Open point.** The closure-mode table below assigns the liver to signal
> decay (venous sprouting, about one postnatal week), while this equation
> replaces that window with capacity. The engine shows the two are not
> interchangeable: see `KNOWN_LIMITATIONS.md`.

**Equation 10 — capacity saturation (liver).**

$$\dot N = \nu\left(\frac{\mathbb{C}}{N} - c_{\text{set}}\right)_{+}$$

A unit nucleates where there is **spare supply per unit**. $\nu$ is the
nucleation efficiency, $c_{\text{set}}$ the per-unit perfusion setpoint, and
$(\cdot)_+$ the positive part.

Two properties follow without being imposed. Multiplying by $N$,

$$\frac{\mathrm{d}(N^{2})}{\mathrm{d}t} = 2\nu\left(\mathbb{C} - c_{\text{set}}N\right)$$

so in the early regime, with $N \ll \mathbb{C}/c_{\text{set}}$ and $\mathbb{C}$
nearly constant, $N \propto t^{1/2}$ and $\dot N \propto t^{-1/2}$: a power law
with exponent $-1/2$, derived rather than assumed. And nucleation stops when
$N$ reaches $\mathbb{C}/c_{\text{set}}$, so $N_\infty$ is finite **by
mechanism**. An imposed window does not converge: with $\sigma \sim t^{-p}$
the integral diverges at $t\to\infty$ for $p<1$ and at $t\to0$ for $p>1$; an
exponential converges by decree.

Writing the drive as $\mathbb{C}/N - c_{\text{set}}$ rather than
$\mathbb{C} - c_{\text{set}}N$ is not cosmetic. The second form is a linear
relaxation whose variance behaviour differs, and the two are distinguished by
the variance test of §7, not by any fit to the mass curve.

**Equation 11 — progenitor exhaustion (kidney, lung).**

$$\dot N_i = \iota P_i, \qquad \dot P_i = -(\iota + \delta - \rho)P_i$$

with $\rho$ self-renewal, $\iota$ induction and $\delta$ loss. The window
closes when $P_i \to 0$. The closure time is an **output**, not a parameter,
and it is asynchronous between individuals because each starts from its own
$P_i(0)$.

**Equation 12 — limiting pool (pancreas).** A fixed budget set aside at the
start is drawn down; the final count is fixed before growth begins.

---

## 6. Local rules

**Equation 13 — spatial exclusion.** A candidate unit at $\mathbf{x}$ is
accepted only if

$$\phi(\mathbf{x}) < 0 \quad\text{and}\quad \min_j \lVert \mathbf{x} - \mathbf{x}_j \rVert > \lambda$$

$\lambda$ is the minimum unit spacing. The first condition confines growth to
the tissue domain, so the shape of the organ influences the architecture that
emerges.

**Equation 14 — unit growth, limited by local supply.**

$$\dot m_j = g(t)\,m_j\left(1 - \frac{m_j}{\kappa\,c_{u,j}}\right), \qquad g(t) = g_0 e^{-t/\tau_g}$$

Each unit is bounded by the supply reaching **it**, not by a global ceiling.

---

## 7. Quenched disorder and what it does to variance

Each individual is one **realisation** of a disordered system. The disorder
lives in the parameters, not in the initial state: in a disordered system what
varies between realisations are the couplings.

Two sources are kept separate, because they behave differently.

| Source | Symbol | What it is | Can feedback correct it? |
| --- | --- | --- | --- |
| Initial condition | $\sigma_{\text{init}}$ | spread of seeded units | yes |
| Parameter | $\sigma_{\text{param}}$ | per-individual $c_{\text{set}}$, $P_0$ | no |

**Equation 15: the dilution law.** If the increment of $N^{2}$ is the same
for every individual, the dispersion of $N^{2}$ is preserved and

$$\frac{\mathrm{CV}_i}{\mathrm{CV}_f} = \left(\frac{N_f}{N_i}\right)^{2}$$

Multiplying the nucleation law by $2N$ gives
$\mathrm{d}(N^2)/\mathrm{d}t = 2\nu(\mathbb{C} - c_{\text{set}}N)$, so the
increment is equal across individuals only while $c_{\text{set}}N$ is
negligible. Once it is not, it acts as a restoring term: an individual ahead
of the others gains less. Contraction then **exceeds** the dilution value; the
open engine contracts three to four times more than it. Equation 15 is the
weak-feedback limit, not an upper bound. It fits the one measurement
available, mouse liver, where $N$ grows 1.82-fold and the CV contracts
3.09-fold against 3.31 predicted. The general form is in §13. With
per-individual setpoints,

**Equation 16 — the floor on terminal variance.**

$$\mathrm{CV}_\infty \gtrsim \mathrm{CV}(c_{\text{set}})$$

Each individual converges to its own ceiling, so terminal variance cannot fall
below the parameter's own spread. This says **when the bound fails**, which the
earlier equality form did not.

The discriminating measurement is the **absolute** dispersion $\mathrm{sd}(N)$,
not the CV. A falling CV is ambiguous: it falls both when variance is destroyed
and when it is merely diluted by a growing $N$. Capacity saturation destroys
$\mathrm{sd}(N)$; exhaustion and limiting pool conserve it exactly. In the
open engine that last statement holds by construction, because those two modes
read nothing that growth changes. The closed model runs all three as one law
(§13), so there the answer is a measurement rather than an input.

---

## 8. Self-averaging: which properties are specifiable

An observable **self-averages** if its dispersion between realisations decays
with system size,

**Equation 17.**

$$\mathrm{CV}(N) \propto N^{-a}, \qquad a > 0$$

Quantities that are sums over the volume self-average. Properties of the
particular realisation do not: in this engine the nearest-vessel identity has
$a \approx 0$, and the reported disagreement between real vascular networks is
81–87 %.

This is the criterion that decides **what to score**, and the reason topology
is excluded from the invariant vector. Comparing topology between a construct
and a native organ measures contingency, not fidelity.

### A caution on self-normalised distances

A distance in invariant space normalised by the cohort's own spread,

$$\varepsilon_{ij} = \sqrt{\left\langle\left(\frac{v_i - v_j}{\sigma}\right)^{2}\right\rangle}$$

has mean square $\langle\varepsilon_{ij}^2\rangle = 2$ for **any**
exchangeable sample with the cohort's own normalisation, whatever the cohort,
and so, by Jensen's inequality, a mean of at most $\sqrt{2}$. Its value is therefore a property of the normalisation, not a
measurement; random vectors give the same number. It is useful only as a test
of exchangeability, which is how `scripts/06_exchangeability.py` reports it.

What does carry information are quantities measured in a **fixed** metric: raw
dispersions, the within-to-between ratio, and the annotator disagreement.

---

## 9. Structural identifiability

Before fitting, ask what is determinable. With $J$ the Jacobian of the
log-observables with respect to the log-parameters,

**Equation 18.**

$$J_{ij} = \frac{\partial \log y_i}{\partial \log \theta_j}, \qquad
\text{rank}(J) \text{ and } \mathrm{cond}(J) = \frac{s_1}{s_{\text{rank}}}$$

A rank below the number of parameters means flat directions: those parameters
are not recoverable from those observables at any sample size. The condition
number sets how much noise is amplified.

The number of replicates needed to pin a parameter to $\pm\epsilon$ follows
from the Fisher information $J^{\mathsf{T}}J/\varsigma^{2}$, with $\varsigma$
the relative noise on each observable.

---

## 10. Three sources of variation in a real organ

Measured, not postulated, in paired human kidneys.

| Source | Governs | How it is seen |
| --- | --- | --- |
| **Systemic** | size | variance appears *between* donors |
| **Positional** | mean shape | a systematic left–right difference *within* a donor |
| **Contingent** | the remainder of shape | the residual within-donor scatter |

The consequence for building organs: mass is specifiable, shape is set by the
surroundings, and the rest is contingent. A construct can reach its mass target
and still not match native proportions, and that is not a process failure.

---

## 11. Symbols at a glance

| Symbol | Name | Where |
| --- | --- | --- |
| $\nu$ | nucleation efficiency | Eq. 6 |
| $c_{\text{set}}$ | per-unit perfusion setpoint | Eq. 6 |
| $\iota,\ \rho,\ \delta$ | induction, self-renewal, loss | Eq. 7 |
| $\lambda$ | minimum unit spacing | Eq. 9 |
| $\beta$ | demand coefficient | Eq. 2 |
| $C_0$ | root conductance at $t=0$ | Eq. 3 |
| $\alpha$ | developmental allometric exponent | Eq. 3 |
| $\tau_b$ | body growth time constant | Eq. 3 |
| $g_0,\ \tau_g$ | unit growth rate and its decay | Eq. 10 |
| $\kappa$ | mass supported per unit of supply | Eq. 10 |
| $\ell_a$ | supply attenuation length | Eq. 8 |
| $\tau_{set}$ | wall-shear setpoint | Eq. 6 |
| $\eta$ | fraction of the calibre gap closed per day | Eq. 6 |
| $r_k,\ \ell_k$ | radius and length of segment $k$ | Eq. 2 |
| $f_k,\ \tau_k$ | flow and wall shear on segment $k$ | Eq. 5 |
| $\gamma$ | exponent in $\mathbb{C}\sim N^{\gamma}$ | §4 |
| $\sigma_{\text{init}},\ \sigma_{\text{param}}$ | the two disorder widths | §7 |
| $a$ | self-averaging exponent | Eq. 15 |

All values are in `config/params.yaml`.

---

## 12. What the framework does not claim

It is not a theory of development. It does not explain organ specification,
the particular branching pattern, or macroscopic form; those are separate
problems that couple to this one rather than follow from it. The mesenchyme
appears as the instructor of the organ-specific pattern but is not formalised.

Flow is solved on a tree, not on a network: anastomoses, pulsatility and
non-Newtonian rheology are absent. Pressure is implicit in the open engine and
explicit in the closed model.

The closed model carries two latent parameters, the resistance of one unit and
the number of units a terminal branch serves. Their defaults reproduce the
reference organ; what they correspond to physically is to be decided by the
fit, and a latent earns an interpretation only if the data constrain it and it
predicts something it was not fitted to.

Further limits, including which parts rest on which evidence, are in
`KNOWN_LIMITATIONS.md`.

---

## 13. What the closed model establishes, and what it withdrew

*This section is the record of how the closed model was built and corrected,
including results later withdrawn; numbers in it come from the version
current when each was written. The numbers the paper uses come only from
`PAPER_NUMBERS.md`, generated by `make paper` with the model version and
commit.*

Sections 4 to 8 describe the open engine that produces the figures. The closed
model (`src/organgrow/closed.py`) replaces its three open loops with one law,
and where the two disagree this section is the current statement.

### One law, two bounds

At every tip, with supply per unit $q_i$ and progenitor pool $P_i$,

$$\lambda_i = \nu\,P_i\,(q_i - q_{\text{set}})_+, \qquad k_i \sim \mathrm{Poisson}(\lambda_i), \qquad \dot P_i = (\rho - \delta)\,P_i - k_i$$

A growing organ meets one of two ceilings: supply per unit falls to
$q_{\text{set}}$, or the pool runs out. Which comes first is set by $\nu$ and
$\rho - \delta$, so organ identity is a position on that plane rather than a
branch in the code. Nucleation is a Poisson draw, so the demographic floor on
the spread is part of the law rather than added afterwards.

### The cost exponent sets three exponents

With a wall cost exponent $b$,

$$r_p^{\,b+2} = \sum_c r_c^{\,b+2}, \qquad \tau_{\text{set}} \propto p^{\,b/2}, \qquad \frac{r_c}{r_p} = 2^{-1/(b+2)}$$

The three follow from one objective per unit length,
$\mathcal{C} = Q^2/r^4 + k\,(p\,r^2)^b$: the first term is viscous dissipation,
the second the metabolic cost of the wall, whose cross-section scales as
$p\,r^2$ by Laplace's law (thickness $\propto p\,r$), raised to $b$.
Minimising in $r$ gives $Q^2 \propto k\,p^b r^{2b+4}$, so $Q \propto r^{b+2}$
at fixed pressure (the branching exponent $b + 2$) and wall shear
$\tau \propto Q/r^3 \propto p^{b/2} r^{b-1}$. The model keeps the pressure
factor and drops $r^{b-1}$, which at $b = 0.9$ is $r^{-0.1}$. In this
notation, cost $\propto r^{\alpha}$ with $\alpha = 2b$. Murray's $b = 1$ is
the case where cost is blood volume. $b = 2/3$ and $3/4$
give branching exponents 2.67 and 2.75; $b \approx 0.9$, appropriate where
flow is pulsatile, gives 2.9. Measured exponents of 2.7 to 3.0 sit in that
range. The engine uses $b = 0.9$.

### A branch buys reach, not conductance

Two daughters at the cost-law ratio conduct $2 \cdot 2^{-4/(b+2)}$ of their
parent, which is 0.77 at $b = 0.9$ and below one for every $b < 2$, which includes the
measured range. Dividing a
tip loses conductance whenever $b < 2$. What it gains is the tissue it can serve, so a
split is worth making when the tissue it relieves exceeds the tissue its
junction occupies. Vessel volume alone is flat in the branching number, and
without the junction term the tree grows until an implementation cap stops
it; with it, the tree stops on its own.

### Contraction of variance, general form

The spread of unit number relaxes from its initial value toward a floor,

$$\sigma^{2}(t) = \sigma_0^{2}\,e^{-2\lambda t} + \sigma_\infty^{2}\left(1 - e^{-2\lambda t}\right)$$

which fits the open engine with $R^2$ between 0.92 and 0.99. $\lambda$ is set
by how strongly the closure reads supply and $\sigma_\infty$ by what keeps
injecting disorder. Equation 15 is its limit when $\lambda$ is small, which is
the regime the mouse liver occupies.

### Where the between-individual spread came from

Two features of the implementation, not of the biology, dominated the spread
in the first versions of the closed model, and both are now fixed.

The vascular root was placed at a random point of the domain. A root on the
boundary leaves most nucleation candidates outside the organ, so the first
branch never pays for its junction and the tree stays at one vessel; a root
deep inside branches explosively. Root depth predicted final vessel count with
Spearman 0.72 ($P = 0.045$, eight seeds). In an organ the root enters at a
fixed site, the hilum, and the model now places it at the deepest point of the
domain.

Sprouting for lack of room counted refused placements, and that count grows
with the progenitor pool, so a pool of a few hundred made hundreds of vessels
in a day and every tree reached the vessel cap. It now responds to the share
of placements refused.

With both fixed, the spread of unit number across replicates fell from 0.255
to 0.124 and that of vessel number from 0.62 to 0.079.

A third error made the tree a star. A unit stayed attached to the tip where it
was born, so when that tip branched to relieve its load, the load never moved:
the parent kept every unit, kept passing the junction test and kept sprouting.
The root ended with 270 daughters and 409 of 492 vessels were tips one
generation from it. Each unit is now fed by the nearest terminal branch, which
is what makes a branch relieve the tissue that justified it. The largest node
now has 20 to 42 daughters, and the trees again reach the vessel cap, so the
cap is again an implementation limit and not physics. An early-history lock-in, of the kind a Polya urn produces above its
threshold, looked like the explanation before the root was fixed; it was the
root.

A fourth defect let a tip make more units than it held progenitors: the
number of units was drawn as Poisson in $\nu P_i (q_i - q_{\text{set}})$
without a cap, and the pool was clamped at zero afterwards. A pool-bound
organ exceeded its own budget by up to 40 %, by an amount that grew with the
drive. A tip now makes at most the progenitors it holds, and the fractional
part of a pool is resolved once, so that it cannot retry until it succeeds.

### The organ does not fill its domain

Scaling the domain with every extensive quantity scaled as $s^3$ gave
$N \propto V^{0.642 \pm 0.028}$ over a 33-fold range in volume, and the
fraction of the progenitor pool left rose from 0.04 to 0.69 with size. Neither
is a scaling law. Measured directly, 61 % of the domain is more than two unit
spacings from any unit in the smallest organ and 87 % in the largest: the
organ grows as a colony around its root and never reaches its walls. Supply is
not what stops it (supply per unit stays near 0.85, above $q_{\text{set}}$),
nor the pool, nor time (doubling the run adds 5 %).

Nothing in the closed model drives growth into unoccupied space. Units appear
only next to tips, and tips sprout toward existing tissue, a short step away
from crowding, or toward tissue short of supply, which nowhere is. This is the
gap §12 names: the mesenchyme instructs the pattern and is not formalised. A
formulation exists outside biology, space colonisation (Runions and
colleagues, for leaf venation): attractor points in unserved tissue guide the
nearest tips and are removed when reached, which fills a domain and gives a
finite tree. The open engine already used it; the closed model dropped it.

It is now in the closed model. Attractors spread at a fixed density through
the domain pull the nearest tip within reach and are removed once a vessel or
a unit comes close; a tip that extends toward them carries its progenitor
niche with it, so the front does not starve. With a pool large enough not to
bind, the unfilled fraction of the smaller domain fell from 54 % to 24 %, the
tree stopped on its own at 1342 vessels against a cap of 5999, and supply per
unit fell to 0.57, below $q_{\text{set}}$: the liver reached its supply bound
for the first time in the closed model. In a domain 4.6 times larger the same
bound binds with 58 % of the domain still empty, which is the expected
behaviour of a supply-limited organ: it grows to what its perfusion sustains,
not to what its container holds.

Two conditions came out of getting there. The liver regime needs a pool that
does not bind, about four progenitors per unit made; with fewer, the
progenitor bound binds first and the organ behaves as a kidney. And the
junction test has to stay on the deficit path: removing it let every
supply-short tip sprout every day, and the tree ran to its cap.

### Self-averaging, and why it is not yet measured

Whether unit number self-averages is decided by how $R = \mathrm{sd}/\mathrm{mean}$
behaves as size grows with nothing else changing: $N^{-1/2}$ for
self-averaging, a constant under quenched disorder. `scripts/14_size_sweep.py`
checks four controls first and reports no exponent if any drifts. In the
current model the supply per unit and the median radius hold across sizes,
but the fraction of the pool left rises from 0.04 to 0.69. That drift comes
from the unfilled domain, not from the pool itself, and no exponent is
reported until organs fill their domains.

### Calibration against the mouse liver

The mouse target is a contraction of 3.09 in the CV of lobule number while
the count grows 1.82-fold. With only demographic noise there is little to
contract: the model starts the interval near its floor. Livers grown at one
rate and counted with an offset in maturity that reproduces the mouse's
initial spread (CV 0.228 against 0.229) close at 0.108, a contraction of 2.1.
The mouse closes at 0.074, below the model's floor, so capacity closure in
this model contracts the spread in the right direction but not enough.

An earlier version reported that capacity erased a spread in maturity but
not one in growth rate, with each individual's final count tracking its rate,
and proposed that as the test between capacity and signal. That difference
came from a defect: a tip could make more units than it had progenitors, by an
amount that grew with its nucleation rate, so fast individuals overshot. With
the defect fixed, a spread in rate ($\sigma(\log\nu) = 1$) contracts from
0.141 to 0.085 and the final count no longer tracks the rate (Spearman 0.33,
$P$ = 0.42). The question of §5, capacity or signal for the liver, is open
again; what the model now says is that capacity alone leaves a floor above
the mouse's.

### Which individual spread reaches the final count

The open engine answered which closure destroys variance by construction. In
the closed model it is measured: each replicate gets its own progenitor pool
or its own supply setpoint, both with a spread of 0.30 in log, and the
transmission is the part of the final CV of $N$ above the demographic
baseline, divided by the input spread. 24 replicates per condition, 95 %
bootstrap intervals.

| Regime | Bound | Pool spread passes | Setpoint spread passes |
|---|---|---|---|
| Liver | supply | 0.23 [0.00, 0.42] | 1.42 [0.79, 2.10] |
| Kidney | progenitors, with loss | 0.59 [0.41, 0.72] | 1.28 [0.69, 1.88] |
| Pancreas | progenitors, no loss | 1.12 [0.85, 1.33] | 1.08 [0.32, 1.69] |

The liver filters the spread in its pool, because once supply binds, more
progenitors cannot make more units; the budget passes it whole. No regime
filters the setpoint. Supply and progenitors are both local to a tip, and a
tip converts its pool only where its supply exceeds the setpoint, so the
setpoint decides how much of the pool becomes tissue in every regime.

The closure mode therefore decides which spread is removed, not which is
kept: a supply-bound organ removes spread in its progenitors, and nothing
removes spread in the supply setpoint.

Two measurement lessons came with this. These numbers are after the fix of
the progenitor overdraft described above, which had raised the setpoint
transmission of the pool-bound regimes. And an input spread of 0.15 in the
setpoint was too small against the liver's demographic baseline (CV 0.091):
its interval reached zero, which looked like filtering and was noise.

### Elasticities, and what sets the setpoint's

A comparison of spreads between cohorts is dominated, for large input spreads,
by the individuals at the tail. The local measure is the elasticity of the
final count, $\varepsilon_x = \partial \ln N / \partial \ln x$, estimated with
the same seed at $x e^{\pm 0.3}$ (`scripts/19_elasticity.py`, six paired
replicates, mean ± s.e.).

| Regime | Setting | $\varepsilon_{P_0}$ | $\varepsilon_{q_{\text{set}}}$ |
|---|---|---|---|
| Liver | reference | +0.38 ± 0.13 | −1.19 ± 0.08 |
| Kidney | $q_{\text{set}}$ 0.5 / 0.6 / 0.7 | +0.59 / +0.55 / +0.51 | −0.40 / −0.91 / −2.95 |
| Kidney | $\nu$ 4 / 16 | +0.62 / +0.53 | −1.20 / −0.75 |
| Pancreas | $q_{\text{set}}$ 0.5 / 0.6 / 0.7 | +0.92 / +0.86 / +0.87 | −0.13 / −0.36 / −2.55 |
| Pancreas | $\nu$ 4 / 16 | +0.68 / +0.83 | −0.53 / −0.29 |

The pool elasticity is stable within a regime and rises from the liver to the
budget. The setpoint elasticity is not a property of the regime: it depends on
where the setpoint sits relative to the supply at the tips.

The liver was then measured at the same five settings (six paired replicates
each): $\varepsilon_{P_0}$ = +0.38, +0.39, +0.35, +0.32, +0.36 at reference,
$\nu$ 4, $\nu$ 16, $q_{\text{set}}$ 0.5 and 0.7, and $\varepsilon_{q_{\text{set}}}$ =
−1.19, −1.42, −0.94, −0.41, −3.03. At every setting the pool elasticity is
ordered liver < kidney < budget.

**One mechanism sets the setpoint elasticity wherever pool is left over.**
Lowering the setpoint by $\mathrm{d}q$ activates the tips whose supply sits
just below it, which convert what they hold, so

$$\varepsilon_{q_{\text{set}}} = -\,q_{\text{set}}\, w(q_{\text{set}}) / N,$$

with $w$ the pool held per unit of supply just below the threshold, read from
the final state ($\Delta q = 0.05$) with no fitted parameter. Across the five
liver settings the prediction follows the measurement with correlation 0.98
(mean absolute error 0.33): −1.04 against −1.19 at reference, −3.13 against
−3.03 at the highest setpoint. In the budget it holds at three of five
settings (−0.44 against −0.36 at reference) and misses where the stranded
pool sits far below the threshold, where a ±30 % step is no longer local. So
the liver does not respond to its setpoint through a global capacity: its
size moves with the setpoint because of the pool waiting tip by tip just below
it. Where the leftover pool sits differs between organs: in the budget
all of it lies below the threshold, 42 % in the slice just under it; in the
liver 85 % sits at tips whose supply is above the threshold, and only 5 % in
that slice. Supply does not hold the liver's remaining pool, and the 5 % at the
threshold is what its elasticity reads. In the kidney the formula gives zero and the measurement does not, as it
should: loss spends the kidney's pool, nothing is left stranded, and the
setpoint acts through the race between converting a progenitor and losing it,
a fraction $\nu D/(\nu D + \delta)$ per tip with $D = q - q_{\text{set}}$; that
derivation is written, not tested.

**This held before the critical radius was corrected, and does not after.**
The formula is the direct, frozen-state part of the response: it counts the
pool released when the threshold drops, with the tree, the pressures and the
assignment of units held fixed. The full derivative adds
$\partial N_f/\partial \mathbf{z}\;\mathrm{d}\mathbf{z}/\mathrm{d}q_{\text{set}}$
for the vascular and tissue state $\mathbf{z}$: units made from released pool
lower the supply per unit and strand part of it again. With the corrected
critical radius the liver holds 37 % of its pool below the threshold, the
feedback dominates, and the direct term predicts −4 to −9 against measured
−1.1 to −3.4 (correlation −0.15 over the five settings). The claim is
withdrawn; what survives is that the setpoint elasticity is not a property of
the regime.

A global derivation for the supply-bound side was tried and failed. With
total root flow $Q \propto N^\gamma$ at closure, $\varepsilon_{q_{\text{set}}}$
would be $-1/(1-\gamma)$. With $\gamma$ read from the last 15 % of growth, a
window fixed before testing, it predicts about −1.3 at every liver setting
against measurements from −0.41 to −3.03 (correlation −0.60). Its agreement at
the reference setting alone was a coincidence.

### The human partition, and what the model can and cannot predict

If each kidney's log size is a systemic deviation shared by the pair plus a
local one and a contingent one, the fraction of size variance between
individuals is
$(\varepsilon_q^2\sigma_{q,s}^2 + \varepsilon_P^2\sigma_{P,s}^2) /
(\varepsilon_q^2\sigma_{q,s}^2 + \varepsilon_P^2\sigma_{P,s}^2 +
\varepsilon_q^2\sigma_{q,\ell}^2 + \varepsilon_P^2\sigma_{P,\ell}^2 + \sigma_c^2)$,
with $s$ systemic, $\ell$ local and $\sigma_c$ the contingent floor. The model
supplies the weights; the input spreads must come from elsewhere. TRUSTED
gives 0.887 for volume, so the systemic part must be 2.8 times the
within-individual part in standard deviation. That is a requirement, not a
prediction, until the input spreads are measured independently: the local
one from paired measures such as split renal function, the systemic ones from
paired nephron counts. The model's contingent floor (0.061 in the kidney
regime) cannot stand in for the human one, because it is measured in organs
of about 130 units and whether it falls with organ size is not known.

### How each organ closes, read from where its pool is

A mean supply per unit hides mixed states. The model now reports, at each day,
the share of the initial pool that is spent, held at tips whose supply is at
or below the threshold, held at tips above it whose placements were all
refused, and still active. With the corrected critical radius, one replicate
each: the kidney spends its pool (1.01); the pancreas spends 82 % and holds
18 % by supply; the liver spends 13 %, holds 37 % by supply and 48 % by lack of
room. The liver's closure is not supply alone: almost half its pool is held
where supply is sufficient and there is no space to place a unit.

### Robustness to the auxiliary rules

The vascular root at the deepest point, at the surface (entering from the
boundary along $-y$, as a hilum) or at random; each sprouting rule removed in
turn (deficit, load, crowding, mesenchymal signal); and units assigned by
vascular territory, the nearest segment and then the nearest tip in its
subtree, instead of the nearest tip in space. Six seeds per regime and
variant (`scripts/21_robustness.py`).

Version 4.5 measures the vascular space with one quantity wherever it acts
(lumens plus junctions, for the sprouting stop and for the calibre feedback;
version 4.4 used lumens alone in the second), takes the volume of a junction
from its geometry, a sphere of its radius ($4\pi r^3/3$, where earlier versions
used $100\,r^3$, chosen to stop branching), and estimates the room left at a
tip from trial placements independent of the dynamics. With the old junction
volume counted as space, 42 % of the budget was filled and two thirds of the
vessels sat at the critical radius, with Murray ratios of 0.5 to 0.7; with the
geometric volume the tree still stops on its own and the ratio returns to 0.9.

The closure modes hold in all eight variants: the kidney spends its pool; the
pancreas spends 55 to 88 % and holds the rest by supply; the liver spends 7 to
15 % and holds the rest mostly by supply (61 to 75 %) and partly by lack of
room (15 to 32 %). The pool elasticity is higher in the pancreas than in the
liver in all eight (z from 2.3 to 8.3 in six; 1.4 with a random root and 1.2
without load-driven sprouting). The kidney lies between them at every setting
of the reference rules (liver +0.13 to +0.26, kidney +0.35 to +0.40, pancreas
+0.61 to +0.70) but is not separated from the liver robustly. Between 19 and
25 % of the liver's units change territory between the two assignment rules
without changing these conclusions.

Two further variants test the progenitor process itself: nucleation drawn as
Binomial$(P_i, 1 - e^{-\nu (q_i - q_{\text{set}})_+ \Delta t})$, each progenitor
differentiating independently, instead of a Poisson draw capped at $P_i$; and
the same with the day split into two steps for nucleation, division and loss.
The closure modes and the pancreas-liver contrast hold in both (z = 6.6 and
3.9). Vascular adaptation and sprouting still run once a day.

### Murray's relation in grown trees

With flow driven by pressure, the branching relation is a result rather than
a restatement of the rule. With the exponent $b + 2 = 2.9$, the ratio
$r_p^{2.9}/\sum_c r_c^{2.9}$ over the bifurcations of grown trees has median
0.95 in the liver regime (313 bifurcations, interquartile 0.94 to 0.96) and
0.96 in the kidney regime (68). The deficit of 3 to 4 % has a sign the theory
gives: the setpoint rises with pressure, so vessels near the root, at higher
pressure, run at higher shear and are slightly narrower than the pure law asks.

### Branching exponent of a human liver: not yet measurable

The liver of donor LADAF-2020-27 was analysed at full resolution (157 µm) in
40 overlapping blocks, with the radius of each branch read just outside its
junction (1.2 to 3 junction radii from the centre), which is what densely
branching vessels allow. Over 112 bifurcations the median exponent is 2.32
(95 % interval 1.99 to 2.42), but that number is not a measurement of the
tree. It moves from 1.6 to 4.4 as the minimum radius required of the smaller
daughter rises from 1.5 to 3 voxels, and it splits by geometry: 3.37 where a
trunk continues and gives off a side branch, 1.92 where a vessel divides. An
exponent below 2 at a division means the daughters carry more cross-section
than their parent, which points to the parent being measured short.

The likely cause is that the vessels of an organ imaged ex vivo are empty and
flatten, the large veins most, and a distance transform reads the minor
half-axis of a flattened lumen rather than its radius. The bias grows with
size and differs between branch geometries, which is the pattern observed.
The earlier figure of 2.87 to 3.02 at 314 µm rests on the same radius and is
withdrawn with it. The exponent becomes measurable with a radius that does not
depend on flattening: the perimeter of a collapsed thin-walled tube is
preserved where its area and minor axis are not, so a perimeter-equivalent
radius on cross-sections normal to the centreline is the next step, with the
portal and hepatic venous trees analysed separately.

### Withdrawn

- Equation 15 as an upper bound on contraction. It is the weak-feedback limit.
- A resonance of vascular disorder in the calibre relaxation rate, an
  artefact of three errors in the open engine.
- An early-history lock-in as the source of spread in the closed model; it
  was the random root.
- Exact conservation of variance under exhaustion and limiting pool as a
  result of the open engine, where it is built in. In the closed model a
  budget without loss does pass its pool's spread whole (1.12, interval
  0.85 to 1.33), while a pool with loss passes 0.59.
- A diagonal transmission in which each closure mode filters the other's
  variable. The setpoint passes in every regime.
- A branching exponent of 2.87 to 3.02 in a human liver; the radius behind
  it reads the minor axis of flattened vessels.
- $N \propto V^{0.64}$ as a law of the closed model; the organs did not fill
  their domains.
- A difference between spread in maturity and in growth rate under capacity
  closure; it came from the progenitor overdraft.

### Open

- Whether the liver window closes by capacity, as §5 has it when it
  eliminates $\sigma(t)$, or by signal decay, as the closure table has it.
  The open engine could not satisfy both the measured contraction and closure
  within the run.
- The size sweep on the closed model (`make closed`), which needs several
  minutes per size and replicate and is for a multicore machine.
- The human partition of Figure 6 as a prediction: two kidneys of one
  individual share its parameters, differ in domain and seeds, and the
  fractions of variance between and within individuals follow from the
  transmission above once the input spreads are fixed independently.
- The density, reach and kill distance of the mesenchymal signal are set,
  not measured.
- The kidney's setpoint elasticity from the race between conversion and
  loss.
- What holds 85 % of the liver's leftover pool at tips whose supply exceeds
  the threshold.
- Independent measures of systemic and local input spreads, which turn the
  partition formula into a prediction.
- Why the model's floor of spread (0.108) lies above the mouse liver's
  final spread (0.074), and whether a signal-timed stop closes the gap.
