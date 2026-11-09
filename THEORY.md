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

**Equation 7 — supply per unit.**

$$c_u(t) = \frac{\mathbb{C}(t)}{N(t)}$$

**Equation 8 — what reaches one unit.** Supply falls with distance $d_j$ to
the nearest vessel and rises with that vessel's calibre,

$$c_{u,j} = c_u \, e^{-d_j/\ell_a} \sqrt{r_{v(j)}/r_{\text{root}}}$$

**Equation 9 — homogenised shear.** For a population of calibres,

$$r_{\text{eff}} = \frac{\langle r^{5}\rangle}{\langle r^{4}\rangle}$$

exact, with no closure assumption, but dependent on the whole distribution and
not its mean. For a log-normal of width $\sigma$, $r_{\text{eff}} =
r_m e^{9\sigma^{2}/2}$: at $\sigma = 1$ the factor is 90, so the width
dominates any exponent.

### How supply scales with unit count

Write $\mathbb{C} \sim N^{\gamma}$. Since $c_u = \mathbb{C}/N \sim
N^{\gamma-1}$, the sign of $\gamma - 1$ decides whether adding units raises or
lowers the supply each one gets, and so whether a capacity-limited window can
close. Measured along grown trajectories:
$\\gamma$ = 0.82 (liver).
$\\gamma$ = 0.42 (kidney).
$\\gamma$ = 0.60 (pancreas).
$\\gamma$ = 0.51 (lung).
All are below 1, so the feedback is stabilising and no runaway occurs. They are
also well below the 4/3 that uniform-shear adaptation alone would give, because
the root is in series and caps the total.

### The tree adds little disorder of its own

With calibres frozen, replicates that reach the same unit count differ in
conductance by CV = 0.0003. Letting the calibres relax raises it to 0.154 at the
fastest rate tested, and the dependence on the relaxation rate is
monotone, not resonant: it rises smoothly from 0.0003 to 0.154 as the rate
goes from 0 to 1. An earlier version of this engine reported a resonant peak near
CV = 0.25; that was an artefact of three errors since fixed (integer
truncation in nucleation, a pressure field clamped at zero, and a circular
per-unit supply term), and it does not survive their correction.

The vasculature is therefore a conduit here, not a source of quenched
disorder, at least at the scale this engine resolves.

## 5. The fabrication window and how it closes

**Organ identity is one parameter: the closure mode.** Everything else in the
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

**Equation 15 — the contraction bound.** Under capacity saturation with
$\mathbb{C}$ fixed, the increment of $N^{2}$ is the same for every individual,
so the dispersion of $N^{2}$ is preserved and

$$\frac{\mathrm{CV}_i}{\mathrm{CV}_f} \le \left(\frac{N_f}{N_i}\right)^{2}$$

It is an **inequality**, not an identity, and the bound is reached only if the
setpoint were identical across individuals. With per-individual setpoints,

**Equation 16 — the floor on terminal variance.**

$$\mathrm{CV}_\infty \gtrsim \mathrm{CV}(c_{\text{set}})$$

Each individual converges to its own ceiling, so terminal variance cannot fall
below the parameter's own spread. This says **when the bound fails**, which the
earlier equality form did not.

The discriminating measurement is the **absolute** dispersion $\mathrm{sd}(N)$,
not the CV. A falling CV is ambiguous: it falls both when variance is destroyed
and when it is merely diluted by a growing $N$. Capacity saturation destroys
$\mathrm{sd}(N)$; exhaustion and limiting pool conserve it exactly.

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

has expectation $\sqrt{2}$ for **any** exchangeable sample, whatever the
cohort. Its value is therefore a property of the normalisation, not a
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
non-Newtonian rheology are absent, and pressure is implicit rather than
carried as a variable.

Further limits, including which parts rest on which evidence, are in
`KNOWN_LIMITATIONS.md`.
