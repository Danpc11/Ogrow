# The paper: framework, claims and evidence

Target: *Frontiers in Cell and Developmental Biology*, Hypothesis and Theory.
Reference model: version 4.7.1 (`src/organgrow/closed.py`). Every number from
the model that the paper uses is in `PAPER_NUMBERS.md`, written by
`make paper` with the version and commit that produced it; this file does not
copy model numbers by hand. Numbers given here are either data (the human
cohort) or analytic. Model numbers depend on the random stream, which numpy
does not keep fixed across versions, so they are produced with the versions
pinned in `requirements.txt` and recorded in the stamp.

## The hypothesis

An organ stops making functional units when local resources run out at the
tips of its growing tree: perfusion per unit, progenitors, or room to place a
unit. Given a geometry and the rules that organise the tree in space, which of
these holds an organ's remaining progenitors at closure distinguishes a liver
from a pancreas. In the runs reported the geometry is held fixed and only the
pool and its turnover change between regimes. Supply is hydraulic: the outlet
pressure at a tip relative to the driving pressure, not oxygen or metabolite
delivery.

## Claims

**1. In human kidneys, volume varies mostly between donors and shape within
a donor.** In 96 kidneys from 48 donors (TRUSTED, clinical CT), 89 % of the
variance in volume lies between donors and 80 % of the variance in sphericity
within a donor; the within-to-between distance ratio is 1.54 (permutation
$P$ = 0.0005). Left and right differ systematically: by paired Wilcoxon with
Bonferroni over six invariants, volume ($P$ = 0.0007), flatness and
surface-per-volume ($P$ < 0.0001) survive. Figure 2.

**2. Organs differ in what holds their progenitors when growth stops, and
the difference does not depend on the auxiliary rules.** One nucleation law,
with supply, progenitors and room local to each tip, gives organs that close
in different ways: the kidney spends its pool; the pancreas spends most of it
and keeps the rest at tips short of supply; the liver keeps nearly all of it,
mostly at tips short of supply and partly at tips with no room. Read from the distribution of the pool over tips at
closure, with the room left at each tip tested by trial placements that do
not depend on the dynamics, not from a mean. The three
closure modes persist in each of twelve variants: the vascular root at the
deepest point, at the surface or at random; each of the four sprouting rules
removed in turn; units assigned by vascular territory instead of by nearest
tip; binomial nucleation with one, two and four steps a day; and division and
loss as competing fates. Closure is stable in all three regimes: forty more
days add little to the liver and the pancreas and nothing to the kidney. It
differs in whether
it can be undone. Relieving perfusion (a 20 % lower threshold) or room (a 20 %
shorter least spacing) on the day of closure restarts the liver's growth;
nothing restarts the kidney's, whose pool is spent; the pancreas gains at most
a few units (`PAPER_NUMBERS.md`, Interventions). The shares of the pool held by each
constraint describe where the pool sits, not how much each constraint weighs:
the liver is held by both. `PAPER_NUMBERS.md`, Regimes and Robustness; Figures 1, 3 and 5A.

**3. Progenitor number reaches size more in a budget organ than in a liver.**
The elasticity of the final unit count to the progenitor pool is higher in the
pancreas than in the liver in every variant, with the z of each in
`PAPER_NUMBERS.md`, Robustness. This
is the pattern found experimentally by Stanger, Tanaka and Melton (*Nature*
2007): halving pancreatic progenitors halved the pancreas, which never
compensated, while the liver did. The model was not fitted to that
experiment. The model does not separate the kidney, which spends its pool
against loss, from the liver robustly, so it does not order all three organs.
`PAPER_NUMBERS.md`, Transmission, Elasticities and Robustness; Figures 4, 5B and 5C.

**4. One wall-cost exponent sets the vascular rules, and a branch buys reach,
not conductance.** Minimising $\mathcal{C} = Q^2/r^4 + k\,(p\,r^2)^b$ per unit
length, with wall cross-section $\propto p\,r^2$ by Laplace's law, gives
branching $r_p^{b+2} = \sum r_c^{b+2}$ and wall shear
$\propto p^{b/2} r^{b-1}$. Two daughters conduct $2\cdot 2^{-4/(b+2)}$ of
their parent, 0.77 at $b = 0.9$ and below one for every $b < 2$, which includes the
measured range. Shear on a
segment in series with its tissue peaks at $r_c = (\ell/3R_t)^{1/4}$, below
which a vessel closes. The Murray ratio in grown trees:
`PAPER_NUMBERS.md`, Murray; Figure 6.

**5. Capacity closure contracts the spread of unit number.** If every
individual gains the same increment of $N^2$, the ratio of initial to final
CV is $(N_f/N_i)^2$; this is the weak-feedback limit of the nucleation law,
derived, and gives 3.31 for the mouse liver against 3.09 measured (Azizoglu et
al. 2025). The relaxation form
$\sigma^2(t) = \sigma_0^2 e^{-2\lambda t} + \sigma_\infty^2(1 - e^{-2\lambda t})$
was fitted to simulations of the earlier open engine ($R^2$ 0.92 to 0.99); it
is not derived, has not been fitted to the closed model, and no mechanism is
claimed for $\lambda$. The model contracts more than the mouse liver, from a
larger initial spread of its own; its final spread is in
`PAPER_NUMBERS.md`, beside the mouse's:
`PAPER_NUMBERS.md`, Contraction; Figure 7.

## Predictions for experiment

- **Closure is reversible in the liver and not in the kidney.** At the end
  of growth the liver keeps progenitors that more perfusion or more room puts
  back to work, while the kidney keeps none. This matches, qualitatively, the
  regenerative growth of the adult liver and the absence of new nephrons after
  birth; it is a correspondence, not a test.
- **A budget organ tracks its progenitors.** Pancreas size should scale with
  progenitor number, as observed, while the liver's should not; a kidney whose
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
contingent spread cannot stand in for a human organ's. The progenitor processes converge with one, two
and four steps a day with the tree held fixed within the day; the coupled
vascular and tissue dynamics run once a day, and their convergence with a
smaller step was not tested. The two interventions are not matched in
strength, so their effects show that both constraints hold the liver, not
which holds it more. The critical radius is the floor of calibre
adaptation, so simulated trees cannot show calibres below it. The pool elasticity of the liver is not
measurable without load-driven sprouting (s.e. 0.55), and the pancreas-liver
contrast is not significant with a random root. The room left at a tip is
estimated from trial placements in the final state, which measures residual
capacity, not whether growth could resume later. The
human data give shape and size but no unit counts, so model and cohort meet
only through the partition formula, which awaits independent input spreads.

## Before submission

1. Commit, then `make test` and `make paper` on that commit, so that
   `PAPER_NUMBERS.md` names version 4.7.1 and a clean commit.
2. Write the manuscript from this file, `PAPER_NUMBERS.md` and `FIGURES.md`.
