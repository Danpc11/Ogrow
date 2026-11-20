# The paper: framework, claims and evidence

Target: *Frontiers in Cell and Developmental Biology*, Hypothesis and Theory.
Reference model: version 4.4 (`src/organgrow/closed.py`). Every number from
the model that the paper uses is in `PAPER_NUMBERS.md`, written by
`make paper` with the version and commit that produced it; this file does not
copy model numbers by hand. Numbers given here are either data (the human
cohort) or analytic.

## The hypothesis

An organ stops making functional units when local resources run out at the
tips of its growing tree: perfusion per unit, progenitors, or room to place a
unit. Given a geometry and the rules that organise the tree in space, which of
these holds an organ's remaining progenitors at closure distinguishes a liver
from a pancreas. In the runs reported the geometry is held fixed and only the
pool and its turnover change between regimes.

## Claims

**1. In human kidneys, size is set by the individual and shape by position.**
In 96 kidneys from 48 donors (TRUSTED, clinical CT), 89 % of the variance in
volume lies between donors and 80 % of the variance in sphericity within a
donor; the within-to-between distance ratio is 1.54 (permutation
$P$ = 0.0005). Left and right differ systematically: by paired Wilcoxon with
Bonferroni over six invariants, volume ($P$ = 0.0007), flatness and
surface-per-volume ($P$ < 0.0001) survive. Figure 2.

**2. Organs differ in what holds their progenitors when growth stops, and
the difference does not depend on the auxiliary rules.** One nucleation law,
with supply, progenitors and room local to each tip, gives organs that close
in different ways: the kidney spends its pool; the pancreas spends most of it
and keeps the rest at tips short of supply; the liver keeps nearly all of it,
held both at tips short of supply and at tips with no room. Read from the
distribution of the pool over tips at closure, not from a mean. The three
closure modes persist in each of eight variants: the vascular root at the
deepest point, at the surface or at random; each of the four sprouting rules
removed in turn; and units assigned by vascular territory instead of by
nearest tip. `PAPER_NUMBERS.md`, Regimes and Robustness; Figures 1, 3 and 5A.

**3. Progenitor number reaches size more in a budget organ than in a liver.**
The elasticity of the final unit count to the progenitor pool is higher in the
pancreas than in the liver in all eight variants, significantly in six. This
is the pattern found experimentally by Stanger, Tanaka and Melton (*Nature*
2007): halving pancreatic progenitors halved the pancreas, which never
compensated, while the liver did. The model was not fitted to that
experiment. The kidney, which spends its pool against loss, transmits about
as little as the liver, so the model does not order all three organs.
`PAPER_NUMBERS.md`, Transmission, Elasticities and Robustness; Figures 4, 5B and 5C.

**4. One wall-cost exponent sets the vascular rules, and a branch buys reach,
not conductance.** Minimising $\mathcal{C} = Q^2/r^4 + k\,(p\,r^2)^b$ per unit
length, with wall cross-section $\propto p\,r^2$ by Laplace's law, gives
branching $r_p^{b+2} = \sum r_c^{b+2}$ and wall shear
$\propto p^{b/2} r^{b-1}$. Two daughters conduct $2\cdot 2^{-4/(b+2)}$ of
their parent, 0.77 at $b = 0.9$ and below one for every $b$. Shear on a
segment in series with its tissue peaks at $r_c = (\ell/3R_t)^{1/4}$, below
which a vessel closes. The Murray ratio in grown trees:
`PAPER_NUMBERS.md`, Murray; Figure 6.

**5. Capacity closure contracts the spread of unit number.** The spread
relaxes as $\sigma^2(t) = \sigma_0^2 e^{-2\lambda t} + \sigma_\infty^2(1 - e^{-2\lambda t})$,
with weak-feedback limit $\mathrm{CV}_i/\mathrm{CV}_f = (N_f/N_i)^2$, which
gives 3.31 for the mouse liver against 3.09 measured (Azizoglu et al. 2025).
The model closes at a spread close to the mouse liver's final one, from a
larger initial spread of its own: `PAPER_NUMBERS.md`, Contraction; Figure 7.

## Predictions for experiment

- **Two holds on the liver's progenitors.** At the end of hepatic growth,
  residual progenitors should split between sites short of perfusion and
  sites crowded by existing units, and relieving either one alone should
  release only part of them.
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
- A difference between spread in maturity and spread in growth rate, which
  came from a defect that let a tip make more units than it had progenitors.
- A branching exponent for human liver vessels, whose radius reads the minor
  axis of flattened ex vivo vessels.
- Glomerular density and spacing in human kidneys, which need a validated
  detector.

## Limitations

The model's parameters are chosen, not fitted; the one calibration is the
mouse contraction. Results come from one domain scale, so the model's
contingent spread cannot stand in for a human organ's. The time step is one
day and has not been varied. The pool elasticity of the liver is not
measurable when the root is placed at random or load-driven sprouting is
removed (s.e. 0.5 to 0.7). The
human data give shape and size but no unit counts, so model and cohort meet
only through the partition formula, which awaits independent input spreads.

## Before submission

1. Commit, then `make test` and `make paper` on that commit, so that
   `PAPER_NUMBERS.md` names version 4.4 and a clean commit.
2. Write the manuscript from this file, `PAPER_NUMBERS.md` and `FIGURES.md`.
