# The paper: framework, claims and evidence

Target: *Frontiers in Cell and Developmental Biology*, article type
Hypothesis and Theory. Each claim gives the number behind it, the figure that
shows it and the command that reproduces it. `THEORY.md` holds the
derivations; nothing in the paper rests on anything not listed here.

## The hypothesis

An organ stops making functional units when one of two local resources runs
out at the tips of its growing tree: perfusion per unit, or progenitors. Both
are local to each tip, and which one binds first is what distinguishes a liver
from a pancreas. Two consequences follow. Organs differ in which individual
variation reaches their size, and an organ's sensitivity to its perfusion
threshold is set by how many progenitors wait at tips whose supply sits just
below that threshold.

## Claims

**1. In human kidneys, size is set by the individual and shape by position.**
In 96 kidneys from 48 donors (TRUSTED, clinical CT), 89 % of the variance in
volume lies between donors and 80 % of the variance in sphericity within a
donor. The within-to-between distance ratio is 1.54 (permutation
$P$ = 0.0005). Left and right differ systematically: by paired Wilcoxon with
Bonferroni over six invariants, volume ($P$ = 0.0007), flatness and
surface-per-volume ($P$ < 0.0001) survive. Figure 6; `make paper` from
`source_data/`, or `scripts/07_trusted.py` and `scripts/08_partition.py` with
the cohort.

**2. Organ identity is which bound binds first.** One nucleation law with a
supply threshold and a local progenitor pool gives three regimes from two
numbers. The liver stops at its supply bound (supply per unit 0.87 of the
threshold, 81 % of its pool unused); the kidney stops when loss and
conversion have spent its pool (supply 1.32 of the threshold); a budget
without loss, the pancreas, stops short of its pool with supply to spare
(1.34). Figure 2; `scripts/17_paper_runs.py regimes`.

**3. Progenitor number reaches size in a budget organ and much less in a
supply-bound one.** The elasticity of the final unit count to the progenitor
pool, measured with paired runs at five settings per regime, is 0.32 to 0.39
in the liver, 0.51 to 0.62 in the kidney and 0.68 to 0.92 in the pancreas,
ordered the same way at every setting. Transmission of an individual spread
(24 replicates, 95 % bootstrap intervals) agrees: 0.23 [0.00, 0.42] in the
liver, 0.59 [0.41, 0.72] in the kidney and 1.12 [0.85, 1.33] in the budget.
This is the pattern found experimentally by Stanger, Tanaka and Melton
(*Nature* 2007): halving pancreatic progenitors halved the pancreas, which
never compensated, while the liver did. The model was not fitted to that
experiment. It filters the liver's pool less than the full compensation they
observed. Figure 3; `scripts/15_transmission.py` with `OG_NREP=24`,
`scripts/19_elasticity.py`.

**4. The sensitivity to the perfusion threshold is set by the pool waiting
just below it.** Lowering the threshold activates the progenitors held at tips
whose supply sits just below it, so $\varepsilon_{q_{\text{set}}} = -q_{\text{set}}\, w(q_{\text{set}})/N$, with
$w$ the pool held per unit of supply just below the threshold, read from the
final state with no fitted parameter. Across five liver settings it predicts
the measured elasticity with correlation 0.98 (mean absolute error 0.33),
from −0.41 to −3.03; in the budget it holds at three of five settings. In
the pancreas all of the leftover pool lies below the threshold, 42 % of it in
the slice just under it. In the liver 85 % of the leftover pool sits at tips
whose supply is above the threshold, where supply does not hold it; only 5 %
lies in that slice, and that 5 % sets the liver's elasticity. In the
kidney, whose pool loss spends, nothing is stranded and the threshold acts
instead through the race between converting a progenitor and losing it. A
global alternative, closure on total root flow, fails a test with its window
fixed in advance. Figure 5; `THEORY.md` §13; `scripts/19_elasticity.py`.

**5. One wall-cost exponent sets the vascular rules, and a branch buys reach,
not conductance.** With wall cost exponent $b$, branching follows
$r_p^{b+2} = \sum r_c^{b+2}$, the shear setpoint rises as $p^{b/2}$ and
daughters are $2^{-1/(b+2)}$ of their parent. Two daughters conduct
$2\cdot 2^{-4/(b+2)}$ of their parent, 0.77 at $b = 0.9$, so a split is
justified by the tissue it reaches. Shear on a segment in series with its
tissue peaks at $r_c = (3\ell/R_t)^{1/4}$, below which a vessel closes. In
grown trees $r_p^{2.9}/\sum r_c^{2.9}$ has median 0.952 (313 bifurcations,
liver) and 0.964 (68, kidney): 3 to 4 % below one, the sign a
pressure-dependent setpoint predicts. Figure 5;
`scripts/17_paper_runs.py regimes`.

**6. Capacity closure contracts the spread of unit number.** The spread
relaxes as $\sigma^2(t) = \sigma_0^2 e^{-2\lambda t} + \sigma_\infty^2(1 - e^{-2\lambda t})$,
with weak-feedback limit $\mathrm{CV}_i/\mathrm{CV}_f = (N_f/N_i)^2$, which
gives 3.31 for the mouse liver against 3.09 measured (Azizoglu et al. 2025).
The closed model contracts the mouse's initial spread (CV 0.228) to 0.108, a
factor of 2.1: the right direction, but its floor lies above the mouse's
final 0.074. Figure 4; `THEORY.md` §7 and §13.

## Predictions for experiment

- **Stranded progenitors in the liver.** Residual hepatic progenitors should
  sit at the vascular front, at tips whose perfusion is just below the
  threshold that stops nucleation, and the liver's response to a shift in
  that threshold (hypoxia or VEGF sensing) should scale with how many sit
  there.
- **A budget organ's perfusion sensitivity is not fixed.** Pancreas size
  should track progenitor number, as observed, while its sensitivity to the
  perfusion threshold rises steeply as the threshold approaches the supply at
  its tips.
- **Paired kidneys.** With the model's elasticities, the fraction of size
  variance between individuals is
  $(\varepsilon_q^2\sigma_{q,s}^2 + \varepsilon_P^2\sigma_{P,s}^2)$ over that
  plus the local and contingent terms. TRUSTED requires the systemic part to
  be 2.8 times the within-individual part in standard deviation. Paired
  nephron counts and split renal function in healthy donors would turn this
  requirement into a test.

## Not claimed

- The earlier figures from the open engine (`make all`), in which exhaustion
  and limiting pool conserved spread by construction.
- A difference between spread in maturity and spread in growth rate under
  capacity closure. It came from a defect, since fixed, that let a tip make
  more units than it had progenitors.
- A branching exponent for human liver vessels: measured at full
  resolution it moves from 1.6 to 4.4 with the selection threshold, because
  the radius reads the minor axis of flattened ex vivo vessels.
- Glomerular density and spacing in human kidneys, which need a validated
  detector.

## Limitations

The model's parameters are chosen, not fitted; the one calibration is the
mouse contraction, which the model undershoots. Results come from one domain
scale, in organs of about 50 to 550 units, so the model's contingent spread
cannot stand in for a human organ's. The kidney's threshold mechanism is
derived, not tested, and what holds most of the liver's leftover pool at tips
above the threshold is not known. The human data give organ shape and size but no unit
counts, so the model and the cohort meet only through the partition formula,
which awaits independent input spreads.
