# The paper: framework, claims and evidence

Target: *Frontiers in Cell and Developmental Biology*, article type
Hypothesis and Theory.

**Status after the external review of 20 November 2026.** Two errors were
confirmed and fixed: the critical radius was $\sqrt{3}$ times too large and
served as the floor of every calibre, and flow was not conserved where tissue
hangs from an internal segment (no run had such tissue, so only the first
error changed results). Every closed-model number below predates the fix and
must be regenerated with `make paper` before writing. Spot checks with the
corrected model, six paired replicates each: the pool elasticity keeps its
order at the reference setting (liver +0.34, kidney +0.42, pancreas +0.77)
and stays low in the liver at all five settings (+0.18 to +0.36); claim 4 is
withdrawn; the liver closes partly by supply and partly by lack of room. Each claim gives the number behind it, the figure that
shows it and the command that reproduces it. `THEORY.md` holds the
derivations; nothing in the paper rests on anything not listed here.

## The hypothesis

An organ stops making functional units when local resources run out at the
tips of its growing tree: perfusion per unit, progenitors, or room to place a
unit. Given a geometry and the rules that organise the tree in space, which of
these binds first is what distinguishes a liver from a pancreas; in the runs
reported, the geometry is held fixed and only the pool and its turnover
change between regimes. Two consequences follow. Organs differ in which individual
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

**4. One wall-cost exponent sets the vascular rules, and a branch buys reach,
not conductance.** With wall cost exponent $b$, branching follows
$r_p^{b+2} = \sum r_c^{b+2}$, the shear setpoint rises as $p^{b/2}$ and
daughters are $2^{-1/(b+2)}$ of their parent. Two daughters conduct
$2\cdot 2^{-4/(b+2)}$ of their parent, 0.77 at $b = 0.9$, so a split is
justified by the tissue it reaches. Shear on a segment in series with its
tissue peaks at $r_c = (\ell/3R_t)^{1/4}$, below which a vessel closes. In
grown trees $r_p^{2.9}/\sum r_c^{2.9}$ has median 0.952 (313 bifurcations,
liver) and 0.964 (68, kidney): 3 to 4 % below one, the sign a
pressure-dependent setpoint predicts. Figure 5;
`scripts/17_paper_runs.py regimes`.

**5. Capacity closure contracts the spread of unit number.** The spread
relaxes as $\sigma^2(t) = \sigma_0^2 e^{-2\lambda t} + \sigma_\infty^2(1 - e^{-2\lambda t})$,
with weak-feedback limit $\mathrm{CV}_i/\mathrm{CV}_f = (N_f/N_i)^2$, which
gives 3.31 for the mouse liver against 3.09 measured (Azizoglu et al. 2025).
The closed model contracts the mouse's initial spread (CV 0.228) to 0.108, a
factor of 2.1: the right direction, but its floor lies above the mouse's
final 0.074. Figure 4; `THEORY.md` §7 and §13.

## Predictions for experiment

- **Two holds on the liver's progenitors.** At the end of hepatic growth,
  residual progenitors should split between sites short of perfusion and
  sites crowded by existing units, and relieving either one alone should
  release only part of them.
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

- That the pool waiting just below the threshold sets the response to the
  threshold. The formula is the direct, frozen-state term; after the critical
  radius was corrected, feedback through the tree dominates and the formula
  fails (correlation −0.15 over five liver settings).

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
cannot stand in for a human organ's. The model's sensitivity to the perfusion
threshold has no derivation that survives the full feedback. The human data give organ shape and size but no unit
counts, so the model and the cohort meet only through the partition formula,
which awaits independent input spreads.
