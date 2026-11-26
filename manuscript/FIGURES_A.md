# Figures of the second manuscript

The figures of `MANUSCRIPT_A.md`, in its order. Citations are numbered as in its
reference list.

| Figure | File in `figures/out/` | Script |
| --- | --- | --- |
| Figure 1 | `Fig1_model.pdf` | `figures/fig_paper.py` |
| Figure 2 | `Fig2_bounds.pdf` | `figures/fig_paper.py` |
| Figure 3 | `Fig10_validation.pdf` | `figures/fig_validation.py` |
| Figure 4 | `Fig11_coupling.pdf` | `figures/fig_coupling.py` |
| Figure 5 | `Fig8_injury.pdf` | `figures/fig_injury.py` |
| Figure 6 | `Fig12_predictions.pdf` | `figures/fig_predictions.py` |

**Figure 1. The model, and a construct grown under it.** Panels A to C and F
are schematic and carry no data; proportions in them are not to scale. Colour
in those panels marks the three resources rather than the three regimes:
supply purple, competent cells teal, units and the room between them ink.
(A) A vascular tip, the tissue it serves, the supply per unit it delivers, the
competent cells it holds and the least spacing between units already built,
with the nucleation law. (B) The three quantities that can run out at a tip,
each with the state it ends in. (C) The one property varied between regimes:
whether the cell that forms a unit is spent in forming it. (D) One renewing
construct at arrest, projected on a plane: vessels with calibre as line width,
units, and the tips that still hold competent cells below the supply
threshold. (E) Disposition of the initial pool at arrest for each regime, mean
over 24 seeds; the annotation splits what the spent-with-loss regime spent
into the share that built units and the share lost to attrition. (F) The bound
that holds each regime at arrest, and whether growth can restart.

**Figure 2. What the binding resource is worth.** 24 seeds throughout, error
bars s.e. (A) Elasticity of the final unit count to the seeded pool, and
(B) to the perfusion setpoint, at five parameter settings per regime; points
are the mean over paired replicates. (C) Paired difference in pool elasticity
between the spent-no-loss and the renewing regime under thirteen variants of
the auxiliary rules, with z. (D) Units added by day 180 relative to day 140
after relieving one bound on day 140, (E) share of the day-140 pool released
by the same relief, left blank for a regime that holds nothing, and
(F) progenitors held on day 140 with the range the seeds span, log scale; the
dashed line is the pool below which a released share is not reported.

**Figure 3. Three numbers of the model beside measurements.** Filled symbols
are the model, open symbols are measurements, and nothing is fitted. (A) Share
of the organ that results against the share of the progenitor pool left.
Model: mass of the pancreas-like and liver-like constructs at day 140 as a
ratio to the reference, 24 paired seeds, bars the 95 % bootstrap interval.
Measured: weight of the mouse pancreas at birth after 8 % and 61 % of its
progenitors were removed, with 95 % intervals, and of the liver four days
after at least 65 % of its progenitors were removed [2]. The dashed line is an
organ in proportion to its pool. (B) The same for the kidney-like construct,
in units, against glomerular number in two mouse studies that each removed
about 40 % of nephron progenitors: through a CITED1 driver, counted three and
ten days after birth [15], and through a Gdnf driver, counted at birth and in
the adult [14]. The two are offset sideways to be told apart. (C) How far one
organ sits from the other grown like it: the kidney-like construct across 24
seeds, as the deviation of log mass from the mean, and the two kidneys of each
of 48 people, as the deviation of log volume from the pair with the mean
difference between sides removed [16]. Bars span one standard deviation
either side of zero.

**Figure 4. An extension in which progenitors make their tip branch.** A
dividing source in which a segment branches only while it holds two caps of
progenitors; on one day each progenitor dies with probability 0.4. All values
are the ablated construct over its control, means of 24 paired seeds. Model
days are not embryonic days, so measurements are drawn at the right-hand edge
and as bands, not along the time axis. (A) With no proliferative response and
the least perfusion-limited setting of the nucleation law: total progenitors,
progenitors per site, and branching sites after the ablation. Open symbols:
cells per tip and tips per kidney one, two and three days after ablation in
the mouse [14]. (B) Units in the same setting, without a response and with the
division rate doubled in a depleted cap. Bands: glomerular number in the study
that found no proliferative response [14] and in the one that found it [15].
(C) Units at the end of the window in the four settings of the nucleation law,
from the reference model's on the left to the least perfusion-limited on the
right; bars are 95 % bootstrap intervals, bands as in B.

**Figure 5. Injury after arrest, and the two things recovery needs.** An
organ loses 30 % of its parenchyma on day 140, after closure, and is read on
day 220; mean over the seeds of `scripts/24_injury.py`. (A) Liver mass against day, relative to the day
before the injury, as grown (the growth term has decayed), with growth
reactivated at a fixed inlet pressure, and with growth reactivated and the
root flow held at its pre-injury value; the dashed line marks the injury.
(B) Liver units, mass, mass per unit and perfused mass on day 220, relative
to the day before the injury, under the three conditions. (C) Kidney after
30 % and 50 % loss, with growth reactivated and the inflow held: nephrons
kept, mass per nephron, and the mass per nephron that full compensation would
require (hatched). (D) Day-220 mass against the flow the
circulation delivers, both relative to the day before the injury, with the
demand per unit mass measured once on the intact ensemble and held fixed; one
point per seed, mean with s.d., the fitted exponent, and the dashed line the
mass the supply budget would give if every unit reached its ceiling.

**Figure 6. Predictions for ablations that have not been done.** Circles
joined by a line are the model at doses already compared with a measurement,
open squares are those measurements, and stars are the model where nothing has
been reported. Bars are 95 % bootstrap intervals over 24 paired seeds. (A)
Share of the organ that results against the share of the progenitor pool left,
for the pancreas-like and liver-like constructs; the dashed line is an organ
in proportion to its pool. Measurements from [2]. (B) Share of the units of
the coupled kidney-like source against the share of its progenitors left, with
no proliferative response (filled, solid) and with the division rate doubled
in a depleted cap (open, dashed). Squares: glomerular number after the two
published ablations of about 40 % [14, 15]. (C) The same 40 % ablation on
three days of a window of 120; day 45 is the one already compared.
