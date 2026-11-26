# Cover letter, *Cell Systems*

Draft. The journal's page asks that the letter state what was known, the
conceptual advance, and the significance, and allows it to name reviewers and
related work. Numbers are those of `PAPER_NUMBERS.md` at version 4.67.0; check
them against the manuscript before sending.

---

Dear Editor,

We submit "Local resource limitation determines how organ size depends on
progenitor number" for consideration as an Article in *Cell Systems*.

**What was known.** Whether the size of an organ depends on the number of its
progenitors differs between organs. Ablating progenitors in the mouse embryo
leaves a pancreas that is smaller in proportion for life and a liver that
returns to normal size within days (Stanger et al., 2007). In the kidney, two
ablations of the same size gave a nephron deficit in proportion to the cells
removed (Cebrian et al., 2014) and one less than half as large (Muthukrishnan
et al., 2018). These results are read as properties of organs, intrinsic or
extrinsic control of size, and no quantitative account derives them from a
common mechanism or explains an organ that behaves both ways.

**The conceptual advance.** We show that one local growth rule produces all of
these outcomes. Functional units form at the tips of a perfused branching tree
at a rate set by the progenitors present and by perfusion above a threshold,
and growth stops at a tip when perfusion, progenitors, or space runs out.
Nothing in the model represents the size of the organ. The resource that is
exhausted first sets the elasticity of organ size to progenitor number. With
no fitted parameter, the model gives 0.44 of the reference mass for the
pancreas-like source after the ablation that left 0.36 in the mouse, and 0.86
for the liver-like source where the liver recovered to about 0.96. The kidney
is the informative failure: the rule reproduces the measured progenitor pool,
cap size, and number of tips only when progenitors also drive branching, and a
proportional nephron deficit only when nucleation is not limited by perfusion.
We then commit the model to the outcome of ablations that have not been done,
each with the result that would refute it.

**Significance.** The work replaces a classification of organs with a
measurable quantity and a mechanism, gives a reading of why nephron endowment
is vulnerable to early loss of progenitors, and states what a tissue grown
from progenitors must be supplied with: cells, perfusion, or room. We expect
it to interest readers working on organ size control, kidney and liver
development, regeneration, and organoid engineering.

**What the paper is not.** It contains no new experiment. Every comparison is
with published measurements, the kidney extension was written after those
measurements were known, and the manuscript says so. All code, and a single
command that regenerates every number and figure, are public.

**Related work.** A second manuscript from the same model, on how unit density
and reproducibility change with the size of an engineered construct, is in
preparation for another journal. It shares the model and its Methods and no
result or figure with this one. We have not discussed this submission with an
editor of the journal.

**Reviewers.** [Suggested reviewers and exclusions, to be completed by the
authors.]

Sincerely,

Daniel Pérez-Calixto, on behalf of all authors
Instituto Nacional de Medicina Genómica, Mexico City
dperez@inmegen.gob.mx
