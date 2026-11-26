# Local resource limitation determines how organ size depends on progenitor number

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

Target journal: *Cell Systems*, Article. Order of sections as Cell Press sets it: Summary, Introduction, Results, Discussion with Limitations of the study, Resource availability, Figure legends, STAR Methods. Highlights, the eTOC paragraph and the graphical abstract are given first because the journal asks for them as separate files.

Reference model version 4.67.0 (`src/organgrow/closed.py`). Every model number below is either held to the result files by `scripts/29_check_manuscript.py` or written between `auto:` markers by the script that computes it, so that no model number in this file is written by hand.

---

## HIGHLIGHTS

- A local growth rule on a perfused tree stops without sensing organ size
- The first resource to run out sets how size scales with progenitor number
- Pancreas and liver sizes after progenitor ablation are reproduced without fitting
- Nephron deficits require progenitor-driven branching; new ablations are predicted

## IN BRIEF

How the size of an organ depends on its progenitors differs between pancreas, liver, and kidney. Using a model in which functional units form at the tips of a perfused branching tree under one local rule, we show that the resource exhausted first determines this dependence, reproduce published progenitor-ablation experiments, and predict the outcome of ablations not yet performed.

## GRAPHICAL ABSTRACT

`figures/out/FigA0_graphical_abstract.pdf`

## SUMMARY

How the size of an organ depends on the progenitors it is built from differs between organs: ablating pancreatic progenitors yields a proportionally smaller pancreas, whereas the liver recovers its size. We show that a single local growth rule accounts for this difference. On a perfused branching tree, functional units form at tips at a rate set by local progenitors and by perfusion above a threshold, until perfusion, progenitors, or space runs out. Which resource runs out first determines the response to progenitor loss. Without fitting, the model reproduces the sizes of mouse pancreas and liver after progenitor ablation. Nephron number is reproduced only when progenitors also drive branching, and then only if nucleation is not perfusion limited. The model predicts the outcome of ablations not yet performed and identifies the measurement that reveals which resource limits a growing tissue.

## KEYWORDS

organ size control; progenitor cells; branching morphogenesis; perfusion; regeneration; nephron endowment; pancreas; liver; kidney; mathematical model

## INTRODUCTION

An organ reaches a characteristic size although none of its cells can observe the whole. Classical work divides the mechanisms into two kinds [1, 2]. In extrinsic control, systemic signals and functional demand set the size, and the tissue grows or regenerates until the demand is met; the liver, which restores its mass after resection, is the standard example [3]. In intrinsic control, the size is fixed from within and does not compensate for a loss.

The clearest experimental statement of the distinction came from ablating progenitors. Stanger and colleagues killed a controlled share of the progenitors of the mouse pancreas early in organogenesis and obtained a pancreas that was smaller in proportion and stayed so for life, whereas a liver deprived of most of its progenitors returned to normal size within days [4]. They concluded that pancreas size is limited by the number of embryonic progenitors and liver size is not. The kidney then gave both answers. Removing about 40% of nephron progenitors left a nephron deficit in proportion to the cells removed in one study [5] and a deficit less than half as large in another, in which the surviving progenitors proliferated faster [6]. Because nephron number is fixed at the end of nephrogenesis [7] and varies about tenfold among people [8], the difference matters beyond the mouse.

These results are usually read as properties of organs: some count their progenitors and some do not. That reading names the difference without deriving it, and it has no place for an organ that behaves both ways. It also leaves open a practical question. Kidney and liver are the two most transplanted organs, 110,467 and 42,497 of the 173,727 solid-organ transplants reported worldwide in 2024, an activity estimated to meet at most a tenth of the need [9], and efforts to grow replacement tissue from progenitors [10] must decide what such a tissue has to be supplied with: more cells, more perfusion, or more room.

Quantitative models have addressed parts of this problem. For the liver, the mass regained after resection has been modeled as the response to a metabolic load shared among the remaining cells [11] or to the portal flow that the remnant receives [12]; both take as given a signal that scales with the size of the organ. For branched organs, the form of the tree and where it stops follow from rules applied locally at growing tips [13], and in the kidney a simple change in the branching rule confines the tips to the surface of the organ [14]. What these models leave open is how the number of progenitors enters the final size, and why it does so differently in different organs.

Here we take a different starting point. Instead of asking whether an organ measures its size, we ask which locally available resource is exhausted first where the tissue is being built. We study a model in which functional units form at the tips of a perfused branching tree under one local rule, and growth ends at a tip when perfusion, progenitors, or space runs out there. The model contains no variable that represents the size of the organ. We show that the resource that runs out determines how the final size responds to the number of progenitors; that this accounts quantitatively, with no fitted parameter, for the pancreas and liver experiments; that the kidney experiments expose a coupling the minimal rule lacks, between progenitors and branching, and indicate that nephron formation in the embryonic kidney is not limited by perfusion; and that the same rule states when arrested growth can be restarted and what sets the size regained after injury. Finally, we commit the model to the outcome of ablations that have not been performed.

## RESULTS

### A local nucleation law stops growth without a size signal

We represent a growing tissue as functional units, such as lobules, nephrons, or acini, that form at the tips of a perfused branching tree inside a fixed domain (Figure 1A; STAR Methods). Each tip holds a pool of progenitors, $P$, and receives a perfusion per unit, $q$, computed from pressure-driven flow through the whole tree. Units form at a tip at the rate $\nu P (q - q^{*})_{+}$, where $q^{*}$ is the perfusion a unit requires. Every new unit draws on the supply of its tip, so nucleation lowers $q$, and the tree sprouts toward tissue that is short of supply. Growth ends at a tip when one of three local resources is exhausted: perfusion per unit has fallen to the threshold, the progenitors are gone, or no free site remains within reach.

Four parameters describe how the progenitor pool is maintained: whether the cell that forms a unit is spent in forming it, the size of the pool seeded, and the rates at which progenitors divide and are lost. We examine three settings, which we name for the tissues whose cell source they resemble (Figure 1; Table 1). In the liver-like source the cell is not spent, as hepatocytes are their own progenitors [15] and are renewed throughout life [16], and division balances loss. In the kidney-like source the cell is spent and the pool decays. In the pancreas-like source a small pool is spent and neither divides nor decays. The domain and every other parameter are the same for all three.

All three constructs stop growing (Figure 1B and 1C), at 208 ± 4, 89 ± 2 and 45 ± 1 units (mean ± SEM, 24 seeds), and each stops against a different bound (Figure 1D and 1E). The liver-like construct keeps its pool: at arrest 0.78 ± 0.03 of it sits at tips where perfusion has fallen to the threshold, and 0.19 ± 0.02 at tips with no free site. The kidney-like construct has no pool left: 0.30 ± 0.01 was consumed into units and 0.70 ± 0.01 lost to attrition. The pancreas-like construct spent 0.75 ± 0.02 of its budget on units and holds the remaining 0.25 ± 0.02 at tips that perfusion no longer reaches. Size is thus the outcome of local exhaustion, and three sources arrive at it by three routes. The ordering of the three holds under 13 variations of the auxiliary rules of the model (STAR Methods).

**Table 1. The three cell sources.** All other parameters are shared.

| Source | Cell spent in forming a unit | Progenitors seeded, $P_0$ | Division, $\rho$ (per day) | Loss, $\delta$ (per day) |
| --- | --- | --- | --- | --- |
| Liver-like | no | 2,600 | 0.10 | 0.10 |
| Kidney-like | yes | 300 | 0 | 0.05 |
| Pancreas-like | yes | 60 | 0 | 0 |

### The turnover of the cell source decides which resource runs out

The three sources differ in several parameters at once, so we asked which of them sets the bound. We varied renewal (two levels), seeded pool (three), and turnover (static, decaying, or balanced) in a full factorial and recorded which bound held the largest share of the pool at arrest (Figure 2A). Changing whether the cell is spent changed the bound in 2 of the 9 settings of the other two. Changing the turnover changed it in 5 of the 6 settings of renewal and the pool, and the pool alone in 2 of 6 (Figure 2B). A pool that decays is exhausted whatever else is true of it. A pool that does not decay outlasts the perfusion that could use it, unless it is both small and spent.

The intuitive distinction, between tissues that keep a self-renewing compartment and tissues that do not, is therefore the weakest of the three axes. What decides whether progenitors or perfusion ends growth is chiefly whether the pool is lost with time.

### The binding resource sets how size responds to progenitor number

If a construct stops because perfusion ran out, adding progenitors should change little; if it stops because progenitors ran out, its size should follow their number. We measured the elasticity of the final unit count to the seeded pool, $\varepsilon_P = \partial \ln N / \partial \ln P_0$, by paired perturbation (Figure 2C). It is +0.07 ± 0.06 in the liver-like construct, +0.37 ± 0.03 in the kidney-like, and +0.71 ± 0.05 in the pancreas-like. The elasticity to the perfusion threshold has the opposite sign and is larger in each case: −1.83 ± 0.22, −1.65 ± 0.06, and −0.77 ± 0.07. A tissue held by perfusion is nearly indifferent to its progenitors and sensitive to what each unit demands. A tissue held by its progenitors is the reverse, though never completely, because progenitors stranded at poorly perfused tips are not converted into units.

When pool and turnover are matched between a spent and a renewing source, the difference in $\varepsilon_P$ that remains is +0.17 at the liver-like setting, +0.04 at the kidney-like, and +0.28 at the pancreas-like. Renewal contributes, but less than half of the range across the three sources.

The elasticity $\varepsilon_P$ is what a progenitor-ablation experiment measures. We therefore compared the model with such experiments.

### The model reproduces the sizes of pancreas and liver after progenitor ablation

Stanger and colleagues ablated Pdx1-expressing pancreatic progenitors in the mouse embryo by tetracycline-regulated expression of diphtheria toxin and weighed the organ at birth [4]. Ablation through E9.5 removed 8% of the progenitor pool and left a pancreas at 77% of control weight (95% CI, 53%–102%); ablation through E10.5 removed 61% and left 36% (19%–53%). The small pancreas grew in parallel with controls for 11 weeks and did not catch up. With the same strategy in the liver, removal of at least 65% of hepatic progenitors left a liver at about one-third of control weight when the ablation ended and at normal weight 4 days later.

We seeded each construct with the share of its pool that these ablations left and changed nothing else (Figure 3A). The pancreas-like construct ended with <!-- auto:v.panc39 -->0.44 [0.37, 0.52]<!-- /auto:v.panc39 --> of its reference mass on 39% of the pool and with <!-- auto:v.panc92 -->0.81 [0.72, 0.90]<!-- /auto:v.panc92 --> on 92% (ratio of means over 24 paired seeds, 95% bootstrap interval), against 0.36 and 0.77 measured (Figure 3B). The liver-like construct seeded with 35% of its pool ended with <!-- auto:v.liv35 -->0.86 [0.84, 0.89]<!-- /auto:v.liv35 -->, against approximately 0.96. The model thus reproduces a pancreas whose size follows its progenitors and a liver whose size largely does not, with parameters that were set before the comparison was made.

Two qualifications apply. The liver-like construct recovers less completely than the liver, and its shortfall grows with deeper cuts: it ends with <!-- auto:v.liv50 -->0.92<!-- /auto:v.liv50 --> of its mass on half the pool and <!-- auto:v.liv25 -->0.77<!-- /auto:v.liv25 --> on a quarter. And in the embryo the liver was already built when its progenitors were removed, so that experiment is closer to an injury than to a smaller seeding; we compare only the outcome.

### Nephron number after ablation exposes a missing coupling

The kidney is the harder test, because two experiments that removed the same share of progenitors disagree. Muthukrishnan and colleagues ablated about 40% of CITED1-expressing nephron progenitors at E15.5 and found 12% fewer glomeruli at P3 and 21% fewer at P10; the surviving progenitors more than doubled their proliferation and repopulated the niche [6]. Cebrian and colleagues ablated 40% of the cap mesenchyme at E12.5 through a Gdnf driver and found 52% fewer glomeruli at E19.5 and 42% fewer in the adult; the survivors did not proliferate faster, ureteric branching slowed, and nephrogenesis ceased on schedule [5].

The kidney-like construct with 40% of its pool removed retained <!-- auto:v.kid60 -->0.80 [0.77, 0.83]<!-- /auto:v.kid60 --> of its units (Figure 3C). This agrees with the first experiment and is less than half the deficit of the second. More importantly, the model contains the mechanism of neither. Its kidney-like progenitors do not divide, so they cannot repopulate a niche, and its tree does not depend on them, so branching cannot slow. We return to this below.

One comparison needs no perturbation. Two constructs grown from identical parameters differ only by chance, as do the two kidneys of one person. Across seeds, the standard deviation of the log mass of the kidney-like construct is <!-- auto:v.sd.model -->0.104<!-- /auto:v.sd.model -->. In CT volumes of <!-- auto:v.pairs.n -->48<!-- /auto:v.pairs.n --> pairs of human kidneys [17], the within-pair standard deviation of log volume, after removal of the mean difference between sides, is <!-- auto:v.sd.pairs -->0.098<!-- /auto:v.sd.pairs --> (Figure 3D). Nothing in the model was set by the cohort. The agreement is nonetheless limited: part of the measured variation is segmentation error; a mass in model units is compared with a volume; and the variation in the model grows with the size of its domain, whereas within-pair variation in the cohort does not rise detectably with kidney volume (Spearman <!-- auto:v.pairs.size -->$\rho = +0.25$, $P = 0.09$<!-- /auto:v.pairs.size -->).

### Progenitor-driven branching accounts for both kidney experiments where perfusion is not limiting

The disagreement between the two kidney experiments points to what the minimal rule lacks. In the embryonic kidney, nephron progenitors renew themselves until a timed cessation [18, 7], and the cap of progenitors around each ureteric tip supplies the signals that make the tip branch [19, 5]. We added these features as an extension that leaves the reference model unchanged (Figure 4A; STAR Methods). Progenitors divide and are not lost to attrition. A tip may branch only while it holds two caps' worth of progenitors, and it passes half to the new branch. Growth ends at a fixed time. On one day, every progenitor is killed with probability 0.4. A single further switch lets a depleted cap divide faster. Because these rules were written with the experiments in hand, the final deficit cannot count as a prediction. The test is whether one coupling yields several observations that were made separately.

It does for the branching side of the experiment of Cebrian and colleagues, in every setting of the nucleation law that we examined (Figure 4B). In the embryo, the total cap mesenchyme of ablated kidneys grew as fast as in controls, so the deficit in cells persisted; the number of progenitors per tip, 0.60 of control one day after tamoxifen, returned to 0.95 and 1.04 on the following two days; and the number of tips fell to 0.71 and then 0.59 of control [5]. In the model, the pool stayed at <!-- auto:k.pool -->0.53 to 0.57<!-- /auto:k.pool --> of control to the end of the window, the cap returned to 0.95 of control <!-- auto:k.capback -->within 3 to 50 days in all 4 settings<!-- /auto:k.capback -->, and the number of branching sites settled at <!-- auto:k.sites -->0.52 to 0.60<!-- /auto:k.sites -->. The explanation is the one those authors proposed. A tip that has lost part of its cap does not branch until the cap is restored, so the size of the cap is corrected at the expense of the tree.

Whether fewer tips produce proportionally fewer nephrons depended on perfusion (Figure 4C and 4D). With the nucleation law set as in the reference model, the units of the ablated construct recovered to <!-- auto:k.units.ref -->0.81 [0.72, 0.94]<!-- /auto:k.units.ref --> of control, because each remaining site is better perfused and forms more units. This compensation is the central mechanism of the model, and the embryonic kidney does not show it. When we lowered the perfusion threshold fourfold and slowed the law, so that perfusion seldom limits nucleation, units ended at <!-- auto:k.units.low -->0.64 [0.57, 0.70]<!-- /auto:k.units.low --> of control, and those formed after the ablation at <!-- auto:k.since.low -->0.52<!-- /auto:k.since.low -->. The measured 0.48–0.58 lies at the lower end of that interval. Read through the model, a nephron deficit in proportion to the progenitors removed implies that nephron formation in that kidney was limited by progenitors and branching and not by perfusion.

Finally, the one known difference between the two experiments is sufficient to separate their outcomes. Doubling the division rate in a depleted cap, which is the response Muthukrishnan and colleagues measured, raised the units of the ablated construct to <!-- auto:k.resp.units.low -->0.86 [0.83, 0.89]<!-- /auto:k.resp.units.low --> of control in the same setting, within the measured 0.79–0.88. The agreement is incomplete: the pool in the model recovered only to <!-- auto:k.resp.pool -->0.65 to 0.68<!-- /auto:k.resp.pool --> of control, whereas the measured pool reached 0.85 within 2 days [6].

### The pool left at arrest decides whether growth can restart, and delivered flow sets the size regained

The bound that ends growth also determines what can restart it (Figure 5A). We relieved each constraint in arrested constructs, by lowering the perfusion threshold or the minimum spacing between units by one-fifth, and followed them for 40 days. Left alone, the liver-like construct gained 2.8 ± 0.3 % in units; with perfusion relieved it gained 22.6 ± 2.7 %, and with room relieved 32.3 ± 1.7 %. The pancreas-like construct, which still holds a quarter of its budget, gained 6.2 ± 1.4 % from the perfusion relief and 0.3 ± 0.2 % from the room relief. The kidney-like construct did not respond to either, because nothing remained to respond. Regenerative capacity is therefore not a separate property of a tissue in this model. It is what the progenitor pool has left when growth stops, which is set mostly by its turnover.

The same logic governs recovery from injury. We removed 30% of an arrested construct and restored the capacity of its units to enlarge. What followed depended on how the circuit was driven (Figure 5B and 5C). With inlet pressure held, each unit added draws more flow, and the liver-like construct overshot to 1.82 [1.54, 2.07] times its mass before the loss and the kidney-like construct to 1.85 [1.53, 2.21]. With the delivered flow held at its value before the injury, as for an organ fed by a circulation that does not rise on demand, the liver-like construct returned to 1.00 [0.95, 1.03] of its mass and the kidney-like construct to 0.97 [0.93, 0.99]. Across a 2.3-fold range of delivered flow, the mass regained rose with flow to an exponent of 0.90 (Figure 5D). No element of the tissue senses its size: the size regained is the flow delivered divided by what a unit of mass requires.

Mass and architecture recover separately (Figure 5E and 5F). The kidney-like construct, which cannot form new units, kept 0.70 [0.69, 0.70] of its units after a 30% loss and regained its mass by enlarging them to 1.39 [1.34, 1.43] times their mass, against 1.43 for full compensation; after a 50% loss the corresponding values are 1.88 [1.78, 1.94] and 2.00. This is compensatory hypertrophy, and the model overstates it: the remnant human kidney gains 30%–50% in mass and recovers 70%–80% of filtration rather than all of it [20]. Nothing in the model limits how large a unit may become.

### The model predicts the outcome of ablations that have not been performed

A model compared only with known results cannot fail. We therefore committed it to the outcome of ablations that have not been reported, in the same systems and by the same method (Figure 6A–6C; STAR Methods).

For the pancreas, removing 25% and 40% of the progenitors should leave an organ at <!-- auto:p.panc75 -->0.70 [0.58, 0.81]<!-- /auto:p.panc75 --> and <!-- auto:p.panc60 -->0.61 [0.51, 0.70]<!-- /auto:p.panc60 --> of control weight. A pancreas of normal weight at either dose would refute the model. For the liver, removing 25% should have no detectable effect (<!-- auto:p.liv75 -->0.97 [0.94, 1.00]<!-- /auto:p.liv75 -->), whereas removing 85% should leave a liver at <!-- auto:p.liv15 -->0.70 [0.68, 0.72]<!-- /auto:p.liv15 -->. In the model the liver is not indifferent to its progenitors; it has a wide margin. A liver of normal weight after an 85% ablation would refute this.

For the kidney, the extended model without a proliferative response predicts <!-- auto:p.kid20 -->0.82 [0.75, 0.89]<!-- /auto:p.kid20 --> of control nephron number after a 20% ablation and <!-- auto:p.kid60 -->0.56 [0.50, 0.63]<!-- /auto:p.kid60 --> after 60%. The second value is well above the 0.40 that proportionality would give, because the sites that remain are each better supplied as the tree thins. With a proliferative response the predictions are <!-- auto:p.kid20r -->0.93 [0.90, 0.97]<!-- /auto:p.kid20r --> and <!-- auto:p.kid60r -->0.79 [0.77, 0.82]<!-- /auto:p.kid60r -->. The model also suggests that the same 40% ablation costs somewhat fewer nephrons late in the window than early when the survivors do not compensate (<!-- auto:p.late -->0.75 [0.71, 0.79]<!-- /auto:p.late --> against <!-- auto:p.early -->0.67 [0.56, 0.79]<!-- /auto:p.early --> of control), although these two intervals overlap.

Taken together, the comparisons place each organ on a single axis, the elasticity of its size to progenitor number (Figure 6D). The liver and the kidney with compensatory proliferation lie near zero, where perfusion or renewed progenitors make up the loss. The pancreas and the kidney without compensation lie near one, where progenitors are limiting. The model places its liver-like and pancreas-like sources close to their organs, and its extended kidney-like source on the correct side in each case.

## DISCUSSION

We have asked what determines whether the size of an organ depends on the number of its progenitors, and found that one local rule is enough to produce both answers. The distinction between organs of intrinsically and extrinsically determined size [1, 2, 4] then becomes a statement about which resource is exhausted first where the tissue is built. Where perfusion runs out, size is nearly independent of progenitor number; where progenitors run out, size follows them. This recasts a classification as a measurable quantity, the elasticity of organ size to progenitor number, and it removes the need for a mechanism that counts cells or senses size.

Two consequences run against the usual reading. First, the liver-like tissue is not indifferent to its progenitors. It has a margin, and the model says where the margin ends. Second, the axis that matters most is not whether a tissue keeps a self-renewing compartment but whether its pool is lost with time. This is consistent with the view that in liver and pancreas no dedicated stem-cell niche carries out repair [10], and it shifts attention from the identity of a progenitor to its turnover.

The kidney was the informative failure. The minimal rule matched one experiment by a route the tissue does not use and missed the other by a factor of two. Adding the coupling that kidney development is known to have, progenitors that drive the branching that makes room for them [19, 5], reproduced the persistence of the cell deficit, the recovery of the cap, and the loss of tips, and it showed that the nephron deficit becomes proportional only when perfusion does not limit nucleation. We read this as a statement about the embryonic kidney: at the stage examined, nephron formation is paced by progenitors and branching. The direction of the coupling, from progenitors to branching and not from nephrons to branching, agrees with the finding that ureteric branching is unaffected when nephron formation is perturbed [21]. It also explains how two ablations of the same size can differ twofold in outcome. The difference is whether the survivors compensate by dividing, which the model needs only as a switch. Nephron endowment, which varies widely among people and is fixed before birth [8, 7], would on this account be sensitive to any insult that reduces progenitors without eliciting that response.

The same rule bears on regeneration. A tissue can resume growth only if it stopped with progenitors in hand, and the size it regains after a loss is set by the flow it is given, not by a memory of its former size. The first statement gives regenerative capacity a mechanistic definition. In models of liver regeneration, a load or a flow shared among the remaining cells plays this role by assumption [11, 12]. Here it follows from the local rule and the boundary condition on the circuit, which offers one reading of how a liver can restore its mass after resection [3] without a signal that reports its size.

For efforts to grow replacement tissue, these results suggest an order of questions. The first is which resource stopped the construct, which the two reliefs of Figure 5A distinguish. A tissue of the kidney kind must be given progenitors in proportion to the units wanted, because better perfusion makes up a shortfall only where perfusion was the limit. In human kidney organoids, the number of nephron progenitors seeded correlates with the size of the organoid [22], as expected of a tissue limited by its progenitors. Lowering what each unit requires does more than seeding more cells in every source we examined. And recovery should be judged by units as well as by mass, since the two return separately.

### Limitations of the study

No parameter of the model is fitted, and none is measured in a tissue; we report ratios, orderings, and contrasts, which is what a chosen parameter set can support. All three sources are grown in one fixed domain, whereas an embryonic organ makes its own, and model time is not calibrated to embryonic days, so that sequences and endpoints are compared and timings are not. The pancreas and liver experiment was known when the three sources were described, although their parameters were not adjusted to it. The liver comparison sets a reduced seeding beside an ablation performed on an organ already formed. The kidney extension was written after the measurements it is compared with, carries six parameters of which none is measured and two were swept, reaches the measured deficit only at the edge of its interval, and recovers the progenitor pool less completely than the embryo does when proliferation responds. In the extension the tree still sprouts toward units that are short of supply, a dependence of branching on nephrons that the embryonic kidney does not show [21]. The model has no readout of function, so every recovery reported is structural, and it places no limit on the size of a unit, so it overstates renal compensation. The human cohort gives CT volumes with segmentation error that we do not separate. Exact values hold on the platform on which they were computed; another platform reproduces them within sampling error (STAR Methods). No prediction has yet been tested.

## RESOURCE AVAILABILITY

### Lead contact

Requests for further information should be directed to and will be fulfilled by the lead contact, Daniel Pérez-Calixto (dperez@inmegen.gob.mx).

### Materials availability

This study did not generate new materials.

### Data and code availability

- This paper analyzes existing, publicly available data. The human kidney volumes are from the TRUSTED dataset [17]; the derived table used here is in the repository under `source_data/`. The values of the mouse experiments are read from the cited publications and listed in `scripts/36_validation.py` and `scripts/37_kidney_coupling.py` with the figure each is taken from.
- All original code is available at https://github.com/Danpc11/Ogrow and will be archived with a DOI at publication. `make paper` regenerates every number, table, and figure, and `PAPER_NUMBERS.md` records them with the version, commit, and platform that produced them.
- Any additional information required to reanalyze the data reported in this paper is available from the lead contact upon request.

## ACKNOWLEDGMENTS

To be completed by the authors.

## AUTHOR CONTRIBUTIONS

To be completed by the authors, using CRediT roles.

## DECLARATION OF INTERESTS

The authors declare no competing interests.

## FIGURE LEGENDS

**Figure 1. A local nucleation law stops growth without a size signal**

(A) The model. Functional units (open circles) form at the tips of a perfused branching tree (purple) that grows in a fixed domain. Each tip holds a pool of progenitors (orange). Units form at the rate shown and stop when perfusion per unit, progenitors, or free sites run out at that tip.
(B) A liver-like construct at arrest, projected onto a plane. Line width follows vessel radius.
(C) Units over time for the three cell sources. Lines, mean of 24 seeds; bands, 10th to 90th percentile.
(D) Progenitors relative to the pool seeded, as in (C).
(E) Where the progenitor pool stands at arrest: spent or lost, held at tips where perfusion has fallen to the threshold, or held at tips with no free site. Mean of 24 seeds.
See also Table 1.

**Figure 2. The turnover of the cell source decides which resource runs out**

(A) Factorial of renewal, seeded pool, and turnover. Each tile is colored by the bound that holds the largest share of the pool at arrest; numbers are units at arrest. Outlines mark the three sources of Figure 1.
(B) Fraction of the settings of the other two factors in which changing one factor changes the bound.
(C) Elasticity of the final unit count to the seeded pool (filled circles) and to the perfusion threshold (open squares; magnitude), by paired perturbation. Error bars, SEM over 24 seeds.

**Figure 3. The model beside progenitor ablations in mouse and paired kidneys in humans**

(A) The experiments and their counterpart in the model.
(B) Share of the organ that results against the share of progenitors left. Filled circles and lines, mass of the pancreas-like and liver-like constructs at day 140 relative to the reference (ratio of means, 24 paired seeds; bars, 95% bootstrap interval). Open squares, weight of the mouse pancreas at birth after ablation of 8% and 61% of its progenitors (bars, 95% CI) and of the liver 4 days after ablation of at least 65% (Stanger et al.). Dashed line, an organ in proportion to its pool.
(C) As in (B) for units of the kidney-like construct, against glomerular number after ablation of about 40% of nephron progenitors: through a CITED1 driver, at P3 and P10 (circles; Muthukrishnan et al.), and through a Gdnf driver, at E19.5 and P50 (diamonds; Cebrian et al.). Measurements are offset horizontally for clarity.
(D) Deviation of one organ from another grown like it: log mass of the kidney-like construct across 24 seeds, and log volume of each kidney relative to its pair in 48 people, with the mean difference between sides removed. Bars span one standard deviation on either side of zero.

**Figure 4. Progenitor-driven branching accounts for both kidney experiments where perfusion is not limiting**

(A) The extension. A tip branches only when its cap of progenitors has doubled, and each branch keeps one cap. After an ablation the tip waits for the cap to be restored, and the delay in branching is not made up.
(B) Total progenitors, progenitors per tip, and number of tips after ablation of 40% of progenitors, relative to the unablated control, with no proliferative response and where perfusion limits nucleation least. Open symbols at right, progenitors per tip and tips per kidney 1, 2, and 3 days after ablation in the mouse (Cebrian et al.); model days are not embryonic days.
(C) Units, ablated over control, in the same setting, with survivors dividing as before (solid) or twice as fast in a depleted cap (dashed). Gray bands, glomerular number in the study without a proliferative response (Cebrian et al.) and in the study with one (Muthukrishnan et al.).
(D) Units at the end of the window in four settings of the nucleation law, from that of the reference model (left) to the least perfusion limited (right). Bars, 95% bootstrap interval over 24 paired seeds.

**Figure 5. The pool left at arrest decides whether growth can restart, and delivered flow sets the size regained**

(A) Units gained in 40 days by arrested constructs left alone, with the perfusion threshold lowered by one-fifth, or with the spacing between units lowered by one-fifth. Error bars, SEM.
(B) Two ways to drive the circuit: inlet pressure held, or delivered flow held.
(C) Mass 80 days after loss of 30% of an arrested construct, relative to mass before the loss. Points, seeds; bars, means.
(D) Mass regained against flow delivered, both relative to their values before the loss, with flow held. Lines, fitted power laws; dashed line, proportionality.
(E) Units and mass after the loss, with flow held, relative to before.
(F) Mass per unit against units kept in the kidney-like construct after losses of 30% and 50%. Line, full compensation.

**Figure 6. Predictions for ablations that have not been performed**

(A) As Figure 3B, with the model at doses not yet tested (stars).
(B) Share of nephron-like units that results against the share of progenitors left, in the extended kidney-like source, with survivors dividing as before (filled, solid) or twice as fast (open, dashed). Squares, the four published measurements of Figure 3C.
(C) The same 40% ablation early, at the time already compared, and late in the window.
(D) Elasticity of organ size to progenitor number, $\log$(share of organ)/$\log$(share of progenitors left), for the published experiments (open squares) and for the model at the same doses (filled circles).

## STAR METHODS

### KEY RESOURCES TABLE

| REAGENT or RESOURCE | SOURCE | IDENTIFIER |
| --- | --- | --- |
| **Deposited data** | | |
| Paired human kidney CT segmentations (TRUSTED) | Ndzimbong et al. [17] | https://doi.org/10.6084/m9.figshare.27981050 |
| Pancreas and liver weight after progenitor ablation | Stanger et al. [4] | Figures 2, 3, and 6 of that study |
| Glomerular number, tips, and cap size after nephron progenitor ablation | Cebrian et al. [5] | Figures 4 and 6 of that study |
| Glomerular number and progenitor number after nephron progenitor ablation | Muthukrishnan et al. [6] | Figures 1 and 2 of that study |
| **Software and algorithms** | | |
| organgrow 4.67 (model, analysis, figures) | This paper | https://github.com/Danpc11/Ogrow; DOI to be added |
| Python 3.9.6 | Python Software Foundation | https://www.python.org |
| NumPy 1.26.4 | NumPy developers | https://numpy.org |
| SciPy 1.13.1 | SciPy developers | https://scipy.org |
| Matplotlib | Matplotlib developers | https://matplotlib.org |

### METHOD DETAILS

#### The nucleation law

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
property of the cell source below are drawn in Figure 1A.

#### The cell source

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
its pool, while neither of the other two divides at all. Results treats
the pair as one axis for that reason, and reports that it is the axis that
most often decides which bound ends the growth. All runs use
$\nu = 8$, $q_{\mathrm{set}} = 0.6$, and the same vascular and geometric
settings.

#### The vascular tree and its boundary condition

The tree grows by sprouting at tips and adapts its calibres under a wall-cost
objective $\mathcal{C} = Q^2/r^4 + k\,(p\,r^2)^b$ minimised per unit length,
in the sense of Murray [23], with wall cross-section proportional to $p\,r^2$
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
every seed, which is why Results sets $\phi$ once on an intact ensemble
and then varies $Q$. With $\xi = 1$ a unit's ceiling moves with its own mass,
so the identity above is approached rather than reached, and how closely is
measured in Results.

#### Readouts, replicates, and robustness

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

#### Published measurements, the extension, and the predictions

Three numbers of the model are set beside published measurements (Section
2.4). The seeded pool of each regime was reduced to a half and a quarter of
its reference value, and to the shares that three mouse experiments left after
removing progenitors with diphtheria toxin: 0.92 and 0.39 of pancreatic
progenitors and at most 0.35 of hepatic ones [4], and 0.6 of nephron
progenitors [5, 6]. Each reduced pool was grown over the same 24 seeds as
the reference, and mass and unit count at day 140 are reported as ratios of
means with a paired bootstrap interval. They are compared with organ weight
for pancreas and liver and with glomerular number for kidney, each as a share
of littermate controls, read from the figures and text of those papers. The
variation across seeds grown from identical parameters is compared with the
left-right variation of kidney volume in 48 people, 96 kidneys segmented from
clinical CT [17]: the standard deviation of the log ratio of the two sides,
divided by $\sqrt{2}$, after the mean difference between sides is removed.

The reference model has no path from the progenitor pool to the tree, and
Results reports an extension that adds one. In it a source divides and is
not lost to attrition, a segment of the tree branches only while it holds two
caps' worth of progenitors and hands one to the new branch, a cap below its
size may divide faster by a set factor, and the run ends on a fixed day. On
one day each progenitor dies with probability 0.4, and the same seed grown
without that is the control. The extension is switched off in every other
result of this paper, and its settings are in `PAPER_NUMBERS.md`, An
extension.

The predictions of Results use the same two procedures at doses and times
that were not published. The pancreas-like and liver-like constructs were
seeded with 0.75 and 0.60, and with 0.75 and 0.15, of their reference pools.
The extended source was run in the setting of the nucleation law in which
perfusion limits least, with 20 %, 40 % and 60 % of its progenitors killed on
day 45 of a window of 120, and with 40 % killed on days 20 and 80, each with
and without a doubled division rate in a depleted cap. They are stamped with
the commit that produced them (`PAPER_NUMBERS.md`, Predictions).

### QUANTIFICATION AND STATISTICAL ANALYSIS

Every reported quantity is computed over 24 independent seeds of the stochastic model. Unit counts and shares of the pool are given as mean ± SEM. Elasticities are central differences in log space for a perturbation of ±0.3 in the logarithm of the parameter, computed seed by seed and reported as mean ± SEM. Comparisons between a perturbed construct and its reference use the same seeds for both and are reported as the ratio of the two means with a 95% interval from 4,000 paired bootstrap resamples of the seeds; for the injury experiments, brackets give the range across seeds around the mean. The within-pair variation of human kidney volume is the standard deviation of the log ratio of left to right volume divided by $\sqrt{2}$, after subtraction of its mean; its association with kidney size is tested by Spearman correlation between the absolute deviation and the mean volume of the pair. No data were excluded. The published values with which the model is compared were read from the text and figures of the cited studies and are listed with their source in the code. Numbers in the text are regenerated from the result files by the pipeline and are exact on the platform recorded in `PAPER_NUMBERS.md` (macOS, arm64); with the same package versions, Linux gives different single realizations and the same means within sampling error, which the test suite checks against a reference ensemble.

## REFERENCES

[1] Bryant, P.J., and Simpson, P. (1984). Intrinsic and extrinsic control of growth in developing organs. Q. Rev. Biol. *59*, 387–415.

[2] Conlon, I., and Raff, M. (1999). Size control in animal development. Cell *96*, 235–244.

[3] Taub, R. (2004). Liver regeneration: from myth to mechanism. Nat. Rev. Mol. Cell Biol. *5*, 836–847.

[4] Stanger, B.Z., Tanaka, A.J., and Melton, D.A. (2007). Organ size is limited by the number of embryonic progenitor cells in the pancreas but not the liver. Nature *445*, 886–891.

[5] Cebrian, C., Asai, N., D'Agati, V., and Costantini, F. (2014). The number of fetal nephron progenitor cells limits ureteric branching and adult nephron endowment. Cell Rep. *7*, 127–137.

[6] Muthukrishnan, S.D., Ryzhov, S., Karolak, M., and Oxburgh, L. (2018). Nephron progenitor cell death elicits a limited compensatory response associated with interstitial expansion in the neonatal kidney. Dis. Model. Mech. *11*, dmm030544.

[7] Hartman, H.A., Lai, H.L., and Patterson, L.T. (2007). Cessation of renal morphogenesis in mice. Dev. Biol. *310*, 379–387.

[8] Bertram, J.F., Douglas-Denton, R.N., Diouf, B., Hughson, M.D., and Hoy, W.E. (2011). Human nephron number: implications for health and disease. Pediatr. Nephrol. *26*, 1529–1533.

[9] Martin, F., Carmona, M., Mahillo, B., Alvarez, M., Luengo, A., Chatzixiros, E., López-Fraga, M., Domínguez-Gil, B., and Tullius, S.G. (2026). Organ donation and transplantation worldwide: the Global Observatory on Donation and Transplantation 2024 report. Transplantation *110*, e655–e669.

[10] Sljukic, A., Jenkinson, J.G., Niksic, A., Prior, N., and Huch, M. (2026). Advances in liver and pancreas organoids: how far we have come and where we go next. Nat. Rev. Gastroenterol. Hepatol. *23*, 44–64.

[11] Furchtgott, L.A., Chow, C.C., and Periwal, V. (2009). A model of liver regeneration. Biophys. J. *96*, 3926–3935.

[12] Hohmann, N., Weiwei, W., Dahmen, U., Dirsch, O., Deutsch, A., and Voss-Böhme, A. (2014). How does a single cell know when the liver has reached its correct size? PLoS One *9*, e93207.

[13] Hannezo, E., Scheele, C.L.G.J., Moad, M., Drogo, N., Heer, R., Sampogna, R.V., van Rheenen, J., and Simons, B.D. (2017). A unifying theory of branching morphogenesis. Cell *171*, 242–255.

[14] Yu, W., Marshall, W.F., Metzger, R.J., Brakeman, P.R., Morsut, L., Lim, W., and Mostov, K.E. (2019). Simple rules determine distinct patterns of branching morphogenesis. Cell Syst. *9*, 221–227.

[15] Wang B, Zhao L, Fish M, Logan C Y and Nusse R 2015 Self-renewing diploid Axin2+ cells fuel homeostatic renewal of the liver *Nature* **524** 180–5

[16] Heinke, P., Rost, F., Rode, J., Trus, P., Simonova, I., et al. (2022). Diploid hepatocytes drive physiological liver renewal in adult humans. Cell Syst. *13*, 499–507.

[17] Ndzimbong, W., Fourniol, C., Themyr, L., Thome, N., Keeza, Y., Sauer, B., Piéchaud, P.-T., Méjean, A., Marescaux, J., George, D., et al. (2025). TRUSTED: the paired 3D transabdominal ultrasound and CT human data for kidney segmentation and registration research. Sci. Data *12*, 615.

[18] Kobayashi, A., Valerius, M.T., Mugford, J.W., Carroll, T.J., Self, M., Oliver, G., and McMahon, A.P. (2008). Six2 defines and regulates a multipotent self-renewing nephron progenitor population throughout mammalian kidney development. Cell Stem Cell *3*, 169–181.

[19] Costantini, F., and Kopan, R. (2010). Patterning a complex organ: branching morphogenesis and nephron segmentation in kidney development. Dev. Cell *18*, 698–712.

[20] Qu X, Shi Y, Li L and Liu T 2026 Compensatory nephron hypertrophy: mechanisms, pathophysiology dynamics, and clinical implications in acute kidney injury and chronic kidney disease *Kidney Res. Clin. Pract.* advance online publication, doi:10.23876/j.krcp.25.405

[21] Short, K.M., Combes, A.N., Lisnyak, V., Lefevre, J.G., Jones, L.K., Little, M.H., Hamilton, N.A., and Smyth, I.M. (2018). Branching morphogenesis in the developing kidney is not impacted by nephron formation or integration. eLife *7*, e38992.

[22] Porter, C.M., Qian, G.C., Grindel, S.H., and Hughes, A.J. (2024). Highly parallel production of designer organoids by mosaic patterning of progenitors. Cell Syst. *15*, 649–661.

[23] Murray C D 1926 The physiological principle of minimum work. I. The vascular system and the cost of blood volume *Proc. Natl Acad. Sci. USA* **12** 207–14

## Notes for the authors, to be removed before submission

**What was checked against the journal, and how.** The journal's page for authors says that an initial submission needs no particular format beyond a title, the authors, a summary, and enough methods to judge the results; that an Article has up to seven figures or tables; that the journal is figure driven, so the argument has to be made through the figures and their legends; that the text should be short enough to read in one sitting; and that STAR Methods and references do not count toward its length. Five research papers from the journal were then read for form (Wang et al. 2022, Kreten et al. 2024, Porter et al. 2024, Heinke et al. 2022, Beck et al. 2022). What they share, and this draft follows: a summary of 130 to 170 words; three or four highlights; result headings written as statements; an introduction of 400 to 1,100 words that ends by saying what the paper shows; and, in the 2024 papers, numbered references with DOIs and the order Discussion, Resource availability, Acknowledgments, Author contributions, Declaration of interests, STAR Methods. None of the five has a section headed "Limitations of the study"; they state limits inside the Discussion. The heading is kept here because a paper with no experiment of its own gains from stating its limits in one place, and it can be folded into the Discussion if an editor prefers. The "In brief" paragraph is written in the first person, as in those papers. Still to check on the journal's page: DOIs in the reference list, and the size of the graphical abstract.

**What an editor will look for.** Whether the biological conclusion needs the model. Here it does: the claim is that one rule yields both answers to the ablation experiments, and no experiment in the paper stands without it. The cover letter should say that in its first two sentences. A presubmission inquiry with the Summary and the six figures is available and worth sending first.

**References.** Every entry added to the full manuscript's list was taken from a paper read in full or from the reference list of one. Three need a look: the author list of Heinke et al. 2022 is abbreviated because the copy read garbled part of it; Hannezo et al. 2017 and Short et al. 2018 are cited from the reference list of Yu et al. 2019 and from what their titles state, and should be read before submission, the second especially, since a sentence of the Discussion and one of the Limitations rest on it. Add DOIs when the list is formatted.

**Numbers.** Reference-model numbers are held to the results by `scripts/29_check_manuscript.py`. Numbers from the comparison with measurements, the kidney extension, and the predictions are written between `auto:` markers by scripts 36 to 38; run those and commit this file before `make paper`. The measured values quoted from the literature are typed here and in the scripts; the two liver values are read from a fitted curve in Fig. 6b of Stanger et al. and should be confirmed.

**Still to do.** Acknowledgments, author contributions, ORCID identifiers, and a DOI for the code. A supplement with the thirteen robustness variants, the involution sweep, and the full table of the factorial. Decide whether the sentence on when the three sources were chosen relative to reading Stanger et al. can be made stronger.
