# Coupled local constraints modulate organ size and growth arrest

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

Reference model version 4.68.0 (`src/organgrow/closed.py`). Every model number below is either held to the result files by `scripts/29_check_manuscript.py` or written between `auto:` markers by the script that computes it, so that no model number in this file is written by hand.

---

## HIGHLIGHTS

- Coupled local constraints arrest tissue growth without a global size variable
- Size follows progenitor number less as a tissue nears its perfusion-set size
- Published pancreas and liver ablation outcomes are recapitulated without optimization
- Progenitor-driven branching reconciles contrasting nephron deficits in the kidney

## IN BRIEF

How the size of an organ depends on its progenitors differs between pancreas, liver, and kidney. Using a model in which perfusion, progenitors, and space jointly limit the formation of functional units at the tips of a branching tree, we show that these coupled constraints modulate both where growth arrests and how size responds to progenitor loss, recapitulate published ablation experiments, and predict ablations not yet performed.

## GRAPHICAL ABSTRACT

`figures/out/FigA0_graphical_abstract.pdf`

## SUMMARY

How the size of an organ depends on the progenitors it is built from differs between organs: ablating pancreatic progenitors yields a proportionally smaller pancreas, whereas the liver recovers its size. We show that coupled local constraints account for this difference. On a perfused branching tree, functional units form at tips at a rate set jointly by local progenitors and by perfusion above a threshold, and growth arrests where perfusion, progenitors, or space runs out. No constraint acts alone: the dependence of size on progenitor number weakens continuously as a tissue nears the size perfusion can support. With parameters not optimized against the measurements, the model recapitulates the sizes of mouse pancreas and liver after progenitor ablation. Contrasting nephron deficits in two kidney experiments are reconciled when progenitors also drive branching, although not all progenitor dynamics are reproduced. The model predicts the outcome of ablations not yet performed.

## KEYWORDS

organ size control; growth arrest; progenitor cells; branching morphogenesis; perfusion; regeneration; nephron endowment; pancreas; liver; kidney; mathematical model

## INTRODUCTION

An organ reaches a characteristic size although none of its cells can observe the whole. Classical work divides the mechanisms into two kinds [1, 2]. In extrinsic control, systemic signals and functional demand set the size, and the tissue grows or regenerates until the demand is met; the liver, which restores its mass after resection, is the standard example [3]. In intrinsic control, the size is fixed from within and does not compensate for a loss.

The clearest experimental statement of the distinction came from ablating progenitors. Stanger and colleagues killed a controlled share of the progenitors of the mouse pancreas early in organogenesis and obtained a pancreas that was smaller in proportion and stayed so for life, whereas a liver deprived of most of its progenitors returned to normal size within days [4]. They concluded that pancreas size is limited by the number of embryonic progenitors and liver size is not. The kidney then gave both answers. Removing about 40% of nephron progenitors left a nephron deficit in proportion to the cells removed in one study [5] and a deficit less than half as large in another, in which the surviving progenitors proliferated faster [6]. Because nephron number is fixed at the end of nephrogenesis [7] and varies about tenfold among people [8], the difference matters beyond the mouse.

These results are usually read as properties of organs: some count their progenitors and some do not. That reading names the difference without deriving it, it has no place for an organ that behaves both ways, and it treats as alternatives two limits, the cells available and the demand to be met, that a growing tissue faces at the same time. It also leaves open a practical question. Kidney and liver are the two most transplanted organs, 110,467 and 42,497 of the 173,727 solid-organ transplants reported worldwide in 2024, an activity estimated to meet at most a tenth of the need [9], and efforts to grow replacement tissue from progenitors [10] must decide what such a tissue has to be supplied with: more cells, more perfusion, or more room.

Quantitative models have addressed parts of this problem. For the liver, the mass regained after resection has been modeled as the response to a metabolic load shared among the remaining cells [11] or to the portal flow that the remnant receives [12]; both take as given a signal that scales with the size of the organ. For branched organs, the form of the tree and where it stops have been derived from rules applied locally at growing tips, which branch at random and become inactive near neighboring ducts [13]. In the kidney the tips are not stopped by the nephrons that form at them [14], and a simple change in the branching rule confines them to the surface of the organ [15]. What these models leave open is how the number of progenitors enters the final size, and why it does so differently in different organs.

Here we take a different starting point. Instead of asking whether an organ measures its size, we ask how the resources available locally, where the tissue is being built, jointly limit its growth. We study a model in which functional units form at the tips of a perfused branching tree under one local rule, and growth arrests at a tip when perfusion, progenitors, or space runs out there. The model contains no variable that represents the size of the organ. We show that these constraints are coupled, so that no cell source is limited by one of them alone, and that they modulate two things together: where growth arrests, and how strongly the final size then depends on the number of progenitors, a dependence that weakens continuously as the tissue nears the size that perfusion can support. The model recapitulates the pancreas and liver experiments quantitatively, with parameters that were not optimized against them. The kidney experiments expose a further coupling that the minimal rule lacks, between progenitors and branching, which reconciles their contrasting nephron deficits but not every measured feature of the progenitor pool. The progenitors left at arrest decide whether growth can be restarted. Finally, we commit the model to the outcome of ablations that have not been performed.

## RESULTS

### A local nucleation law arrests growth without a global size variable

We represent a growing tissue as functional units, such as lobules, nephrons, or acini, that form at the tips of a perfused branching tree inside a fixed domain (Figure 1A; STAR Methods). Each tip holds a pool of progenitors, $P$. Flow through the tree is solved for a fixed inlet pressure, with every unit a hydraulic load at its tip, and $q$ is the pressure that reaches a tip as a fraction of the inlet pressure. Because the flow through a unit is its conductance times that pressure, $q$ is equally the perfusion each unit at the tip receives, relative to what it would receive at the inlet. It is a hydraulic quantity: the model does not represent oxygen or any transported solute, and the threshold stands for whatever a unit requires of its perfusion. Units form at a tip at the rate $\nu P (q - q^{*})_{+}$, where $q^{*}$ is the perfusion a unit requires. Every new unit draws on the supply of its tip, so nucleation lowers $q$, and the tree sprouts toward tissue that is short of supply. Growth ends at a tip when one of three local resources is exhausted: perfusion per unit has fallen to the threshold, the progenitors are gone, or no free site remains within reach.

Four parameters describe how the progenitor pool is maintained: whether the cell that forms a unit is spent in forming it, the size of the pool seeded, and the rates at which progenitors divide and are lost. We examine three settings, which we name for the tissues whose cell source they resemble (Figure 1; Table 1). In the liver-like source the cell is not spent, as hepatocytes are their own progenitors [16] and are renewed throughout life [17], and division balances loss. In the kidney-like source the cell is spent and the pool decays. In the pancreas-like source a small pool is spent and neither divides nor decays. The domain and every other parameter are the same for all three.

All three constructs stop growing (Figure 1B, a single realization, and Figure 1C), at 208 ± 4, 89 ± 2 and 45 ± 1 units (mean ± SEM, 24 seeds), and each stops against a different bound (Figure 1D and 1E). The liver-like construct keeps its pool: at arrest 0.78 ± 0.03 of it sits at tips where perfusion has fallen to the threshold, and 0.19 ± 0.02 at tips with no free site. The kidney-like construct has no pool left: 0.30 ± 0.01 was consumed into units and 0.70 ± 0.01 lost to attrition. The pancreas-like construct spent 0.75 ± 0.02 of its budget on units and holds the remaining 0.25 ± 0.02 at tips that perfusion no longer reaches. Size is thus the outcome of local exhaustion, and three sources arrive at it by three routes. The ordering of the three holds under 13 variations of the auxiliary rules of the model (STAR Methods).

**Table 1. The three cell sources.** All other parameters are shared.

| Source | Cell spent in forming a unit | Progenitors seeded, $P_0$ | Division, $\rho$ (per day) | Loss, $\delta$ (per day) |
| --- | --- | --- | --- | --- |
| Liver-like | no | 2,600 | 0.10 | 0.10 |
| Kidney-like | yes | 300 | 0 | 0.05 |
| Pancreas-like | yes | 60 | 0 | 0 |

### The turnover of the cell source decides which resource runs out

The three sources differ in several parameters at once, so we asked which of them sets the bound. We varied renewal (two levels), seeded pool (three), and turnover (static, decaying, or balanced) in a full factorial and recorded which bound held the largest share of the pool at arrest (Figure 2A). Changing whether the cell is spent changed the bound in 2 of the 9 settings of the other two. Changing the turnover changed it in 5 of the 6 settings of renewal and the pool, and the pool alone in 2 of 6 (Figure 2B). A pool that decays is exhausted whatever else is true of it. A pool that does not decay outlasts the perfusion that could use it, unless it is both small and spent.

The intuitive distinction, between tissues that keep a self-renewing compartment and tissues that do not, is therefore the weakest of the three axes. What decides whether progenitors or perfusion ends growth is chiefly whether the pool is lost with time.

### The dependence of size on progenitor number weakens as a tissue nears its perfusion-set size

If a construct stops because perfusion ran out, adding progenitors should change little; if it stops because progenitors ran out, its size should follow their number. We measured the elasticity of the final unit count to the seeded pool, $\varepsilon_P = \partial \ln N / \partial \ln P_0$, by paired perturbation (Figure 2C). It is +0.07 ± 0.06 in the liver-like construct, +0.37 ± 0.03 in the kidney-like, and +0.71 ± 0.05 in the pancreas-like. The elasticity to the perfusion threshold has the opposite sign and is substantial in every source: −1.83 ± 0.22, −1.65 ± 0.06, and −0.77 ± 0.07. The constraints are therefore coupled. No source is limited by one resource alone, and even the pancreas-like construct, which spends most of its pool, remains sensitive to what each unit demands, because progenitors stranded at poorly perfused tips are not converted into units.

Three sources are three points, so we measured $\varepsilon_P$ in all <!-- auto:f.cells -->18<!-- /auto:f.cells --> sources of the factorial (Figure 2D). The label of the bound sorts them only coarsely. $\varepsilon_P$ ranges from <!-- auto:f.range.spent -->+0.22 to +0.73<!-- /auto:f.range.spent --> in the <!-- auto:f.n.spent -->8<!-- /auto:f.n.spent --> sources whose pool is gone at arrest and from <!-- auto:f.range.supply -->+0.07 to +0.62<!-- /auto:f.range.supply --> in the <!-- auto:f.n.supply -->10<!-- /auto:f.n.supply --> whose pool is held by perfusion; the ranges overlap, and the correlation of $\varepsilon_P$ with the share of the pool spent is weak (Spearman <!-- auto:f.rho.spent -->+0.52<!-- /auto:f.rho.spent -->). What orders the sources is the size they reach: $\varepsilon_P$ falls steadily with the number of units at arrest (Spearman <!-- auto:f.rho.units -->-0.92<!-- /auto:f.rho.units -->). A source that stops far below the size reached by the largest pools in this domain answers to its pool, whatever stopped it, and one that approaches that size does not. The dependence of size on progenitor number is thus graded. It is set by how near the tissue comes to the size perfusion can support, and the bound that holds the pool at arrest is a coarse summary of that distance.

When pool and turnover are matched between a spent and a renewing source, the difference in $\varepsilon_P$ that remains is +0.17 at the liver-like setting, +0.04 at the kidney-like, and +0.28 at the pancreas-like. Renewal contributes, but less than half of the range across the three sources.

The elasticity $\varepsilon_P$ is what a progenitor-ablation experiment measures. We therefore compared the model with such experiments.

### Progenitors remaining at arrest predict whether growth can restart

The bound that ends growth also determines what can restart it (Figure 2E). We relieved each constraint in arrested constructs, by lowering the perfusion threshold or the minimum spacing between units by one-fifth, and followed them for 40 days. Left alone, the liver-like construct gained 2.8 ± 0.3 % in units; with perfusion relieved it gained 22.6 ± 2.7 %, and with room relieved 32.3 ± 1.7 %. The pancreas-like construct, which still holds a quarter of its budget, gained 6.2 ± 1.4 % from the perfusion relief and 0.3 ± 0.2 % from the room relief. The kidney-like construct did not respond to either, because nothing remained to respond. Regenerative capacity is therefore not a separate property of a tissue in this model. It is what the progenitor pool has left when growth stops, which is set mostly by its turnover.

### The model recapitulates the sizes of pancreas and liver after progenitor ablation

Stanger and colleagues ablated Pdx1-expressing pancreatic progenitors in the mouse embryo by tetracycline-regulated expression of diphtheria toxin and weighed the organ at birth [4]. Ablation through E9.5 removed 8% of the progenitor pool and left a pancreas at 77% of control weight (95% CI, 53%–102%); ablation through E10.5 removed 61% and left 36% (19%–53%). The small pancreas grew in parallel with controls for 11 weeks and did not catch up. With the same strategy in the liver, removal of at least 65% of hepatic progenitors left a liver at about one-third of control weight when the ablation ended (35% of littermate controls at E13.5) and at normal weight 4 days later (96% at E17.5).

We seeded each construct with the share of its pool that these ablations left and changed nothing else (Figure 3A). The pancreas-like construct ended with <!-- auto:v.panc39 -->0.44 [0.37, 0.52]<!-- /auto:v.panc39 --> of its reference mass on 39% of the pool and with <!-- auto:v.panc92 -->0.81 [0.72, 0.90]<!-- /auto:v.panc92 --> on 92% (ratio of means over 24 paired seeds, 95% bootstrap interval), against 0.36 and 0.77 measured (Figure 3B). The liver-like construct seeded with 35% of its pool ended with <!-- auto:v.liv35 -->0.86 [0.84, 0.89]<!-- /auto:v.liv35 -->, against approximately 0.96. The model thus recapitulates a pancreas whose size follows its progenitors and a liver whose size largely does not, with parameters that were set before the comparison was made.

Two qualifications apply. The liver-like construct recovers less completely than the liver, and its shortfall grows with deeper cuts: it ends with <!-- auto:v.liv50 -->0.92<!-- /auto:v.liv50 --> of its mass on half the pool and <!-- auto:v.liv25 -->0.77<!-- /auto:v.liv25 --> on a quarter. And in the embryo the liver was already built when its progenitors were removed, so that experiment is closer to an injury than to a smaller seeding; we compare only the outcome.

### Nephron number after ablation exposes a missing coupling

The kidney is the harder test, because two experiments that removed the same share of progenitors disagree. Muthukrishnan and colleagues ablated about 40% of CITED1-expressing nephron progenitors at E15.5 and found 12% fewer glomeruli at P3 and 21% fewer at P10; the surviving progenitors more than doubled their proliferation and repopulated the niche [6]. Cebrian and colleagues ablated 40% of the cap mesenchyme at E12.5 through a Gdnf driver and found 52% fewer glomeruli at E19.5 and 42% fewer in the adult; the survivors did not proliferate faster, ureteric branching slowed, and nephrogenesis ceased on schedule [5].

The kidney-like construct with 40% of its pool removed retained <!-- auto:v.kid60 -->0.80 [0.77, 0.83]<!-- /auto:v.kid60 --> of its units (Figure 3C). This agrees with the first experiment and is less than half the deficit of the second. More importantly, the model contains the mechanism of neither. Its kidney-like progenitors do not divide, so they cannot repopulate a niche, and its tree does not depend on them, so branching cannot slow. We return to this below.

### Progenitor-dependent branching reconciles contrasting nephron deficits but not all progenitor dynamics

The disagreement between the two kidney experiments points to what the minimal rule lacks. In the embryonic kidney, nephron progenitors renew themselves until a timed cessation [18, 7, 14], and the cap of progenitors around each ureteric tip supplies the signals that make the tip branch [19, 5]. We added these features as an extension that leaves the reference model unchanged (Figure 4A; STAR Methods). Progenitors divide and are not lost to attrition. A tip may branch only while it holds two caps' worth of progenitors, and it passes half to the new branch. Growth ends at a fixed time. On one day, every progenitor is killed with probability 0.4. A single further switch lets a depleted cap divide faster. Because these rules were written with the experiments in hand, the final deficit cannot count as a prediction. The test is whether one coupling yields several observations that were made separately.

It does for the branching side of the experiment of Cebrian and colleagues, in every setting of the nucleation law that we examined (Figure 4B). In the embryo, the total cap mesenchyme of ablated kidneys grew as fast as in controls, so the deficit in cells persisted; the number of progenitors per tip, 0.60 of control one day after tamoxifen, returned to 0.95 and 1.04 on the following two days; and the number of tips fell to 0.71 and then 0.59 of control [5]. In the model, the pool stayed at <!-- auto:k.pool -->0.53 to 0.57<!-- /auto:k.pool --> of control to the end of the window, the cap returned to 0.95 of control <!-- auto:k.capback -->within 3 to 50 days in all 4 settings<!-- /auto:k.capback -->, and the number of branching sites settled at <!-- auto:k.sites -->0.52 to 0.60<!-- /auto:k.sites -->. The explanation is the one those authors proposed. A tip that has lost part of its cap does not branch until the cap is restored, so the size of the cap is corrected at the expense of the tree.

Whether fewer tips produce proportionally fewer nephrons depended on perfusion (Figure 4C and 4D). With the nucleation law set as in the reference model, the units of the ablated construct recovered to <!-- auto:k.units.ref -->0.81 [0.72, 0.94]<!-- /auto:k.units.ref --> of control, because each remaining site is better perfused and forms more units. This compensation is the central mechanism of the model, and the embryonic kidney does not show it. When we lowered the perfusion threshold fourfold and slowed the law, so that perfusion seldom limits nucleation, units ended at <!-- auto:k.units.low -->0.64 [0.57, 0.70]<!-- /auto:k.units.low --> of control, and those formed after the ablation at <!-- auto:k.since.low -->0.52<!-- /auto:k.since.low -->. The measured 0.48–0.58 touches only the lower end of that interval, so the agreement is marginal. In the model, a nephron deficit that approaches proportion requires that perfusion limit nucleation weakly. This suggests, and does not show, that nephron formation in the embryonic kidney at that stage is paced by progenitors and branching more than by perfusion; the direct test is to vary perfusion in the embryo or in culture.

Finally, the one known difference between the two experiments is sufficient to separate their outcomes. Doubling the division rate in a depleted cap, which is the response Muthukrishnan and colleagues measured, raised the units of the ablated construct to <!-- auto:k.resp.units.low -->0.86 [0.83, 0.89]<!-- /auto:k.resp.units.low --> of control in the same setting, within the measured 0.79–0.88. The agreement is incomplete (Table 2): the pool in the model recovered only to <!-- auto:k.resp.pool -->0.65 to 0.68<!-- /auto:k.resp.pool --> of control, whereas the measured pool reached 0.85 within 2 days [6].

**Table 2. What the extension reproduces of the two kidney experiments.** Values are ablated over control. Model days are not embryonic days, so timings are not compared.

| Observable | Measured | Model | Outcome |
| --- | --- | --- | --- |
| Total progenitors, no proliferative response [5] | grow as fast as in control, so the deficit persists | <!-- auto:k.pool -->0.53 to 0.57<!-- /auto:k.pool --> to the end of the window | reproduced |
| Tips [5] | 0.71, then 0.59 | <!-- auto:k.sites -->0.52 to 0.60<!-- /auto:k.sites --> | reproduced |
| Progenitors per tip [5] | 0.60, then 0.95 and 1.04 | back to 0.95 <!-- auto:k.capback -->within 3 to 50 days in all 4 settings<!-- /auto:k.capback -->, after an overshoot | approximated |
| Nephrons, no proliferative response [5] | 0.48–0.58 | <!-- auto:k.units.low -->0.64 [0.57, 0.70]<!-- /auto:k.units.low --> | approximated: intervals touch |
| Nephrons, with a proliferative response [6] | 0.79–0.88 | <!-- auto:k.resp.units.low -->0.86 [0.83, 0.89]<!-- /auto:k.resp.units.low --> | reproduced |
| Progenitor pool, with a proliferative response [6] | 0.85 after 2 days | <!-- auto:k.resp.pool -->0.65 to 0.68<!-- /auto:k.resp.pool --> | not reproduced |

### The model predicts the outcome of ablations that have not been performed

A model compared only with known results cannot fail. We therefore committed it to the outcome of ablations that have not been reported, in the same systems and by the same method (Figure 5A–5C; STAR Methods).

For the pancreas, removing 25% and 40% of the progenitors should leave an organ at <!-- auto:p.panc75 -->0.70 [0.58, 0.81]<!-- /auto:p.panc75 --> and <!-- auto:p.panc60 -->0.61 [0.51, 0.70]<!-- /auto:p.panc60 --> of control weight. A pancreas of normal weight at either dose would refute the model. For the liver, removing 25% should have little effect (<!-- auto:p.liv75 -->0.97 [0.94, 1.00]<!-- /auto:p.liv75 -->), whereas removing 85% should leave a liver at <!-- auto:p.liv15 -->0.70 [0.68, 0.72]<!-- /auto:p.liv15 -->. In the model the liver is not indifferent to its progenitors; it has a wide margin. A liver of normal weight after an 85% ablation would refute this. These dose predictions move little when the rate of the law is halved or doubled and when the perfusion threshold is lowered or raised: across those four settings the pancreas-like construct ends at <!-- auto:p.panc75.alt -->0.68 to 0.74<!-- /auto:p.panc75.alt --> and <!-- auto:p.panc60.alt -->0.61 to 0.66<!-- /auto:p.panc60.alt --> of its reference, and the liver-like construct at <!-- auto:p.liv75.alt -->0.93 to 0.95<!-- /auto:p.liv75.alt --> and <!-- auto:p.liv15.alt -->0.69 to 0.73<!-- /auto:p.liv15.alt --> (Figure 5A).

For the kidney, the extended model without a proliferative response predicts <!-- auto:p.kid20 -->0.82 [0.75, 0.89]<!-- /auto:p.kid20 --> of control nephron number after a 20% ablation and <!-- auto:p.kid60 -->0.56 [0.50, 0.63]<!-- /auto:p.kid60 --> after 60%. The second value is well above the 0.40 that proportionality would give, because the sites that remain are each better supplied as the tree thins. With a proliferative response the predictions are <!-- auto:p.kid20r -->0.93 [0.90, 0.97]<!-- /auto:p.kid20r --> and <!-- auto:p.kid60r -->0.79 [0.77, 0.82]<!-- /auto:p.kid60r -->. The model also suggests that the same 40% ablation costs somewhat fewer nephrons late in the window than early when the survivors do not compensate (<!-- auto:p.late -->0.75 [0.71, 0.79]<!-- /auto:p.late --> against <!-- auto:p.early -->0.67 [0.56, 0.79]<!-- /auto:p.early --> of control), although these two intervals overlap.

Taken together, the comparisons place each organ on a single axis (Figure 5D). The axis is an effective elasticity, $\log$(share of organ)/$\log$(share of progenitors left), computed over one finite ablation; it is what an ablation yields, and it is not the differential elasticity of Figure 2. The liver and the kidney with compensatory proliferation lie near zero, where perfusion or renewed progenitors make up the loss. The pancreas and the kidney without compensation lie near one, where progenitors are limiting. The model places its liver-like and pancreas-like sources close to their organs, and its extended kidney-like source on the correct side in each case.

## DISCUSSION

We have asked what determines where the growth of an organ arrests and whether its size then depends on the number of its progenitors, and found that one local rule with three coupled constraints is enough to produce the range of answers that experiments give. The distinction between organs of intrinsically and extrinsically determined size [1, 2, 4] becomes a matter of degree: how near a tissue comes, where it is built, to the size that perfusion can support. Far from that size, the tissue follows its progenitors; near it, size is almost independent of them. Because the constraints are coupled, every source we examined remains sensitive to the perfusion threshold, and no organ is of one kind only. This recasts a classification as a measurable quantity, the elasticity of organ size to progenitor number, and it removes the need for a mechanism that counts cells or senses size.

Two consequences run against the usual reading. First, the liver-like tissue is not indifferent to its progenitors. It has a margin, and the model says where the margin ends. Second, the axis that matters most is not whether a tissue keeps a self-renewing compartment but whether its pool is lost with time. This is consistent with the view that in liver and pancreas no dedicated stem-cell niche carries out repair [10], and it shifts attention from the identity of a progenitor to its turnover.

The kidney was the informative failure. The minimal rule matched one experiment by a route the tissue does not use and missed the other by a factor of two. Adding the coupling that kidney development is known to have, progenitors that drive the branching that makes room for them [19, 5], reproduced the persistence of the cell deficit, the recovery of the cap, and the loss of tips, and it brought the nephron deficit toward proportion only where perfusion limits nucleation weakly. It does not reproduce how far the pool recovers when the survivors divide faster, and its agreement with the proportional deficit is marginal (Table 2). We therefore take the requirement for weak perfusion limitation as a hypothesis about the embryonic kidney, to be tested, and not as a result. The direction of the coupling, from progenitors to branching and not from nephrons to branching, agrees with the finding that neither the number of nephrons attached to a ureteric tip nor their formation alters its branching or the proliferation of its tip and cap cells [14]. It also explains how two ablations of the same size can differ twofold in outcome. The difference is whether the survivors compensate by dividing, which the model needs only as a switch. Nephron endowment, which varies widely among people and is fixed before birth [8, 7], would on this account be sensitive to any insult that reduces progenitors without eliciting that response.

The same rule bears on regeneration. A tissue can resume growth only if it stopped with progenitors in hand, which gives regenerative capacity a mechanistic definition: it is what the pool has left at arrest, and that is set mostly by its turnover. How an arrested tissue regains mass after injury depends in addition on how its circulation is driven, a question that models of liver regeneration address through a load or a flow shared among the remaining cells [11, 12] and that we leave to a separate study.

For efforts to grow replacement tissue, these results suggest an order of questions. The first is which resource stopped the construct, which the two reliefs of Figure 2E distinguish. A tissue of the kidney kind must be given progenitors in proportion to the units wanted, because better perfusion makes up a shortfall only where perfusion was the limit. In human kidney organoids, the number of nephron progenitors seeded correlates with the size of the organoid [20], as expected of a tissue limited by its progenitors. Lowering what each unit requires does more than seeding more cells in every source we examined.

In sum, neither the size at which a tissue stops nor its response to a loss of progenitors is set by a single limiting factor. Perfusion, progenitors, and space act together at each growing tip, and in the kidney the progenitors act on the tree itself. Organ size and growth arrest are modulated by these coupled local constraints, and the differences among pancreas, liver, and kidney follow from where each tissue sits among them. The model states that position in terms that can be measured, and it says which ablations would show it to be wrong.

### Limitations of the study

No parameter of the model was optimized against the measurements it is compared with, and none is measured in a tissue; we report ratios, orderings, and contrasts, which is what a chosen parameter set can support. All three sources are grown in one fixed domain, whereas an embryonic organ makes its own, and model time is not calibrated to embryonic days, so that sequences and endpoints are compared and timings are not. The pancreas and liver experiment was known when the three sources were defined, so those comparisons are retrospective: they show that the model recapitulates the experiment, and they are not blind predictions. The liver comparison sets a reduced seeding beside an ablation performed on an organ already formed. The kidney extension was written after the measurements it is compared with, carries six parameters of which none is measured and two were swept, reaches the measured deficit only at the edge of its interval, and recovers the progenitor pool less completely than the embryo does when proliferation responds. In the extension the tree still sprouts toward units that are short of supply, a dependence of branching on nephrons that the embryonic kidney does not show [14]. The model has no readout of function, so every recovery reported is structural. The variable $q$ is hydraulic; whether a threshold in it stands for oxygen, nutrients, or a flow-dependent signal is not resolved by the model. Exact values hold on the platform on which they were computed; another platform reproduces them within sampling error (STAR Methods). No prediction has yet been tested.

## RESOURCE AVAILABILITY

### Lead contact

Requests for further information should be directed to and will be fulfilled by the lead contact, Daniel Pérez-Calixto (dperez@inmegen.gob.mx).

### Materials availability

This study did not generate new materials.

### Data and code availability

- This paper analyzes existing, publicly available data. The human kidney volumes are from the TRUSTED dataset [21]; the derived table used here is in the repository under `source_data/`. The values of the mouse experiments are read from the cited publications and listed in `scripts/36_validation.py` and `scripts/37_kidney_coupling.py` with the figure each is taken from.
- All original code is available at https://github.com/Danpc11/Ogrow and will be archived with a DOI at publication. `make paper` regenerates every number, table, and figure, and `PAPER_NUMBERS.md` records them with the version, commit, and platform that produced them.
- Any additional information required to reanalyze the data reported in this paper is available from the lead contact upon request.

## ACKNOWLEDGMENTS

To be completed by the authors.

## AUTHOR CONTRIBUTIONS

To be completed by the authors, using CRediT roles.

## DECLARATION OF INTERESTS

The authors declare no competing interests.

## FIGURE LEGENDS

**Figure 1. A local nucleation law arrests growth without a global size variable**

(A) The model. Functional units (open circles) form at the tips of a perfused branching tree (purple) that grows in a fixed domain. Each tip holds a pool of progenitors (orange). Units form at the rate shown and stop when perfusion per unit, progenitors, or free sites run out at that tip.
(B) One liver-like construct at arrest, projected onto a plane; a single realization, whereas the text gives the mean of 24. Line width follows vessel radius.
(C) Units over time for the three cell sources. Lines, mean of 24 seeds; bands, 10th to 90th percentile.
(D) Progenitors relative to the pool seeded, as in (C).
(E) Where the progenitor pool stands at arrest: spent or lost, held at tips where perfusion has fallen to the threshold, or held at tips with no free site. Mean of 24 seeds.
See also Table 1.

**Figure 2. The dependence of size on progenitor number weakens as a tissue nears its perfusion-set size**

(A) Factorial of renewal, seeded pool, and turnover. Each tile is colored by the bound that holds the largest share of the pool at arrest; numbers are units at arrest. Outlines mark the three sources of Figure 1.
(B) Fraction of the settings of the other two factors in which changing one factor changes the bound.
(C) Elasticity of the final unit count to the seeded pool (filled circles) and to the perfusion threshold (open squares; magnitude) in the three sources, by paired perturbation. Error bars, SEM over 24 seeds.
(D) Elasticity to the seeded pool in all 18 sources of (A), against units at arrest. Fill, the bound that holds the pool at arrest; outlines, the three sources of Figure 1. Error bars, SEM over 24 paired seeds.
(E) Units gained in 40 days by arrested constructs left alone, with the perfusion threshold lowered by one-fifth, or with the minimum spacing between units lowered by one-fifth. Error bars, SEM.

**Figure 3. The model beside published progenitor ablations in the mouse**

(A) The experiments and their counterpart in the model.
(B) Share of the organ that results against the share of progenitors left. Filled circles and lines, mass of the pancreas-like and liver-like constructs at day 140 relative to the reference (ratio of means, 24 paired seeds; bars, 95% bootstrap interval). Open squares, weight of the mouse pancreas at birth after ablation of 8% and 61% of its progenitors (bars, 95% CI) and of the liver 4 days after ablation of at least 65% [4]. Dashed line, an organ in proportion to its pool.
(C) As in (B) for units of the kidney-like construct, against glomerular number after ablation of about 40% of nephron progenitors: through a CITED1 driver, at P3 and P10 (circles) [6], and through a Gdnf driver, at E19.5 and P50 (diamonds) [5]. Measurements are offset horizontally for clarity.
See also Figure S1.

**Figure 4. Progenitor-dependent branching reconciles contrasting nephron deficits but not all progenitor dynamics**

(A) The extension. A tip branches only when its cap of progenitors has doubled, and each branch keeps one cap. After an ablation the tip waits for the cap to be restored, and the delay in branching is not made up.
(B) Total progenitors, progenitors per tip, and number of tips after ablation of 40% of progenitors, relative to the unablated control, with no proliferative response and where perfusion limits nucleation least. The horizontal axis is in model days. Open symbols at right, progenitors per tip and tips per kidney 1, 2, and 3 days after ablation in the mouse [5], drawn beside the axis and not on it because the two time scales are not calibrated to each other.
(C) Units, ablated over control, in the same setting, with survivors dividing as before (solid) or twice as fast in a depleted cap (dashed). Gray bands, glomerular number in the study without a proliferative response [5] and in the study with one [6].
(D) Units at the end of the window in four settings of the nucleation law, from that of the reference model (left) to the least perfusion limited (right). Bars, 95% bootstrap interval over 24 paired seeds.
See also Table 2.

**Figure 5. Predictions for ablations that have not been performed**

(A) As Figure 3B, with the model at doses not yet tested (stars). Pale bars beside the stars, the range of the same prediction across four other settings of the nucleation law.
(B) Share of nephron-like units that results against the share of progenitors left, in the extended kidney-like source, with survivors dividing as before (filled, solid) or twice as fast (open, dashed). Squares, the four published measurements of Figure 3C.
(C) The same 40% ablation early, at the time already compared, and late in the window. The intervals for early and late overlap.
(D) Effective elasticity of organ size to progenitor number over one finite ablation, $\log$(share of organ)/$\log$(share of progenitors left), for the published experiments (open squares) and for the model at the same doses (filled circles). This is not the differential elasticity of Figure 2.

## SUPPLEMENTAL INFORMATION

**Figure S1. Two organs grown alike, related to Figure 3**

Deviation of one organ from another grown like it: log mass of the kidney-like construct across 24 seeds, and log volume of each kidney relative to its pair in 48 people [21], with the mean difference between sides removed. The standard deviation of each group is given above it.

**Note S1. Variation between two organs grown alike**

Two constructs grown from identical parameters differ only by chance, as do the two kidneys of one person. Across seeds, the standard deviation of the log mass of the kidney-like construct is <!-- auto:v.sd.model -->0.104<!-- /auto:v.sd.model -->. In CT volumes of <!-- auto:v.pairs.n -->48<!-- /auto:v.pairs.n --> pairs of human kidneys [21], the within-pair standard deviation of log volume, after removal of the mean difference between sides, is <!-- auto:v.sd.pairs -->0.098<!-- /auto:v.sd.pairs --> (Figure S1). Nothing in the model was set by the cohort. The agreement is nonetheless limited: part of the measured variation is segmentation error; a mass in model units is compared with a volume; and the variation in the model grows with the size of its domain, whereas within-pair variation in the cohort does not rise detectably with kidney volume (Spearman <!-- auto:v.pairs.size -->$\rho = +0.25$, $P = 0.09$<!-- /auto:v.pairs.size -->). We report the comparison for completeness and draw no conclusion from it: the two kidneys of one person share a body, and no hypothesis here connects the variation of the model across seeds to bilateral variation in people.

## STAR METHODS

### KEY RESOURCES TABLE

| REAGENT or RESOURCE | SOURCE | IDENTIFIER |
| --- | --- | --- |
| **Deposited data** | | |
| Paired human kidney CT segmentations (TRUSTED) | Ndzimbong et al. [21] | https://doi.org/10.6084/m9.figshare.27981050 |
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

where $P_i$ is the number of competent cells held at that tip, $q_i$ is the pressure delivered at tip $i$ as a fraction of the inlet pressure, which is also the flow through each unit at that tip relative to the flow the same unit would draw at the inlet (defined below), $q_{\mathrm{set}}$ is a
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
in the sense of Murray [22], with wall cross-section proportional to $p\,r^2$
by Laplace's law. This gives the branching relation
$r_p^{\,b+2} = \sum r_c^{\,b+2}$ and wall shear
proportional to $p^{b/2} r^{\,b-1}$. Shear on a segment in series with its
tissue peaks at $r_c = (\ell / 3R_t)^{1/4}$, the radius below which a vessel
regresses.

Flow is solved for a fixed inlet pressure $\Delta P$. Each segment has a
Poiseuille conductance proportional to $r^4/\ell$, each unit is a terminal load
of resistance $R_u$ at its tip, and a tip that serves no unit is a dead end.
Conductances are reduced from the leaves to the root, the root draws the flow
its total conductance allows, and that flow divides among daughters in
proportion to their conductances. The pressure at the distal end of tip $i$,
$p_i$, is the pressure across the units it serves, and

$$q_i = p_i / \Delta P .$$

The flow through one unit at that tip is $p_i/R_u$, so $q_i$ is also that flow
in units of $\Delta P/R_u$, the flow a unit would draw if it were connected
directly to the inlet. The variable lies between 0 and 1 and is purely
hydraulic: no solute is transported, and the threshold $q_{\mathrm{set}}$
stands for whatever a unit requires of its perfusion. The tree can also be
driven with a fixed flow; results that depend on that choice are not part of
this paper.

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
clinical CT [21]: the standard deviation of the log ratio of the two sides,
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

Every reported quantity is computed over 24 independent seeds of the stochastic model. Unit counts and shares of the pool are given as mean ± SEM. Elasticities are central differences in log space for a perturbation of ±0.3 in the logarithm of the parameter, computed seed by seed and reported as mean ± SEM. Comparisons between a perturbed construct and its reference use the same seeds for both and are reported as the ratio of the two means with a 95% interval from paired bootstrap resamples of the seeds. Associations between the elasticity and properties of a source across the 18 sources of the factorial are Spearman rank correlations. The within-pair variation of human kidney volume is the standard deviation of the log ratio of left to right volume divided by $\sqrt{2}$, after subtraction of its mean; its association with kidney size is tested by Spearman correlation between the absolute deviation and the mean volume of the pair. No data were excluded. The published values with which the model is compared were read from the text and figures of the cited studies and are listed with their source in the code. Numbers in the text are regenerated from the result files by the pipeline and are exact on the platform recorded in `PAPER_NUMBERS.md` (macOS, arm64); with the same package versions, Linux gives different single realizations and the same means within sampling error, which the test suite checks against a reference ensemble.

## REFERENCES

[1] Bryant, P.J., and Simpson, P. (1984). Intrinsic and extrinsic control of growth in developing organs. Q. Rev. Biol. *59*, 387–415.

[2] Conlon, I., and Raff, M. (1999). Size control in animal development. Cell *96*, 235–244.

[3] Taub, R. (2004). Liver regeneration: from myth to mechanism. Nat. Rev. Mol. Cell Biol. *5*, 836–847.

[4] Stanger, B.Z., Tanaka, A.J., and Melton, D.A. (2007). Organ size is limited by the number of embryonic progenitor cells in the pancreas but not the liver. Nature *445*, 886–891. https://doi.org/10.1038/nature05537.

[5] Cebrian, C., Asai, N., D'Agati, V., and Costantini, F. (2014). The number of fetal nephron progenitor cells limits ureteric branching and adult nephron endowment. Cell Rep. *7*, 127–137. https://doi.org/10.1016/j.celrep.2014.02.033.

[6] Muthukrishnan, S.D., Ryzhov, S., Karolak, M., and Oxburgh, L. (2018). Nephron progenitor cell death elicits a limited compensatory response associated with interstitial expansion in the neonatal kidney. Dis. Model. Mech. *11*, dmm030544. https://doi.org/10.1242/dmm.030544.

[7] Hartman, H.A., Lai, H.L., and Patterson, L.T. (2007). Cessation of renal morphogenesis in mice. Dev. Biol. *310*, 379–387.

[8] Bertram, J.F., Douglas-Denton, R.N., Diouf, B., Hughson, M.D., and Hoy, W.E. (2011). Human nephron number: implications for health and disease. Pediatr. Nephrol. *26*, 1529–1533.

[9] Martin, F., Carmona, M., Mahillo, B., Alvarez, M., Luengo, A., Chatzixiros, E., López-Fraga, M., Domínguez-Gil, B., and Tullius, S.G. (2026). Organ donation and transplantation worldwide: the Global Observatory on Donation and Transplantation 2024 report. Transplantation *110*, e655–e669.

[10] Sljukic, A., Jenkinson, J.G., Niksic, A., Prior, N., and Huch, M. (2026). Advances in liver and pancreas organoids: how far we have come and where we go next. Nat. Rev. Gastroenterol. Hepatol. *23*, 44–64. https://doi.org/10.1038/s41575-025-01116-1.

[11] Furchtgott, L.A., Chow, C.C., and Periwal, V. (2009). A model of liver regeneration. Biophys. J. *96*, 3926–3935. https://doi.org/10.1016/j.bpj.2009.01.061.

[12] Hohmann, N., Weiwei, W., Dahmen, U., Dirsch, O., Deutsch, A., and Voss-Böhme, A. (2014). How does a single cell know when the liver has reached its correct size? PLoS One *9*, e93207. https://doi.org/10.1371/journal.pone.0093207.

[13] Hannezo, E., Scheele, C.L.G.J., Moad, M., Drogo, N., Heer, R., Sampogna, R.V., van Rheenen, J., and Simons, B.D. (2017). A unifying theory of branching morphogenesis. Cell *171*, 242–255. https://doi.org/10.1016/j.cell.2017.08.026.

[14] Short, K.M., Combes, A.N., Lisnyak, V., Lefevre, J.G., Jones, L.K., Little, M.H., Hamilton, N.A., and Smyth, I.M. (2018). Branching morphogenesis in the developing kidney is not impacted by nephron formation or integration. eLife *7*, e38992. https://doi.org/10.7554/eLife.38992.

[15] Yu, W., Marshall, W.F., Metzger, R.J., Brakeman, P.R., Morsut, L., Lim, W., and Mostov, K.E. (2019). Simple rules determine distinct patterns of branching morphogenesis. Cell Syst. *9*, 221–227. https://doi.org/10.1016/j.cels.2019.08.001.

[16] Wang, B., Zhao, L., Fish, M., Logan, C.Y., and Nusse, R. (2015). Self-renewing diploid Axin2+ cells fuel homeostatic renewal of the liver. Nature *524*, 180–185.

[17] Heinke, P., Rost, F., Rode, J., Trus, P., Simonova, I., Lázár, E., Feddema, J., Welsch, T., Alkass, K., Salehpour, M., et al. (2022). Diploid hepatocytes drive physiological liver renewal in adult humans. Cell Syst. *13*, 499–507. https://doi.org/10.1016/j.cels.2022.05.001.

[18] Kobayashi, A., Valerius, M.T., Mugford, J.W., Carroll, T.J., Self, M., Oliver, G., and McMahon, A.P. (2008). Six2 defines and regulates a multipotent self-renewing nephron progenitor population throughout mammalian kidney development. Cell Stem Cell *3*, 169–181.

[19] Costantini, F., and Kopan, R. (2010). Patterning a complex organ: branching morphogenesis and nephron segmentation in kidney development. Dev. Cell *18*, 698–712.

[20] Porter, C.M., Qian, G.C., Grindel, S.H., and Hughes, A.J. (2024). Highly parallel production of designer organoids by mosaic patterning of progenitors. Cell Syst. *15*, 649–661. https://doi.org/10.1016/j.cels.2024.06.004.

[21] Ndzimbong, W., Fourniol, C., Themyr, L., Thome, N., Keeza, Y., Sauer, B., Piéchaud, P.-T., Méjean, A., Marescaux, J., George, D., et al. (2025). TRUSTED: the paired 3D transabdominal ultrasound and CT human data for kidney segmentation and registration research. Sci. Data *12*, 615. https://doi.org/10.1038/s41597-025-04467-1.

[22] Murray, C.D. (1926). The physiological principle of minimum work. I. The vascular system and the cost of blood volume. Proc. Natl. Acad. Sci. USA *12*, 207–214.

## Notes for the authors, to be removed before submission

**What was checked against the journal, and how.** The journal's page for authors says that an initial submission needs no particular format beyond a title, the authors, a summary, and enough methods to judge the results; that an Article has up to seven figures or tables; that the journal is figure driven, so the argument has to be made through the figures and their legends; that the text should be short enough to read in one sitting; and that STAR Methods and references do not count toward its length. Five research papers from the journal were then read for form (Wang et al. 2022, Kreten et al. 2024, Porter et al. 2024, Heinke et al. 2022, Beck et al. 2022). What they share, and this draft follows: a summary of 130 to 170 words; three or four highlights; result headings written as statements; an introduction of 400 to 1,100 words that ends by saying what the paper shows; and, in the 2024 papers, numbered references with DOIs and the order Discussion, Resource availability, Acknowledgments, Author contributions, Declaration of interests, STAR Methods. None of the five has a section headed "Limitations of the study"; they state limits inside the Discussion. The heading is kept here because a paper with no experiment of its own gains from stating its limits in one place, and it can be folded into the Discussion if an editor prefers. The "In brief" paragraph is written in the first person, as in those papers. Still to check on the journal's page: DOIs in the reference list, and the size of the graphical abstract.

**Changes made in answer to the internal review (version 4.68).** (1) "Without fitting" is replaced throughout by "parameters not optimized against the measurements", and the Limitations say the pancreas and liver comparisons are retrospective. (2) The elasticity was measured in all 18 sources of the factorial. The label of the dominant bound predicts it poorly; the size reached predicts it well. The central claim is restated accordingly as a graded one, in the Summary, the third Results section, Figure 2D and the Discussion. (3) The variable q is defined: tip pressure over inlet pressure, equal to flow per unit in units of the flow at the inlet, and hydraulic only. (4) The kidney section no longer says the extension "accounts for both experiments"; Table 2 lists what is reproduced, approximated and not reproduced, and the statement about perfusion in the embryonic kidney is a hypothesis. (5) The comparison with paired human kidneys is moved to Figure S1 and Note S1. (6) The section on injury, delivered flow and compensatory enlargement is removed, with its figure; it belongs to the engineering manuscript, and the identity M* = Q/phi is no longer presented here. The restart result stays, as Figure 2E. (7) The dose predictions are given with their range across four other settings of the law. (8) Figure 5D is labelled as an effective elasticity over a finite ablation.

**The title.** Changed to "Coupled local constraints modulate organ size and growth arrest". The Summary, the last two paragraphs of the Introduction, and the first and last paragraphs of the Discussion were rewritten to carry it: the constraints act together, they modulate where growth arrests and how size answers the pool, and no organ is of one kind only. The graded result of Figure 2D is what supports "modulate" in place of "determine".

**References, as checked.** DOIs are given for the twelve entries whose DOI was read off the paper or its record; add the rest when the list is formatted. Four papers that sentences here depend on were read in full for this version. Hannezo et al. 2017: tips branch at random and become inactive near neighbouring ducts, in mammary gland, kidney and prostate; for the kidney it proposes that branching stops through nephron differentiation. Short et al. 2018 tested that proposal and found against it: nearly all tips carry a nephron while they go on branching, and nephron number per tip relates neither to branching history nor to the proliferation of tip or cap cells. The study is observational, so the Discussion no longer says nephron formation was "perturbed". The Introduction now cites the two together. Heinke et al. 2022: author list completed. Stanger et al. 2007: the liver values were confirmed from the fitted curve of its Fig. 6b, y = 1.839x² + 4.3992x + 28.407 at x = 1, 3, 5, which gives 35 %, 58 % and 96 % of littermate control at E13.5, E15.5 and E17.5; the supplementary figures that give the 65 % reduction of the progenitor pool were not available.

**What an editor will look for.** Whether the biological conclusion needs the model. Here it does: the claim is that one rule yields both answers to the ablation experiments, and no experiment in the paper stands without it. The cover letter should say that in its first two sentences. A presubmission inquiry with the Summary and the six figures is available and worth sending first.

**Numbers.** Reference-model numbers are held to the results by `scripts/29_check_manuscript.py`. Numbers from the comparison with measurements, the kidney extension, and the predictions are written between `auto:` markers by scripts 36 to 38; run those and commit this file before `make paper`. The measured values quoted from the literature are typed here and in the scripts; the two liver values were read from the fitted curve in Fig. 6b of Stanger et al. and confirmed against its equation.

**Still to do.** Acknowledgments, author contributions, ORCID identifiers, and a DOI for the code. A supplement with the thirteen robustness variants, the involution sweep, and the full table of the factorial. Decide whether the sentence on when the three sources were chosen relative to reading Stanger et al. can be made stronger.
