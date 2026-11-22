# The closed model (version 4.14.0)

`src/organgrow/closed.py` is the model the paper rests on; `make paper`
builds its results. Its parameters are in `config/params.yaml`, section
`closed`; the three regimes are set in the scripts that run them (below).
How the model was built, what was corrected and what was withdrawn is in
`THEORY.md` §13; this file states what the model is now.

## State

- A tree of straight segments with one radius each. The root enters at the
  deepest point of the domain (`root_site`: also `"surface"`, `"random"`).
- Units, with a position and a mass. Each unit is fed by the tip nearest to
  it (`assign_mode`: also `"segment"`, the nearest segment by exact
  projection and then the nearest tip in its subtree).
- Progenitors, a whole number per tip.

## One day

1. **Perfusion** (`hemodynamics.perfuse`). Flow is solved by series-parallel
   reduction with dead-end tips and the tissue of each tip as a load
   (conductance $m^{\xi}/R_{\text{unit}}$ per unit). Flow is conserved at
   every node, including nodes that carry tissue. The supply of a tip is
   hydraulic, $q_i = p_{\text{out},i}/\Delta P$.
2. **Calibres** (`relax`). Each radius moves toward the shear setpoint, which
   rises with pressure as $p^{b/2}$ ($b = 0.9$). The critical radius
   $r_c = (\ell/3R_{\text{below}})^{1/4}$, where shear peaks, is the floor of
   the adaptation. As vascular space fills its budget the setpoint rises.
   Vascular space is lumens plus junctions, a junction counting as a sphere
   of its radius, $4\pi r^3/3$; the budget is a fraction of the organ's
   volume. The same space and budget stop sprouting.
3. **Nucleation and progenitors.** At tip $i$ the number of units attempted
   is $\min(\text{Poisson}(\nu P_i (q_i - q_{\text{set}})_+), P_i)$, or
   Binomial$(P_i, 1 - e^{-\nu (q_i - q_{\text{set}})_+})$ with
   `nucleation="binomial"`. A unit lands within `nucleation_spread` λ of the
   tip and is refused outside the organ or closer than λ to another unit; each
   unit placed takes one cell from the competent pool. Whether that cell is
   spent is the organ's defining property, set by `renewal`: the cell stays
   competent with probability `renewal` and is spent otherwise. The liver uses
   renewal = 1, because its differentiated cells are its own progenitors;
   mononuclear hepatocytes divide and restore the compartment, so forming a
   lobule does not draw it down. The kidney and the pancreas use renewal = 0,
   because the cap mesenchyme and the multipotent progenitors are separate,
   consumable populations. The survivors then divide with probability ρ and
   are lost with probability δ, independently or as competing fates
   (`turnover="competing"`). `substeps` splits this step within the day.
4. **Units** grow logistically toward $\kappa q_i$ at a rate that decays with
   time.
5. **Sprouting**, by four rules, each gated by the junction test (a split must
   relieve more tissue than its junction occupies) and each stopped by the
   vascular budget or the vessel cap: a tip whose supply is below the
   threshold calls a vessel (deficit); a tip that feeds more units than one
   branch serves divides (load); a tip whose placements are refused sprouts in
   proportion to the share refused (crowding); and a tip extends toward
   unserved tissue, the mesenchymal signal, formalised as space colonisation
   with attractors that are removed once reached. A daughter's radius is
   $2^{-1/(b+2)}$ of its parent's, and not below $1.05\,r_c$. A tip that
   extends carries its progenitors with it.

## What the model reports

The trace of units, pool, supply, vessels and progenitors born and lost per
day. At the end, `tr.closure`: the shares of the initial pool spent, held at
tips short of supply, held for lack of room (from trial placements, 40 per
tip and 400 if none succeeds) and still active, and the residual nucleation
rate. `intervene_day` and `intervene` change the state during a run, lowering
the threshold (`"supply"`) or the least spacing (`"room"`). `injure_day` and
`injury` destroy a part of the organ past a plane cut along `injury_axis`,
holding `injury_fraction` of the units: `"parenchyma"` removes those units,
`"vascular"` occludes the vessels there so they carry no flow, and
`"combined"` does both. `injury_reactivates` is the share of the day-zero
growth rate the injury restores to the units, 0 for none and 1 for full, and
`inflow_limited` holds the flow entering the organ at what it was receiving
before the injury, in place of a fixed pressure at the inlet. The network is
linear, so the perfusion solve is rescaled rather than redone: the flow
distribution is unchanged and the driving pressure falls as the tree grows its
conductance, so the supply the model reports and the shear the calibres chase
obey the same boundary condition as the growth ceiling. Each unit's share is
its own conductance times the pressure its tip delivers, and the mass it can
support is set by the flow that reaches it. It takes effect from the injury
onward, since what it holds the inflow to is the pre-injury value. Each day the trace also records `perfused`, the mass
whose supply meets the threshold, counting a unit short of supply in
proportion to what reaches it.

## The regimes

| Regime | Initial pool | ρ | δ | renewal |
| --- | --- | --- | --- | --- |
| Liver | 2600 | 0.10 | 0.10 | 1 |
| Kidney | 300 | 0 | 0.05 | 0 |
| Pancreas | 60 | 0 | 0 | 0 |

All with ν = 8, $q_{\text{set}}$ = 0.6, domain scale 0.6, 140 days, a vessel
cap of 6000 that no run reaches, $R_{\text{unit}}$ = 650 and 650/600 units per
terminal branch (`scripts/15`, `17`, `19`, `21`, `22`).

## Tests

`tests/test_smoke.py` checks, among others, that $r_c$ is the shear peak of
`perfuse` and satisfies $\tau'(r_c) = 0$, $\tau''(r_c) < 0$; that flow is
conserved at a node carrying tissue and the root flow equals the flow the
tissue takes; and that progenitors balance exactly in every realisation, with
division and loss, binomial nucleation with substeps and competing fates.
