# Figure provenance

Every panel traces to a command, a seed and an input. Run `make all` to
regenerate from scratch; the intermediate files are listed so panels can be
rebuilt without repeating the ensembles.

| Figure | Command | Seeds | Inputs | Outputs |
| --- | --- | --- | --- | --- |
| 1 | `python figures/fig1.py` | none (timings median of 7 after warm-up) | `results/cost.json` | `figures/out/Fig1.pdf` |
| 2 | `python scripts/01_ensembles.py && python figures/fig2.py` | representative 0; replicates 1-10 | `config/params.yaml` | `results/ens_main.npz` |
| 3 | `python scripts/02_closure.py && python scripts/11_vascular_disorder.py && python figures/fig3.py` | 300-315; sweep 600-611 | `config/params.yaml` | `results/f3_main.npz`, `results/f3_sweep_*.npz` |
| 4 | `python scripts/03_selfaveraging.py && python figures/fig4.py` | 900-911 | `config/params.yaml` | `results/scaling_*.json` |
| 5 | `python scripts/04_identifiability.py && python figures/fig5.py` | 2000-2004 | `config/params.yaml` | `results/ident.json`, `results/sens.npz` |
| 6 | `python scripts/07_trusted.py && python scripts/08_partition.py && python figures/fig6.py` | none (no simulation) | `data/trusted/CT_masks/` | `results/partition.json` |

## Timings, Figure 1e

Measured with `scripts/00_cost.py`: three warm-up calls are discarded, then the
median of seven repeats is reported, so the first-call cost of the sparse
factorisation does not enter. Rerunning reproduces the ordering and the slope,
not the absolute milliseconds, which depend on the machine.

## Datasets excluded

`results/hoa_rejected.csv` records every dataset the quality gate refused and
why, plus the outlier rule applied afterwards (|z| > 3 on volume or
sphericity, `config/params.yaml: hoa.outlier_z`).

## Audit

`scripts/09_audit.py` prints, for every figure, the value shown beside the
value computed, and exits non-zero on any mismatch. Current state:

| Figure | Checked |
| --- | --- |
| 1 | tree-traversal cost is linear; sparse solve is 9.8 to 12.5x slower |
| 3 | sd(N) 5.01 to 1.05 under saturation; 7.24 unchanged under the other three; sigma_init matches config |
| 4 | exponents 0.96, 1.48, 1.26 additive; -0.01 and -0.22 topological; sigma_param matches config |
| 5 | rank 1/5 from final N; condition 705 to 110; 287 replicates for lambda, 41 for beta |
| 6 | 96 kidneys, 48 donors; ratio 1.541, P = 0.0005; all six invariants systematic |
