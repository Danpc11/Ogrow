# Figure provenance

Every panel traces to a command, a seed and an input. Run `make all` to
regenerate from scratch; the intermediate files are listed so panels can be
rebuilt without repeating the ensembles.

| Figure | Command | Seeds | Inputs | Outputs |
| --- | --- | --- | --- | --- |
| 1 | `python figures/fig1.py` | none (timings median of 7 after warm-up) | `results/cost.json` | `figures/out/Fig1.pdf` |
| 2 | `python scripts/01_ensembles.py && python figures/fig2.py` | representative 0; replicates 1-9; reference cohort 3001-3009 | `config/params.yaml` | `results/ens_main.npz`, `results/rep_main.npz`, `results/ens_snapshots.npy` |
| 3 | `python scripts/02_closure.py && python scripts/11_vascular_disorder.py && python figures/fig3.py` | 300-315; sweep 600-611; vascular 200-205, 500-511, 300-313 | `config/params.yaml` | `results/f3_main.npz`, `results/f3_sweep_*.npz`, `results/vascular_disorder.json` |
| 4 | `python scripts/03_selfaveraging.py && python figures/fig4.py` | 900-911 | `config/params.yaml` | `results/scaling_*.json` |
| 5 | `python scripts/04_identifiability.py && python figures/fig5.py` | 2000-2004 | `config/params.yaml` | `results/ident.json`, `results/sens.npz` |
| 6 | `python scripts/07_trusted.py && python scripts/08_partition.py && python figures/fig6.py` | none (no simulation) | `data/trusted/CT_masks/` | `results/partition.json` |

## Timings, Figure 1e

Measured with `scripts/00_cost.py`: three warm-up calls are discarded, then the
median of seven repeats is reported, so the first-call cost of the sparse
factorisation does not enter. Rerunning reproduces the ordering and the slope,
not the absolute milliseconds, which depend on the machine.

## Datasets excluded

`results/hoa_rejected_<organ>.csv` records every dataset the quality gate refused and
why, plus the outlier rule applied afterwards (|z| > 3 on volume or
sphericity, `config/params.yaml: hoa.outlier_z`).

Figures read only from `results/`; none of them runs the engine, so a panel
cannot differ from the source data written by `scripts/10_source_data.py`.

## Tests

`make test` runs `tests/test_smoke.py`: the engine grows every organ with the
spacing rule intact and no unit at the mass floor; the tree conserves flow,
gives dead-end tips no flow and never exceeds the tissue's conductance; the
closed model grows without a calibre runaway. Each is a property an earlier
version violated.

## Audit

`scripts/09_audit.py` prints, for every figure, the value shown beside the
value computed, and exits non-zero on any mismatch. It runs before the figures
in `make all`, so a caption cannot drift from the data. Run it to see the
current numbers; they are not reproduced here, because a table in a document
is exactly what goes stale.
