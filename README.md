# Growth of exchange organs, and the limit of what can be specified

Code, data tables and figures for the manuscript.

A generative model in which liver, kidney, pancreas and lung emerge from the
same local rules, with organ identity entering as a single parameter — how the
fabrication window closes. Paired human kidneys then separate the variation in
a real organ into what the individual fixes, what its surroundings fix, and
what is contingent.

`THEORY.md` states the model formally, with every symbol defined.
`CLOSED_MODEL.md` describes a rewrite of the engine that is under way in
`src/organgrow/closed.py`; it runs and passes its smoke tests, but is not yet
calibrated, and the figures and the audit still run on `engine.py`.

## The paper

`PAPER.md` lists each claim with the number behind it, the figure that shows
it and the command that produces it. `make paper` builds those figures from
the closed model and the human cohort tables:

| Figure | Content | Built from |
| --- | --- | --- |
| 2 | which bound stops each regime | `scripts/17_paper_runs.py` |
| 3 | which individual spread reaches organ size | `scripts/15_transmission.py` |
| 4 | calibration against the mouse liver | `scripts/17_paper_runs.py` |
| 5 | vascular rules: Murray, critical radius, reach | `scripts/17_paper_runs.py` |
| 6 | size and shape partition in 48 human donors | `source_data/` or the cohort |

`make all` still builds the earlier figures from the open engine; they are
kept for the record and are not the paper's.

## Install and reproduce

    python3 -m venv .venv && source .venv/bin/activate
    python3 -m pip install -r requirements.txt
    make test      # ten seconds: engine, tree and closed model, no data
    make all

No install step is needed: `make` puts `src` on `PYTHONPATH` and picks up
`python3` when `python` is not on the path. To force an interpreter, use
`make all PY=/path/to/python`. `pip install -e .` also works if you prefer it.

Python 3.10+, Linux or macOS. A full run takes roughly 1.5 h on one core; the vascular sweep in
`scripts/11_vascular_disorder.py` is the longest step.
`make all` runs the analyses, checks every number the figures state against
the pipeline, and writes `figures/out/Fig1.pdf` to `Fig6.pdf`.

The human cohorts are not redistributed; `data/README.md` gives the DOI and the
layout for each. Figures 1 to 5 need no external data; Figure 6 is skipped with
a message when the cohort is absent, and the build continues.

## Layout

    THEORY.md            the model, formally
    config/params.yaml   every parameter and seed
    src/organgrow/       engine, invariants, image handling, metric, style
    scripts/             analyses, numbered in dependency order
    figures/             one script per figure, drawn from results/ only
    tests/               smoke tests (make test)
    source_data/         the numbers behind every panel, one CSV per table
    results/             generated intermediates
    MANIFEST.md          figure -> command -> seed -> input
    KNOWN_LIMITATIONS.md what the work does not establish

## Source data

`source_data/` holds one CSV per table, named for the panels it supports, so a
reader can redraw any panel without running the engine. Regenerate with
`python scripts/10_source_data.py`.

## Numbers in figures

`scripts/09_audit.py` compares every number a figure states against what the
pipeline produces and exits non-zero on any mismatch. It runs before the
figures in `make all`, so a caption cannot drift from the data without the
build breaking.

## Modules

`engine.py` grows an organ in 3D from local rules: branching confined to a
signed-distance domain, unit nucleation limited by the closure mode, a second
tree that either pursues demand or runs bundled with the primary tree, and unit
growth limited by local supply.

`metrics.py` computes distances in invariant space. Note the caution in
`THEORY.md` §8: a self-normalised distance has expectation sqrt(2) for any
exchangeable sample and is a test, not a measurement.

`hoa_io.py`, `hoa_invariants.py`, `mask_invariants.py` and `qc.py` load the
human cohorts, segment them, and gate the result; every rejected dataset is
recorded with its reason in `results/hoa_rejected_<organ>.csv`.
