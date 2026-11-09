# Growth of exchange organs, and the limit of what can be specified

Code, data tables and figures for the manuscript.

A generative model in which liver, kidney, pancreas and lung emerge from the
same local rules, with organ identity entering as a single parameter — how the
fabrication window closes. Paired human kidneys then separate the variation in
a real organ into what the individual fixes, what its surroundings fix, and
what is contingent.

`THEORY.md` states the model formally, with every symbol defined.

## Install and reproduce

    python -m venv .venv && source .venv/bin/activate
    pip install -r requirements.txt && pip install -e .
    make all

Python 3.10+, Linux or macOS. A full run takes about 40 min on one core.
`make all` runs the analyses, checks every number the figures state against
the pipeline, and writes `figures/out/Fig1.pdf` to `Fig6.pdf`.

The human cohorts are not redistributed; `data/README.md` gives the DOI and the
layout for each. Figures 1 to 5 need no external data.

## Layout

    THEORY.md            the model, formally
    config/params.yaml   every parameter and seed
    src/organgrow/       engine, invariants, image handling, metric, style
    scripts/             analyses, numbered in dependency order
    figures/             one script per figure
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
recorded with its reason in `results/hoa_rejected.csv`.
