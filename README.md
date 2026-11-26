# Ogrow: perfusion-limited growth of vascularised tissue from one nucleation law

[![tests](https://github.com/Danpc11/Ogrow/actions/workflows/tests.yml/badge.svg)](https://github.com/Danpc11/Ogrow/actions/workflows/tests.yml)

Code, results and the recipe that reproduces them.

A generative model of how a perfused tissue construct stops growing. A
vascular tree whose tips hold supply, progenitors and room, and one nucleation
law: a construct arrests when one of those three runs out locally, and which
one it is decides how much of the variation in what is seeded reaches the
final size. The manuscript is a tissue-engineering paper and asks what that
costs a specification. Its three cell sources are named liver-like,
kidney-like and pancreas-like after the published turnover rates they take,
and all three grow in the same domain and at the same geometric settings, so
the contrast between them is the cell source and nothing else.

The repository is wider than the manuscript. It also holds the earlier work
this grew out of, on liver, kidney and pancreas as organs and on paired human
kidneys; those tables and figures are kept for the record, and the manuscript
does not rest on them.

`PAPER_NUMBERS.md` holds every model number, stamped with the version, the
commit and the platform that produced it; `source_data/` holds the numbers
behind every panel. `THEORY.md` states the model formally; its section 13 records how the closed
model (`src/organgrow/closed.py`) was built and corrected.

**Two pipelines, not one.** `make paper` is the manuscript: closed model,
human tables, the paper's figures and `PAPER_NUMBERS.md`. `make all` builds
the earlier figures from the open engine (`src/organgrow/engine.py`); they are
kept for the record, are not the paper's, and should not be used to check it.

The manuscript's own text is under `manuscript/`: the draft, the ledger of
claims and evidence it is written from, and the figure legends. It is
versioned here because `make check` holds its numbers to the results, and it
is left out of release archives (`.gitattributes`), so a release is the code,
the results and how to reproduce them. Without that folder `make check` says
so and checks the rest.

| Manuscript figure | File in `figures/out/` | Script |
| --- | --- | --- |
| Figure 1 | `Fig1_model.pdf` | `figures/fig_paper.py` |
| Figure 2 | `Fig2_bounds.pdf` | `figures/fig_paper.py` |
| Figure 3 | `Fig8_injury.pdf` | `figures/fig_injury.py` |
| Figure 4 | `Fig9_transcriptome.pdf` | `figures/fig_transcriptome.py` |
| Figure 5 | `Fig6_scaling.pdf` | `figures/fig_design.py` |

File names are historical and do not follow the manuscript's numbering.

## Install and reproduce

    python3 -m venv .venv && source .venv/bin/activate
    python3 -m pip install -r requirements.txt -r requirements-dev.txt
    make test      # under a minute: derivations, conservation, closed model
    make check     # static analysis; version and stamped numbers agree

When the model's version changes, run `make paper` before committing: the
regenerated `PAPER_NUMBERS.md` belongs in the same commit as the change, so
that no commit carries numbers from two versions and the check stays green.
    make paper     # the manuscript

The three heaviest steps, the elasticities, the robustness variants and the
injury experiment, run their replicates in parallel and use every core by
default. Set `OGROW_JOBS` to cap them, for instance `OGROW_JOBS=4 make paper`.
Each replicate is cached on disk, so an interrupted run resumes where it
stopped and only `make clean` forces the work again.

No install step is needed: `make` puts `src` on `PYTHONPATH` and picks up
`python3` when `python` is not on the path. To force an interpreter, use
`make all PY=/path/to/python`. `pip install -e .` also works if you prefer it.

Python 3.9 or later, Linux or macOS, as `pyproject.toml` requires and as CI
tests on 3.9 and 3.12. The stamped numbers in `PAPER_NUMBERS.md` were produced
on Python 3.9.6 on macOS (arm64), and they are exact only there: with the same
package versions Linux gives different single realisations, in unit count and
not in rounding, and the same means within sampling error. `make test`
therefore holds another platform to a reference ensemble rather than to a
digest, and `scripts/29_check_manuscript.py`, which requires the manuscript's
numbers digit for digit, is a check for the stamping machine. A full run takes a few hours on one core and
well under that on several; the vascular sweep in
`scripts/11_vascular_disorder.py` is the longest step that does not
parallelise.
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
    PAPER_NUMBERS.md     every model number, with the stamp of what produced it
    manuscript/          the draft, its ledger and legends; not in a release

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
`THEORY.md` §8: a self-normalised distance has mean square 2 (and so a mean of
at most sqrt(2)) for any exchangeable sample, and is a test, not a
measurement.

`hoa_io.py`, `hoa_invariants.py`, `mask_invariants.py` and `qc.py` load the
human cohorts, segment them, and gate the result; every rejected dataset is
recorded with its reason in `results/hoa_rejected_<organ>.csv`.
