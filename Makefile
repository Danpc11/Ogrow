.PHONY: all clean distclean figures analyses audit source_data test closed paper check
PY ?= $(shell command -v python3 2>/dev/null || command -v python)
# the package is importable from src without installing anything
export PYTHONPATH := $(CURDIR)/src:$(PYTHONPATH)

all: figures

analyses:
	$(PY) scripts/00_cost.py
	$(PY) scripts/01_ensembles.py
	$(PY) scripts/02_closure.py
	$(PY) scripts/03_selfaveraging.py
	$(PY) scripts/04_identifiability.py
	$(PY) scripts/06_exchangeability.py
	$(PY) scripts/11_vascular_disorder.py
	$(PY) scripts/05_hoa.py kidney
	$(PY) scripts/05_hoa.py lung
	$(PY) scripts/07_trusted.py
	$(PY) scripts/08_partition.py

audit: analyses
	$(PY) scripts/09_audit.py

source_data: audit
	$(PY) scripts/10_source_data.py

figures: source_data
	mkdir -p figures/out
	@for f in figures/fig[1-6].py; do \
	  $(PY) $$f || exit 1; \
	done

test:
	$(PY) tests/test_smoke.py

# what CI runs besides the tests: static analysis, and that the version and
# the stamped numbers of the paper agree
check:
	$(PY) -m pyflakes src scripts figures tests
	$(PY) scripts/check_consistency.py
	$(PY) scripts/29_check_manuscript.py

# Exploratory runs from the development of the closed model, with parameters
# from earlier versions. The paper does not use them; see `make paper`.
closed:
	$(PY) scripts/12_calibrate_closed.py
	$(PY) scripts/13_phase_diagram.py
	$(PY) scripts/14_size_sweep.py
	for r in liver kidney pancreas; do for c in baseline pool setpoint; do \
	  OG_NREP=24 $(PY) scripts/15_transmission.py $$r $$c 0 24; done; done
	OG_NREP=24 $(PY) scripts/15_transmission.py

# clean keeps source_data: the Figure 6 tables come from the human cohort and
# cannot be regenerated without it
clean:
	rm -rf results/* figures/out/*

distclean: clean
	rm -rf source_data/*

# The paper: closed model and human cohort. The regime replicates and the
# transmission are the long parts, tens of minutes on one core each and a few
# on a multicore machine; every run is cached, so an interrupted build resumes.
paper:
	for r in liver kidney pancreas; do for c in baseline pool setpoint; do \
	  OG_NREP=24 $(PY) scripts/15_transmission.py $$r $$c 0 24; done; done
	OG_NREP=24 $(PY) scripts/15_transmission.py
	$(PY) scripts/17_paper_runs.py
	for r in liver kidney pancreas; do for s in ref nu_lo nu_hi q_lo q_hi; do \
	  $(PY) scripts/19_elasticity.py $$r $$s 0 24; done; done
	$(PY) scripts/19_elasticity.py
	$(PY) scripts/21_robustness.py reference root_surface root_random no_deficit no_load no_crowding no_mesenchyme territory binomial binomial_half_day binomial_quarter_day competing_fates no_self_renewal
	$(PY) scripts/21_robustness.py
	$(PY) scripts/22_interventions.py
	$(PY) scripts/23_scaling.py
	$(PY) scripts/24_injury.py
	$(PY) scripts/26_transcriptome.py
	# figures and tables last: they read every result above
	$(PY) figures/fig_framework.py
	$(PY) figures/fig_injury.py
	$(PY) figures/fig_transcriptome.py
	$(PY) figures/fig_design.py
	$(PY) figures/fig_closed_paper.py
	$(PY) figures/fig3_closed.py
	$(PY) figures/fig6.py
	$(PY) scripts/18_paper_tables.py
	$(PY) scripts/20_paper_numbers.py
