.PHONY: all clean distclean figures analyses audit source_data test
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
	@for f in figures/fig*.py; do \
	  $(PY) $$f || echo "  $$f skipped"; \
	done

test:
	$(PY) tests/test_smoke.py

# clean keeps source_data: the Figure 6 tables come from the human cohort and
# cannot be regenerated without it
clean:
	rm -rf results/* figures/out/*

distclean: clean
	rm -rf source_data/*
