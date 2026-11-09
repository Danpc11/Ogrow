.PHONY: all clean figures analyses audit source_data
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
	for f in figures/fig*.py; do $(PY) $$f; done

clean:
	rm -rf results/* figures/out/* source_data/*
