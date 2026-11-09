#!/usr/bin/env bash
# Full reproduction. See MANIFEST.md for what each step produces.
set -euo pipefail
export PYTHONPATH="${PYTHONPATH:-}:$(pwd)/src"
mkdir -p results figures/out

echo "== simulation =="
python scripts/00_cost.py
python scripts/01_ensembles.py
python scripts/02_closure.py
python scripts/03_selfaveraging.py
python scripts/04_identifiability.py

echo "== human cohorts (needs data/hoa and data/trusted; see data/README.md) =="
for organ in kidney lung; do
  [ -d data/hoa ] && python scripts/05_hoa.py "$organ" || echo "  skipping $organ: no data/hoa"
done
[ -d data/trusted ] && python scripts/07_trusted.py || echo "  skipping TRUSTED: no data/trusted"
python scripts/06_exchangeability.py || true
python scripts/08_partition.py || true

echo "== audit: every stated number must match the pipeline =="
python scripts/09_audit.py

python scripts/10_source_data.py

echo "== figures =="
for f in figures/fig*.py; do python "$f"; done
echo "done; figures in figures/out/"
