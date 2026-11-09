#!/usr/bin/env bash
# Full reproduction. See MANIFEST.md for what each step produces.
set -euo pipefail
PY="${PY:-$(command -v python3 || command -v python)}"
export PYTHONPATH="$(cd "$(dirname "$0")" && pwd)/src:${PYTHONPATH:-}"
mkdir -p results figures/out

echo "== simulation =="
"$PY" scripts/00_cost.py
"$PY" scripts/01_ensembles.py
"$PY" scripts/02_closure.py
"$PY" scripts/03_selfaveraging.py
"$PY" scripts/04_identifiability.py
"$PY" scripts/11_vascular_disorder.py

echo "== human cohorts (needs data/hoa and data/trusted; see data/README.md) =="
for organ in kidney lung; do
  [ -d data/hoa ] && "$PY" scripts/05_hoa.py "$organ" || echo "  skipping $organ: no data/hoa"
done
[ -d data/trusted ] && "$PY" scripts/07_trusted.py || echo "  skipping TRUSTED: no data/trusted"
"$PY" scripts/06_exchangeability.py || true
"$PY" scripts/08_partition.py || true

echo "== audit: every stated number must match the pipeline =="
"$PY" scripts/09_audit.py

"$PY" scripts/10_source_data.py

echo "== figures =="
for f in figures/fig*.py; do "$PY" "$f" || echo "  $f skipped"; done
echo "done; figures in figures/out/"
