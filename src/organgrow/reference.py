"""One short run, used to test that two environments agree on a number.

`PAPER_NUMBERS.md` is produced on one machine, with one Python and one numpy,
and the paper's claim is that the repository reproduces it. Nothing checked
that across platforms: a difference in BLAS or in the sparse solver moves a
pressure in the last bits, and a nucleation that sits on the supply threshold
then falls the other way, which changes a unit count and not a rounding.

This is the smallest run that exercises the whole loop -- perfusion, the three
bounds, nucleation, turnover, growth and remodelling -- and gives one digest
over the geometry it produces. `scripts/20_paper_numbers.py` records the
digest with the environment that produced it, and `make test` compares, so
that CI on Linux checks the numbers stamped on macOS.

A mismatch is a result, not a broken test: it says the stamped numbers are not
portable to that environment, and the paper has to say so.
"""
import hashlib

import numpy as np

from . import sdf
from .closed import grow

#: short, but long enough that all three bounds act
SPEC = dict(S=0.6, days=40, seed=7, nu=8.0, q_set=0.6, rho=0.10, delta=0.10,
            P0=600.0, renewal=1.0, R_unit=650.0, units_per_branch=650 / 600.0,
            max_vessels=3000)


def digest():
    """A digest of the geometry one reference run produces."""
    p = dict(SPEC)
    S, seed = p.pop("S"), p.pop("seed")
    U, m, vt, tr = grow(lambda x, y, z: sdf("liver", x / S, y / S, z / S),
                        190 * S, rng=np.random.default_rng(seed), **p)
    h = hashlib.sha256()
    # rounded to nine places, so that a difference in the last bits of a float
    # repr is not reported as a difference in the model
    for a in (np.asarray(U, float), np.asarray(m, float),
              np.asarray(vt.radius, float)):
        h.update(np.round(a, 9).tobytes())
    h.update(f"N={len(U)} vessels={len(vt)}".encode())
    return h.hexdigest()[:16], len(U), len(vt)
