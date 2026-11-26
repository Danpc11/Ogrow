"""Generative growth of exchange organs, and the invariant-space metric.

The `grow` exported here is the open engine of `engine.py`, which belongs to
the earlier work this repository grew out of. The manuscript does not use it.
The model the paper runs is `organgrow.closed.grow`, which has a different
signature and is imported explicitly everywhere it is used. The export is kept
because the legacy scripts rely on it; `README.md` separates the two.
"""
__version__ = "4.67.0"
from .engine import grow, sdf, IDENTITY
from .sim_invariants import invariants as sim_invariants
from .metrics import epsilon, floor, effective_axes

__all__ = ["grow", "sdf", "IDENTITY", "sim_invariants", "epsilon", "floor",
           "effective_axes", "__version__"]
