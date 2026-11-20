"""Generative growth of exchange organs, and the invariant-space metric."""
__version__ = "4.1.0"
from .engine import grow, sdf, IDENTITY
from .sim_invariants import invariants as sim_invariants
from .metrics import epsilon, floor, effective_axes
