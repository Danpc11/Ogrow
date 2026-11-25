"""Fingerprints for cached simulations.

A cached run is only the answer to the question that produced it. File names
carry a regime and a seed, never the parameters, so changing a parameter and
rerunning silently reuses runs that answered a different question. Every cached
file therefore also carries a fingerprint of what determined it, and is
re-simulated when that fingerprint does not match.

`OGROW_ADOPT_CACHE=1` stamps whatever is on disk as current instead. It is an
assertion by the person running it that the existing files came from this code
and these parameters; nothing checks it.
"""
import hashlib
import inspect
import json
import os

from . import closed as _closed

#: fingerprint of the engine that produces the runs
CODE = hashlib.sha256(inspect.getsource(_closed).encode()).hexdigest()[:12]


def adopting():
    """Is the caller asserting that the cache on disk is current?"""
    return os.environ.get("OGROW_ADOPT_CACHE") == "1"


def fingerprint(**parts):
    """A short digest of everything that decides what a simulation returns.

    Pass whatever identifies the run: the regime, the seed, the parameters, the
    horizon. Floats are formatted rather than repr'd so that the digest does
    not move with a platform's float repr.
    """
    def norm(v):
        if isinstance(v, float):
            return f"{v:.12e}"
        if isinstance(v, dict):
            return {k: norm(x) for k, x in sorted(v.items())}
        if isinstance(v, (list, tuple)):
            return [norm(x) for x in v]
        return v

    payload = dict(code=CODE, **{k: norm(v) for k, v in sorted(parts.items())})
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()[:16]


def cached(path, fp):
    """The stored run if it answers the same question, otherwise None.

    Under OGROW_ADOPT_CACHE the stored run is stamped with `fp` and returned.
    """
    if not path.exists():
        return None
    r = json.load(open(path))
    if r.get("fingerprint") == fp:
        return r
    if adopting():
        r["fingerprint"] = fp
        json.dump(r, open(path, "w"))
        return r
    return None
