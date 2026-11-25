"""Fingerprints for cached simulations.

A cached run is only the answer to the question that produced it. File names
carry a regime and a seed, never the parameters, so changing a parameter and
rerunning silently reuses runs that answered a different question. Every cached
file therefore also carries a fingerprint of what determined it, and is
re-simulated when that fingerprint does not match.

The fingerprint covers the parameters the caller passes and, through `CODE`,
everything else the answer depends on: the source of `closed.py` and of
`hemodynamics.py`, the bytes of `config/params.yaml`, and the versions of
numpy and scipy. The whole config file is hashed rather than the sections the
model reads, because the scripts read their own keys from it and because the
stamp in `PAPER_NUMBERS.md` is a hash of the same bytes: the two then move
together. The cost of that choice is that editing a comment in `params.yaml`
invalidates every cached run.

`OGROW_ADOPT_CACHE=1` stamps whatever is on disk as current instead. It is an
assertion by the person running it that the existing files came from this code
and these parameters; nothing checks it. It is the way out of a recompute that
is known to be unnecessary: after a change proved not to move any number, or
after an edit to a comment.

Writes go through `save_json` and `save_npz`, which write a temporary file and
rename it over the target. A killed process therefore leaves either the old
file or the new one, never a truncated file that `exists()` would read as a
finished run.
"""
import hashlib
import importlib
import inspect
import json
import os

import numpy as np

from . import closed as _closed
from . import hemodynamics as _hemo
from .config import ROOT


def _version(name):
    try:
        return importlib.import_module(name).__version__
    except Exception:                                    # pragma: no cover
        return "absent"


def _code_digest():
    """A digest of everything outside the caller that decides a run."""
    h = hashlib.sha256()
    for mod in (_closed, _hemo):
        h.update(inspect.getsource(mod).encode())
    h.update((ROOT / "config" / "params.yaml").read_bytes())
    for name in ("numpy", "scipy"):
        h.update(f"{name}={_version(name)}".encode())
    return h.hexdigest()[:12]


#: fingerprint of the engine, the configuration and the libraries
CODE = _code_digest()


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


def _atomic(path, write):
    """Write through a temporary file, so an interrupted run leaves no stub."""
    tmp = path.with_name(f"{path.name}.{os.getpid()}.tmp")
    try:
        write(tmp)
        os.replace(tmp, path)
    except BaseException:
        try:
            tmp.unlink()
        except OSError:
            pass
        raise


def save_json(path, obj, fp):
    """Store a run and its fingerprint. Returns what was stored."""
    r = dict(obj, fingerprint=fp)
    _atomic(path, lambda t: t.write_text(json.dumps(r)))
    return r


def save_npz(path, fp, **arrays):
    """Store arrays and their fingerprint."""
    def write(t):
        with open(t, "wb") as fh:
            np.savez(fh, fingerprint=np.array(fp), **arrays)
    _atomic(path, write)


def cached(path, fp):
    """The stored run if it answers the same question, otherwise None.

    Under OGROW_ADOPT_CACHE the stored run is stamped with `fp` and returned.
    A file that cannot be parsed is treated as absent: a run interrupted before
    this module wrote atomically could have left one.
    """
    if not path.exists():
        return None
    try:
        with open(path) as fh:
            r = json.load(fh)
    except Exception:
        return None
    if not isinstance(r, dict):
        return None
    if r.get("fingerprint") == fp:
        return r
    if adopting():
        return save_json(path, {k: v for k, v in r.items() if k != "fingerprint"}, fp)
    return None


def cached_npz(path, fp):
    """Does the stored .npz answer the same question?

    Unlike `cached` this does not return the arrays: the callers load the
    entries they want, and holding every array of every replicate open would
    cost more than loading one twice.
    """
    if not path.exists():
        return False
    try:
        with np.load(path, allow_pickle=False) as z:
            if "fingerprint" in z.files and str(z["fingerprint"]) == fp:
                return True
            if not adopting():
                return False
            keep = {k: z[k] for k in z.files if k != "fingerprint"}
    except Exception:
        # However it fails to read, the answer is the same: simulate it again.
        # A truncated .npz raises zipfile.BadZipFile, which is not an OSError.
        return False
    save_npz(path, fp, **keep)
    return True


def stale(path, fp):
    """Does this run have to be simulated? The one test the scripts call."""
    if str(path).endswith(".npz"):
        return not cached_npz(path, fp)
    return cached(path, fp) is None
