"""Quality gate: accept a segmentation only if it is plausibly an organ.

Thresholds live in config/params.yaml (hoa.qc); nothing is set here.
"""
import numpy as np
from scipy import ndimage as ndi

from .config import CFG

_Q = CFG["hoa"]["qc"]


def qc(msk, vol_mm3, dxy, dz):
    """Return (ok, reason). Rejects shells, slabs and runaway masks."""
    frac = msk.mean()
    lo, hi = _Q["mask_fraction"]
    if not (lo < frac < hi):
        return False, "mask fraction %.2f out of range" % frac
    # solidity: an organ is compact, a shell or a leak is not
    fill = ndi.binary_fill_holes(msk)
    solid = msk.sum() / max(fill.sum(), 1)
    if solid < _Q["min_solidity"]:
        return False, "hollow mask (solidity %.2f)" % solid
    # the mask must not hug the field border
    edge = np.zeros_like(msk)
    edge[0] = edge[-1] = True
    edge[:, 0] = edge[:, -1] = True
    edge[:, :, 0] = edge[:, :, -1] = True
    touch = (msk & edge).sum() / max(edge.sum(), 1)
    if touch > _Q["max_border_touch"]:
        return False, "mask touches field border (%.2f)" % touch
    vlo, vhi = _Q["volume_mm3"]
    if not (vlo < vol_mm3 < vhi):
        return False, "volume %.0f cm3 implausible" % (vol_mm3 / 1000)
    return True, "ok"
