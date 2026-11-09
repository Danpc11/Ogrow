"""Quality gate: accept a segmentation only if it is plausibly an organ."""
import numpy as np
from scipy import ndimage as ndi

def qc(msk, vol_mm3, dxy, dz):
    """Return (ok, reason). Rejects shells, slabs and runaway masks."""
    frac = msk.mean()
    if not (0.03 < frac < 0.60):
        return False, "mask fraction %.2f out of range" % frac
    # solidity: an organ is compact, a shell or a leak is not
    fill = ndi.binary_fill_holes(msk)
    solid = msk.sum() / max(fill.sum(), 1)
    if solid < 0.80:
        return False, "hollow mask (solidity %.2f)" % solid
    # the mask must not hug the field border
    edge = np.zeros_like(msk)
    edge[0] = edge[-1] = True
    edge[:, 0] = edge[:, -1] = True
    edge[:, :, 0] = edge[:, :, -1] = True
    touch = (msk & edge).sum() / max(edge.sum(), 1)
    if touch > 0.12:
        return False, "mask touches field border (%.2f)" % touch
    if not (30e3 < vol_mm3 < 400e3):
        return False, "volume %.0f cm3 implausible" % (vol_mm3 / 1000)
    return True, "ok"
