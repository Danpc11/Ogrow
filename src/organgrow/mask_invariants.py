"""Shape invariants from a binary organ mask, in physical units."""
import numpy as np
from scipy import ndimage as ndi

def invariants(msk, zoom):
    dx, dy, dz = [float(z) for z in zoom]
    vox = dx*dy*dz
    V = msk.sum()*vox / 1000.0                       # cm3
    S = 0.0
    for ax, area in ((0, dy*dz), (1, dx*dz), (2, dx*dy)):
        S += np.abs(np.diff(msk.astype(np.int8), axis=ax)).sum()*area
    S /= 100.0                                       # cm2
    idx = np.array(np.nonzero(msk), float)
    idx[0] *= dx; idx[1] *= dy; idx[2] *= dz
    c = idx.mean(1, keepdims=True)
    ev = np.linalg.eigvalsh(np.cov(idx - c))[::-1]
    ax_len = 2*np.sqrt(np.maximum(ev, 1e-9))/10.0    # cm
    fill = ndi.binary_fill_holes(msk)
    return dict(
        volume_cm3=V,
        surface_cm2=S,
        sphericity=(np.pi**(1/3)*(6*V)**(2/3))/max(S, 1e-9),
        elong=ax_len[0]/max(ax_len[1], 1e-9),
        flat=ax_len[1]/max(ax_len[2], 1e-9),
        solidity=float(msk.sum()/max(fill.sum(), 1)),
        axes_cm3=float(np.prod(ax_len)),
        surf_vol=S/max(V, 1e-9),
    )

def dice(a, b):
    return 2.0*np.logical_and(a, b).sum()/max(a.sum()+b.sum(), 1)
