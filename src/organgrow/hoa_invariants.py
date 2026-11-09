"""Segment a kidney volume and compute coarse-scale invariants."""
import numpy as np
from scipy import ndimage as ndi

def segment(vol, debug=False):
    """Organ mask from the boundary gradient, not from intensity level.

    The organ/medium interface is a strong closed step whatever its polarity,
    so the background is found by flooding inward from the field border
    through low-gradient voxels. Everything the flood cannot reach is the
    organ and its internal cavities.
    """
    v = ndi.gaussian_filter(vol.astype(np.float32), 1.5)
    lo, hi = np.percentile(v, [0.5, 99.5])
    x = np.clip((v - lo) / (hi - lo), 0, 1)
    g = ndi.gaussian_gradient_magnitude(x, 1.5)
    thr = np.percentile(g, 88)
    ridge = g >= thr

    # the organ outline closes within each slice even when the 3D shell
    # leaks along the coarse axial direction, so flood slice by slice
    msk = np.zeros_like(ridge)
    st2 = np.ones((3, 3), bool)
    for k in range(ridge.shape[0]):
        r = ndi.binary_dilation(ridge[k], st2)
        free = np.pad(~r, 2, constant_values=True)
        lab, n = ndi.label(free)
        if n == 0:
            continue
        border = np.zeros_like(free)
        border[0] = border[-1] = True
        border[:, 0] = border[:, -1] = True
        ids = np.unique(lab[border & free])
        ids = ids[ids > 0]
        inside = ~np.isin(lab, ids)
        msk[k] = inside[2:-2, 2:-2]

    msk = ndi.binary_closing(msk, np.ones((3, 5, 5)))
    msk = ndi.binary_opening(msk, np.ones((3, 3, 3)))
    lab2, n2 = ndi.label(msk)
    if n2 > 1:
        sz = ndi.sum(msk, lab2, range(1, n2 + 1))
        msk = lab2 == (np.argmax(sz) + 1)
    msk = ndi.binary_fill_holes(msk)
    return (msk, x, (thr, float(x[~msk].mean()) if (~msk).any() else float('nan'))) if debug else (msk, x)


def invariants(vol, dxy, dz, seg=None):
    """Coarse-scale invariants of one volume. `seg` may pass a segmentation
    computed already (the pair `segment` returns) so it is not redone."""
    msk, x = segment(vol) if seg is None else seg
    vox = dxy * dxy * dz * 1e-9          # mm^3 per voxel
    V = msk.sum() * vox
    # surface by counting mask boundary faces
    s = 0
    for ax, d in ((0, dz), (1, dxy), (2, dxy)):
        df = np.diff(msk.astype(np.int8), axis=ax)
        area = (dxy * dz if ax != 0 else dxy * dxy) * 1e-6
        s += np.abs(df).sum() * area
    S = s                                  # mm^2
    # principal axes of the organ
    idx = np.array(np.nonzero(msk), float)
    idx[0] *= dz; idx[1] *= dxy; idx[2] *= dxy
    c = idx.mean(1, keepdims=True)
    ev = np.linalg.eigvalsh(np.cov(idx - c))[::-1]
    ax_len = 2 * np.sqrt(np.maximum(ev, 1e-9)) / 1000.0     # mm
    # interior intensity structure: cortex vs medulla contrast
    inside = x[msk]
    # vessels/sinus: darkest fraction inside the organ
    thr = np.percentile(inside, 12)
    dark = (x < thr) & msk
    lab, nlab = ndi.label(dark)
    sz = ndi.sum(dark, lab, range(1, nlab + 1)) if nlab else np.array([0.0])
    return dict(
        volume_mm3=V,
        surface_mm2=S,
        sphericity=(np.pi ** (1 / 3) * (6 * V) ** (2 / 3)) / max(S, 1e-9),
        elong=ax_len[0] / max(ax_len[1], 1e-9),
        flat=ax_len[1] / max(ax_len[2], 1e-9),
        dark_frac=float(dark.sum() / max(msk.sum(), 1)),
        dark_components=int(nlab),
        dark_largest_frac=float(sz.max() / max(dark.sum(), 1)) if nlab else 0.0,
        intensity_cv=float(inside.std(ddof=1) / max(inside.mean(), 1e-9)),
        intensity_skew=float(((inside - inside.mean()) ** 3).mean() / max(inside.std() ** 3, 1e-12)),
    )
