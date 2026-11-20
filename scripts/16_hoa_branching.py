#!/usr/bin/env python
"""Branching exponent of the large vessels of a human liver, from HiP-CT.

Kept for the record: the radius here is the median distance transform along
each segment, which reads the minor axis of flattened ex vivo vessels, and the
exponent it gives is not a measurement (THEORY.md section 13).

The theory predicts a branching exponent b + 2 for vessels adapted to a shear
setpoint, with b near 0.9 above the Womersley crossover (pulsatile flow), so
about 2.9, against 2 for area-preserving branching. This script measures it in
a complete-organ overview from the Human Organ Atlas.

Steps: vessels are the dark lumens inside the filled parenchyma; components
that are sheets rather than tubes (volume less than 0.75 of the sum of discs
along their skeleton: gallbladder, hilum, fissures) are dropped; each
remaining tree is skeletonised, spurs shorter than twice the local radius are
pruned, and a bifurcation is kept only if all three of its branches are tubes
(radius at least RMIN voxels, length at least 3 radii, radius CV at most 0.3).
For each, the exponent x solves r1^x + r2^x = rp^x.

Usage (one zip of JP2 slices; BIN = 1 for full resolution, needs ~20 GB):
    python scripts/16_hoa_branching.py ZIP BIN
"""
import io
import json
import sys
import zipfile

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from scipy.optimize import brentq
from skimage.morphology import skeletonize

from organgrow.config import ROOT

RMIN, L_OVER_R, CV_MAX, Q_MIN = 2.0, 3.0, 0.30, 0.75
K26 = np.ones((3, 3, 3), np.uint8)


def load(zpath, b):
    zf = zipfile.ZipFile(zpath)
    names = sorted(n for n in zf.namelist() if n.endswith(".jp2"))
    names = names[:len(names) // b * b]
    first = np.asarray(Image.open(io.BytesIO(zf.read(names[0]))))
    h, w = (first.shape[0] // b) * b, (first.shape[1] // b) * b
    vol = np.zeros((len(names) // b, h // b, w // b), np.float32)
    for i, n in enumerate(names):
        a = np.asarray(Image.open(io.BytesIO(zf.read(n))))[:h, :w].astype(np.float32)
        vol[i // b] += a.reshape(h // b, b, w // b, b).mean((1, 3)) / b
    return vol


def segment(v):
    """Parenchyma and lumen levels are read from the histogram's two upper modes."""
    inside = v[v > np.percentile(v, 5)]
    par = np.percentile(inside, 90)               # parenchyma, the bright mode
    lum = np.percentile(inside, 30)               # lumens and mounting medium
    thr_p, thr_v = lum + 0.6 * (par - lum), lum + 0.45 * (par - lum)
    P = ndi.median_filter(v, 3) > thr_p
    P = ndi.binary_opening(P)
    lab, _ = ndi.label(P); sz = np.bincount(lab.ravel()); sz[0] = 0
    P = lab == np.argmax(sz)
    O = ndi.binary_closing(P, structure=ndi.generate_binary_structure(3, 1), iterations=3)
    for ax in range(3):
        O |= np.stack([ndi.binary_fill_holes(np.take(O, i, axis=ax))
                       for i in range(O.shape[ax])], axis=ax)
    V = O & ~P & (v < thr_v)
    return V, O


def prune(sk, dt, rounds=3):
    for _ in range(rounds):
        nb = ndi.convolve(sk.astype(np.uint8), K26, mode="constant") - 1
        J = sk & (nb >= 3); S = sk & ~J
        sl, ns = ndi.label(S, structure=K26)
        if ns == 0:
            break
        ends = np.unique(sl[sk & (nb == 1)]); ends = ends[ends > 0]
        touch = set(np.unique(sl[ndi.binary_dilation(J, structure=K26) & S]).tolist())
        idx = np.arange(1, ns + 1)
        cnt = ndi.sum(np.ones_like(sl), sl, idx); rmax = ndi.maximum(dt, sl, idx)
        drop = [L for L in ends if L in touch and cnt[L - 1] < max(3, 2 * rmax[L - 1])]
        if not drop:
            break
        sk &= ~np.isin(sl, drop)
    return sk


def bifurcations(m):
    m = np.pad(m, 2)
    dt = ndi.distance_transform_edt(m)
    sk = prune(skeletonize(m) > 0, dt)
    nb = ndi.convolve(sk.astype(np.uint8), K26, mode="constant") - 1
    J = sk & (nb >= 3)
    jl, nj = ndi.label(J, structure=K26)
    S = sk & ~J
    sl, ns = ndi.label(S, structure=K26)
    if not ns or not nj:
        return []
    lab, d = sl[S], dt[S]; o = np.argsort(lab, kind="stable"); lab, d = lab[o], d[o]
    idx = np.arange(1, ns + 1)
    b0, b1 = np.searchsorted(lab, idx), np.searchsorted(lab, idx, side="right")
    med = np.array([np.median(d[a:b]) if b > a else 0 for a, b in zip(b0, b1)])
    sd = np.array([d[a:b].std() if b > a else 0 for a, b in zip(b0, b1)])
    cnt = b1 - b0
    ok = (med >= RMIN) & (cnt * 1.15 >= L_OVER_R * med) & (sd / np.maximum(med, 1e-9) <= CV_MAX)
    jd = ndi.grey_dilation(jl, footprint=K26)
    touch = {}
    for j, s in set(zip(jd[S].tolist(), sl[S].tolist())):
        if j > 0:
            touch.setdefault(j, set()).add(s)
    return [sorted([float(med[s - 1]) for s in segs], reverse=True)
            for segs in touch.values() if len(segs) == 3 and all(ok[s - 1] for s in segs)]


def exponent(rp, r1, r2):
    a, b = r1 / rp, r2 / rp
    if max(a, b) >= 0.999:
        return np.nan
    return brentq(lambda x: a ** x + b ** x - 1, 0.3, 30)


def tube_ratio(m):
    m = np.pad(m, 2)
    r = ndi.distance_transform_edt(m)[skeletonize(m) > 0]
    return m.sum() / max(np.sum(np.pi * r ** 2), 1)


if __name__ == "__main__":
    zpath, b = sys.argv[1], int(sys.argv[2])
    V, O = segment(load(zpath, b))
    lab, n = ndi.label(V); sz = np.bincount(lab.ravel()); sz[0] = 0
    sl = ndi.find_objects(lab)
    B, kept, dropped = [], 0, 0
    for k in np.argsort(sz)[::-1]:
        if sz[k] < 500:
            break
        m = lab[sl[k - 1]] == k
        if tube_ratio(m) < Q_MIN:
            dropped += 1
            continue
        kept += 1
        B += bifurcations(m)
    x = np.array([exponent(*r) for r in B]); x = x[np.isfinite(x)]
    rng = np.random.default_rng(0)
    boot = [np.median(rng.choice(x, len(x))) for _ in range(2000)] if len(x) else [np.nan]
    res = dict(source=zpath.split("/")[-1], bin=b, trees=kept, sheets_dropped=dropped,
               n=int(len(x)), median=float(np.median(x)) if len(x) else None,
               ci95=[float(v) for v in np.percentile(boot, [2.5, 97.5])],
               exponents=sorted(map(float, x)), radii=B)
    out = ROOT / "results" / "hoa_branching"
    out.mkdir(parents=True, exist_ok=True)
    json.dump(res, open(out / (res["source"].replace(".zip", "") + f"_bin{b}.json"), "w"), indent=1)
    print(f"{kept} trees ({dropped} sheets dropped), {len(x)} bifurcations, "
          f"median exponent {res['median']}, 95% CI {res['ci95']}")
