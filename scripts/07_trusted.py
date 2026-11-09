#!/usr/bin/env python
"""TRUSTED: split each fused kidney mask into left and right, then measure.

The published masks contain both kidneys in one volume. Shape invariants are
only comparable with the synchrotron cohort once they are separated, so the
two largest connected components are taken and ordered by position.
"""
import glob
import json
import os
import re
import numpy as np
import nibabel as nib
from scipy import ndimage as ndi
from organgrow.config import ROOT
from organgrow.mask_invariants import invariants

OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
SRC = ROOT / "data" / "trusted" / "CT_masks" / "GT_estimated_masksCT"
BIN = 2                      # shape invariants are robust to 2x binning
MIN_VOX = 150

if not SRC.exists():
    print(f"  skipped: no masks at {SRC} (see data/README.md)")
    raise SystemExit(0)

outp = OUT / "trusted_split.json"
out = json.load(open(outp)) if outp.exists() else {}
for f in sorted(glob.glob(str(SRC / "*.nii.gz"))):
    pid = re.match(r"(\d+)", os.path.basename(f)).group(1)
    if f"{pid}L" in out:
        continue
    im = nib.load(f)
    a = np.asanyarray(im.dataobj)[::BIN, ::BIN, ::BIN] > 0
    z = [BIN * float(v) for v in im.header.get_zooms()]
    if a.sum() < MIN_VOX * 4:
        continue
    lab, n = ndi.label(a)
    if n >= 2:
        sz = ndi.sum(a, lab, range(1, n + 1))
        parts = [lab == t for t in np.argsort(sz)[::-1][:2] + 1]
    else:                                   # fused: split at the midline
        xs = np.nonzero(a.any(axis=(1, 2)))[0]
        cut = int(xs.mean())
        parts = [a.copy(), a.copy()]
        parts[0][cut:] = False
        parts[1][:cut] = False
    cx = [float(np.nonzero(p.any(axis=(1, 2)))[0].mean()) for p in parts]
    for side, k in zip("LR", np.argsort(cx)):
        p = parts[k]
        if p.sum() < MIN_VOX:
            continue
        out[f"{pid}{side}"] = dict(invariants(p, z), pid=pid, side=side)
    json.dump(out, open(outp, "w"))

donors = len({v["pid"] for v in out.values()})
print(f"{len(out)} kidneys from {donors} donors -> results/trusted_split.json")
