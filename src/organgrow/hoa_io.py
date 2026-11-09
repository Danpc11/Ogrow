"""Load an HOA complete-organ stack from its ZIP at a common working scale."""
import io, re, zipfile
import numpy as np
from PIL import Image

def stack(zpath, target_um=200.0, src_um=None, max_slices=260):
    z = zipfile.ZipFile(zpath)
    names = sorted(n for n in z.namelist() if n.lower().endswith('.jp2'))
    if not names:
        raise RuntimeError('empty archive: %s' % zpath)
    if src_um is None:
        m = re.search(r'(\d+(?:\.\d+)?)um', names[0])
        src_um = float(m.group(1))
    step_z = max(1, int(round(target_um / src_um)))
    sel = names[::step_z][:max_slices]
    # in-plane binning to the same physical voxel
    probe = np.array(Image.open(io.BytesIO(z.read(sel[0]))))
    binxy = max(1, int(round(target_um / src_um)))
    h, w = probe.shape
    H, W = h // binxy, w // binxy
    vol = np.empty((len(sel), H, W), np.float32)
    for i, n in enumerate(sel):
        a = np.array(Image.open(io.BytesIO(z.read(n))), np.float32)
        a = a[:H * binxy, :W * binxy].reshape(H, binxy, W, binxy).mean((1, 3))
        vol[i] = a
    eff = src_um * binxy
    return vol, eff, eff * step_z, len(names)
