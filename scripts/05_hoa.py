#!/usr/bin/env python
"""Human Organ Atlas: load, segment, gate, and tabulate invariants."""
import csv, glob, json, os, re, sys
from organgrow.config import CFG, ROOT
from organgrow.hoa_io import stack
from organgrow.hoa_invariants import invariants, segment
from organgrow.qc import qc

organ = sys.argv[1] if len(sys.argv) > 1 else "kidney"
OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
tgt = CFG["hoa"]["target_um"][organ]
mx = CFG["hoa"]["max_slices"][organ]

keep, rejected, seen = [], [], set()
for f in sorted(glob.glob(str(ROOT / "data" / "hoa" / f"*{organ}*complete*.zip"))):
    b = os.path.basename(f)
    m = re.match(r'_(\d+_?\d*)um_([A-Za-z0-9-]+(?:_\d+)?)_%s-?(\w*)_' % organ, b)
    if not m:
        continue
    key = (m.group(2), m.group(3) or "-")
    if key in seen:
        continue
    seen.add(key)
    try:
        v, xy, zz, n = stack(f, target_um=tgt, max_slices=mx)
        inv = invariants(v, xy, zz)
        msk, _ = segment(v)
        ok, why = qc(msk, inv["volume_mm3"], xy, zz)
        inv.update(donor=m.group(2), part=key[1], organ=organ,
                   src_um=float(m.group(1).replace("_", ".")), file=b)
        (keep.append(inv) if ok else rejected.append(dict(donor=key[0], part=key[1], reason=why)))
        print(f"  {'keep' if ok else 'drop'} {key[0]:>14} {key[1]:>8} "
              f"V={inv['volume_mm3']/1000:7.1f} cm3  {'' if ok else why}", flush=True)
    except Exception as exc:
        rejected.append(dict(donor=key[0], part=key[1], reason=str(exc)))
        print(f"  drop {key[0]:>14} {key[1]:>8}  {exc}", flush=True)

json.dump(keep, open(OUT / f"hoa_{organ}.json", "w"), indent=1)
with open(OUT / "hoa_rejected.csv", "a", newline="") as fh:
    w = csv.writer(fh)
    for r in rejected:
        w.writerow([organ, r.get("donor"), r.get("part"), r.get("reason")])
print(f"kept {len(keep)}, rejected {len(rejected)} -> results/hoa_{organ}.json")
