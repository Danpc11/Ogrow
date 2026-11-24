#!/usr/bin/env python
"""Three developing organs, read against the law that is supposed to govern them.

A pool that is a budget empties as its organ finishes; a pool that renews is
not what limits its organ at all, so it can go early and leave the organ
building. That is what these three series are read for, and it is a question
the last measured stage answers.

The nucleation law also says that while the supply term stays open, induction
is proportional to the pool, so the number induced per progenitor still held
should hold steady as the pool falls. That one is reported and not concluded
from: the series carries a single sample per stage, and across the five the
ratio does not order with stage, so a change over development and a
difference between samples cannot be told apart here.

    kidney    GSE114530, human fetal, weeks 9 to 18. Cells of the nephron
              lineage are put in one of three states by marker score: cap
              mesenchyme still uncommitted, a pretubular aggregate or renal
              vesicle being induced, and a unit already formed.
    pancreas  GSE132188, mouse, E12.5 to E15.5, the authors' own annotation.
              Strict progenitors are the multipotent pool, broad adds tip and
              trunk, which are upstream of a fate but no longer multipotent.
    liver     GSE90047, mouse, E10.5 to E17.5, cells sorted on DLK1, so the
              readout is the share of the lineage still running the
              hepatoblast programme. A two-component mixture on the maturation
              axis sets the boundary rather than a chosen cutoff; where the
              fit places it moves the middle stages by a few points and
              leaves the two ends, which is what the reading rests on.

Raw matrices are not kept here. The derived tables are committed under
source_data/, so the figure and every number stand without the download; with
the files in place under data/transcriptome/ the script rebuilds them. See
data/transcriptome/README.md.

Usage:  python scripts/26_transcriptome.py
"""
import csv
import gzip
import json

import numpy as np

from organgrow.config import ROOT

DATA = ROOT / "data" / "transcriptome"
OUT = ROOT / "results"
SD = ROOT / "source_data"
RNG = np.random.default_rng(0)

KIDNEY_SAMPLES = {9: "GSM3509837_1_wk09", 11: "GSM3509838_2_wk11",
                  13: "GSM3509839_3_wk13", 16: "GSM3143601_4_w16",
                  18: "GSM3509840_5_wk18"}
NPC = ["SIX2", "CITED1", "MEOX1", "EYA1", "SALL1", "CRYM", "DAPL1", "TMEM100"]
INDUCED = [g for g in ["WNT4", "LHX1", "PAX8", "JAG1", "DLL1", "HEYL", "CCND1"]
           if g not in NPC]
NEPHRON = ["NPHS1", "NPHS2", "PODXL", "LRP2", "SLC34A1", "UMOD", "GATA3", "WT1"]
OFF = ["PDGFRB", "COL3A1", "FOXD1", "PECAM1", "PTPRC", "HBB", "HBA1"]
# Dlk1 is the sorting marker, so it is not in the hepatoblast signature
HEPATOBLAST = ["Afp", "H19", "Igf2", "Lin28b", "Mest", "Peg3", "Cd24a",
               "Dlk1os", "Sox9", "Lgals1", "Mdk", "Npm1", "Hmga2", "Cnn3", "Tnc"]
HEPATOCYTE = ["Alb", "Ttr", "Apoa1", "Apoa2", "Apob", "Serpina1a", "Serpina1b",
              "Serpina1c", "Tat", "G6pc", "Cps1", "Fabp1", "Trf", "Ahsg",
              "Hamp", "Cyp2e1", "Cyp3a11", "Itih4", "Apoc3", "Hpx"]

# Where each organ's unit-forming window ends, from the literature, in the same
# units as its stages. Mouse gestation is taken as 19.5 days, so a postnatal
# day P is E(19.5 + P). The liver's window is the one that matters: hepatocyte
# number keeps rising to about P20, long after the hepatoblasts are gone.
WINDOW = {"kidney": dict(unit="gestational week", start=5.0, end=35.0),
          "pancreas": dict(unit="embryonic day", start=9.5, end=15.5),
          "liver": dict(unit="embryonic day", start=10.5, end=39.5)}


def boot_ci(v, n=2000):
    """Percentile interval of a proportion, by resampling the cells."""
    v = np.asarray(v, float)
    if len(v) < 2:
        return float("nan"), float("nan")
    m = RNG.choice(v, size=(n, len(v)), replace=True).mean(1)
    return tuple(float(x) for x in np.percentile(m, [2.5, 97.5]))


def write_table(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def read_table(path):
    with open(path) as f:
        r = list(csv.DictReader(f))
    return {k: [float(d[k]) if _num(d[k]) else d[k] for d in r] for k in r[0]}


def _num(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def zscore_mean(v):
    """Mean of each row z-scored across cells; v is genes by cells."""
    z = (v - v.mean(1, keepdims=True)) / (v.std(1, keepdims=True) + 1e-9)
    return np.asarray(z.mean(0)).ravel()


# ----------------------------------------------------------------- kidney
def kidney_states(d):
    from scipy.io import mmread
    from scipy import sparse

    rows = []
    for wk, pre in KIDNEY_SAMPLES.items():
        with gzip.open(d / f"{pre}_genes.tsv.gz", "rt") as f:
            sym = np.array([ln.split("\t")[1].strip() for ln in f])
        M = mmread(str(d / f"{pre}_matrix.mtx.gz")).tocsc()
        tot = np.asarray(M.sum(0)).ravel()
        keep = (tot >= 1000) & (np.asarray((M > 0).sum(0)).ravel() >= 500)
        X = M[:, keep].multiply(sparse.csr_matrix(1e4 / tot[keep])).tocsr()

        def sc(genes):
            r = [i for g in genes for i in np.flatnonzero(sym == g)]
            return (zscore_mean(np.log1p(np.asarray(X[r].todense())))
                    if r else np.zeros(X.shape[1]))

        s = np.vstack([sc(NPC), sc(INDUCED), sc(NEPHRON)])
        off = sc(OFF)
        lin = (off < 0.25) & (s.max(0) > 0.15)
        state = s[:, lin].argmax(0)
        rows.append([wk, int(lin.sum()), int((state == 0).sum()),
                     int((state == 1).sum()), int((state == 2).sum())])
    return rows


# --------------------------------------------------------------- pancreas
def pancreas_states(d):
    import h5py

    f = h5py.File(d / "GSE132188_adata.h5", "r")
    obs = f["obs"][:]

    def cat(k):
        v = f["uns"][k][:]
        return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v])

    day = cat("day_categories")[obs["day"]]
    cl = cat("clusters_fig2_final_categories")[obs["clusters_fig2_final"]]
    strict = np.isin(cl, ["Multipotent"])
    broad = np.isin(cl, ["Multipotent", "Tip", "Trunk"])
    rows = []
    for s in sorted(set(day)):
        m = day == s
        lo1, hi1 = boot_ci(strict[m])
        lo2, hi2 = boot_ci(broad[m])
        rows.append([f"E{s}", int(m.sum()), float(strict[m].mean()), lo1, hi1,
                     float(broad[m].mean()), lo2, hi2])
    return rows


# ------------------------------------------------------------------ liver
def _mixture_2(x, iters=200):
    """Two-component one-dimensional Gaussian mixture, fitted by EM.

    A chosen cutoff on the maturation axis would be a free parameter; this
    lets the data place the boundary. Written out rather than imported to
    keep the dependency list as it is.
    """
    x = np.asarray(x, float)
    mu = np.array([np.quantile(x, 0.2), np.quantile(x, 0.8)])
    sd = np.full(2, x.std() + 1e-9)
    w = np.array([0.5, 0.5])
    for _ in range(iters):
        p = w * np.exp(-0.5 * ((x[:, None] - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))
        r = p / np.maximum(p.sum(1, keepdims=True), 1e-300)
        n = r.sum(0) + 1e-9
        w, mu = n / len(x), (r * x[:, None]).sum(0) / n
        sd = np.sqrt((r * (x[:, None] - mu) ** 2).sum(0) / n) + 1e-9
    hi = int(np.argmax(mu))
    return r[:, hi]


def liver_states(d):
    with gzip.open(d / "GSE90047_sc.txt.gz", "rt") as f:
        head = f.readline().rstrip("\n").split("\t")
        sym, cols, vals = [], head[3:], []
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            sym.append(p[1])
            vals.append([float(x) for x in p[3:]])
    sym = np.array(sym)
    X = np.array(vals, float)
    keep = [i for i, c in enumerate(cols)
            if c.startswith("E") and c.split("_")[0].endswith("D")]
    X, cols = X[:, keep], [cols[i] for i in keep]
    stage = np.array([c.split("_")[0][:-1] for c in cols])
    cpm = np.log1p(X / np.maximum(X.sum(0), 1) * 1e6)

    def sc(genes):
        r = [i for g in genes for i in np.flatnonzero(sym == g)]
        return zscore_mean(cpm[r]) if r else np.zeros(cpm.shape[1])

    maturation = sc(HEPATOCYTE) - sc(HEPATOBLAST)
    p_mature = _mixture_2(maturation)
    prog = p_mature < 0.5
    rows = []
    for s in sorted(set(stage), key=lambda z: float(z[1:])):
        m = stage == s
        lo, hi = boot_ci(prog[m])
        rows.append([s, int(m.sum()), float(prog[m].mean()), lo, hi,
                     float(np.median(maturation[m]))])
    return rows


# ------------------------------------------------------------------- main
def build(name, folder, fn, header, path):
    d = DATA / folder
    if d.exists():
        rows = fn(d)
        write_table(path, header, rows)
        return read_table(path)
    if path.exists():
        print(f"{name}: using the committed table; put the files under {d} to redo it")
        return read_table(path)
    raise SystemExit(f"{name}: nothing under {d} and no table at {path}")


def main():
    SD.mkdir(parents=True, exist_ok=True)
    k = build("kidney", "GSE114530", kidney_states,
              ["week", "lineage", "npc", "induced", "nephron"],
              SD / "Fig9_kidney_states.csv")
    p = build("pancreas", "GSE132188", pancreas_states,
              ["stage", "n", "frac_strict", "s_lo", "s_hi",
               "frac_broad", "b_lo", "b_hi"],
              SD / "Fig9_pancreas_states.csv")
    lv = build("liver", "GSE90047", liver_states,
               ["stage", "n", "frac_prog", "lo", "hi", "maturation"],
               SD / "Fig9_liver_states.csv")

    w = np.array(k["week"], float)
    npc, ind, lin = (np.array(k[x], float) for x in ("npc", "induced", "lineage"))
    share, per = npc / np.maximum(lin, 1), ind / np.maximum(npc, 1)
    a = float(np.polyfit(w, np.log(share), 1)[0])
    b = float(np.polyfit(w, np.log(per), 1)[0])
    m = w > w.min()
    a2 = float(np.polyfit(w[m], np.log(share[m]), 1)[0])
    b2 = float(np.polyfit(w[m], np.log(per[m]), 1)[0])
    # How far either ratio orders with stage at all. Spearman over the stages,
    # which are the units here: the cells inside one sample are not replicates
    # of each other, so a test that counts them would answer a question nobody
    # asked.
    from scipy.stats import spearmanr
    rs_share, p_share = spearmanr(w, share)
    rs_per, p_per = spearmanr(w, per)

    print(f"\n{'week':>5} {'lineage':>8} {'NPC':>6} {'induced':>8} {'nephron':>8} "
          f"{'NPC share':>10} {'induced per NPC':>16}")
    for i in range(len(w)):
        print(f"{w[i]:5.0f} {lin[i]:8.0f} {npc[i]:6.0f} {ind[i]:8.0f} "
              f"{k['nephron'][i]:8.0f} {share[i]:10.3f} {per[i]:16.2f}")
    print(f"\nslope of log(NPC share) per week        {a:+.3f}   "
          f"without week {w.min():.0f}: {a2:+.3f}")
    print(f"slope of log(induced per NPC) per week  {b:+.3f}   "
          f"without week {w.min():.0f}: {b2:+.3f}")
    print(f"ordering with stage, Spearman over the {len(w)} samples: "
          f"NPC share rho = {rs_share:+.2f} (P = {p_share:.2f}), "
          f"induced per NPC rho = {rs_per:+.2f} (P = {p_per:.2f})")
    print("one sample per stage, so neither ratio separates development from"
          "\nthe difference between samples")

    print(f"\n{'organ':9s} {'last stage':>11} {'progenitor share':>17} "
          f"{'of its window':>14}")
    done = {}
    for name, stages, frac, key in (
            ("kidney", w, share, "kidney"),
            ("pancreas", [float(s[1:]) for s in p["stage"]], p["frac_strict"], "pancreas"),
            ("liver", [float(s[1:]) for s in lv["stage"]], lv["frac_prog"], "liver")):
        win = WINDOW[key]
        f = (float(stages[-1]) - win["start"]) / (win["end"] - win["start"])
        done[key] = dict(last_stage=float(stages[-1]), share=float(frac[-1]),
                         window_fraction=f, **win)
        print(f"{name:9s} {float(stages[-1]):11.1f} {float(frac[-1]):17.4f} {f:13.0%}")
    print("the liver loses its named pool a quarter of the way through its window"
          "\nand goes on building; the budget organs lose theirs as they finish")

    OUT.mkdir(parents=True, exist_ok=True)
    json.dump(dict(kidney=k, pancreas=p, liver=lv, slope_log_npc_share=a,
                   slope_log_induced_per_npc=b, slope_log_npc_share_trimmed=a2,
                   slope_log_induced_per_npc_trimmed=b2,
                   spearman=dict(npc_share=[float(rs_share), float(p_share)],
                                 induced_per_npc=[float(rs_per), float(p_per)],
                                 n_samples=int(len(w))),
                   windows=done),
              open(OUT / "transcriptome.json", "w"), indent=1)


if __name__ == "__main__":
    main()
