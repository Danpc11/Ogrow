#!/usr/bin/env python
"""What a batch correction does to a progenitor share, and what it cannot do.

Two studies of human fetal kidney at nearby stages disagree fourfold on the
share of the nephron lineage that is progenitor. The obvious move is to
correct for batch and compare the corrected shares. This script runs that
move, on the same cells and the same markers, under three scalings.

    per sample   each gene z-scored within its own sample, which is what
                 scripts/26_transcriptome.py does. The threshold then moves
                 with the composition of the sample it is applied to.
    pooled       each gene z-scored over every kidney cell at once, so one
                 threshold serves all samples.
    centred      pooled, and each gene additionally centred within its sample,
                 which is the location adjustment a batch correction makes.

The result is that centring does not bring the two studies together. An
expression-level correction aligns shared states; it does not change which
cells a preparation contains, and composition is what these two studies
differ in. That is the reason only presence and absence are read across
studies in the manuscript, and it is measured here rather than asserted.

This is a diagnostic, not part of the manuscript pipeline, so `make paper`
does not run it. It needs the raw matrices; see data/transcriptome/README.md.

Usage:  python scripts/28_batch.py
"""
import csv
import gzip

import numpy as np

from organgrow.config import ROOT

DATA = ROOT / "data" / "transcriptome"
OUT = ROOT / "source_data" / "Fig9_batch.csv"
# week, folder, file stem, whether the matrices are gzipped
SETS = [(9, "GSE114530", "GSM3509837_1_wk09", True),
        (11, "GSE114530", "GSM3509838_2_wk11", True),
        (13, "GSE114530", "GSM3509839_3_wk13", True),
        (14, "GSE139280", "GSM4135994_count_table_H14w_0h", False),
        (16, "GSE114530", "GSM3143601_4_w16", True),
        (18, "GSE114530", "GSM3509840_5_wk18", True)]


def load(tx):
    """Marker-gene expression of every kidney cell, with its sample."""
    from scipy.io import mmread
    from scipy import sparse

    genes = sorted(set(tx.NPC + tx.INDUCED + tx.NEPHRON + tx.OFF))
    blocks, batch = [], []
    for wk, folder, stem, gz in SETS:
        d = DATA / folder
        g = d / (stem + ("_genes.tsv.gz" if gz else "_genes.tsv"))
        with (gzip.open(g, "rt") if gz else open(g)) as f:
            sym = np.array([ln.split("\t")[1].strip() for ln in f])
        M = mmread(str(d / (stem + ("_matrix.mtx.gz" if gz else "_matrix.mtx")))).tocsc()
        tot = np.asarray(M.sum(0)).ravel()
        keep = (tot >= 1000) & (np.asarray((M > 0).sum(0)).ravel() >= 500)
        X = M[:, keep].multiply(sparse.csr_matrix(1e4 / tot[keep])).tocsr()
        rows = {x: np.flatnonzero(sym == x) for x in genes}
        blocks.append(np.vstack(
            [np.log1p(np.asarray(X[rows[x]].todense())).mean(0) if len(rows[x])
             else np.zeros(X.shape[1]) for x in genes]))
        batch += [wk] * blocks[-1].shape[1]
        print(f"  wk{wk:>3} {folder}: {blocks[-1].shape[1]} cells", flush=True)
    return np.hstack(blocks), np.array(batch), {g: i for i, g in enumerate(genes)}


def z(V, within=None):
    """Z-score each gene, over everything or within each group."""
    if within is None:
        return (V - V.mean(1, keepdims=True)) / (V.std(1, keepdims=True) + 1e-9)
    Z = np.zeros_like(V)
    for b in sorted(set(within)):
        m = within == b
        Z[:, m] = ((V[:, m] - V[:, m].mean(1, keepdims=True))
                   / (V[:, m].std(1, keepdims=True) + 1e-9))
    return Z


def shares(Z, batch, idx, tx):
    """Progenitor share of the nephron lineage, per sample."""
    out = {}
    for b in sorted(set(batch)):
        m = batch == b

        def sc(gs):
            return Z[[idx[g] for g in gs]][:, m].mean(0)

        s = np.vstack([sc(tx.NPC), sc(tx.INDUCED), sc(tx.NEPHRON)])
        lin = (sc(tx.OFF) < 0.25) & (s.max(0) > 0.15)
        out[b] = float((s[:, lin].argmax(0) == 0).mean())
    return out


def main():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "tx", str(ROOT / "scripts" / "26_transcriptome.py"))
    tx = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tx)

    missing = sorted({f for _, f, stem, gz in SETS
                      if not (DATA / f / (stem + ("_matrix.mtx.gz" if gz
                                                  else "_matrix.mtx"))).exists()})
    if missing:
        raise SystemExit(f"needs the raw matrices under {DATA}: {missing}")

    V, batch, idx = load(tx)
    per = shares(z(V, within=batch), batch, idx, tx)
    pooled = shares(z(V), batch, idx, tx)
    Vc = V.copy()
    for b in sorted(set(batch)):
        m = batch == b
        Vc[:, m] -= Vc[:, m].mean(1, keepdims=True)
    centred = shares(z(Vc), batch, idx, tx)

    rows = []
    print(f"\n{'week':>5} {'study':>10} {'per sample':>11} {'pooled':>8} {'centred':>9}")
    for wk, folder, _, _ in sorted(SETS):
        rows.append([wk, folder, per[wk], pooled[wk], centred[wk]])
        print(f"{wk:>5} {folder:>10} {per[wk]:11.3f} {pooled[wk]:8.3f} {centred[wk]:9.3f}")

    def gap(d):
        """The second study against the mean of the stages that bracket it."""
        return d[14] / max(np.mean([d[13], d[16]]), 1e-9)

    print(f"\nthe week-14 study against the mean of weeks 13 and 16 that bracket "
          f"it:\n  per sample {gap(per):.1f}x, pooled {gap(pooled):.1f}x, "
          f"centred {gap(centred):.1f}x")
    print("centring each gene within its sample is the location adjustment a batch"
          "\ncorrection makes, and it leaves the two studies apart: what separates"
          "\nthem is which cells the preparation holds, which no expression-level"
          "\ncorrection changes")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["week", "study", "per_sample", "pooled", "centred_within_sample"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
