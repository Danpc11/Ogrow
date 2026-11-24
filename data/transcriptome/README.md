# Transcriptome data

`scripts/26_transcriptome.py` reads three public series. The matrices run to
tens of megabytes and are not kept here. The tables the script derives from
them are committed under `source_data/`, so the figure and every number in the
manuscript stand without the download; with the files in place the script
rebuilds those tables instead of reading them.

Each series goes in its own folder, with the file names the Gene Expression
Omnibus gives them.

## GSE114530, human fetal kidney, weeks 9 to 18

    GSE114530/GSM3509837_1_wk09_genes.tsv.gz
    GSE114530/GSM3509837_1_wk09_matrix.mtx.gz
    ...and the same for GSM3509838 (wk11), GSM3509839 (wk13), GSM3509840 (wk18)

Each cell is scored against marker sets for cap mesenchyme, for a pretubular
aggregate or renal vesicle, and for a formed nephron, and takes the highest
scoring state. Anything scoring as stroma, endothelium or blood is dropped.

## GSE132188, mouse pancreas, E12.5 to E15.5

    GSE132188/GSE132188_adata.h5

The cell types are the authors' own annotation, read from the `obs` table.
Strict progenitors are the multipotent pool; broad adds tip and trunk, which
are past multipotency but upstream of a fate.

## GSE90047, mouse liver, E10.5 to E17.5

    GSE90047/GSE90047_sc.txt.gz

Cells are sorted on DLK1, so the readout is the share of the lineage still
running the hepatoblast programme rather than the share of the organ. DLK1
itself is left out of the signature, being the sorting marker. A
two-component mixture on the maturation axis places the boundary.
