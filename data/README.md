# Human Organ Atlas datasets used

Open access, CC-BY-4.0, from https://human-organ-atlas.esrf.eu
Publications using these data must cite the DOI of each dataset.

Download the complete-organ archives at the coarsest available resolution and
place them, unrenamed, in `data/hoa/`. The scripts read the voxel size from
the file names, so do not rename them.

## Kidney

| Donor | Side | DOI |
| --- | --- | --- |
| LADAF-2020-27 | right | 10.15151/ESRF-DC-2191663251 |
| LADAF-2021-17 | right | 10.15151/ESRF-DC-1773966439 |
| LADAF-2021-17 | left | 10.15151/ESRF-DC-1773965419 |
| LADAF-2022-27 | bottom | 10.15151/ESRF-DC-2222696470 |
| S-20-28 | - | 10.15151/ESRF-DC-2008212555 |
| S-20-29 | - | 10.15151/ESRF-DC-1773965411 |

## Heart and lung

The heart and lung cohorts are listed in `data/datasets_heart.csv` and
`data/datasets_lung.csv` with the same fields. Donor metadata, including the
clinical context that qualifies each cohort, come from the JSON file that
accompanies every dataset on the portal.

## TRUSTED (kidney, clinical CT, 48 donors)

96 kidneys from 48 patients, CT at 0.81 x 0.81 x 1.05 mm, segmented
independently by two radiologists with the reference produced by STAPLE.

    https://doi.org/10.6084/m9.figshare.27981050

Place the extracted tree at `data/trusted/`, so that the masks are at
`data/trusted/CT_masks/GT_estimated_masksCT/`. `Annotator1/` and
`Annotator2/` beside it give the annotation error, which `08_partition.py`
reports separately from biological variability.

## Cohort composition

These are cadaveric donors, not a healthy population. The kidney cohort is
elderly with hypertension; the heart cohort is congenital disease (tetralogy
of Fallot); the lung cohort is COVID-19. The floor measured here is therefore
an upper bound on developmental variability, and the manuscript says so.
