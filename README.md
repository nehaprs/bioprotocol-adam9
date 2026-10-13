# O2PLS integration of RNA-seq and cell-surface proteomics to identify ADAM9 targets

This R script integrates a transcriptome (RNA-seq) dataset and a cell-surface proteome dataset from control and ADAM9-knockdown colon cancer cells using two-way orthogonal partial least squares (O2PLS; R package `OmicsPLS`). It then uses the fitted model to sort proteins into candidate **transcriptional** and **post-transcriptional** ADAM9 targets.

## Requirements

- R (≥ 4.0)
- R packages: `dplyr`, `OmicsPLS`, `readxl`, `writexl`. Any missing package is installed automatically when the script starts.

## Input files

Put both files in the R working directory:

| File | Contents |
|---|---|
| `rna.xlsx` | RNA-seq expression: one row per gene |
| `cp.xlsx` | Cell-surface proteomics: one row per protein |

Both files must have the same layout:

| Column 1 | Columns 2–6 | Columns 7–11 |
|---|---|---|
| Gene/protein name | `control1` … `control5` | `ADAM9KD1` … `ADAM9KD5` |

- **The sample columns must be in the same order in both files.** O2PLS pairs the two datasets sample by sample.
- RNA genes and proteins are matched on the name in column 1, so the proteomics names need to be gene symbols (e.g. `ADAM9`).
- **Optional:** add a column named `FDR` to `cp.xlsx` with differential-expression FDRs. When it is there, the `ttDEP` and `pttDEP` outputs keep only proteins with FDR < 0.05. Without it, those outputs are not filtered and the script prints a warning.

## How to run

In RStudio, set the working directory to the folder holding the input files, open the script, and click **Source**. From a terminal:

```
Rscript o2pls_adam9_pipeline.R
```

To rerun cross-validation for choosing the number of components, set `run_cv <- TRUE` at the top of the script. This step is slow.

## Method

1. **Preprocessing.** Each dataset is transposed to samples × features. Every feature is then centered and scaled to unit variance.
2. **O2PLS model.** `o2m()` fits a model with **2 joint components** and **1 orthogonal component per dataset**. These numbers were chosen by 5-fold cross-validation (`crossval_o2m_adjR2`, `crossval_o2m`).
3. **Gene–protein joint correlation.** Genes and proteins are matched by name. For each matched pair, the script rebuilds the joint (shared) part of the gene and of the protein across samples from the joint scores and loadings. It then computes their correlation (`Corr_xyJnt`).
4. **Transcriptional targets.** Each correlation gets a one-sided (positive) p-value from Fisher's z-transformation with n − 3 degrees of freedom, where n = 10 samples. P-values are then Benjamini–Hochberg (BH) adjusted. Pairs with BH < 0.05 are potential transcriptional targets.
5. **Post-transcriptional targets.** Pairs with a negative joint correlation, or with BH > 0.1, are potential post-transcriptional targets.
6. **Correlation with ADAM9.** The model-predicted proteome is computed as the joint part plus the orthogonal part. Each protein's predicted profile is then correlated with ADAM9's predicted profile (Pearson, two-sided). Within each target class, proteins are split by the sign of this correlation and BH-adjusted separately.

## Output files

All outputs go to the working directory, except `ptt.xlsx`, which goes to `analysis after integration/`. That folder is created if it does not exist.

**Model outputs**

| File | Contents |
|---|---|
| `lv_rnajnt.xlsx` | RNA joint loadings (components 1–2) |
| `lv_cpjnt.xlsx` | Protein joint loadings (components 1–2) |
| `lv_rnaorth.xlsx` | RNA orthogonal loadings |
| `lv_cpOrth.xlsx` | Protein orthogonal loadings |
| `scores_rnajnt.xlsx` | RNA joint scores per sample (t1, t2) |
| `scores_cpjnt.xlsx` | Protein joint scores per sample (u1, u2) |
| `scores_rnaorth.xlsx` | RNA orthogonal scores per sample |
| `scores_cporth.xlsx` | Protein orthogonal scores per sample |

**Target identification**

| File | Contents |
|---|---|
| `alltargets_bh1xy.xlsx` | All matched gene–protein pairs, with joint correlation and BH-adjusted p-value |
| `tt.xlsx` | Potential transcriptional targets (BH < 0.05) |
| `ttDEP.xlsx` | Transcriptional targets with expression values (FDR < 0.05 if an `FDR` column is supplied) |
| `analysis after integration/ptt.xlsx` | Potential post-transcriptional targets |
| `pttDEP.xlsx` | Post-transcriptional targets with expression values (FDR < 0.05 if an `FDR` column is supplied) |
| `ttpossigBH.xlsx`, `ttnegsigBH.xlsx` | Transcriptional targets positively / negatively correlated with ADAM9 (BH < 0.05) |
| `pttpossigBH.xlsx`, `pttnegsigBH.xlsx` | Post-transcriptional targets positively / negatively correlated with ADAM9 (BH < 0.05) |

## Notes

- If a gene name appears more than once in either file, the gene–protein join may warn about "many-to-many" matches. In that case, every combination of the duplicated entries is kept.
- Thresholds (0.05, 0.1), the gene of interest (`"ADAM9"`), and the numbers of components (`o2m(rnaTS, cpTS, 2, 1, 1)`) are set directly in the script.

## Citation

If you use this code, please cite [Citation for this paper] and the OmicsPLS package:
Bouhaddani, S. el, et al. (2018). Integrating omics datasets with the OmicsPLS package. *BMC Bioinformatics*, 19, 371.
