# 22123 Computational Precision Medicine — group project (Fall 2026)

Gene-expression biomarkers of immunotherapy response in IMvigor210
(urothelial cancer, anti-PD-L1; Mariathasan et al., 2018, Nature).

## Research questions and scripts

| RQ | Topic | Script |
|----|-------|--------|
| RQ1 (part A) | Data loading, QC, filtering, VST + PCA, DESeq2 (R vs NR), volcano plot | `scripts/RQ1.1.qmd` |
| RQ1 (part B) | Heatmaps of DEGs, GSEA (Hallmark, GO:BP) and leading-edge genes | `scripts/RQ1.2.qmd` |
| RQ2 | T cell–inflamed GEP and response | `scripts/RQ2.qmd` |
| RQ3 | Biomarker comparison: TMB, PD-1 (PDCD1) and GEP | `scripts/RQ3.qmd` |
| RQ4 | Tumor purity (ESTIMATE) and immune-related biomarkers | `scripts/RQ4.qmd` |

## Structure
```
data/              IMvigor210_plain.rds  (course data: counts, tpm, clinical)
scripts/           RQ1.1.qmd, RQ1.2.qmd, RQ2.qmd, RQ3.qmd, RQ4.qmd
results/           tables (.csv), shared processed objects (.rds), sessionInfo
results/figures/   figures (.png)
slides/      slides per RQ (merged into the group presentation)
```

## Run
Render RQ1.1.qmd first, the rest of the Quarto files can be in ran in any order:

1. `scripts/RQ1.1.qmd` — writes the shared objects used by RQ1 part B and RQ4
2. `scripts/RQ1.2.qmd`
3. `scripts/RQ2.qmd`, `scripts/RQ3.qmd` — independent (read the data file directly)
4. `scripts/RQ4.qmd`

R >= 4.3. Packages:
```r
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "apeglm", "fgsea", "limma"))
install.packages(c("ggplot2", "ggrepel", "msigdbr", "pheatmap", "data.table",
                   "patchwork", "tidyverse", "remotes"))
remotes::install_github("omnideconv/immunedeconv")   # RQ4 (ESTIMATE)
```

## Shared objects (written by RQ1.1)
- `results/processed_shared.rds` — list with `vst` (VST matrix, 17,470 filtered genes x 192
  samples) and `clinical` (pat_id, BOR factor NR/R, TMB, sample_id, lib_size)
- `results/dds.rds` — fitted DESeq2 object (design `~ BOR`, NR = reference)
- TPM for GEP / PDCD1 / ESTIMATE: `readRDS("data/IMvigor210_plain.rds")$tpm`

```r
p   <- readRDS("results/processed_shared.rds")
dds <- readRDS("results/dds.rds")
```

## Key numbers (RQ1 part A)
- 192 patients: 25 R (complete response) vs 167 NR (progressive disease); TMB missing for 42
- 17,470 / 31,087 genes kept (>= 10 counts in >= 25 samples)
- PCA: no separation by response (PC1 p = 0.45, PC2 p = 0.14), no depth effect, no outliers
- 277 DEGs (FDR < 0.05, |log2FC| >= 1): 73 up, 204 down in R
- 50 robust DEGs (survive apeglm shrinkage), 49 up in R — IFN-γ / cytotoxic T & NK genes
  (CXCL9, CXCL13, IFNG, CCL5, GBP5, NCR1, KLRC3, SLAMF7 ...)
- Most "down in R" genes are liver/gut/skin genes (ALB, FGB, MUC2, KRT1) driven by a few
  samples — likely biopsy site; shrunk to ~0
- Tables: `results/DE_R_vs_NR_all.csv` (all genes; column `robust`),
  `results/DE_R_vs_NR_significant.csv`
  
## Key numbers (RQ1 part B)

## Key numbers (RQ2)
- GEP score higher in R than NR (Wilcoxon p = 0.00018)
- Logistic regression: OR 3.1 per unit of GEP (95% CI 1.7–5.8)
- Response rate above vs below median GEP: 21% vs 5% (Fisher p = 0.002)
- 15 of 18 GEP genes higher in R (FDR < 0.05); not CD276, HLA-DQA1, HLA-DRB1
- High GEP + low TMB: 2 of 35 respond; high GEP + high TMB: 39% respond

## Key numbers (RQ3)
- All three biomarkers are higher in R than NR (Wilcoxon, BH-adjusted p < 0.001):
  - TMB: median 14 vs 7 mut/Mb (p = 1.7e-5; n = 150, 42 missing)
  - PD-1 (PDCD1): median 1.29 vs 0.51 TPM (p = 0.0004)
  - GEP: median 0.60 vs -0.11 (p = 0.0002)

## Key numbers (RQ4)
- Tumor purity (ESTIMATE): 0.30–0.96, median 0.65
- Purity vs biomarkers (Spearman): PD-1 rho = -0.69, GEP rho = -0.77 (both p < 0.001);
  TMB rho = 0.12 (p = 0.15)
- Purity R vs NR: median 0.63 vs 0.66 (Wilcoxon p = 0.17)
- After regressing purity out of the expression data (`limma::removeBatchEffect`),
  PD-1 and GEP remain higher in R (Wilcoxon p = 0.002 and p < 0.001)
- Tables: `results/RQ4_purity_correlations.csv`, `results/RQ4_before_after_adjustment.csv`

## Notes
- GEP gene list: the course (Exercise 6) version, used in RQ2–RQ4. It differs from
  Ayers et al. (2017) by containing PDCD1 instead of CCL5.