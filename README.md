# 22123 Computational Precision Medicine — group project (Fall 2026)

Gene-expression biomarkers of immunotherapy response in IMvigor210
(urothelial cancer, anti-PD-L1; Mariathasan et al., 2018, Nature).

## Structure
```
data/        IMvigor210_plain.rds   (not tracked by git — put the course file here)
scripts/     01_load_qc_deseq2.R    RQ1 part A: QC, filtering, VST + PCA, DESeq2, volcano
results/     DE tables, shared processed object, sessionInfo
results/figures/
slides/      slides per RQ (merged into the group presentation)
```

## Run
From the repository root (R >= 4.3):
```r
# install.packages("BiocManager"); BiocManager::install(c("DESeq2", "apeglm"))
# install.packages(c("ggplot2", "ggrepel"))
source("scripts/01_load_qc_deseq2.R")   # ~1-2 min
```

## Shared objects for RQ1 part B and RQ2-4
Running the script writes (not tracked by git, regenerate locally):
- `results/processed_shared.rds` — list with `vst` (VST matrix, 17,470 filtered genes x 192
  samples) and `clinical` (pat_id, BOR factor NR/R, TMB, sample_id, lib_size)
- `results/dds.rds` — fitted DESeq2 object (design `~ BOR`, NR = reference)
- TPM for GEP / PDCD1: `readRDS("data/IMvigor210_plain.rds")$tpm`

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
