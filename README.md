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
Open `22123_CPM_gr9.Rproj` (or work from the repository root) and render the
Quarto files in this order:

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
remotes::install_github("omnideconv/immunedeconv")   
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

## Notes
- GEP gene list: the course (Exercise 6) version, used in RQ2–RQ4. It differs from
  Ayers et al. (2017) by containing PDCD1 instead of CCL5.