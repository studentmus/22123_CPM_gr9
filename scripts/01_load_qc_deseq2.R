# =============================================================================
# 22123 CPM group project — RQ1 (part A): data loading, QC, filtering,
# VST + PCA, differential expression (DESeq2) R vs NR, volcano plot.
#
# Input : data/IMvigor210_plain.rds  (list: counts, tpm, clinical)
# Output: results/  (figures, DE tables, shared processed objects)
#
# Run from the repository root:  Rscript scripts/01_load_qc_deseq2.R
# Packages: DESeq2, apeglm, ggplot2, ggrepel
#   BiocManager::install(c("DESeq2", "apeglm")); install.packages("ggrepel")
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(ggrepel)
})

set.seed(22123)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 12))
pal <- c(NR = "#7F7F7F", R = "#D55E00")

# -----------------------------------------------------------------------------
# 1. Load data and sanity checks
# -----------------------------------------------------------------------------
d        <- readRDS("data/IMvigor210_plain.rds")
counts   <- as.matrix(d$counts)
tpm      <- as.matrix(d$tpm)
clinical <- as.data.frame(d$clinical)

stopifnot(identical(colnames(counts), clinical$sample_id),
          identical(colnames(tpm),    clinical$sample_id),
          identical(rownames(counts), rownames(tpm)))

# Response as factor with NR as reference -> log2FC > 0 means higher in R
clinical$BOR <- factor(clinical$BOR, levels = c("NR", "R"))
rownames(clinical) <- clinical$sample_id

cat("Genes x samples:", dim(counts), "\n")
print(table(clinical$BOR))
cat("TMB missing:", sum(is.na(clinical$TMB)), "\n")

# -----------------------------------------------------------------------------
# 2. QC: library sizes and gene filtering
# -----------------------------------------------------------------------------
clinical$lib_size <- colSums(counts)
cat("Library size (M reads): ",
    paste(round(summary(clinical$lib_size) / 1e6, 1), collapse = " / "), "\n")

p_lib <- ggplot(clinical, aes(BOR, lib_size / 1e6, fill = BOR)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(x = NULL, y = "Library size (million reads)",
       title = "Sequencing depth per sample")
ggsave("results/figures/qc_library_size.png", p_lib, width = 4, height = 4, dpi = 300)

# Keep genes with >= 10 counts in at least as many samples as the smallest
# group (25 responders): a gene expressed only in responders can still pass.
min_group <- min(table(clinical$BOR))
keep      <- rowSums(counts >= 10) >= min_group
cat("Genes all-zero:", sum(rowSums(counts) == 0),
    "| genes kept after filtering:", sum(keep), "of", nrow(counts), "\n")
counts_f <- counts[keep, ]

# -----------------------------------------------------------------------------
# 3. DESeq2 object, VST and PCA
# -----------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts_f,
                              colData   = clinical,
                              design    = ~ BOR)

# Variance-stabilising transform for exploration (blind to design)
vsd <- vst(dds, blind = TRUE)

# PCA on the 500 most variable genes
rv   <- apply(assay(vsd), 1, var)
top  <- order(rv, decreasing = TRUE)[1:500]
pca  <- prcomp(t(assay(vsd)[top, ]), center = TRUE, scale. = FALSE)
pve  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pc_df <- cbind(clinical, PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3])

p_pca <- ggplot(pc_df, aes(PC1, PC2, colour = BOR)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_manual(values = pal) +
  labs(x = paste0("PC1 (", pve[1], "%)"), y = paste0("PC2 (", pve[2], "%)"),
       colour = "Response", title = "PCA of VST expression (top 500 variable genes)")
ggsave("results/figures/pca_response.png", p_pca, width = 6, height = 4.5, dpi = 300)

p_pca_lib <- ggplot(pc_df, aes(PC1, PC2, colour = lib_size / 1e6)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_viridis_c() +
  labs(x = paste0("PC1 (", pve[1], "%)"), y = paste0("PC2 (", pve[2], "%)"),
       colour = "Library\nsize (M)", title = "PCA coloured by sequencing depth")
ggsave("results/figures/pca_libsize.png", p_pca_lib, width = 6, height = 4.5, dpi = 300)

# Is response or depth associated with the main axes of variation?
cat("\nPC1 vs response (Wilcoxon p):",
    signif(wilcox.test(PC1 ~ BOR, data = pc_df)$p.value, 3), "\n")
cat("PC2 vs response (Wilcoxon p):",
    signif(wilcox.test(PC2 ~ BOR, data = pc_df)$p.value, 3), "\n")
cat("PC1 vs library size (Spearman rho):",
    round(cor(pc_df$PC1, pc_df$lib_size, method = "spearman"), 2), "\n")

# Outlier check: samples > 4 SD from the centre on PC1 or PC2
z <- scale(pca$x[, 1:2])
outl <- rownames(z)[apply(abs(z) > 4, 1, any)]
cat("Potential PCA outliers (>4 SD):",
    if (length(outl)) paste(outl, collapse = ", ") else "none", "\n")

# -----------------------------------------------------------------------------
# 4. Differential expression: R vs NR
# -----------------------------------------------------------------------------
dds <- DESeq(dds, quiet = TRUE)
print(resultsNames(dds))

res <- results(dds, contrast = c("BOR", "R", "NR"), alpha = 0.05)
summary(res)

# Shrunken log2 fold changes (apeglm) for ranking and plotting
res_shr <- lfcShrink(dds, coef = "BOR_R_vs_NR", type = "apeglm", quiet = TRUE)

de <- data.frame(gene           = rownames(res),
                 baseMean       = res$baseMean,
                 log2FC         = res$log2FoldChange,
                 log2FC_shrunk  = res_shr$log2FoldChange,
                 lfcSE          = res$lfcSE,
                 stat           = res$stat,
                 pvalue         = res$pvalue,
                 padj           = res$padj)
de <- de[order(de$pvalue), ]

padj_cut <- 0.05; lfc_cut <- 1
is_sig <- !is.na(de$padj) & de$padj < padj_cut
# DEGs: adjusted p < 0.05 and |MLE log2FC| >= 1 (standard definition)
de$direction <- "not significant"
de$direction[is_sig & de$log2FC >=  lfc_cut] <- "up in R"
de$direction[is_sig & de$log2FC <= -lfc_cut] <- "down in R"
# Robust DEGs: effect survives apeglm shrinkage with the same sign.
# Shrinkage pulls fold changes driven by a few extreme samples / high
# dispersion towards 0, so these are the most reliable hits.
de$robust <- de$direction != "not significant" &
             abs(de$log2FC_shrunk) >= lfc_cut &
             sign(de$log2FC_shrunk) == sign(de$log2FC)

sig <- subset(de, direction != "not significant")
cat("\nDEGs (padj <", padj_cut, ", |log2FC| >=", lfc_cut, "):", nrow(sig),
    "| up in R:", sum(sig$direction == "up in R"),
    "| down in R:", sum(sig$direction == "down in R"), "\n")
cat("Robust DEGs (also |shrunken log2FC| >= 1):", sum(de$robust),
    "| up:", sum(de$robust & de$log2FC > 0),
    "| down:", sum(de$robust & de$log2FC < 0), "\n")
cat("Genes with padj < 0.05 (any fold change):", sum(is_sig), "\n")

write.csv(de,  "results/DE_R_vs_NR_all.csv",         row.names = FALSE)
write.csv(sig, "results/DE_R_vs_NR_significant.csv", row.names = FALSE)

cat("\nTop 20 genes by p-value:\n")
print(head(de[, c("gene", "baseMean", "log2FC", "log2FC_shrunk", "padj")], 20), row.names = FALSE)
cat("\nRobust DEGs:\n")
print(de[de$robust, c("gene", "baseMean", "log2FC", "log2FC_shrunk", "padj")], row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Volcano plot
# -----------------------------------------------------------------------------
# Label: top DEGs by adjusted p-value + robust DEGs with the strongest effect
lab <- unique(c(head(sig$gene[order(sig$padj)], 12),
                head(de$gene[de$robust][order(-abs(de$log2FC_shrunk[de$robust]))], 12)))
de$label <- ifelse(de$gene %in% lab, de$gene, NA)

p_vol <- ggplot(subset(de, !is.na(padj)),
                aes(log2FC, -log10(padj), colour = direction)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", linewidth = 0.3) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 40,
                  colour = "black", na.rm = TRUE, min.segment.length = 0) +
  scale_colour_manual(values = c("down in R" = "#0072B2", "up in R" = "#D55E00",
                                 "not significant" = "grey80")) +
  labs(x = "log2 fold change (R vs NR)", y = "-log10 adjusted p-value",
       colour = NULL,
       title = "Complete response (n=25) vs progressive disease (n=167)") +
  theme(legend.position = "top")
ggsave("results/figures/volcano_R_vs_NR.png", p_vol, width = 7, height = 6, dpi = 300)

# MLE vs shrunken fold change: which DEGs are driven by a few samples?
p_shr <- ggplot(subset(de, direction != "not significant"),
                aes(log2FC, log2FC_shrunk, colour = robust)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_point(size = 1.6, alpha = 0.8) +
  geom_text_repel(aes(label = ifelse(abs(log2FC) > 3.5 | (robust & abs(log2FC_shrunk) > 1.4),
                                     gene, NA)),
                  size = 2.8, colour = "black", na.rm = TRUE, max.overlaps = 30) +
  scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey55"),
                      labels = c(`TRUE` = "robust", `FALSE` = "shrunk to ~0")) +
  labs(x = "MLE log2 fold change", y = "apeglm-shrunken log2 fold change", colour = NULL,
       title = "Effect of fold-change shrinkage on DEGs") +
  theme(legend.position = "top")
ggsave("results/figures/lfc_shrinkage.png", p_shr, width = 6, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 6. Shared objects for the rest of the group (RQ1 part B, RQ2-4)
# -----------------------------------------------------------------------------
# processed_shared.rds : VST matrix of filtered genes + clinical table
#                        (BOR as factor, lib_size) -> heatmaps, signature scores
# dds.rds              : fitted DESeq2 object -> enrichment / more contrasts
# TPM stays in data/IMvigor210_plain.rds (GEP score, PDCD1).
# Both .rds files are regenerated by this script and are not tracked by git.
saveRDS(list(vst = assay(vsd), clinical = clinical), "results/processed_shared.rds")
saveRDS(dds, "results/dds.rds")

writeLines(capture.output(sessionInfo()), "results/sessionInfo_01.txt")
cat("\nDone. Outputs written to results/\n")
