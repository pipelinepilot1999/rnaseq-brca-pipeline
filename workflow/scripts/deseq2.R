# ============================================================
# DESeq2 Differential Expression Analysis
# Project: RNA-seq Pipeline — GSE60450 Phase 1
# Contrast: basal vs luminal (mouse mammary)
# ============================================================

# Block 1 — Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
counts_file   <- args[1]
metadata_file <- args[2]
output_file   <- args[3]

cat("Input counts:", counts_file, "\n")
cat("Metadata:", metadata_file, "\n")
cat("Output:", output_file, "\n")

# Block 2 — Load libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)

# Block 3 — Read input files
counts <- read.table(counts_file,
                     header = TRUE,
                     skip = 1,
                     row.names = 1)

# Keep only sample columns (remove annotation columns)
counts <- counts[, 6:ncol(counts)]

# Fix column names — strip full path and lwts keep sample ID only
colnames(counts) <- gsub("results.star.(.*?)_Aligned.*", "\\1", colnames(counts))
metadata <- read.csv(metadata_file,
                     row.names = 1)

# Block 4 — Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = metadata,
    design    = ~ condition
)

# Set reference level (basal is baseline)
dds$condition <- relevel(dds$condition, ref = "basal")

# Run DESeq2
dds <- DESeq(dds)

# Block 5 — Extract results
res <- results(dds, contrast = c("condition", "luminal", "basal"))

# Order by adjusted p-value
res <- res[order(res$padj), ]

# Save results to CSV
write.csv(as.data.frame(res), file = output_file)

cat("DESeq2 analysis complete\n")
cat("Results saved to:", output_file, "\n")

# Block 6 — Plots
output_dir <- dirname(output_file)

# PCA plot
vsd <- vst(dds, blind = FALSE)
pca_plot <- plotPCA(vsd, intgroup = "condition")
ggsave(file.path(output_dir, "pca_plot.png"), pca_plot)

# Volcano plot
res_df <- as.data.frame(res)
res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal() +
    labs(title = "Volcano Plot: Luminal vs Basal",
         x = "log2 Fold Change",
         y = "-log10 adjusted p-value")

ggsave(file.path(output_dir, "volcano_plot.png"), volcano)

# Heatmap — top 50 DEGs
top50 <- head(rownames(res), 50)
mat <- assay(vsd)[top50, ]
pheatmap(mat,
         annotation_col = as.data.frame(colData(dds)["condition"]),
         filename = file.path(output_dir, "heatmap.png"))

cat("Plots saved\n")
