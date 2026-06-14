# ============================================================
# DESeq2 Differential Expression Analysis
# Project: RNA-seq Pipeline — GSE58135 Phase 3
# Contrast: TNBC vs ER+ (human breast cancer)
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
counts_file   <- args[1]
metadata_file <- args[2]
output_file   <- args[3]

cat("Input counts:", counts_file, "\n")
cat("Metadata:", metadata_file, "\n")
cat("Output:", output_file, "\n")

library(DESeq2)
library(ggplot2)
library(pheatmap)

counts <- read.table(counts_file, header = TRUE, skip = 1, row.names = 1)
counts <- counts[, 6:ncol(counts)]
colnames(counts) <- gsub(".*star\\.(SRR[0-9]+)_Aligned.*", "\\1", colnames(counts))

metadata <- read.csv(metadata_file, row.names = 1)

stopifnot(all(colnames(counts) == rownames(metadata)))

dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = metadata,
    design    = ~ condition
)

dds$condition <- relevel(dds$condition, ref = "ERpos")
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "TNBC", "ERpos"))
res <- res[order(res$padj), ]
write.csv(as.data.frame(res), file = output_file)

cat("DESeq2 analysis complete\n")

output_dir <- dirname(output_file)

vsd <- vst(dds, blind = FALSE)
pca_plot <- plotPCA(vsd, intgroup = "condition")
ggsave(file.path(output_dir, "pca_plot.png"), pca_plot)

res_df <- as.data.frame(res)
res_df$significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1

volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("grey", "red")) +
    theme_minimal() +
    labs(title = "Volcano Plot: TNBC vs ER+", x = "log2 Fold Change", y = "-log10 adjusted p-value")
ggsave(file.path(output_dir, "volcano_plot.png"), volcano)

top50 <- head(rownames(res), 50)
mat <- assay(vsd)[top50, ]
pheatmap(mat,
         annotation_col = as.data.frame(colData(dds)["condition"]),
         scale = "row",
         main = "Top 50 DEGs: TNBC vs ER+",
         filename = file.path(output_dir, "heatmap.png"))

cat("Plots saved\n")
