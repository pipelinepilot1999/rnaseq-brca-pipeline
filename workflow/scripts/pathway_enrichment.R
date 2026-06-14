args <- commandArgs(trailingOnly = TRUE)
results_file <- args[1]
output_dir   <- args[2]

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

res <- read.csv(results_file, row.names = 1)
cat("Total genes:", nrow(res), "\n")

up_tnbc  <- rownames(res[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange >  1, ])
up_erpos <- rownames(res[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -1, ])

cat("Genes UP in TNBC:", length(up_tnbc), "\n")
cat("Genes UP in ER+:", length(up_erpos), "\n")

to_entrez <- function(genes) {
    converted <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    return(converted$ENTREZID)
}

tnbc_entrez  <- to_entrez(up_tnbc)
erpos_entrez <- to_entrez(up_erpos)

go_tnbc <- enrichGO(gene = tnbc_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                    pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
go_erpos <- enrichGO(gene = erpos_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                     pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

kegg_tnbc <- enrichKEGG(gene = tnbc_entrez, organism = "hsa", pvalueCutoff = 0.05)
kegg_erpos <- enrichKEGG(gene = erpos_entrez, organism = "hsa", pvalueCutoff = 0.05)

ggsave(file.path(output_dir, "GO_TNBC_dotplot.png"),
       dotplot(go_tnbc, showCategory = 20) + ggtitle("GO BP: UP in TNBC"), width = 10, height = 8)
ggsave(file.path(output_dir, "GO_ERpos_dotplot.png"),
       dotplot(go_erpos, showCategory = 20) + ggtitle("GO BP: UP in ER+"), width = 10, height = 8)
ggsave(file.path(output_dir, "KEGG_TNBC_dotplot.png"),
       dotplot(kegg_tnbc, showCategory = 20) + ggtitle("KEGG: UP in TNBC"), width = 10, height = 8)
ggsave(file.path(output_dir, "KEGG_ERpos_dotplot.png"),
       dotplot(kegg_erpos, showCategory = 20) + ggtitle("KEGG: UP in ER+"), width = 10, height = 8)

write.csv(as.data.frame(go_tnbc),   file.path(output_dir, "GO_TNBC_results.csv"))
write.csv(as.data.frame(go_erpos),  file.path(output_dir, "GO_ERpos_results.csv"))
write.csv(as.data.frame(kegg_tnbc), file.path(output_dir, "KEGG_TNBC_results.csv"))
write.csv(as.data.frame(kegg_erpos),file.path(output_dir, "KEGG_ERpos_results.csv"))

cat("Pathway enrichment complete\n")
