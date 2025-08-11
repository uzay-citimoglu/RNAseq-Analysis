# --- USER INPUTS (EDIT THESE AS NEEDED) ---------------------------------------
# BAM paths (aligned, sorted BAMs)  # EDIT
# GTF path (annotation)             # EDIT
# Output paths (counts/fpkm TSVs)   # EDIT
# Threads for featureCounts         # EDIT (nthreads)
# Condition labels / group factor   # EDIT (group <- factor(...))
# Filtering thresholds              # EDIT (FPKM > 1 in >= 3 samples; FDR/logFC cutoffs)
# PCA/plot labeling and colors      # EDIT (sample names, color map)
# -----------------------------------------------------------------------------

library(Rsubread)
library(BiocParallel)
library(limma)
library(edgeR)
library(rtracklayer)

bams <- c(                                                            # EDIT
  "/PATH/TO/Alignment/yourdata_1.sorted.bam",
  "/PATH/TO/Alignment/yourdata_2.sorted.bam",
  "/PATH/TO/Alignment/yourdata2_1.sorted.bam",
  "/PATH/TO/Alignment/yourdata2_2.sorted.bam"
)

counts <- featureCounts(
  files = bams,
  annot.ext = "/PATH/TO/ANNOTATION/gencode.v48.primary_assembly.basic.annotation.gtf", # EDIT
  isGTFAnnotationFile = TRUE,
  isPairedEnd = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  countMultiMappingReads = TRUE,
  nthreads = 20                                                     # EDIT
)

write.table(counts$counts,
            "/PATH/TO/Alignment/counts.tsv",  # EDIT
            quote = FALSE, sep = "\t", col.names = NA, row.names = TRUE)
write.table(counts$counts, "counts.tsv",
            quote = FALSE, sep = "\t", col.names = NA, row.names = TRUE)

z1 <- DGEList(counts = counts$counts, genes = counts$annotation[, c("GeneID","Length")])
y1 <- calcNormFactors(z1)
FPKM1 <- rpkm(y1)

write.table(FPKM1,
            "/PATH/TO/Alignment/fpkm_values.tsv",  # EDIT
            quote = FALSE, sep = "\t", col.names = NA, row.names = TRUE)
write.table(FPKM1, "fpkm_values.tsv",
            quote = FALSE, sep = "\t", col.names = NA, row.names = TRUE)

counts <- read.delim("/PATH/TO/Alignment/counts.tsv",   # EDIT
                     sep = "\t", row.names = 1, check.names = FALSE)
fpkm   <- read.delim("/PATH/TO/Alignment/fpkm_values.tsv", # EDIT
                     sep = "\t", row.names = 1, check.names = FALSE)

counts$GeneId <- sub("\\.\\d+$", "", rownames(counts))
fpkm$GeneId   <- sub("\\.\\d+$", "", rownames(fpkm))

gtf <- import("/PATH/TO/ANNOTATION/gencode.v48.primary_assembly.basic.annotation.gtf")  # EDIT
gtf <- gtf[gtf$type == "gene"]
gtf <- data.frame(
  seqid     = as.character(seqnames(gtf)),
  gene_id   = gtf$gene_id,
  gene_type = gtf$gene_type,
  gene_name = gtf$gene_name,
  stringsAsFactors = FALSE
)
gtf$gene_id <- sub("\\.\\d+$", "", gtf$gene_id)

annotated_counts <- merge(counts, gtf, by.x = "GeneId", by.y = "gene_id", all.x = TRUE)
annotated_fpkm  <- merge(fpkm,   gtf, by.x = "GeneId", by.y = "gene_id", all.x = TRUE)

fpkm_filt <- annotated_fpkm[ rowSums(annotated_fpkm[, 2:5] > 1) >= 3, ]     # EDIT thresholds
counts_filt <- annotated_counts[ annotated_counts$GeneId %in% fpkm_filt$GeneId, ]

countsnum  <- as.matrix(counts_filt[, 2:5])
logcounts  <- log2(countsnum + 1)
logcounts  <- logcounts[ apply(logcounts, 1, function(x) var(x) > 0), ]
pca        <- prcomp(t(logcounts), scale. = TRUE)
pca_df     <- data.frame(pca$x[, 1:2])
pca_df$Sample <- colnames(logcounts)

library(ggplot2)
variance   <- pca$sdev^2
percentVar <- round(variance / sum(variance) * 100, 1)
pca_df$Sample <- sub("_Aligned\\.sortedByCoord\\.out\\.bam", "", pca_df$Sample)

ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample, color = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3) +
  theme_minimal() +
  labs(
    title = "PCA of the Samples",
    x = paste0("PC1 (", percentVar[1], "%)"),
    y = paste0("PC2 (", percentVar[2], "%)")
  ) +
  scale_color_manual(values = c(                                   # EDIT
    "MV_1"        = "red",
    "MV_2"        = "red",
    "MV_SORRES_1" = "green",
    "MV_SORRES_2" = "green"
  ))
ggsave("pca.png", width = 8, height = 6, dpi = 300)

## Differential expression
library(edgeR)
group <- factor(c("yourdata", "yourdata", "yourdata2", "yourdata2"))                 # EDIT groups (order matches BAMs)
dge <- DGEList(counts = counts_filt[, 2:5], group = group)
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge)
fit <- glmFit(dge)
lrt <- glmLRT(fit)
deg_results <- topTags(lrt, n = Inf)$table
sig_genes <- deg_results[ deg_results$FDR < 0.05 & abs(deg_results$logFC) > 1, ]  # EDIT cutoffs
deg_results$significant <- with(deg_results, FDR < 0.05 & abs(logFC) > 1)         # EDIT cutoffs

## TMM normalization & exact test
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge, verbose = TRUE)
dge <- estimateTagwiseDisp(dge)
test <- exactTest(dge)
deg_results2 <- topTags(test, n = Inf)$table
deg_results2$significant <- with(deg_results2, FDR < 0.05 & abs(logFC) > 1)      # EDIT cutoffs

## Save DEG tables
write.table(deg_results,  "edgeR_glm_DEG.tsv",   sep = "\t", quote = FALSE, row.names = TRUE)
write.table(deg_results2, "edgeR_exact_DEG.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

## Volcano plots
ggplot(deg_results, aes(x = logFC, y = -log10(FDR), color = significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +             # EDIT colors
  theme_minimal() +
  labs(
    title = "Volcano Plot of Differentially Expressed Genes",
    x = "log2 Fold Change",
    y = "-log10 FDR",
    color = "Significant"
  )
ggsave("volcano_glm.png", width = 8, height = 6, dpi = 300)

ggplot(deg_results2, aes(x = logFC, y = -log10(FDR), color = significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "green")) +           # EDIT colors
  theme_minimal() +
  labs(
    title = "Volcano Plot of Differentially Expressed Genes",
    x = "log2 Fold Change",
    y = "-log10 FDR",
    color = "Significant"
  )
ggsave("volcano_exact.png", width = 8, height = 6, dpi = 300)

## Heatmap
library(pheatmap)
logcpm <- cpm(dge, log = TRUE)
sig_logcpm <- logcpm[ rownames(logcpm) %in% rownames(sig_genes), ]
scaled_matrix <- t(scale(t(sig_logcpm)))
pheatmap(scaled_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         filename = "heatmap.png",
         width = 8, height = 10)

## KEGG and GO analysis
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

entrez_ids <- unique(rownames(sig_genes))                              # EDIT keyType if needed

# KEGG
kegg_result <- enrichKEGG(gene = entrez_ids, organism = "hsa",
                          pAdjustMethod = "BH", pvalueCutoff = 0.05)    # EDIT thresholds
if (nrow(as.data.frame(kegg_result)) > 0) {
  kegg_result <- enrichplot::pairwise_termsim(kegg_result)
  pdf("kegg_plots.pdf", width = 10, height = 8)                         # EDIT filenames if desired
  print(dotplot(kegg_result, showCategory = 20) + ggtitle("KEGG Pathway Enrichment"))
  print(barplot(kegg_result, showCategory = 20))
  try(print(emapplot(kegg_result, showCategory = 20)), silent = TRUE)
  dev.off()
  ggsave("kegg_dot.png", dotplot(kegg_result, showCategory = 20) + ggtitle("KEGG Pathway Enrichment"),
         width = 8, height = 6, dpi = 300)
  ggsave("kegg_bar.png", barplot(kegg_result, showCategory = 20),
         width = 8, height = 6, dpi = 300)
  tmp_kegg_emap <- try(emapplot(kegg_result, showCategory = 20), silent = TRUE)
  if (!inherits(tmp_kegg_emap, "try-error"))
    ggsave("kegg_emap.png", tmp_kegg_emap, width = 10, height = 8, dpi = 300)

  write.table(as.data.frame(kegg_result), "kegg_enrichment.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# GO
ego <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", # EDIT keyType if needed
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)          # EDIT thresholds
if (nrow(as.data.frame(ego)) > 0) {
  ego <- enrichplot::pairwise_termsim(ego)
  pdf("go_plots.pdf", width = 10, height = 8)                                   # EDIT filenames if desired
  print(dotplot(ego, showCategory = 20) + ggtitle("GO Biological Process Enrichment"))
  print(barplot(ego, showCategory = 20))
  try(print(emapplot(ego, showCategory = 20)), silent = TRUE)
  dev.off()
  ggsave("go_dot.png", dotplot(ego, showCategory = 20) + ggtitle("GO Biological Process Enrichment"),
         width = 8, height = 6, dpi = 300)
  ggsave("go_bar.png", barplot(ego, showCategory = 20),
         width = 8, height = 6, dpi = 300)
  tmp_go_emap <- try(emapplot(ego, showCategory = 20), silent = TRUE)
  if (!inherits(tmp_go_emap, "try-error"))
    ggsave("go_emap.png", tmp_go_emap, width = 10, height = 8, dpi = 300)

  write.table(as.data.frame(ego), "go_enrichment.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE)
}
