#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(paste(
    "Usage: differential_expression.R <count_matrix.tsv> <metadata.tsv> <output_dir>",
    "[contrast_column] [numerator] [denominator]"
  ))
}

count_path <- args[[1]]
metadata_path <- args[[2]]
output_dir <- args[[3]]
contrast_column <- ifelse(length(args) >= 4, args[[4]], "analysis_group")
numerator <- ifelse(length(args) >= 5, args[[5]], "NAFLD_F3_F4")
denominator <- ifelse(length(args) >= 6, args[[6]], "NAFLD_F0_F1")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read.delim(count_path, row.names = 1, check.names = FALSE)
metadata <- read.delim(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)

if (!("sample_id" %in% colnames(metadata))) {
  stop("metadata must contain sample_id")
}
if (!(contrast_column %in% colnames(metadata))) {
  stop(paste("metadata is missing contrast column:", contrast_column))
}
if (anyDuplicated(metadata$sample_id)) {
  stop("metadata contains duplicate sample_id values")
}

rownames(metadata) <- metadata$sample_id
missing_metadata <- setdiff(colnames(counts), rownames(metadata))
if (length(missing_metadata) > 0) {
  stop(paste("metadata missing count samples:", paste(missing_metadata, collapse = ", ")))
}

metadata <- metadata[colnames(counts), , drop = FALSE]
selected <- metadata[[contrast_column]] %in% c(numerator, denominator)
metadata_sub <- metadata[selected, , drop = FALSE]
counts_sub <- counts[, rownames(metadata_sub), drop = FALSE]

if (nrow(metadata_sub) < 4) {
  stop("contrast contains too few samples")
}
if (!all(c(numerator, denominator) %in% unique(metadata_sub[[contrast_column]]))) {
  stop("both numerator and denominator levels must be present")
}

if (any(is.na(counts_sub)) || any(counts_sub < 0) || any(counts_sub %% 1 != 0)) {
  stop("count matrix must contain non-negative integer counts with no missing values")
}

# Conservative expression filter before DESeq2 model fitting.
minimum_samples <- min(10L, ncol(counts_sub))
keep_gene <- rowSums(counts_sub >= 10) >= minimum_samples
counts_sub <- counts_sub[keep_gene, , drop = FALSE]
if (nrow(counts_sub) == 0) {
  stop("expression filter removed all genes")
}

coldata <- data.frame(
  group = factor(metadata_sub[[contrast_column]]),
  row.names = rownames(metadata_sub)
)
coldata$group <- relevel(coldata$group, ref = denominator)

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_sub)),
  colData = coldata,
  design = ~ group
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", numerator, denominator), alpha = 0.05)

result_table <- as.data.frame(res)
result_table$gene_id <- rownames(result_table)
result_table <- result_table[, c(
  "gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
)]
result_table <- result_table[order(result_table$padj, na.last = TRUE), ]
write.table(
  result_table,
  file.path(output_dir, "deseq2_all_results.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

significant <- subset(
  result_table,
  !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1
)
write.table(
  significant,
  file.path(output_dir, "deseq2_significant_fdr05_lfc1.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

normalized <- counts(dds, normalized = TRUE)
write.table(
  data.frame(gene_id = rownames(normalized), normalized, check.names = FALSE),
  file.path(output_dir, "deseq2_normalized_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

plot_table <- result_table[!is.na(result_table$padj), ]
plot_table$neg_log10_padj <- -log10(pmax(plot_table$padj, .Machine$double.xmin))
plot_table$significant <- with(
  plot_table,
  padj < 0.05 & abs(log2FoldChange) >= 1
)

p <- ggplot(plot_table, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(shape = significant), alpha = 0.5, size = 1.2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = paste(numerator, "vs", denominator),
    subtitle = "DESeq2 on public GSE135251 gene counts",
    x = "log2 fold change",
    y = "-log10 adjusted p-value"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "volcano_deseq2.png"), p, width = 7, height = 5, dpi = 160)

summary_table <- data.frame(
  metric = c(
    "contrast_column", "numerator", "denominator", "samples_in_contrast",
    "genes_after_filter", "significant_fdr05_abs_lfc1"
  ),
  value = c(
    contrast_column, numerator, denominator, ncol(counts_sub),
    nrow(counts_sub), nrow(significant)
  )
)
write.table(
  summary_table,
  file.path(output_dir, "deseq2_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
