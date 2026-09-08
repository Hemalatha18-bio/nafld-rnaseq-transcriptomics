#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: pathway_enrichment.R <all_deseq2.tsv> <significant_deseq2.tsv> <output_dir>")
}

all_path <- args[[1]]
sig_path <- args[[2]]
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

all_results <- read.delim(all_path, stringsAsFactors = FALSE, check.names = FALSE)
significant <- read.delim(sig_path, stringsAsFactors = FALSE, check.names = FALSE)

required <- c("gene_id", "padj", "log2FoldChange")
if (!all(required %in% colnames(all_results))) {
  stop("all-results table is missing required DESeq2 columns")
}
if (!all(required %in% colnames(significant))) {
  stop("significant-results table is missing required DESeq2 columns")
}

strip_version <- function(ids) sub("\\..*$", "", as.character(ids))
all_ensembl <- unique(strip_version(all_results$gene_id))
sig_ensembl <- unique(strip_version(significant$gene_id))

all_entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = all_ensembl,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)
sig_entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = sig_ensembl,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

mapping <- data.frame(
  ensembl_gene_id = names(all_entrez),
  entrez_id = unname(all_entrez),
  stringsAsFactors = FALSE
)
write.table(
  mapping,
  file.path(output_dir, "ensembl_to_entrez_mapping.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

universe <- unique(na.omit(unname(all_entrez)))
genes <- unique(na.omit(unname(sig_entrez)))

if (length(genes) == 0 || length(universe) == 0) {
  enrichment_table <- data.frame()
} else {
  enrichment <- enrichGO(
    gene = genes,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  enrichment_table <- as.data.frame(enrichment)
}

write.table(
  enrichment_table,
  file.path(output_dir, "go_bp_enrichment.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (nrow(enrichment_table) > 0) {
  plot_data <- head(enrichment_table[order(enrichment_table$p.adjust), ], 20)
  plot_data$Description <- factor(
    plot_data$Description,
    levels = rev(plot_data$Description)
  )
  p <- ggplot(plot_data, aes(x = -log10(p.adjust), y = Description)) +
    geom_point(aes(size = Count), alpha = 0.75) +
    labs(
      title = "GO Biological Process enrichment",
      subtitle = "FDR-significant DESeq2 genes; tested genes used as background",
      x = "-log10 adjusted p-value",
      y = NULL,
      size = "Gene count"
    ) +
    theme_bw(base_size = 11)
} else {
  p <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No GO BP terms passed the enrichment threshold") +
    xlim(-1, 1) + ylim(-1, 1) +
    theme_void() +
    labs(title = "GO Biological Process enrichment")
}

ggsave(file.path(output_dir, "go_bp_enrichment.png"), p, width = 8, height = 6, dpi = 160)

summary_table <- data.frame(
  metric = c(
    "tested_ensembl_genes",
    "mapped_background_entrez_genes",
    "significant_ensembl_genes",
    "mapped_significant_entrez_genes",
    "enriched_go_bp_terms_fdr05"
  ),
  value = c(
    length(all_ensembl),
    length(universe),
    length(sig_ensembl),
    length(genes),
    nrow(enrichment_table)
  )
)
write.table(
  summary_table,
  file.path(output_dir, "go_bp_enrichment_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
