# ==============================================================================
# 01 - Gene filtering and differential expression analysis (DESeq2)
#
# Input : annotated count matrix (Ensembl IDs x samples, last row = group labels)
# Output: filtered count matrix, expression distribution plot, DESeq2 results,
#         and an Excel workbook with significant DEGs for every pairwise
#         comparison between conditions (Treg / Tfr / Tfh).
#
# Note: run this script with the project root as the working directory
# (e.g. open the .Rproj file, or setwd() to the repo root manually).
# ==============================================================================

library(dplyr)
library(readxl)
library(ggplot2)
library(biomaRt)
library(DESeq2)
library(openxlsx)
library(here)

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
alpha_thresh <- 0.1     # FDR threshold for DESeq2 significance
lfc_thresh   <- 0.5     # log2 fold-change threshold
expr_thresh  <- 5       # minimum count to consider a gene "expressed" in a sample
min_frac     <- 0.75    # fraction of samples (within a group) that must clear expr_thresh

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

setwd('/Users/amin/Library/Mobile Documents/com~apple~CloudDocs/PhD Project/Simulation/TGFb-IL2-Tfr-RNAseq')

annotated_data_path <- "data/raw_counts/annotated_count_data.xlsx"

results_tag <- paste0(
  "Thr", expr_thresh, "_Frac", min_frac,
  "_FDR", alpha_thresh, "_LFC", lfc_thresh
)
tables_dir  <- file.path("results", "tables", results_tag)
figures_dir <- file.path("results", "figures", results_tag)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load annotated count matrix
# ---------------------------------------------------------------------------
# Read everything as text first so the group-label row (mixed types) doesn't
# get coerced/mangled by read_excel's type guessing.
annotated_raw <- read_excel(annotated_data_path, col_types = "text")
annotated_raw <- as.data.frame(annotated_raw)

# Group labels live in the last row, between the ID column and the gene-symbol column
group_labels <- as.character(unlist(annotated_raw[nrow(annotated_raw), -c(1, ncol(annotated_raw))]))

# Drop that row now that the labels are stored separately
annotated_data <- annotated_raw[-nrow(annotated_raw), ]

# Gene symbols / Ensembl IDs are in the first and last column respectively
gene_symbols <- annotated_data[, ncol(annotated_data)]
ensembl_ids  <- annotated_data[, 1]

# Everything in between is the actual count data
data_numeric <- annotated_data[, -c(1, ncol(annotated_data))]
expression_values <- data.frame(lapply(data_numeric, as.numeric))
rownames(expression_values) <- ensembl_ids

saveRDS(list(ensembl_ids = ensembl_ids, gene_symbols = gene_symbols),
        file.path(tables_dir, "gene_id_map.rds"))

# ---------------------------------------------------------------------------
# Restrict to protein-coding genes (via Ensembl BioMart)
# ---------------------------------------------------------------------------
connect_to_biomart <- function() {
  tryCatch({
    useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "https://www.ensembl.org")
  }, error = function(e1) {
    cat("Main Ensembl server unavailable, trying mirror...\n")
    tryCatch({
      useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", mirror = "useast")
    }, error = function(e2) {
      stop("Failed to connect to Ensembl. Please check your internet connection or try again later.")
    })
  })
}
mart <- connect_to_biomart()

protein_coding_genes <- getBM(
  attributes = c("ensembl_gene_id", "transcript_biotype"),
  filters = "transcript_biotype",
  values = "protein_coding",
  mart = mart
)
if (nrow(protein_coding_genes) == 0) stop("No protein-coding genes retrieved.")
protein_coding_gene_ids <- unique(protein_coding_genes$ensembl_gene_id)

expression_values <- expression_values[rownames(expression_values) %in% protein_coding_gene_ids, ]

initial_gene_count     <- nrow(annotated_data)
protein_coding_count   <- nrow(expression_values)
message("Initial genes: ", initial_gene_count)
message("Protein-coding genes: ", protein_coding_count)

# ---------------------------------------------------------------------------
# Expression-based gene filtering
# ---------------------------------------------------------------------------
# Keep a gene if, in at least one condition, it clears expr_thresh in at
# least min_frac of that condition's samples.
filter_genes <- function(data, groups, threshold = expr_thresh, min_fraction = min_frac) {
  unique_groups <- unique(groups)
  keep_genes <- rep(FALSE, nrow(data))
  for (group in unique_groups) {
    group_indices <- groups %in% group
    count_above_threshold <- rowSums(data[, group_indices] > threshold)
    keep_genes <- keep_genes | (count_above_threshold >= ceiling(sum(group_indices) * min_fraction))
  }
  return(data[keep_genes, ])
}

filtered_expression_values <- filter_genes(expression_values, group_labels)
expression_filtered_count  <- nrow(filtered_expression_values)
message("After filtering: ", expression_filtered_count)

# ---------------------------------------------------------------------------
# Expression distribution before/after filtering (sanity check plot)
# ---------------------------------------------------------------------------
expression_long_before <- data.frame(expression = as.vector(as.matrix(expression_values)), state = "Before Filtering")
expression_long_after  <- data.frame(expression = as.vector(as.matrix(filtered_expression_values)), state = "After Filtering")

expression_long_before$expression <- as.numeric(expression_long_before$expression)
expression_long_after$expression  <- as.numeric(expression_long_after$expression)

expression_long_before <- na.omit(expression_long_before)
expression_long_after  <- na.omit(expression_long_after)

expression_long_before$log_expression <- log1p(expression_long_before$expression)
expression_long_after$log_expression  <- log1p(expression_long_after$expression)

expression_long <- rbind(expression_long_before, expression_long_after)
expression_long$state <- factor(expression_long$state, levels = c("Before Filtering", "After Filtering"))

png(file.path(figures_dir, "Expression_Distribution_Filtering.png"), width = 1600, height = 800)
print(
  ggplot(expression_long, aes(x = log_expression, fill = state)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ state, ncol = 2, scales = "free_x") +
    labs(
      title = "Expression Distribution Before and After Filtering",
      x = "Log(Expression + 1)",
      y = "Density"
    ) +
    theme_minimal() +
    geom_vline(
      data = subset(expression_long, state == "Before Filtering"),
      aes(xintercept = log1p(expr_thresh)), linetype = "dashed", color = "red", size = 1
    ) +
    geom_text(
      data = subset(expression_long, state == "Before Filtering"),
      aes(x = log1p(expr_thresh), y = 0.55, label = paste0("Cut-off = ", expr_thresh)),
      color = "black", angle = 90, vjust = -0.5, size = 4
    ) +
    theme(
      plot.title = element_text(size = 22),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      strip.text = element_text(size = 18),
      legend.position = "top",
      legend.justification = "center"
    )
)
dev.off()

# ---------------------------------------------------------------------------
# DESeq2: build dataset and fit the model
# ---------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = filtered_expression_values,
  colData = data.frame(condition = factor(group_labels)),
  design = ~ condition
)
dds <- DESeq(dds)

saveRDS(dds, file.path(tables_dir, "dds.rds"))

# ---------------------------------------------------------------------------
# DEG extraction helper
# ---------------------------------------------------------------------------
get_deg_results <- function(dds, contrast_group, alpha = alpha_thresh, lfc_threshold = lfc_thresh) {
  res <- results(dds, contrast = c("condition", contrast_group[1], contrast_group[2]), alpha = alpha)
  res <- res[!is.na(res$padj), ]

  cat("\nComparison:", paste(contrast_group, collapse = " vs "), "\n")
  cat("Total Genes in Results:", nrow(res), "\n")

  sig_genes    <- res[which(res$padj < alpha & abs(res$log2FoldChange) >= lfc_threshold), ]
  upregulated  <- sig_genes[sig_genes$log2FoldChange >= lfc_threshold, ]
  downregulated <- sig_genes[sig_genes$log2FoldChange <= -lfc_threshold, ]

  cat("Significant Genes (FDR <", alpha, ", |log2FC| >", lfc_threshold, "):", nrow(sig_genes), "\n")
  cat("Upregulated Genes:", nrow(upregulated), "\n")
  cat("Downregulated Genes:", nrow(downregulated), "\n")

  return(list(up = rownames(upregulated), down = rownames(downregulated), all = res))
}

# Readable condition names, used for plot labels and comparison titles downstream
group_labels_map <- list(
  "Treg" = "Treg",
  "Tfr"  = "Tfr",
  "Tfh"  = "Tfh"
)

get_nice_comparison_label <- function(comp_name, label_map) {
  groups <- unlist(strsplit(comp_name, "_vs_"))
  nice_names <- sapply(groups, function(g) label_map[[g]])
  paste(nice_names, collapse = " vs ")
}

# ---------------------------------------------------------------------------
# Run every pairwise comparison between conditions
# ---------------------------------------------------------------------------
unique_groups <- unique(group_labels)
deg_results <- list()

for (i in 1:(length(unique_groups) - 1)) {
  for (j in (i + 1):length(unique_groups)) {
    contrast_group <- c(unique_groups[i], unique_groups[j])
    comparison_name <- paste(contrast_group, collapse = "_vs_")
    deg_results[[comparison_name]] <- get_deg_results(dds, contrast_group)
  }
}

de_gene_lists <- deg_results

# ---------------------------------------------------------------------------
# Export significant DEGs to Excel (one up/down sheet pair per comparison)
# ---------------------------------------------------------------------------
wb <- createWorkbook()

for (comp_name in names(deg_results)) {
  res_table <- deg_results[[comp_name]]$all
  sig_res <- res_table[which(res_table$padj < alpha_thresh & abs(res_table$log2FoldChange) >= lfc_thresh), ]
  sig_res$qvalue <- sig_res$padj

  # Upregulated
  up_res <- sig_res[sig_res$log2FoldChange >= lfc_thresh, ]
  if (nrow(up_res) > 0) {
    up_res <- up_res[order(-up_res$log2FoldChange), ]
    up_genes <- rownames(up_res)
    up_symbols <- gene_symbols[match(up_genes, ensembl_ids)]

    df_up <- data.frame(
      Rank = seq_len(nrow(up_res)),
      GeneSymbol = up_symbols,
      log2FoldChange = up_res$log2FoldChange,
      pvalue = up_res$pvalue,
      padj = up_res$padj,
      qvalue = up_res$qvalue
    )
    addWorksheet(wb, paste0(comp_name, "_up"))
    writeData(wb, sheet = paste0(comp_name, "_up"), df_up)
  }

  # Downregulated
  down_res <- sig_res[sig_res$log2FoldChange <= -lfc_thresh, ]
  if (nrow(down_res) > 0) {
    down_res <- down_res[order(down_res$log2FoldChange), ]
    down_genes <- rownames(down_res)
    down_symbols <- gene_symbols[match(down_genes, ensembl_ids)]

    df_down <- data.frame(
      Rank = seq_len(nrow(down_res)),
      GeneSymbol = down_symbols,
      log2FoldChange = down_res$log2FoldChange,
      pvalue = down_res$pvalue,
      padj = down_res$padj,
      qvalue = down_res$qvalue
    )
    addWorksheet(wb, paste0(comp_name, "_down"))
    writeData(wb, sheet = paste0(comp_name, "_down"), df_down)
  }
}

saveWorkbook(wb, file = file.path(tables_dir, "DEG_Significant_Results.xlsx"), overwrite = TRUE)
