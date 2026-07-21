# ==============================================================================
# 03 - Volcano plots
#
# Input : deg_results + gene_id_map saved by script 01
# Output: three sets of volcano plots per comparison, saved as SVG:
#           1) plain volcano (colored by significance only)
#           2) volcano labeled with a curated gene list
#           3) volcano labeled with a separate, stricter FDR/log2FC cutoff
#
# Note: run this script with the project root as the working directory.
# ==============================================================================

library(ggplot2)
library(ggrepel)
library(readxl)
library(dplyr)

# ---------------------------------------------------------------------------
# Parameters (must match script 01, since they define the results folder tag)
# ---------------------------------------------------------------------------
alpha_thresh <- 0.1
lfc_thresh   <- 0.5
expr_thresh  <- 5
min_frac     <- 0.75

results_tag <- paste0(
  "Thr", expr_thresh, "_Frac", min_frac,
  "_FDR", alpha_thresh, "_LFC", lfc_thresh
)
tables_dir  <- file.path("results", "tables", results_tag)
figures_dir <- file.path("results", "figures", results_tag, "VolcanoPlots")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

deg_results <- readRDS(file.path(tables_dir, "deg_results.rds"))
gene_id_map <- readRDS(file.path(tables_dir, "gene_id_map.rds"))
ensembl_ids  <- gene_id_map$ensembl_ids
gene_symbols <- gene_id_map$gene_symbols

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

# ==============================================================================
# 1) Plain volcano plot (no gene labels)
# ==============================================================================
plot_volcano <- function(deg_results, comparison, label_map) {
  res <- deg_results[[comparison]]$all
  res$Significance <- "Not Significant"
  res$Significance[res$padj < alpha_thresh & res$log2FoldChange >= lfc_thresh] <- "Positive"
  res$Significance[res$padj < alpha_thresh & res$log2FoldChange <= -lfc_thresh] <- "Negative"

  nice_label <- get_nice_comparison_label(comparison, label_map)

  ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(size = 3, alpha = 0.5, aes(color = Significance)) +
    geom_point(
      data = subset(res, Significance != "Not Significant"),
      aes(x = log2FoldChange, y = -log10(padj), color = Significance), alpha = 0.7, size = 3
    ) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "black", alpha = 0.3, linewidth = 0.7) +
    geom_hline(yintercept = -log10(alpha_thresh), linetype = "dashed", color = "black", alpha = 0.3, linewidth = 0.7) +
    scale_color_manual(values = c("Not Significant" = "gray", "Positive" = "red", "Negative" = "blue")) +
    labs(
      title = paste("DEGs: ", nice_label),
      x = "Log2 Fold Change",
      y = "-Log10 (Adj. p-value)",
      color = "Regulation"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.position = "top",
      legend.justification = "center",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8)
    )
}

for (comparison in names(deg_results)) {
  p <- plot_volcano(deg_results, comparison, group_labels_map)
  out <- file.path(figures_dir, paste0("Volcano_", comparison, ".svg"))
  svglite::svglite(filename = out, width = 10, height = 8, bg = "transparent")
  print(p)
  dev.off()
}

# ==============================================================================
# 2) Volcano labeled with a curated gene symbol list
# ==============================================================================
curated_list_path <- "data/metadata/Curated_List.xlsx"
my_gene_symbol_list <- as.character(read_excel(curated_list_path)[[1]])
my_gene_symbol_list <- trimws(my_gene_symbol_list)

plot_volcano_labelled_by_symbol <- function(deg_results, comparison, label_map, gene_symbol_list, ensembl_ids, gene_symbols) {
  res <- deg_results[[comparison]]$all

  res$Significance <- "Not Significant"
  res$Significance[res$padj < alpha_thresh & res$log2FoldChange >= lfc_thresh] <- "Positive"
  res$Significance[res$padj < alpha_thresh & res$log2FoldChange <= -lfc_thresh] <- "Negative"

  nice_label <- get_nice_comparison_label(comparison, label_map)

  res$GeneSymbol <- gene_symbols[match(rownames(res), ensembl_ids)]
  res$GeneSymbol <- trimws(res$GeneSymbol)

  # Only significant genes that also appear in the curated list get a label
  sig_indices <- which(res$Significance != "Not Significant")
  sig_symbols <- res$GeneSymbol[sig_indices]

  label_indices <- sig_indices[sig_symbols %in% gene_symbol_list]
  label_df <- res[label_indices, , drop = FALSE]

  res$Highlight <- "Other"
  res$Highlight[label_indices] <- res$Significance[label_indices]

  p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(data = subset(res, Highlight == "Other"), aes(color = Significance), alpha = 0.7, size = 3) +
    geom_point(data = label_df, aes(color = Significance), size = 5, stroke = 1.1, shape = 21, fill = "white", show.legend = FALSE) +
    geom_point(data = label_df, aes(color = Significance), size = 5, stroke = 2, show.legend = FALSE) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "black", alpha = 0.3, linewidth = 0.7) +
    geom_hline(yintercept = -log10(alpha_thresh), linetype = "dashed", color = "black", alpha = 0.3, linewidth = 0.7) +
    scale_color_manual(values = c(
      "Not Significant" = "gray",
      "Positive" = "#C62828",
      "Negative" = "#1565C0"
    )) +
    labs(
      title = paste("DEGs: ", nice_label),
      x = "Log2 Fold Change",
      y = "-Log10 (Adj. p-value)",
      color = "Regulation"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.position = "top",
      legend.justification = "center",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8)
    ) +
    geom_text_repel(
      data = label_df,
      aes(label = GeneSymbol),
      color = "black",
      size = 6,
      max.overlaps = 20,
      box.padding = 1.2,
      point.padding = 0.6,
      segment.color = "black",
      segment.size = 1.2,
      fontface = "bold"
    )

  return(p)
}

for (comparison in names(deg_results)) {
  p <- plot_volcano_labelled_by_symbol(
    deg_results,
    comparison,
    group_labels_map,
    gene_symbol_list = my_gene_symbol_list,
    ensembl_ids = ensembl_ids,
    gene_symbols = gene_symbols
  )
  out <- file.path(figures_dir, paste0("Volcano_", comparison, "_Labeled_OldList.svg"))
  svglite::svglite(filename = out, width = 10, height = 8, bg = "transparent")
  print(p)
  dev.off()
}

# ==============================================================================
# 3) Volcano labeled by a separate, stricter FDR/log2FC cutoff
# ==============================================================================
# Note: the fdr_cutoff / log2fc_cutoff below are independent of alpha_thresh /
# lfc_thresh above - they only decide which points get a text label, not
# which points count as significant (that's still alpha_thresh / lfc_thresh).
plot_volcano_labelled_by_threshold <- function(deg_results, comparison, label_map, ensembl_ids, gene_symbols,
                                                fdr_cutoff = 0.05, log2fc_cutoff = 1.5) {
  res <- as.data.frame(deg_results[[comparison]]$all)

  res$Significance <- "Not Significant"
  res$Significance[res$padj < alpha_thresh & res$log2FoldChange >= lfc_thresh] <- "Positive"
  res$Significance[res$padj < alpha_thresh & res$log2FoldChange <= -lfc_thresh] <- "Negative"

  nice_label <- get_nice_comparison_label(comparison, label_map)

  res$GeneSymbol <- gene_symbols[match(rownames(res), ensembl_ids)]
  res$GeneSymbol <- trimws(res$GeneSymbol)

  label_df <- res %>%
    dplyr::filter(
      padj < fdr_cutoff,
      abs(log2FoldChange) > log2fc_cutoff,
      Significance != "Not Significant"
    )

  res$Highlight <- "Other"
  if (nrow(label_df) > 0) {
    label_indices <- which(rownames(res) %in% rownames(label_df))
    res$Highlight[label_indices] <- res$Significance[label_indices]
  }

  p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(data = subset(res, Highlight == "Other"), aes(color = Significance), alpha = 0.7, size = 3) +
    geom_point(data = label_df, aes(color = Significance), size = 3, stroke = 1.1, shape = 21, fill = "white", show.legend = FALSE) +
    geom_point(data = label_df, aes(color = Significance), size = 3, stroke = 2, show.legend = FALSE) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "black", alpha = 0.3, linewidth = 0.7) +
    geom_hline(yintercept = -log10(alpha_thresh), linetype = "dashed", color = "black", alpha = 0.3, linewidth = 0.7) +
    scale_color_manual(values = c(
      "Not Significant" = "gray",
      "Positive" = "#C62828",
      "Negative" = "#1565C0"
    )) +
    labs(
      title = paste("DEGs: ", nice_label),
      x = "Log2 Fold Change",
      y = "-Log10 (Adj. p-value)",
      color = "Regulation"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 16),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.position = "top",
      legend.justification = "center",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8)
    ) +
    geom_text_repel(
      data = label_df,
      aes(label = GeneSymbol),
      color = "black",
      size = 6,
      max.overlaps = 20,
      box.padding = 1.2,
      point.padding = 0.6,
      segment.color = "black",
      segment.size = 1.2,
      fontface = "bold"
    )

  return(p)
}

for (comparison in names(deg_results)) {
  p <- plot_volcano_labelled_by_threshold(
    deg_results,
    comparison,
    group_labels_map,
    ensembl_ids = ensembl_ids,
    gene_symbols = gene_symbols,
    fdr_cutoff = 0.009,
    log2fc_cutoff = 0.5
  )
  out <- file.path(figures_dir, paste0("Volcano_", comparison, "_Labeled_byThreshold.svg"))
  svglite::svglite(filename = out, width = 10, height = 8, bg = "transparent")
  print(p)
  dev.off()
}
