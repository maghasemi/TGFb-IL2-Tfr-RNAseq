# ==============================================================================
# 02 - PCA and expression QC plots
#
# Input : dds object saved by script 01 (fitted DESeq2 dataset)
# Output: PCA plots (PC1 vs PC2, PC1 vs PC3), a scree plot, and a boxplot of
#         VST-normalized expression per condition.
#
# Note: run this script with the project root as the working directory.
# ==============================================================================

library(DESeq2)
library(ggplot2)
library(ggfortify)

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
figures_dir <- file.path("results", "figures", results_tag)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

dds <- readRDS(file.path(tables_dir, "dds.rds"))

# Readable condition names, same mapping used across all scripts
group_labels_map <- list(
  "Treg" = "Treg",
  "Tfr"  = "Tfr",
  "Tfh"  = "Tfh"
)
group_labels <- as.character(colData(dds)$condition)

# ---------------------------------------------------------------------------
# Variance-stabilizing transformation + PCA
# ---------------------------------------------------------------------------
vst_data <- vst(dds)
pca_data <- prcomp(t(assay(vst_data)))
readable_conditions <- unlist(group_labels_map[group_labels])

explained_var <- summary(pca_data)$importance[2, ] * 100
pc_labels_pc1_pc2 <- paste0(c("PC1", "PC2"), " (", round(explained_var[1:2], 1), "%)")
pc_labels_pc1_pc3 <- paste0(c("PC1", "PC3"), " (", round(explained_var[c(1, 3)], 1), "%)")

# PC1 vs PC2
png(file.path(figures_dir, "PCA_PC1_vs_PC2.png"), width = 1000, height = 800)
print(
  autoplot(pca_data,
    data = data.frame(condition = factor(readable_conditions)),
    colour = "condition", size = 6
  ) +
    labs(
      title = "PCA Plot of RNA-seq Data (PC1 vs PC2)",
      x = pc_labels_pc1_pc2[1],
      y = pc_labels_pc1_pc2[2]
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
)
dev.off()

# PC1 vs PC3
png(file.path(figures_dir, "PCA_PC1_vs_PC3.png"), width = 1000, height = 800)
print(
  autoplot(pca_data,
    data = data.frame(condition = factor(readable_conditions)),
    colour = "condition", size = 4, x = 1, y = 3
  ) +
    labs(
      title = "PCA Plot of RNA-seq Data (PC1 vs PC3)",
      x = pc_labels_pc1_pc3[1],
      y = pc_labels_pc1_pc3[2]
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
)
dev.off()

# ---------------------------------------------------------------------------
# Scree plot
# ---------------------------------------------------------------------------
scree_data <- data.frame(PC = paste0("PC", 1:length(explained_var)), Variance = explained_var)

png(file.path(figures_dir, "Scree_Plot.png"), width = 1000, height = 600)
print(
  ggplot(scree_data, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(
      title = "Scree Plot of Principal Components",
      x = "Principal Component",
      y = "Percentage of Explained Variance"
    ) +
    theme_minimal()
)
dev.off()

# ---------------------------------------------------------------------------
# Expression boxplot (VST-normalized), per condition
# ---------------------------------------------------------------------------
expression_long <- data.frame(
  expression = as.vector(assay(vst_data)),
  sample = rep(colnames(vst_data), each = nrow(vst_data)),
  condition = rep(group_labels, each = nrow(vst_data))
)

expression_long$condition <- unlist(group_labels_map[expression_long$condition])
expression_long$condition <- factor(expression_long$condition, levels = c("Treg", "Tfr", "Tfh"))

png(file.path(figures_dir, "Boxplot_VST_Expression.png"), width = 1000, height = 800)
print(
  ggplot(expression_long, aes(x = condition, y = expression, fill = condition)) +
    geom_boxplot() +
    labs(
      title = "Boxplot of Gene Expression per Condition",
      x = "Condition",
      y = "Expression (VST-normalized)"
    ) +
    theme_minimal()
)
dev.off()
