# ==============================================================================
# 04 - DEG overlap: Euler diagram, Venn diagram, and overlapping gene lists
#
# Input : deg_results + gene_id_map saved by script 01
# Output: Euler diagram (SVG), Venn diagram (SVG), and an Excel workbook
#         listing which genes fall into each overlap region across the
#         three pairwise comparisons.
#
# Note: run this script with the project root as the working directory.
# ==============================================================================

library(ggplot2)
library(ggvenn)
library(openxlsx)

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

deg_results <- readRDS(file.path(tables_dir, "deg_results.rds"))
gene_id_map <- readRDS(file.path(tables_dir, "gene_id_map.rds"))
ensembl_ids  <- gene_id_map$ensembl_ids
gene_symbols <- gene_id_map$gene_symbols

# ---------------------------------------------------------------------------
# DEG sets per comparison (up + down combined)
# ---------------------------------------------------------------------------
tfr_treg <- unique(c(deg_results[["Tfr_vs_Treg"]]$up, deg_results[["Tfr_vs_Treg"]]$down))
tfr_tfh  <- unique(c(deg_results[["Tfr_vs_Tfh"]]$up,  deg_results[["Tfr_vs_Tfh"]]$down))
tfh_treg <- unique(c(deg_results[["Tfh_vs_Treg"]]$up, deg_results[["Tfh_vs_Treg"]]$down))

# ---------------------------------------------------------------------------
# Extract overlapping genes for each region of the Euler/Venn diagram
# ---------------------------------------------------------------------------
get_symbols <- function(ensembl_id_vec, ensembl_ids, gene_symbols) {
  syms <- gene_symbols[match(ensembl_id_vec, ensembl_ids)]
  data.frame(
    EnsemblID = ensembl_id_vec,
    GeneSymbol = syms,
    stringsAsFactors = FALSE
  )
}

triple_overlap <- Reduce(intersect, list(tfr_treg, tfr_tfh, tfh_treg))

only_tfr_treg_and_tfr_tfh  <- setdiff(intersect(tfr_treg, tfr_tfh),  triple_overlap)
only_tfr_treg_and_tfh_treg <- setdiff(intersect(tfr_treg, tfh_treg), triple_overlap)
only_tfr_tfh_and_tfh_treg  <- setdiff(intersect(tfr_tfh,  tfh_treg), triple_overlap)

only_tfr_treg <- setdiff(tfr_treg, union(tfr_tfh, tfh_treg))
only_tfr_tfh  <- setdiff(tfr_tfh,  union(tfr_treg, tfh_treg))
only_tfh_treg <- setdiff(tfh_treg, union(tfr_treg, tfr_tfh))

wb_overlap <- createWorkbook()

sections <- list(
  "Triple_Overlap"      = triple_overlap,
  "TfrTreg_AND_TfrTfh"  = only_tfr_treg_and_tfr_tfh,
  "TfrTreg_AND_TfhTreg" = only_tfr_treg_and_tfh_treg,
  "TfrTfh_AND_TfhTreg"  = only_tfr_tfh_and_tfh_treg,
  "Exclusive_TfrTreg"   = only_tfr_treg,
  "Exclusive_TfrTfh"    = only_tfr_tfh,
  "Exclusive_TfhTreg"   = only_tfh_treg
)

for (sheet_name in names(sections)) {
  genes <- sections[[sheet_name]]

  if (length(genes) == 0) {
    df <- data.frame(EnsemblID = character(), GeneSymbol = character())
  } else {
    df <- get_symbols(genes, ensembl_ids, gene_symbols)
  }

  addWorksheet(wb_overlap, sheet_name)
  writeData(wb_overlap, sheet = sheet_name, df)

  cat(sheet_name, ":", nrow(df), "genes\n")
}

saveWorkbook(
  wb_overlap,
  file = file.path(tables_dir, "Euler_Overlapping_Genes.xlsx"),
  overwrite = TRUE
)

cat("\nSaved: Euler_Overlapping_Genes.xlsx\n")

# ---------------------------------------------------------------------------
# Venn diagram
# ---------------------------------------------------------------------------
deg_gene_sets_venn <- list(
  "Tfr vs Treg" = tfr_treg,
  "Tfh vs Treg" = tfh_treg,
  "Tfr vs Tfh"  = tfr_tfh
)

p_venn <- ggvenn(
  deg_gene_sets_venn,
  fill_color = c("tomato", "palegreen", "skyblue"),
  stroke_size = 0.8,
  set_name_size = 6,
  text_size = 5,
  show_percentage = FALSE
) +
  labs(title = "Venn Diagram of Differentially Expressed Genes (DEGs)") +
  theme(plot.title = element_text(size = 16, hjust = 0.5))

svglite::svglite(
  file.path(figures_dir, "Venn_DEG_Overlap.svg"),
  width = 8, height = 6, bg = "transparent"
)
print(p_venn)
dev.off()
