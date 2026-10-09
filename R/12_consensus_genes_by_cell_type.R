# ============================================================
# 12_consensus_genes_by_cell_type.R
# ============================================================
# Purpose:
#   Several consensus genes are usually associated with other cell
#   types (FLT1 endothelial; DOC2A, NEXMIF, BDNF-AS, TSHZ3 neuronal).
#   This script measures how much each of the 48 consensus genes is
#   expressed in microglia compared with other brain cell types, using
#   the whole-tissue clustered objects already built for 03D.
#
# Method:
#   1. Load each whole-tissue clustered Seurat object (all cell types).
#   2. Label every cluster with the cell type whose marker panel (the
#      03D panels) is expressed in the highest percentage of its nuclei.
#      FLT1 is removed from the endothelial panel for labelling, because
#      FLT1 is one of the genes being tested. The microglial clusters
#      used in the paper are labelled "microglia".
#   3. For each consensus gene and cell type: percent of nuclei
#      expressing it and mean normalized expression.
#   4. Flag genes expressed in a higher percentage of nuclei in another
#      cell type than in microglia.
#
# Inputs:  RDS/[GSE157827|GSE160936|GSE188545]_clustered_seurat.rds
#          Final tables and figures/Supplementary_Table_2_DE_Consensus_Signatures.xlsx
# Outputs (Output/CellType_check/):
#   cluster_labels.csv
#   consensus48_expression_by_cell_type.csv
#   consensus48_cell_type_flags.csv
#   summary.txt
#
# HOW TO RUN: open in RStudio and click Source. Loading the three
# whole-tissue objects is slow (10-30 minutes, needs plenty of RAM).
# Results are written to Output/CellType_check/.
# ============================================================

# ------------------------------------------------------------
# Project folder. Set the AD_PROJECT_ROOT environment variable to the
# folder that contains Data/, RDS/, Output/ and Supplementary_Tables/,
# or start R with that folder as the working directory.
# ------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("AD_PROJECT_ROOT", unset = getwd())

ROOT_DIR <- PROJECT_ROOT
supp2_file <- file.path(PROJECT_ROOT, "Supplementary_Tables", "Supplementary_Table_2_DE_Consensus_Signatures.xlsx")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readxl)
})

rds_dir <- file.path(ROOT_DIR, "RDS")
out_dir <- file.path(ROOT_DIR, "Output", "CellType_check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_lines <- c()
say <- function(...) { s <- paste0(...); cat(s, "\n"); log_lines <<- c(log_lines, s) }

genes48 <- read_excel(supp2_file, sheet = "Consensus48")$gene
stopifnot(length(genes48) == 48)

# Marker panels from 03D (FLT1 removed from endothelial for labelling)
MARKERS <- list(
  microglia = c("P2RY12", "CX3CR1", "CSF1R", "TMEM119", "AIF1", "CTSS", "C3", "APBB1IP", "SALL1"),
  neuron = c("RBFOX3", "SNAP25", "SLC17A7", "GAD1", "GAD2"),
  astrocyte = c("AQP4", "GFAP", "SLC1A2", "ALDH1L1"),
  oligodendrocyte = c("MBP", "MOG", "PLP1", "MOBP"),
  opc = c("PDGFRA", "VCAN", "CSPG4"),
  endothelial_pericyte = c("CLDN5", "PECAM1", "COL1A1", "DCN"),
  macrophage_pvm = c("CD163", "MRC1", "MSR1", "LYVE1")
)
MICROGLIA_CLUSTERS <- list(
  GSE157827 = c("7"),
  GSE160936 = c("1", "2", "6", "9"),
  GSE188545 = c("7", "18")
)

all_expr <- list(); all_labels <- list()

for (ds in names(MICROGLIA_CLUSTERS)) {
  say("\n=== ", ds, " ===")
  seu <- readRDS(file.path(rds_dir, paste0(ds, "_clustered_seurat.rds")))
  DefaultAssay(seu) <- "RNA"
  seu <- tryCatch(JoinLayers(seu, assay = "RNA"), error = function(e) seu)
  expr <- tryCatch(GetAssayData(seu, assay = "RNA", layer = "data"), error = function(e) NULL)
  if (is.null(expr) || max(expr[seq_len(min(100, nrow(expr))), seq_len(min(100, ncol(expr)))]) == 0) {
    say("  Normalizing (no data layer found)")
    seu <- NormalizeData(seu, verbose = FALSE)
    expr <- GetAssayData(seu, assay = "RNA", layer = "data")
  }
  keep_genes <- intersect(c(genes48, unlist(MARKERS)), rownames(expr))
  expr <- expr[keep_genes, , drop = FALSE]
  cl <- as.character(seu$seurat_clusters); names(cl) <- colnames(seu)
  rm(seu); gc(verbose = FALSE)
  say("  Nuclei: ", ncol(expr), "  clusters: ", length(unique(cl)))

  # Label clusters by marker panel
  labels <- bind_rows(lapply(sort(unique(cl)), function(k) {
    cells <- names(cl)[cl == k]
    sc <- sapply(MARKERS, function(m) {
      m <- intersect(m, rownames(expr))
      if (!length(m)) return(NA)
      mean(Matrix::rowMeans(expr[m, cells, drop = FALSE] > 0)) * 100
    })
    best <- names(which.max(sc))
    lab <- if (k %in% MICROGLIA_CLUSTERS[[ds]]) "microglia"
           else if (max(sc, na.rm = TRUE) < 20) "unassigned"
           else if (best == "microglia") "microglia_like_not_used"
           else best
    data.frame(dataset = ds, cluster_id = k, n_nuclei = length(cells), cell_type = lab,
               t(round(sc, 1)), check.names = FALSE)
  }))
  all_labels[[ds]] <- labels
  say("  Nuclei by cell type: ",
      paste(tapply(labels$n_nuclei, labels$cell_type, sum) %>% { paste(names(.), ., sep = "=") }, collapse = ", "))

  ct <- labels$cell_type[match(cl, labels$cluster_id)]
  g_here <- intersect(genes48, rownames(expr))
  res <- bind_rows(lapply(setdiff(unique(ct), "unassigned"), function(t) {
    cells <- names(cl)[ct == t]
    sub <- expr[g_here, cells, drop = FALSE]
    data.frame(dataset = ds, gene = g_here, cell_type = t, n_nuclei = length(cells),
               pct_expressing = as.numeric(Matrix::rowMeans(sub > 0) * 100),
               mean_expression = as.numeric(Matrix::rowMeans(sub)))
  }))
  all_expr[[ds]] <- res
  rm(expr); gc(verbose = FALSE)
}

labels_df <- bind_rows(all_labels)
expr_df <- bind_rows(all_expr)
write.csv(labels_df, file.path(out_dir, "cluster_labels.csv"), row.names = FALSE)
write.csv(expr_df, file.path(out_dir, "consensus48_expression_by_cell_type.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Flags: compare microglia with the highest other major cell type
# ------------------------------------------------------------
major <- c("neuron", "astrocyte", "oligodendrocyte", "opc", "endothelial_pericyte", "macrophage_pvm")
per_ds <- expr_df %>%
  group_by(dataset, gene) %>%
  summarise(
    microglia_pct = pct_expressing[cell_type == "microglia"][1],
    top_other_type = cell_type[cell_type %in% major][which.max(pct_expressing[cell_type %in% major])],
    top_other_pct = max(pct_expressing[cell_type %in% major]),
    .groups = "drop") %>%
  mutate(ratio_microglia_to_top_other = microglia_pct / pmax(top_other_pct, 0.1),
         higher_elsewhere = top_other_pct > microglia_pct)

flags <- per_ds %>%
  group_by(gene) %>%
  summarise(
    n_datasets = n(),
    n_higher_elsewhere = sum(higher_elsewhere),
    median_microglia_pct = round(median(microglia_pct), 1),
    median_top_other_pct = round(median(top_other_pct), 1),
    most_common_top_other = names(sort(table(top_other_type), decreasing = TRUE))[1],
    median_ratio = round(median(ratio_microglia_to_top_other), 2),
    .groups = "drop") %>%
  mutate(category = case_when(
    n_higher_elsewhere == 0 & median_ratio >= 2 ~ "Microglia-enriched",
    n_higher_elsewhere == 0 ~ "Highest in microglia",
    n_higher_elsewhere < n_datasets / 2 ~ "Shared",
    TRUE ~ paste0("Higher in ", most_common_top_other)
  )) %>%
  arrange(desc(n_higher_elsewhere), median_ratio)
write.csv(flags, file.path(out_dir, "consensus48_cell_type_flags.csv"), row.names = FALSE)

say("\n=== Summary across datasets ===")
say("Category counts:")
for (k in names(table(flags$category))) say("  ", k, ": ", table(flags$category)[[k]])
say("\nGenes expressed in more nuclei of another cell type than microglia (>= half of datasets):")
fl <- flags %>% filter(n_higher_elsewhere >= n_datasets / 2)
for (i in seq_len(nrow(fl))) say("  ", fl$gene[i], ": microglia ", fl$median_microglia_pct[i], "% vs ",
                                 fl$most_common_top_other[i], " ", fl$median_top_other_pct[i], "%")
writeLines(log_lines, file.path(out_dir, "summary.txt"))
cat("\nDone. Results written to:", out_dir, "\n")
