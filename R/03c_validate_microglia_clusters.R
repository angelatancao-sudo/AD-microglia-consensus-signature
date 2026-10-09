# ============================================================
# 03c_validate_microglia_clusters.R
# ============================================================
# Project: AD Microglia Multi-Dataset Transcriptomic Analysis
#
# Purpose:
#   Validate candidate microglia cluster annotations before
#   extracting microglia and creating donor-level pseudobulk counts.
#
# Fix in this version:
#   The previous version failed in Seurat v5 because merged objects
#   can contain multiple assay layers. GetAssayData(slot=...) is now
#   defunct in SeuratObject 5. This version uses JoinLayers() first
#   and then reads the normalized "data" layer using layer = "data".
#
# This script does NOT change the saved clustered Seurat objects.
# It only loads each object, joins layers in memory, and writes
# validation tables.
# ============================================================

# ============================================================
# SECTION 1: Setup
# ============================================================

# ------------------------------------------------------------
# Project folder. Set the AD_PROJECT_ROOT environment variable to the
# folder that contains Data/, RDS/, Output/ and Supplementary_Tables/,
# or start R with that folder as the working directory.
# ------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("AD_PROJECT_ROOT", unset = getwd())

setup_path <- file.path(PROJECT_ROOT, "RDS/project_setup.rds")

if (!file.exists(setup_path)) {
  stop("project_setup.rds not found. Run 01_Setup_AD61026.R first.")
}

setup <- readRDS(setup_path)
list2env(setup, envir = .GlobalEnv)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

cat("============================================================\n")
cat("  03D FIXED: Validate candidate microglia cluster calls\n")
cat("============================================================\n\n")

VALIDATION_DIR <- file.path(TAB_DIR, "Microglia_validation")
dir.create(VALIDATION_DIR, recursive = TRUE, showWarnings = FALSE)

# Set to FALSE if you want the fast marker-expression tables only.
# TRUE gives stronger annotation evidence, but takes longer.
RUN_TOP_MARKERS <- TRUE

# ============================================================
# SECTION 2: Candidate cluster calls for review
# ============================================================

PRIMARY_CLUSTERS <- list(
  GSE157827 = c("7"),
  GSE160936 = c("1", "2", "6", "9"),
  GSE188545 = c("7", "18")
)

EXPANDED_CLUSTERS <- list(
  GSE157827 = c("7", "28", "29"),
  GSE160936 = c("1", "2", "6", "9", "27"),
  GSE188545 = c("7", "18", "23", "29", "32")
)

DATASET_IDS <- names(PRIMARY_CLUSTERS)

# ============================================================
# SECTION 3: Marker panels
# ============================================================

MARKERS <- list(
  microglia = c("P2RY12", "CX3CR1", "CSF1R", "TMEM119",
                "AIF1", "CTSS", "C3", "APBB1IP", "SALL1"),
  neuron = c("RBFOX3", "SNAP25", "SLC17A7", "GAD1", "GAD2"),
  astrocyte = c("AQP4", "GFAP", "SLC1A2", "ALDH1L1"),
  oligodendrocyte = c("MBP", "MOG", "PLP1", "MOBP"),
  opc = c("PDGFRA", "VCAN", "CSPG4"),
  endothelial_pericyte = c("CLDN5", "FLT1", "PECAM1",
                           "COL1A1", "DCN"),
  macrophage_pvm = c("CD163", "MRC1", "MSR1", "LYVE1")
)

ALL_MARKER_GENES <- unique(unlist(MARKERS))

# ============================================================
# SECTION 4: Helper functions
# ============================================================

# ------------------------------------------------------------
# Function: prepare_seurat_for_validation
# Purpose:
#   Handle Seurat v5 assay layers safely.
# ------------------------------------------------------------
prepare_seurat_for_validation <- function(seu, dataset_id) {

  DefaultAssay(seu) <- "RNA"

  cat("  Preparing Seurat object for validation:", dataset_id, "\n")

  # JoinLayers is needed for merged Seurat v5 objects where each sample
  # may be stored as a separate layer. This is done only in memory.
  seu <- tryCatch({
    JoinLayers(seu, assay = "RNA")
  }, error = function(e) {
    cat("  JoinLayers note:", e$message, "\n")
    seu
  })

  # Confirm that normalized data are available.
  data_ok <- tryCatch({
    mat <- GetAssayData(seu, assay = "RNA", layer = "data")
    nrow(mat) > 0 && ncol(mat) > 0
  }, error = function(e) {
    FALSE
  })

  # If the data layer is missing, recreate it by log-normalization.
  # This should usually not be needed because Stage 1 already ran NormalizeData.
  if (!data_ok) {
    cat("  Normalized data layer not found. Running NormalizeData in memory.\n")
    seu <- NormalizeData(
      seu,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
  }

  return(seu)
}

# ------------------------------------------------------------
# Function: get normalized expression matrix
# ------------------------------------------------------------
get_norm_data <- function(seu) {
  DefaultAssay(seu) <- "RNA"

  mat <- GetAssayData(seu, assay = "RNA", layer = "data")

  return(mat)
}

# ------------------------------------------------------------
# Function: summarize marker expression by cluster
# ------------------------------------------------------------
summarize_marker_expression <- function(seu, dataset_id) {

  expr <- get_norm_data(seu)
  genes_present <- intersect(ALL_MARKER_GENES, rownames(expr))

  if (length(genes_present) == 0) {
    stop("No marker genes found in expression matrix for ", dataset_id)
  }

  meta <- seu@meta.data %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    mutate(cluster_id = as.character(seurat_clusters))

  cluster_ids <- sort(unique(meta$cluster_id))

  out_list <- list()

  for (cl in cluster_ids) {

    cells <- meta$cell_id[meta$cluster_id == cl]

    if (length(cells) == 0) next

    submat <- expr[genes_present, cells, drop = FALSE]

    avg_expr <- Matrix::rowMeans(submat)
    pct_expr <- Matrix::rowMeans(submat > 0) * 100

    tmp <- data.frame(
      dataset = dataset_id,
      cluster_id = cl,
      gene = genes_present,
      avg_expression = as.numeric(avg_expr),
      pct_expressed = as.numeric(pct_expr),
      stringsAsFactors = FALSE
    )

    out_list[[cl]] <- tmp
  }

  marker_long <- bind_rows(out_list)

  gene_to_panel <- data.frame(
    gene = unlist(MARKERS),
    panel = rep(names(MARKERS), lengths(MARKERS)),
    stringsAsFactors = FALSE
  )

  marker_long <- marker_long %>%
    left_join(gene_to_panel, by = "gene") %>%
    select(dataset, cluster_id, panel, gene, avg_expression, pct_expressed)

  return(marker_long)
}

# ------------------------------------------------------------
# Function: summarize composition and module scores
# ------------------------------------------------------------
summarize_cluster_composition <- function(seu, dataset_id) {

  meta <- seu@meta.data %>%
    as.data.frame() %>%
    mutate(cluster_id = as.character(seurat_clusters))

  score_cols <- grep("_score$", colnames(meta), value = TRUE)

  composition <- meta %>%
    group_by(cluster_id) %>%
    summarise(
      dataset = dataset_id,
      n_cells = n(),
      n_donors = n_distinct(donor_id),
      n_samples = n_distinct(sample_id),
      n_AD = sum(diagnosis_std == "AD", na.rm = TRUE),
      n_Control = sum(diagnosis_std == "Control", na.rm = TRUE),
      pct_AD = round(100 * n_AD / n_cells, 2),
      pct_Control = round(100 * n_Control / n_cells, 2),
      across(all_of(score_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    relocate(dataset, .before = cluster_id) %>%
    arrange(as.numeric(cluster_id))

  return(composition)
}

# ------------------------------------------------------------
# Function: build candidate validation summary
# ------------------------------------------------------------
build_candidate_summary <- function(composition, marker_long, dataset_id) {

  primary <- PRIMARY_CLUSTERS[[dataset_id]]
  expanded <- EXPANDED_CLUSTERS[[dataset_id]]

  panel_summary <- marker_long %>%
    group_by(dataset, cluster_id, panel) %>%
    summarise(
      mean_avg_expression = mean(avg_expression, na.rm = TRUE),
      mean_pct_expressed = mean(pct_expressed, na.rm = TRUE),
      max_avg_expression = max(avg_expression, na.rm = TRUE),
      max_pct_expressed = max(pct_expressed, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = panel,
      values_from = c(mean_avg_expression, mean_pct_expressed,
                      max_avg_expression, max_pct_expressed),
      names_sep = "__"
    )

  micro_detail <- marker_long %>%
    filter(panel == "microglia") %>%
    mutate(detail = paste0(gene, ":avg=", round(avg_expression, 3),
                           ",pct=", round(pct_expressed, 1))) %>%
    group_by(dataset, cluster_id) %>%
    summarise(
      microglia_marker_detail = paste(detail, collapse = "; "),
      .groups = "drop"
    )

  summary <- composition %>%
    left_join(panel_summary, by = c("dataset", "cluster_id")) %>%
    left_join(micro_detail, by = c("dataset", "cluster_id")) %>%
    mutate(
      primary_candidate = cluster_id %in% primary,
      expanded_candidate = cluster_id %in% expanded,
      proposed_use = case_when(
        cluster_id %in% primary ~ "Primary",
        cluster_id %in% expanded ~ "Expanded_sensitivity_only",
        TRUE ~ "Exclude"
      )
    ) %>%
    arrange(as.numeric(cluster_id))

  return(summary)
}

# ------------------------------------------------------------
# Function: run top marker detection
# ------------------------------------------------------------
run_top_markers <- function(seu, dataset_id) {

  cat("  Running FindAllMarkers for", dataset_id, "\n")
  cat("  This is for cluster annotation review only.\n")

  DefaultAssay(seu) <- "RNA"
  Idents(seu) <- "seurat_clusters"

  set.seed(42)
  markers <- FindAllMarkers(
    object = seu,
    only.pos = TRUE,
    min.pct = 0.10,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    max.cells.per.ident = 1000,
    verbose = FALSE
  )

  markers <- markers %>%
    as.data.frame() %>%
    mutate(dataset = dataset_id) %>%
    relocate(dataset)

  fc_col <- if ("avg_log2FC" %in% colnames(markers)) {
    "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(markers)) {
    "avg_logFC"
  } else {
    NULL
  }

  if (!is.null(fc_col)) {
    markers_top <- markers %>%
      group_by(cluster) %>%
      arrange(p_val_adj, desc(.data[[fc_col]]), .by_group = TRUE) %>%
      slice_head(n = 30) %>%
      ungroup()
  } else {
    markers_top <- markers %>%
      group_by(cluster) %>%
      arrange(p_val_adj, .by_group = TRUE) %>%
      slice_head(n = 30) %>%
      ungroup()
  }

  return(markers_top)
}

# ============================================================
# SECTION 5: Run validation for each dataset
# ============================================================

all_candidate_summaries <- list()

for (dataset_id in DATASET_IDS) {

  cat("============================================================\n")
  cat("Validating", dataset_id, "\n")
  cat("============================================================\n")

  rds_path <- file.path(RDS_DIR, paste0(dataset_id, "_clustered_seurat.rds"))

  if (!file.exists(rds_path)) {
    stop("Missing clustered Seurat object: ", rds_path)
  }

  seu <- readRDS(rds_path)
  seu <- prepare_seurat_for_validation(seu, dataset_id)

  seu$cluster_id <- as.character(seu$seurat_clusters)

  composition <- summarize_cluster_composition(seu, dataset_id)
  marker_long <- summarize_marker_expression(seu, dataset_id)
  candidate_summary <- build_candidate_summary(composition, marker_long, dataset_id)

  write.csv(
    composition,
    file.path(VALIDATION_DIR, paste0(dataset_id, "_cluster_composition.csv")),
    row.names = FALSE
  )

  write.csv(
    marker_long,
    file.path(VALIDATION_DIR, paste0(dataset_id, "_marker_expression_long.csv")),
    row.names = FALSE
  )

  write.csv(
    candidate_summary,
    file.path(VALIDATION_DIR, paste0(dataset_id, "_candidate_validation_summary.csv")),
    row.names = FALSE
  )

  if (RUN_TOP_MARKERS) {
    top_markers <- run_top_markers(seu, dataset_id)

    write.csv(
      top_markers,
      file.path(VALIDATION_DIR, paste0(dataset_id, "_top_cluster_markers.csv")),
      row.names = FALSE
    )

    candidate_clusters <- EXPANDED_CLUSTERS[[dataset_id]]
    top_candidate_markers <- top_markers %>%
      filter(as.character(cluster) %in% candidate_clusters)

    write.csv(
      top_candidate_markers,
      file.path(VALIDATION_DIR, paste0(dataset_id, "_top_candidate_cluster_markers.csv")),
      row.names = FALSE
    )
  }

  all_candidate_summaries[[dataset_id]] <- candidate_summary

  rm(seu, composition, marker_long, candidate_summary)
  if (exists("top_markers")) rm(top_markers)
  if (exists("top_candidate_markers")) rm(top_candidate_markers)
  gc(verbose = FALSE)
}

all_summary <- bind_rows(all_candidate_summaries)

write.csv(
  all_summary,
  file.path(VALIDATION_DIR, "all_datasets_candidate_validation_summary.csv"),
  row.names = FALSE
)

cluster_call_table <- bind_rows(lapply(DATASET_IDS, function(dataset_id) {
  clusters <- sort(unique(EXPANDED_CLUSTERS[[dataset_id]]))
  data.frame(
    dataset = dataset_id,
    cluster_id = clusters,
    primary_candidate = clusters %in% PRIMARY_CLUSTERS[[dataset_id]],
    expanded_candidate = TRUE,
    proposed_use = ifelse(
      clusters %in% PRIMARY_CLUSTERS[[dataset_id]],
      "Primary",
      "Expanded_sensitivity_only"
    ),
    stringsAsFactors = FALSE
  )
}))

write.csv(
  cluster_call_table,
  file.path(VALIDATION_DIR, "proposed_microglia_cluster_calls.csv"),
  row.names = FALSE
)

cat("============================================================\n")
cat("Validation complete.\n\n")
cat("Review these files:\n")
cat("  ", file.path(VALIDATION_DIR, "all_datasets_candidate_validation_summary.csv"), "\n")
cat("  ", file.path(VALIDATION_DIR, "proposed_microglia_cluster_calls.csv"), "\n")
cat("  plus dataset-specific top marker files if RUN_TOP_MARKERS = TRUE.\n\n")
cat("Next:\n")
cat("  Upload all_datasets_candidate_validation_summary.csv and the\n")
cat("  top_candidate_cluster_markers files for review.\n")
cat("============================================================\n")
