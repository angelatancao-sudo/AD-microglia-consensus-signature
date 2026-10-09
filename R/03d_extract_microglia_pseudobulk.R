# ============================================================
# 03d_extract_microglia_pseudobulk.R
# ============================================================
# Project: AD Microglia Multi-Dataset Transcriptomic Analysis
#
# Purpose:
#   Extract validated microglia clusters from the three raw GEO
#   datasets and create donor-level pseudobulk count matrices.
#
# This script uses the validated cluster calls from 03D:
#
#   Primary analysis:
#     GSE157827: cluster 7
#     GSE160936: clusters 1, 2, 6, 9
#     GSE188545: clusters 7, 18
#
#   Sensitivity analysis:
#     GSE157827: cluster 7
#     GSE160936: clusters 1, 2, 6, 9, 27
#     GSE188545: clusters 7, 18, 29
#
# Why this is important:
#   The old pipeline selected cells by simple per-cell marker positivity.
#   This script replaces that with a field-aligned approach:
#     1. cluster raw nuclei
#     2. validate microglia clusters by marker enrichment
#     3. extract validated clusters
#     4. aggregate to donor-level pseudobulk counts
#
# GSE160936 note:
#   GSE160936 contains 24 region-level samples from 12 donors.
#   For the main analysis, this script aggregates EC and SSC together
#   by donor_id, so the DE model treats the dataset as 12 donors,
#   not 24 independent samples.
#
# Inputs:
#   RDS/GSE157827_clustered_seurat.rds
#   RDS/GSE160936_clustered_seurat.rds
#   RDS/GSE188545_clustered_seurat.rds
#
# Outputs:
#   RDS/[dataset]_microglia_PRIMARY_seurat.rds
#   RDS/[dataset]_microglia_PRIMARY_pseudobulk_counts.rds
#   RDS/[dataset]_microglia_PRIMARY_pseudobulk_meta.rds
#   Output/Tables/[dataset]_microglia_PRIMARY_pseudobulk_meta.csv
#   Output/Tables/[dataset]_microglia_PRIMARY_cell_counts_by_donor.csv
#
#   Same outputs are also created for SENSITIVITY.
#
# Run after:
#   03D_validate_microglia_cluster_calls_AD61026_FIXED.R
# ============================================================

# ============================================================
# SECTION 1: Load setup and packages
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
  library(tibble)
  library(readr)
})

cat("============================================================\n")
cat("  03E: Extract validated microglia and pseudobulk\n")
cat("============================================================\n\n")

dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# SECTION 2: Validated cluster calls
# ============================================================
# Keep these hard-coded so the extraction is reproducible and does
# not depend on a manually edited spreadsheet.

CLUSTER_CALLS <- list(
  PRIMARY = list(
    GSE157827 = c("7"),
    GSE160936 = c("1", "2", "6", "9"),
    GSE188545 = c("7", "18")
  ),
  SENSITIVITY = list(
    GSE157827 = c("7"),
    GSE160936 = c("1", "2", "6", "9", "27"),
    GSE188545 = c("7", "18", "29")
  )
)

DATASETS_TO_PROCESS <- c("GSE157827", "GSE160936", "GSE188545")

# ============================================================
# SECTION 3: Helper functions
# ============================================================

# ------------------------------------------------------------
# Function: prepare Seurat v5 object for count extraction
# ------------------------------------------------------------
prepare_for_counts <- function(seu, dataset_id) {

  DefaultAssay(seu) <- "RNA"

  # Merged Seurat v5 objects may store each sample as a separate layer.
  # JoinLayers creates a single counts layer in memory.
  seu <- tryCatch({
    JoinLayers(seu, assay = "RNA")
  }, error = function(e) {
    cat("  JoinLayers note for ", dataset_id, ": ", e$message, "\n", sep = "")
    seu
  })

  # Confirm counts are available.
  counts_ok <- tryCatch({
    mat <- GetAssayData(seu, assay = "RNA", layer = "counts")
    nrow(mat) > 0 && ncol(mat) > 0
  }, error = function(e) {
    FALSE
  })

  if (!counts_ok) {
    stop("Could not access RNA counts layer for ", dataset_id)
  }

  return(seu)
}

# ------------------------------------------------------------
# Function: create donor-level pseudobulk counts
# ------------------------------------------------------------
make_pseudobulk_by_donor <- function(seu, dataset_id, analysis_label) {

  DefaultAssay(seu) <- "RNA"

  counts <- GetAssayData(seu, assay = "RNA", layer = "counts")

  meta <- seu@meta.data %>%
    as.data.frame() %>%
    rownames_to_column("cell_id")

  # Make sure donor_id exists.
  if (!"donor_id" %in% colnames(meta)) {
    stop("donor_id column missing from metadata for ", dataset_id)
  }

  # Use donor_id as the biological replicate.
  donor_ids <- sort(unique(meta$donor_id))

  cat("  Building pseudobulk matrix for ", dataset_id,
      " / ", analysis_label, "\n", sep = "")
  cat("    Donors:", length(donor_ids), "\n")

  pb_list <- list()
  pb_meta_list <- list()
  cell_count_list <- list()

  for (donor in donor_ids) {

    donor_cells <- meta$cell_id[meta$donor_id == donor]

    if (length(donor_cells) == 0) next

    # Sum raw counts across all selected microglia nuclei from the donor.
    # This preserves integer count structure for edgeR/limma-voom.
    donor_counts <- Matrix::rowSums(counts[, donor_cells, drop = FALSE])
    pb_list[[donor]] <- donor_counts

    donor_meta <- meta %>%
      filter(donor_id == donor)

    # Diagnosis should be constant within donor.
    diagnosis_values <- unique(as.character(donor_meta$diagnosis_std))
    diagnosis_values <- diagnosis_values[!is.na(diagnosis_values)]

    if (length(diagnosis_values) != 1) {
      warning("Donor ", donor, " in ", dataset_id,
              " has more than one diagnosis label: ",
              paste(diagnosis_values, collapse = ", "))
    }

    # For GSE160936, each donor has EC and SSC region samples.
    # We combine regions into one donor-level pseudobulk sample,
    # but keep a record of which regions contributed.
    region_values <- if ("region" %in% colnames(donor_meta)) {
      paste(sort(unique(as.character(donor_meta$region))), collapse = ";")
    } else {
      NA_character_
    }

    sample_values <- paste(sort(unique(as.character(donor_meta$sample_id))), collapse = ";")

    pb_meta_list[[donor]] <- data.frame(
      dataset = dataset_id,
      analysis = analysis_label,
      pseudobulk_id = paste(dataset_id, donor, analysis_label, sep = "_"),
      donor_id = donor,
      diagnosis_std = diagnosis_values[1],
      n_microglia_cells = nrow(donor_meta),
      n_source_samples = dplyr::n_distinct(donor_meta$sample_id),
      source_samples = sample_values,
      source_regions = region_values,
      stringsAsFactors = FALSE
    )

    # Cluster composition within each donor.
    tmp_counts <- donor_meta %>%
      mutate(cluster_id = as.character(seurat_clusters)) %>%
      count(dataset = dataset_id,
            analysis = analysis_label,
            donor_id,
            diagnosis_std,
            cluster_id,
            name = "n_cells")

    cell_count_list[[donor]] <- tmp_counts
  }

  # Combine donor count vectors into a genes x donors matrix.
  pb_counts <- do.call(cbind, pb_list)
  colnames(pb_counts) <- names(pb_list)

  # Convert dense matrix from rowSums to sparse matrix for storage efficiency.
  pb_counts <- Matrix(pb_counts, sparse = TRUE)

  pb_meta <- bind_rows(pb_meta_list)
  cell_counts <- bind_rows(cell_count_list)

  # Align metadata to count matrix columns.
  pb_meta <- pb_meta %>%
    mutate(donor_id = as.character(donor_id)) %>%
    arrange(match(donor_id, colnames(pb_counts)))

  if (!all(pb_meta$donor_id == colnames(pb_counts))) {
    stop("Pseudobulk metadata and count matrix columns are not aligned for ",
         dataset_id, " / ", analysis_label)
  }

  # Use globally unique pseudobulk IDs as column names.
  colnames(pb_counts) <- pb_meta$pseudobulk_id
  rownames(pb_meta) <- pb_meta$pseudobulk_id

  return(list(
    counts = pb_counts,
    meta = pb_meta,
    cell_counts = cell_counts
  ))
}

# ------------------------------------------------------------
# Function: extract selected clusters and save outputs
# ------------------------------------------------------------
extract_and_save_microglia <- function(dataset_id, analysis_label, selected_clusters) {

  cat("============================================================\n")
  cat("Processing ", dataset_id, " / ", analysis_label, "\n", sep = "")
  cat("Selected clusters: ", paste(selected_clusters, collapse = ", "), "\n", sep = "")
  cat("============================================================\n")

  rds_path <- file.path(RDS_DIR, paste0(dataset_id, "_clustered_seurat.rds"))

  if (!file.exists(rds_path)) {
    stop("Missing clustered Seurat object: ", rds_path)
  }

  seu <- readRDS(rds_path)
  seu <- prepare_for_counts(seu, dataset_id)

  seu$cluster_id <- as.character(seu$seurat_clusters)

  n_before <- ncol(seu)

  keep_cells <- colnames(seu)[seu$cluster_id %in% selected_clusters]

  if (length(keep_cells) == 0) {
    stop("No cells found for selected clusters in ", dataset_id,
         " / ", analysis_label)
  }

  seu_micro <- subset(seu, cells = keep_cells)

  n_after <- ncol(seu_micro)

  cat("  Cells before extraction:", n_before, "\n")
  cat("  Cells after extraction :", n_after, "\n")
  cat("  Donors represented    :", length(unique(seu_micro$donor_id)), "\n")
  cat("  Diagnosis table:\n")
  print(table(seu_micro$diagnosis_std, useNA = "ifany"))
  cat("\n")

  pb <- make_pseudobulk_by_donor(
    seu = seu_micro,
    dataset_id = dataset_id,
    analysis_label = analysis_label
  )

  cat("  Pseudobulk counts dimension:",
      nrow(pb$counts), "genes x", ncol(pb$counts), "donors\n\n")

  # Save selected microglia object.
  saveRDS(
    seu_micro,
    file.path(RDS_DIR, paste0(dataset_id, "_microglia_", analysis_label, "_seurat.rds"))
  )

  # Save pseudobulk counts and metadata.
  saveRDS(
    pb$counts,
    file.path(RDS_DIR, paste0(dataset_id, "_microglia_", analysis_label, "_pseudobulk_counts.rds"))
  )

  saveRDS(
    pb$meta,
    file.path(RDS_DIR, paste0(dataset_id, "_microglia_", analysis_label, "_pseudobulk_meta.rds"))
  )

  write.csv(
    pb$meta,
    file.path(TAB_DIR, paste0(dataset_id, "_microglia_", analysis_label, "_pseudobulk_meta.csv")),
    row.names = FALSE
  )

  write.csv(
    pb$cell_counts,
    file.path(TAB_DIR, paste0(dataset_id, "_microglia_", analysis_label, "_cell_counts_by_donor_cluster.csv")),
    row.names = FALSE
  )

  # Save a compact summary for quick inspection.
  summary_df <- data.frame(
    dataset = dataset_id,
    analysis = analysis_label,
    selected_clusters = paste(selected_clusters, collapse = ";"),
    n_cells_before_extraction = n_before,
    n_microglia_cells = n_after,
    n_pseudobulk_donors = ncol(pb$counts),
    n_genes = nrow(pb$counts),
    n_AD_donors = sum(pb$meta$diagnosis_std == "AD"),
    n_Control_donors = sum(pb$meta$diagnosis_std == "Control"),
    stringsAsFactors = FALSE
  )

  write.csv(
    summary_df,
    file.path(TAB_DIR, paste0(dataset_id, "_microglia_", analysis_label, "_summary.csv")),
    row.names = FALSE
  )

  rm(seu, seu_micro, pb)
  gc(verbose = FALSE)

  return(summary_df)
}

# ============================================================
# SECTION 4: Run primary and sensitivity extraction
# ============================================================

all_summaries <- list()

for (analysis_label in names(CLUSTER_CALLS)) {

  cat("\n############################################################\n")
  cat("Running analysis set:", analysis_label, "\n")
  cat("############################################################\n\n")

  for (dataset_id in DATASETS_TO_PROCESS) {

    selected_clusters <- CLUSTER_CALLS[[analysis_label]][[dataset_id]]

    summary_df <- extract_and_save_microglia(
      dataset_id = dataset_id,
      analysis_label = analysis_label,
      selected_clusters = selected_clusters
    )

    all_summaries[[paste(dataset_id, analysis_label, sep = "_")]] <- summary_df
  }
}

combined_summary <- bind_rows(all_summaries)

write.csv(
  combined_summary,
  file.path(TAB_DIR, "03E_microglia_extraction_pseudobulk_summary.csv"),
  row.names = FALSE
)

saveRDS(
  combined_summary,
  file.path(RDS_DIR, "03E_microglia_extraction_pseudobulk_summary.rds")
)

cat("============================================================\n")
cat("03E extraction complete.\n\n")
cat("Summary saved to:\n")
cat("  ", file.path(TAB_DIR, "03E_microglia_extraction_pseudobulk_summary.csv"), "\n\n")
cat("Next step:\n")
cat("  Review the summary and pseudobulk metadata.\n")
cat("  Then run the DE/meta-analysis pipeline using PRIMARY first,\n")
cat("  followed by SENSITIVITY as robustness analysis.\n")
cat("============================================================\n")
