# ============================================================
# 03e_build_preannotated_GEO_pseudobulk.R
# ============================================================
# Project: AD Microglia Multi-Dataset Transcriptomic Analysis
#
# Purpose:
#   Build donor-level pseudobulk inputs for the two pre-annotated
#   GEO datasets that were missing in the 03F check:
#
#     GSE174367  author cell-type labels, microglia label = MG
#     GSE243292  microglia-specific/pre-annotated H5AD
#
# This script does NOT touch the three raw datasets that were rebuilt
# by cluster-based annotation:
#
#     GSE157827
#     GSE160936
#     GSE188545
#
# Why this script is separate:
#   The raw datasets needed cluster-based microglia annotation.
#   These two datasets already have author-provided microglia information,
#   so they can be converted directly to donor-level pseudobulk counts.
#
# Outputs:
#   RDS/GSE174367_pseudobulk_counts.rds
#   RDS/GSE174367_pseudobulk_meta.rds
#   Output/Tables/GSE174367_pseudobulk_meta.csv
#   Output/Tables/GSE174367_cell_counts_by_donor.csv
#
#   RDS/GSE243292_pseudobulk_counts.rds
#   RDS/GSE243292_pseudobulk_meta.rds
#   Output/Tables/GSE243292_pseudobulk_meta.csv
#   Output/Tables/GSE243292_cell_counts_by_donor.csv
#
# Run after:
#   01_Setup_AD61026.R
#   02_load_SEAAD_AD61026.R
#   03E_extract_validated_microglia_pseudobulk_AD61026.R
#
# Then rerun:
#   03F_check_pseudobulk_inputs_before_DE_AD61026.R
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
  library(readr)
  library(tibble)
  library(hdf5r)
  library(R.utils)
})

cat("============================================================\n")
cat("  03G: Build pre-annotated GEO pseudobulk inputs\n")
cat("============================================================\n\n")

dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR,  recursive = TRUE, showWarnings = FALSE)

# ============================================================
# SECTION 2: Shared pseudobulk helper
# ============================================================

make_donor_pseudobulk <- function(counts,
                                  meta,
                                  dataset_id,
                                  donor_col,
                                  diagnosis_col,
                                  ad_labels,
                                  control_labels,
                                  celltype_col = NULL,
                                  microglia_labels = NULL,
                                  analysis_label = "PRIMARY",
                                  min_cells = PARAMS$min_cells_per_donor) {

  cat("------------------------------------------------------------\n")
  cat("Creating pseudobulk for ", dataset_id, "\n", sep = "")
  cat("------------------------------------------------------------\n")

  # Confirm input alignment.
  if (ncol(counts) != nrow(meta)) {
    stop(dataset_id, ": counts columns and metadata rows do not match.")
  }

  if (is.null(rownames(meta))) {
    rownames(meta) <- colnames(counts)
  }

  # Optional author cell-type filter.
  if (!is.null(celltype_col)) {
    if (!celltype_col %in% colnames(meta)) {
      stop(dataset_id, ": cell type column not found: ", celltype_col)
    }

    if (is.null(microglia_labels)) {
      stop(dataset_id, ": microglia_labels must be provided when celltype_col is used.")
    }

    keep_mg <- as.character(meta[[celltype_col]]) %in% microglia_labels

    cat("  Cells before microglia filter:", ncol(counts), "\n")
    cat("  Cells passing microglia filter:", sum(keep_mg), "\n")

    counts <- counts[, keep_mg, drop = FALSE]
    meta <- meta[keep_mg, , drop = FALSE]
  } else {
    cat("  No cell-type filter used because dataset is already microglia-specific.\n")
    cat("  Cells:", ncol(counts), "\n")
  }

  # Diagnosis filter.
  dx <- as.character(meta[[diagnosis_col]])
  keep_dx <- dx %in% c(ad_labels, control_labels)

  cat("  Cells before diagnosis filter:", ncol(counts), "\n")
  cat("  Cells after diagnosis filter :", sum(keep_dx), "\n")

  counts <- counts[, keep_dx, drop = FALSE]
  meta <- meta[keep_dx, , drop = FALSE]

  meta$diagnosis_std <- ifelse(
    as.character(meta[[diagnosis_col]]) %in% ad_labels,
    "AD",
    "Control"
  )
  meta$diagnosis_std <- factor(meta$diagnosis_std, levels = c("Control", "AD"))

  cat("  Diagnosis table after filtering:\n")
  print(table(meta$diagnosis_std, useNA = "ifany"))

  # Donor-level cell count QC.
  donor_cell_counts <- table(as.character(meta[[donor_col]]))
  donors_keep <- names(donor_cell_counts)[donor_cell_counts >= min_cells]
  donors_fail <- names(donor_cell_counts)[donor_cell_counts < min_cells]

  cat("  Donors before min-cell filter:", length(donor_cell_counts), "\n")
  cat("  Donors passing >=", min_cells, "cells:", length(donors_keep), "\n")

  if (length(donors_fail) > 0) {
    cat("  Donors removed for low microglia count:",
        paste(donors_fail, collapse = ", "), "\n")
  }

  keep_donors <- as.character(meta[[donor_col]]) %in% donors_keep
  counts <- counts[, keep_donors, drop = FALSE]
  meta <- meta[keep_donors, , drop = FALSE]

  # Build donor-level count matrix by summing raw counts per donor.
  donors <- sort(unique(as.character(meta[[donor_col]])))

  pb_list <- list()
  meta_list <- list()
  cell_count_list <- list()

  for (donor in donors) {

    donor_cells <- rownames(meta)[as.character(meta[[donor_col]]) == donor]

    donor_counts <- Matrix::rowSums(counts[, donor_cells, drop = FALSE])
    pb_list[[donor]] <- donor_counts

    donor_meta <- meta[as.character(meta[[donor_col]]) == donor, , drop = FALSE]

    dx_values <- unique(as.character(donor_meta$diagnosis_std))
    dx_values <- dx_values[!is.na(dx_values)]

    if (length(dx_values) != 1) {
      warning(dataset_id, " donor ", donor,
              " has more than one diagnosis_std label: ",
              paste(dx_values, collapse = ", "))
    }

    meta_list[[donor]] <- data.frame(
      dataset = dataset_id,
      analysis = analysis_label,
      pseudobulk_id = paste(dataset_id, donor, analysis_label, sep = "_"),
      donor_id = donor,
      diagnosis_std = dx_values[1],
      n_microglia_cells = nrow(donor_meta),
      n_source_samples = 1,
      source_samples = donor,
      source_regions = NA_character_,
      stringsAsFactors = FALSE
    )

    cell_count_list[[donor]] <- data.frame(
      dataset = dataset_id,
      analysis = analysis_label,
      donor_id = donor,
      diagnosis_std = dx_values[1],
      n_cells = nrow(donor_meta),
      stringsAsFactors = FALSE
    )
  }

  pb_counts <- do.call(cbind, pb_list)
  colnames(pb_counts) <- names(pb_list)
  pb_counts <- Matrix(pb_counts, sparse = TRUE)

  pb_meta <- bind_rows(meta_list) %>%
    arrange(match(donor_id, colnames(pb_counts)))

  if (!all(pb_meta$donor_id == colnames(pb_counts))) {
    stop(dataset_id, ": donor metadata not aligned to pseudobulk count matrix.")
  }

  colnames(pb_counts) <- pb_meta$pseudobulk_id
  rownames(pb_meta) <- pb_meta$pseudobulk_id

  cell_counts <- bind_rows(cell_count_list)

  cat("  Final pseudobulk:", nrow(pb_counts), "genes x",
      ncol(pb_counts), "donors\n")
  cat("  Donor diagnosis table:\n")
  print(table(pb_meta$diagnosis_std, useNA = "ifany"))
  cat("\n")

  return(list(
    counts = pb_counts,
    meta = pb_meta,
    cell_counts = cell_counts
  ))
}

# ============================================================
# SECTION 3: GSE174367 author-labeled microglia
# ============================================================

process_GSE174367 <- function() {

  dataset_id <- "GSE174367"

  cat("============================================================\n")
  cat("Processing ", dataset_id, " author-labeled microglia\n", sep = "")
  cat("============================================================\n")

  h5_file <- DATASETS$GSE174367$h5_file
  meta_file <- DATASETS$GSE174367$meta_file

  if (!file.exists(h5_file)) {
    stop("Missing GSE174367 count H5 file: ", h5_file)
  }
  if (!file.exists(meta_file)) {
    stop("Missing GSE174367 metadata file: ", meta_file)
  }

  counts <- Read10X_h5(h5_file)

  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }

  meta <- readr::read_csv(meta_file, show_col_types = FALSE) %>%
    as.data.frame()

  cat("  Count matrix:", nrow(counts), "genes x", ncol(counts), "cells\n")
  cat("  Metadata rows:", nrow(meta), "\n")

  required_cols <- c("Barcode", "Cell.Type", "SampleID", "Diagnosis")
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("GSE174367 metadata missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  # Match count matrix columns to metadata barcodes.
  common_barcodes <- intersect(colnames(counts), meta$Barcode)

  cat("  Matched barcodes:", length(common_barcodes), "\n")
  cat("  Count-only barcodes removed:", ncol(counts) - length(common_barcodes), "\n")

  counts <- counts[, common_barcodes, drop = FALSE]
  meta <- meta[match(common_barcodes, meta$Barcode), , drop = FALSE]
  rownames(meta) <- common_barcodes

  cat("  Cell type table:\n")
  print(table(meta$Cell.Type, useNA = "ifany"))
  cat("  Diagnosis table:\n")
  print(table(meta$Diagnosis, useNA = "ifany"))

  pb <- make_donor_pseudobulk(
    counts = counts,
    meta = meta,
    dataset_id = dataset_id,
    donor_col = "SampleID",
    diagnosis_col = "Diagnosis",
    ad_labels = c("AD"),
    control_labels = c("Control"),
    celltype_col = "Cell.Type",
    microglia_labels = c("MG"),
    analysis_label = "PRIMARY",
    min_cells = PARAMS$min_cells_per_donor
  )

  saveRDS(pb$counts, file.path(RDS_DIR, "GSE174367_pseudobulk_counts.rds"))
  saveRDS(pb$meta,   file.path(RDS_DIR, "GSE174367_pseudobulk_meta.rds"))

  write.csv(pb$meta,
            file.path(TAB_DIR, "GSE174367_pseudobulk_meta.csv"),
            row.names = FALSE)
  write.csv(pb$cell_counts,
            file.path(TAB_DIR, "GSE174367_cell_counts_by_donor.csv"),
            row.names = FALSE)

  rm(counts, meta, pb)
  gc(verbose = FALSE)

  cat("Saved GSE174367 pseudobulk files.\n\n")
}

# ============================================================
# SECTION 4: GSE243292 microglia-specific H5AD
# ============================================================

# Helper: find an H5AD obs or var column robustly.
find_h5_name <- function(group, candidates) {
  available <- names(group)
  found <- intersect(candidates, available)
  if (length(found) == 0) {
    return(NA_character_)
  }
  found[1]
}

# Helper: decode H5AD categorical columns when needed.
decode_obs_column <- function(obs_group, col_name) {

  if (!col_name %in% names(obs_group)) {
    stop("Column not found in obs: ", col_name)
  }

  # AnnData categorical columns often store integer codes in obs
  # and labels in obs/__categories.
  if ("__categories" %in% names(obs_group) &&
      col_name %in% names(obs_group[["__categories"]])) {

    categories <- obs_group[["__categories"]][[col_name]][]
    codes <- obs_group[[col_name]][]
    codes_int <- as.integer(codes)

    # Handle signed 8-bit overflow if present.
    codes_int[codes_int < 0] <- codes_int[codes_int < 0] + 256L

    decoded <- categories[codes_int + 1L]
    return(as.character(decoded))
  }

  # Non-categorical column.
  values <- obs_group[[col_name]][]
  return(as.character(values))
}

# Helper: read H5AD X matrix as genes x cells sparse matrix.
read_h5ad_X_as_gene_by_cell <- function(h5, dataset_id) {

  X <- h5[["X"]]

  # Most GEO H5AD files store X as sparse CSR:
  #   shape = cells x genes
  #   indptr length = n_cells + 1
  #   indices are gene indices
  # We convert to genes x cells dgCMatrix for pseudobulk.
  if (all(c("data", "indices", "indptr") %in% names(X))) {

    shape <- X$attr_open("shape")$read()
    n_cells <- shape[1]
    n_genes <- shape[2]

    mat <- sparseMatrix(
      i = X[["indices"]][] + 1L,
      p = X[["indptr"]][],
      x = as.numeric(X[["data"]][]),
      dims = c(n_genes, n_cells),
      repr = "C"
    )

    return(mat)
  }

  stop(dataset_id, ": unsupported H5AD X format. X is not sparse CSR.")
}

process_GSE243292 <- function() {

  dataset_id <- "GSE243292"

  cat("============================================================\n")
  cat("Processing ", dataset_id, " microglia-specific H5AD\n", sep = "")
  cat("============================================================\n")

  h5ad_file <- DATASETS$GSE243292$h5ad_file

  # Support either .h5ad or .h5ad.gz, depending on what was copied.
  if (!file.exists(h5ad_file)) {
    gz_file <- paste0(h5ad_file, ".gz")
    if (file.exists(gz_file)) {
      cat("  Decompressing H5AD gz file...\n")
      R.utils::gunzip(gz_file, destname = h5ad_file, remove = FALSE, overwrite = TRUE)
    } else {
      stop("Missing GSE243292 H5AD file: ", h5ad_file)
    }
  }

  h5 <- H5File$new(h5ad_file, mode = "r")
  on.exit(h5$close_all(), add = TRUE)

  counts <- read_h5ad_X_as_gene_by_cell(h5, dataset_id)

  # Gene names are usually stored in var/_index or var/index.
  var_group <- h5[["var"]]
  gene_col <- find_h5_name(var_group, c("_index", "index", "gene_ids", "features"))

  if (is.na(gene_col)) {
    stop("Could not identify gene-name column in GSE243292 var group.")
  }

  gene_names <- as.character(var_group[[gene_col]][])
  rownames(counts) <- gene_names

  obs_group <- h5[["obs"]]

  sample_col <- find_h5_name(obs_group, c("sampleID", "sample_id", "SampleID", "donor_id"))
  atscore_col <- find_h5_name(obs_group, c("atscore", "ATscore", "AT_score"))

  if (is.na(sample_col)) {
    stop("Could not identify sample/donor column in GSE243292 obs group.")
  }
  if (is.na(atscore_col)) {
    stop("Could not identify atscore column in GSE243292 obs group.")
  }

  sample_id <- decode_obs_column(obs_group, sample_col)
  atscore <- decode_obs_column(obs_group, atscore_col)

  meta <- data.frame(
    sampleID = sample_id,
    atscore = atscore,
    stringsAsFactors = FALSE
  )

  # Create cell IDs so count columns and metadata rows align.
  cell_ids <- paste0("cell_", seq_len(ncol(counts)))
  colnames(counts) <- cell_ids
  rownames(meta) <- cell_ids

  cat("  Count matrix:", nrow(counts), "genes x", ncol(counts), "cells\n")
  cat("  atscore table before filtering:\n")
  print(table(meta$atscore, useNA = "ifany"))

  # A+T+ = AD, A-T- = Control.
  # A+T- is excluded because it is not clean AD or clean control.
  meta$diagnosis_for_filter <- ifelse(
    meta$atscore == "A+T+",
    "AD",
    ifelse(meta$atscore == "A-T-", "Control", "Exclude")
  )

  pb <- make_donor_pseudobulk(
    counts = counts,
    meta = meta,
    dataset_id = dataset_id,
    donor_col = "sampleID",
    diagnosis_col = "diagnosis_for_filter",
    ad_labels = c("AD"),
    control_labels = c("Control"),
    celltype_col = NULL,
    microglia_labels = NULL,
    analysis_label = "SUPPLEMENTARY",
    min_cells = PARAMS$min_cells_per_donor
  )

  saveRDS(pb$counts, file.path(RDS_DIR, "GSE243292_pseudobulk_counts.rds"))
  saveRDS(pb$meta,   file.path(RDS_DIR, "GSE243292_pseudobulk_meta.rds"))

  write.csv(pb$meta,
            file.path(TAB_DIR, "GSE243292_pseudobulk_meta.csv"),
            row.names = FALSE)
  write.csv(pb$cell_counts,
            file.path(TAB_DIR, "GSE243292_cell_counts_by_donor.csv"),
            row.names = FALSE)

  rm(counts, meta, pb)
  gc(verbose = FALSE)

  cat("Saved GSE243292 pseudobulk files.\n\n")
}

# ============================================================
# SECTION 5: Run both datasets
# ============================================================

process_GSE174367()
process_GSE243292()

cat("============================================================\n")
cat("03G complete.\n\n")
cat("Next step:\n")
cat("  Rerun 03F_check_pseudobulk_inputs_before_DE_AD61026.R\n")
cat("  and upload 03F_pseudobulk_input_check_summary.csv.\n")
cat("============================================================\n")
