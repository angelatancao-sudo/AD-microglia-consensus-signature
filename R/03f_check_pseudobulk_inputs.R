# ============================================================
# 03f_check_pseudobulk_inputs.R
# ============================================================
# Project: AD Microglia Multi-Dataset Transcriptomic Analysis
#
# Purpose:
#   Check that all pseudobulk count and metadata files needed for
#   the next differential expression step are present and aligned.
#
# This script does NOT modify any data.
#
# It checks:
#   1. SEA-AD discovery pseudobulk files
#   2. GSE174367 author-annotated microglia pseudobulk files
#   3. Newly rebuilt validated raw-GEO microglia pseudobulk files
#      for GSE157827, GSE160936, and GSE188545
#   4. GSE243292 author microglia pseudobulk files, if present
#   5. Optional sensitivity pseudobulk files from 03E
#
# Why this matters:
#   Before running 04_DE_analysis, we need to make sure that every
#   count matrix has matching metadata and that the AD/control donor
#   counts are correct.
#
# Run after:
#   03E_extract_validated_microglia_pseudobulk_AD61026.R
# ============================================================

# ============================================================
# SECTION 1: Load setup
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
  library(Matrix)
  library(dplyr)
})

cat("============================================================\n")
cat("  03F: Check pseudobulk inputs before DE\n")
cat("============================================================\n\n")

dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# SECTION 2: Define expected files
# ============================================================

# Main datasets expected for primary DE/meta-analysis.
PRIMARY_DATASETS <- data.frame(
  dataset = c(
    "SEAAD",
    "GSE174367",
    "GSE157827",
    "GSE160936",
    "GSE188545",
    "GSE243292"
  ),
  analysis = c(
    "DISCOVERY",
    "PRIMARY",
    "PRIMARY",
    "PRIMARY",
    "PRIMARY",
    "SUPPLEMENTARY"
  ),
  counts_file = c(
    file.path(RDS_DIR, "SEAAD_pseudobulk_counts.rds"),
    file.path(RDS_DIR, "GSE174367_pseudobulk_counts.rds"),
    file.path(RDS_DIR, "GSE157827_microglia_PRIMARY_pseudobulk_counts.rds"),
    file.path(RDS_DIR, "GSE160936_microglia_PRIMARY_pseudobulk_counts.rds"),
    file.path(RDS_DIR, "GSE188545_microglia_PRIMARY_pseudobulk_counts.rds"),
    file.path(RDS_DIR, "GSE243292_pseudobulk_counts.rds")
  ),
  meta_file = c(
    file.path(RDS_DIR, "SEAAD_pseudobulk_meta.rds"),
    file.path(RDS_DIR, "GSE174367_pseudobulk_meta.rds"),
    file.path(RDS_DIR, "GSE157827_microglia_PRIMARY_pseudobulk_meta.rds"),
    file.path(RDS_DIR, "GSE160936_microglia_PRIMARY_pseudobulk_meta.rds"),
    file.path(RDS_DIR, "GSE188545_microglia_PRIMARY_pseudobulk_meta.rds"),
    file.path(RDS_DIR, "GSE243292_pseudobulk_meta.rds")
  ),
  stringsAsFactors = FALSE
)

# Optional sensitivity datasets for raw-GEO cluster definition checks.
SENSITIVITY_DATASETS <- data.frame(
  dataset = c("GSE157827", "GSE160936", "GSE188545"),
  analysis = c("SENSITIVITY", "SENSITIVITY", "SENSITIVITY"),
  counts_file = c(
    file.path(RDS_DIR, "GSE157827_microglia_SENSITIVITY_pseudobulk_counts.rds"),
    file.path(RDS_DIR, "GSE160936_microglia_SENSITIVITY_pseudobulk_counts.rds"),
    file.path(RDS_DIR, "GSE188545_microglia_SENSITIVITY_pseudobulk_counts.rds")
  ),
  meta_file = c(
    file.path(RDS_DIR, "GSE157827_microglia_SENSITIVITY_pseudobulk_meta.rds"),
    file.path(RDS_DIR, "GSE160936_microglia_SENSITIVITY_pseudobulk_meta.rds"),
    file.path(RDS_DIR, "GSE188545_microglia_SENSITIVITY_pseudobulk_meta.rds")
  ),
  stringsAsFactors = FALSE
)

FILES_TO_CHECK <- bind_rows(PRIMARY_DATASETS, SENSITIVITY_DATASETS)

# ============================================================
# SECTION 3: Helper function to check one dataset
# ============================================================

check_one_dataset <- function(row) {

  dataset <- row$dataset
  analysis <- row$analysis
  counts_file <- row$counts_file
  meta_file <- row$meta_file

  cat("------------------------------------------------------------\n")
  cat(dataset, "/", analysis, "\n")
  cat("------------------------------------------------------------\n")

  counts_exists <- file.exists(counts_file)
  meta_exists <- file.exists(meta_file)

  if (!counts_exists || !meta_exists) {
    cat("  MISSING FILE(S)\n")
    cat("  counts exists:", counts_exists, "\n")
    cat("  meta exists  :", meta_exists, "\n\n")

    return(data.frame(
      dataset = dataset,
      analysis = analysis,
      counts_exists = counts_exists,
      meta_exists = meta_exists,
      status = "MISSING",
      n_genes = NA_integer_,
      n_pseudobulk_samples = NA_integer_,
      n_meta_rows = NA_integer_,
      metadata_aligned = NA,
      n_AD = NA_integer_,
      n_Control = NA_integer_,
      min_microglia_cells = NA_integer_,
      median_microglia_cells = NA_real_,
      max_microglia_cells = NA_integer_,
      notes = "One or both RDS files are missing.",
      stringsAsFactors = FALSE
    ))
  }

  counts <- readRDS(counts_file)
  meta <- readRDS(meta_file)

  # Convert dense matrix to sparse summary-compatible matrix if needed.
  n_genes <- nrow(counts)
  n_samples <- ncol(counts)
  n_meta_rows <- nrow(meta)

  # Determine diagnosis column.
  dx_col <- if ("diagnosis_std" %in% colnames(meta)) {
    "diagnosis_std"
  } else if ("Diagnosis" %in% colnames(meta)) {
    "Diagnosis"
  } else if ("diagnosis" %in% colnames(meta)) {
    "diagnosis"
  } else {
    NA_character_
  }

  if (is.na(dx_col)) {
    n_AD <- NA_integer_
    n_Control <- NA_integer_
    dx_note <- "No diagnosis column found."
  } else {
    dx_values <- as.character(meta[[dx_col]])
    n_AD <- sum(dx_values == "AD", na.rm = TRUE)
    n_Control <- sum(dx_values == "Control", na.rm = TRUE)
    dx_note <- paste0("Diagnosis column: ", dx_col)
  }

  # Check metadata alignment.
  # For the new 03E pseudobulk files, rownames(meta) should match colnames(counts).
  # For older files, rownames may be donor IDs. We check exact alignment first.
  metadata_aligned <- identical(rownames(meta), colnames(counts))

  # If exact rowname alignment fails, check whether metadata has pseudobulk_id
  # or donor/sample columns matching counts.
  alt_alignment_note <- ""
  if (!metadata_aligned) {
    possible_id_cols <- intersect(
      c("pseudobulk_id", "donor_id", "SampleID", "sampleID", "sample_id"),
      colnames(meta)
    )

    for (id_col in possible_id_cols) {
      if (identical(as.character(meta[[id_col]]), colnames(counts))) {
        metadata_aligned <- TRUE
        alt_alignment_note <- paste0("Aligned by column ", id_col, ".")
        break
      }
    }
  }

  # Capture n_microglia_cells if available.
  if ("n_microglia_cells" %in% colnames(meta)) {
    cell_counts <- meta$n_microglia_cells
    min_cells <- min(cell_counts, na.rm = TRUE)
    median_cells <- median(cell_counts, na.rm = TRUE)
    max_cells <- max(cell_counts, na.rm = TRUE)
  } else if ("n_cells" %in% colnames(meta)) {
    cell_counts <- meta$n_cells
    min_cells <- min(cell_counts, na.rm = TRUE)
    median_cells <- median(cell_counts, na.rm = TRUE)
    max_cells <- max(cell_counts, na.rm = TRUE)
  } else {
    min_cells <- NA_integer_
    median_cells <- NA_real_
    max_cells <- NA_integer_
  }

  notes <- paste(dx_note, alt_alignment_note)

  cat("  Genes:", n_genes, "\n")
  cat("  Pseudobulk samples:", n_samples, "\n")
  cat("  Metadata rows:", n_meta_rows, "\n")
  cat("  Metadata aligned:", metadata_aligned, "\n")
  cat("  AD:", n_AD, " Control:", n_Control, "\n")
  cat("  Min/median/max microglia cells:",
      min_cells, "/", median_cells, "/", max_cells, "\n\n")

  status <- if (metadata_aligned && n_samples == n_meta_rows) {
    "OK"
  } else {
    "CHECK_ALIGNMENT"
  }

  return(data.frame(
    dataset = dataset,
    analysis = analysis,
    counts_exists = counts_exists,
    meta_exists = meta_exists,
    status = status,
    n_genes = n_genes,
    n_pseudobulk_samples = n_samples,
    n_meta_rows = n_meta_rows,
    metadata_aligned = metadata_aligned,
    n_AD = n_AD,
    n_Control = n_Control,
    min_microglia_cells = min_cells,
    median_microglia_cells = median_cells,
    max_microglia_cells = max_cells,
    notes = notes,
    stringsAsFactors = FALSE
  ))
}

# ============================================================
# SECTION 4: Run checks
# ============================================================

check_results <- bind_rows(lapply(seq_len(nrow(FILES_TO_CHECK)), function(i) {
  check_one_dataset(FILES_TO_CHECK[i, ])
}))

out_path <- file.path(TAB_DIR, "03F_pseudobulk_input_check_summary.csv")

write.csv(check_results, out_path, row.names = FALSE)

cat("============================================================\n")
cat("03F check complete.\n\n")
cat("Summary saved to:\n")
cat("  ", out_path, "\n\n")
cat("Rows with status OK are ready for DE input.\n")
cat("Rows with MISSING usually mean that those datasets still need to be loaded.\n")
cat("Rows with CHECK_ALIGNMENT should be reviewed before DE.\n")
cat("============================================================\n")
