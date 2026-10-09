# ============================================================
# 03b_load_GEO_and_cluster.R
# ============================================================
# Project: AD Microglia Multi-Dataset Transcriptomic Analysis
#
# Purpose:
#   Start the field-aligned GEO reprocessing step.
#
#   This script does NOT use the old per-cell marker-positive rule.
#   Instead, for raw unannotated datasets it performs:
#
#     raw 10x loading
#     per-nucleus QC
#     normalization
#     highly variable feature selection
#     PCA
#     graph-based clustering
#     UMAP visualization
#     cluster-level marker summary
#     dot plots for canonical cell-type markers
#
#   The goal is to generate the evidence needed to annotate
#   microglial clusters by cluster-level marker enrichment.
#
# Why this is staged:
#   Cluster annotation should be reviewed before extracting
#   microglia. Automatically selecting clusters using an arbitrary
#   score cutoff would create another reviewer-sensitive criterion.
#   This script creates cluster review tables and plots first.
#
# Datasets handled here:
#   GSE157827  raw 10x, 21 donors, PFC
#   GSE160936  raw 10x tar.gz, 24 region samples from 12 donors, EC + SSC
#   GSE188545  raw 10x, 12 donors, MTG
#
# Output:
#   RDS/GSE157827_clustered_seurat.rds
#   RDS/GSE160936_clustered_seurat.rds
#   RDS/GSE188545_clustered_seurat.rds
#
#   Output/QC/[dataset]_QC_violin.pdf
#   Output/QC/[dataset]_UMAP_by_cluster.pdf
#   Output/QC/[dataset]_UMAP_by_diagnosis.pdf
#   Output/QC/[dataset]_DotPlot_celltype_markers.pdf
#   Output/Tables/[dataset]_cluster_marker_summary.csv
#   Output/Tables/[dataset]_microglia_cluster_review_template.csv
#
# Run after:
#   01_Setup_AD61026.R
#   02_load_SEAAD_AD61026.R
#   03B_parse_GSE160936_metadata_AD61026_fixed.R
#
# Important:
#   Restart R/RStudio before running this script.
#   This script can use substantial memory.
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
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

cat("============================================================\n")
cat("  03 Stage 1: Cluster-review workflow for raw GEO datasets\n")
cat("============================================================\n\n")

# Create output folders if missing.
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# SECTION 2: Field-aligned marker panels
# ============================================================
# These marker panels are used for cluster-level annotation.
# They are not used to select single cells by marker positivity.

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

ALL_MARKERS <- unique(unlist(MARKERS))

# ============================================================
# SECTION 3: Helper functions for reading raw 10x files
# ============================================================

safe_unlink <- function(path) {
  if (dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE)
  }
}

# ------------------------------------------------------------
# Function: read_flat_10x_sample
# Purpose:
#   Read one sample from a flat GEO folder where each sample has
#   separate barcode/features/matrix files.
#
# Handles both:
#   features.tsv.gz
#   genes.tsv.gz
# ------------------------------------------------------------
read_flat_10x_sample <- function(raw_dir, sample_name, dataset_id) {

  temp_dir <- file.path(tempdir(), paste0("read10x_", dataset_id, "_", sample_name))
  safe_unlink(temp_dir)
  dir.create(temp_dir, recursive = TRUE)

  sample_files <- list.files(
    raw_dir,
    pattern = paste0("_", sample_name, "_"),
    full.names = TRUE
  )

  if (length(sample_files) == 0) {
    stop("No files found for sample ", sample_name, " in ", raw_dir)
  }

  copied_barcodes <- FALSE
  copied_features <- FALSE
  copied_matrix   <- FALSE

  for (f in sample_files) {
    bn <- basename(f)

    if (grepl("barcodes", bn)) {
      file.copy(f, file.path(temp_dir, "barcodes.tsv.gz"), overwrite = TRUE)
      copied_barcodes <- TRUE

    } else if (grepl("features", bn) || grepl("genes", bn)) {
      # Read10X expects features.tsv.gz.
      # Some older GEO files use genes.tsv.gz.
      file.copy(f, file.path(temp_dir, "features.tsv.gz"), overwrite = TRUE)
      copied_features <- TRUE

    } else if (grepl("matrix", bn)) {
      file.copy(f, file.path(temp_dir, "matrix.mtx.gz"), overwrite = TRUE)
      copied_matrix <- TRUE
    }
  }

  if (!all(c(copied_barcodes, copied_features, copied_matrix))) {
    print(sample_files)
    stop("Missing one or more 10x files for sample ", sample_name)
  }

  counts <- Read10X(data.dir = temp_dir)
  safe_unlink(temp_dir)

  # In case Read10X returns a list because multiple assays exist,
  # keep Gene Expression.
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }

  return(counts)
}

# ------------------------------------------------------------
# Function: read_tarred_10x_sample
# Purpose:
#   Read one 10x sample from a .tar.gz file.
# ------------------------------------------------------------
read_tarred_10x_sample <- function(tar_path, dataset_id, sample_id) {

  if (!file.exists(tar_path)) {
    stop("Missing tar.gz file: ", tar_path)
  }

  temp_dir <- file.path(tempdir(), paste0("untar_", dataset_id, "_", sample_id))
  safe_unlink(temp_dir)
  dir.create(temp_dir, recursive = TRUE)

  utils::untar(tar_path, exdir = temp_dir)

  matrix_files <- list.files(
    temp_dir,
    pattern = "matrix\\.mtx(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(matrix_files) == 0) {
    cat("Unpacked files:\n")
    print(list.files(temp_dir, recursive = TRUE, full.names = FALSE))
    stop("No matrix.mtx or matrix.mtx.gz found in ", tar_path)
  }

  tenx_dir <- dirname(matrix_files[1])
  counts <- Read10X(data.dir = tenx_dir)

  safe_unlink(temp_dir)

  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }

  return(counts)
}

# ============================================================
# SECTION 4: Helper functions for Seurat processing and review
# ============================================================

# ------------------------------------------------------------
# Function: basic_qc_filter
# Purpose:
#   Apply transparent per-nucleus QC before clustering.
#
# Notes:
#   - pct mitochondrial is calculated from genes starting with MT-.
#   - GSE160936 uses the original Smith-style thresholds:
#       nFeature_RNA >= 200
#       nFeature_RNA <= 6000
#       nCount_RNA   <= 25000
#       percent.mt   <= 5
#   - For other raw datasets, we use the project defaults:
#       min features 200
#       mito <= 5%
#     without imposing an invented upper feature/count threshold
#     unless specified in qc_params.
# ------------------------------------------------------------
basic_qc_filter <- function(seu, dataset_id, qc_params) {

  # Calculate mitochondrial percentage.
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

  n_before <- ncol(seu)

  keep <- rep(TRUE, ncol(seu))

  if (!is.null(qc_params$min_features) && is.finite(qc_params$min_features)) {
    keep <- keep & seu$nFeature_RNA >= qc_params$min_features
  }

  if (!is.null(qc_params$max_features) && is.finite(qc_params$max_features)) {
    keep <- keep & seu$nFeature_RNA <= qc_params$max_features
  }

  if (!is.null(qc_params$max_counts) && is.finite(qc_params$max_counts)) {
    keep <- keep & seu$nCount_RNA <= qc_params$max_counts
  }

  if (!is.null(qc_params$max_mito_pct) && is.finite(qc_params$max_mito_pct)) {
    keep <- keep & seu$percent.mt <= qc_params$max_mito_pct
  }

  seu <- subset(seu, cells = colnames(seu)[keep])

  n_after <- ncol(seu)

  cat("  QC:", dataset_id, "\n")
  cat("    Cells before QC:", n_before, "\n")
  cat("    Cells after QC :", n_after, "\n")
  cat("    Removed        :", n_before - n_after, "\n\n")

  return(seu)
}

# ------------------------------------------------------------
# Function: run_seurat_clustering
# Purpose:
#   Run standard Seurat normalization, HVG selection, scaling,
#   PCA, nearest-neighbor graph, clustering, and UMAP.
# ------------------------------------------------------------
run_seurat_clustering <- function(seu, dataset_id) {

  cat("  Running Seurat workflow for", dataset_id, "\n")

  DefaultAssay(seu) <- "RNA"

  seu <- NormalizeData(
    seu,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  seu <- FindVariableFeatures(
    seu,
    selection.method = "vst",
    nfeatures = PARAMS$seurat_nfeatures,
    verbose = FALSE
  )

  seu <- ScaleData(
    seu,
    features = VariableFeatures(seu),
    verbose = FALSE
  )

  seu <- RunPCA(
    seu,
    features = VariableFeatures(seu),
    npcs = PARAMS$seurat_npcs,
    verbose = FALSE
  )

  seu <- FindNeighbors(
    seu,
    dims = 1:PARAMS$seurat_npcs,
    verbose = FALSE
  )

  seu <- FindClusters(
    seu,
    resolution = PARAMS$seurat_resolution,
    verbose = FALSE
  )

  seu <- RunUMAP(
    seu,
    dims = 1:PARAMS$seurat_npcs,
    verbose = FALSE
  )

  cat("    Clusters found:", length(unique(seu$seurat_clusters)), "\n\n")

  return(seu)
}

# ------------------------------------------------------------
# Function: add_module_scores_by_panel
# Purpose:
#   Add module scores for each marker panel.
#
# Note:
#   Module scores are used only to summarize clusters for review.
#   We do not use a fixed score cutoff to define microglia.
# ------------------------------------------------------------
add_module_scores_by_panel <- function(seu) {

  for (panel_name in names(MARKERS)) {
    genes <- intersect(MARKERS[[panel_name]], rownames(seu))

    if (length(genes) >= 2) {
      seu <- AddModuleScore(
        object = seu,
        features = list(genes),
        name = paste0(panel_name, "_score"),
        assay = "RNA",
        search = FALSE
      )
      # AddModuleScore appends "1" to the name.
      old_col <- paste0(panel_name, "_score1")
      new_col <- paste0(panel_name, "_score")
      colnames(seu@meta.data)[colnames(seu@meta.data) == old_col] <- new_col
    } else {
      seu@meta.data[[paste0(panel_name, "_score")]] <- NA_real_
    }
  }

  return(seu)
}

# ------------------------------------------------------------
# Function: make_cluster_review_outputs
# Purpose:
#   Create dot plots, UMAPs, and a cluster-level summary table
#   to support manual annotation of microglial clusters.
# ------------------------------------------------------------
make_cluster_review_outputs <- function(seu, dataset_id) {

  cat("  Creating cluster review outputs for", dataset_id, "\n")

  # Ensure clusters are characters for cleaner tables.
  seu$cluster_id <- as.character(seu$seurat_clusters)

  # ----------------------------------------------------------
  # QC violin plot
  # ----------------------------------------------------------
  p_qc <- VlnPlot(
    seu,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = "sample_id",
    pt.size = 0,
    ncol = 3
  ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

  ggsave(
    file.path(QC_DIR, paste0(dataset_id, "_QC_violin.pdf")),
    p_qc,
    width = 13,
    height = 5
  )

  # ----------------------------------------------------------
  # UMAP by cluster
  # ----------------------------------------------------------
  p_umap_cluster <- DimPlot(
    seu,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle(paste0(dataset_id, " clusters"))

  ggsave(
    file.path(QC_DIR, paste0(dataset_id, "_UMAP_by_cluster.pdf")),
    p_umap_cluster,
    width = 7,
    height = 6
  )

  # ----------------------------------------------------------
  # UMAP by diagnosis, if available
  # ----------------------------------------------------------
  if ("diagnosis_std" %in% colnames(seu@meta.data)) {
    p_umap_dx <- DimPlot(
      seu,
      reduction = "umap",
      group.by = "diagnosis_std"
    ) +
      ggtitle(paste0(dataset_id, " diagnosis"))

    ggsave(
      file.path(QC_DIR, paste0(dataset_id, "_UMAP_by_diagnosis.pdf")),
      p_umap_dx,
      width = 7,
      height = 6
    )
  }

  # ----------------------------------------------------------
  # Dot plot of marker panels
  # ----------------------------------------------------------
  marker_genes_present <- intersect(ALL_MARKERS, rownames(seu))

  if (length(marker_genes_present) >= 2) {
    p_dot <- DotPlot(
      seu,
      features = marker_genes_present,
      group.by = "seurat_clusters"
    ) +
      RotatedAxis() +
      ggtitle(paste0(dataset_id, " cluster marker dot plot"))

    ggsave(
      file.path(QC_DIR, paste0(dataset_id, "_DotPlot_celltype_markers.pdf")),
      p_dot,
      width = 14,
      height = 7
    )
  }

  # ----------------------------------------------------------
  # Cluster-level module score summary
  # ----------------------------------------------------------
  score_cols <- paste0(names(MARKERS), "_score")
  score_cols <- score_cols[score_cols %in% colnames(seu@meta.data)]

  cluster_summary <- seu@meta.data %>%
    as.data.frame() %>%
    group_by(cluster_id) %>%
    summarise(
      n_cells = n(),
      n_donors = dplyr::n_distinct(donor_id),
      n_samples = dplyr::n_distinct(sample_id),
      diagnosis_mix = paste(names(table(diagnosis_std)),
                            as.integer(table(diagnosis_std)),
                            collapse = "; "),
      across(all_of(score_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(as.numeric(cluster_id))

  write.csv(
    cluster_summary,
    file.path(TAB_DIR, paste0(dataset_id, "_cluster_marker_summary.csv")),
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # Template for manual microglia cluster review
  # ----------------------------------------------------------
  review_template <- cluster_summary %>%
    mutate(
      likely_microglia_candidate = NA,
      keep_as_microglia = NA,
      reviewer_note = NA
    )

  write.csv(
    review_template,
    file.path(TAB_DIR, paste0(dataset_id, "_microglia_cluster_review_template.csv")),
    row.names = FALSE
  )

  cat("    Saved review table and plots for", dataset_id, "\n\n")

  return(seu)
}

# ============================================================
# SECTION 5: Dataset-specific loaders
# ============================================================

# ------------------------------------------------------------
# Function: load_GSE157827_raw
# ------------------------------------------------------------
load_GSE157827_raw <- function() {

  dataset_id <- "GSE157827"
  raw_dir <- DATASETS$GSE157827$raw_dir

  cat("============================================================\n")
  cat("Loading", dataset_id, "\n")
  cat("============================================================\n")

  files <- list.files(raw_dir, full.names = FALSE)

  sample_names <- unique(
    sub(".*_(AD[0-9]+|NC[0-9]+)_.*", "\\1", files)
  )
  sample_names <- sample_names[grepl("^(AD|NC)[0-9]+$", sample_names)]
  sample_names <- sort(sample_names)

  cat("Samples parsed:", length(sample_names), "\n")
  cat("AD samples:", sum(grepl("^AD", sample_names)), "\n")
  cat("Control samples:", sum(grepl("^NC", sample_names)), "\n\n")

  seu_list <- list()

  for (sname in sample_names) {

    cat("  Reading", sname, "\n")

    counts <- read_flat_10x_sample(
      raw_dir = raw_dir,
      sample_name = sname,
      dataset_id = dataset_id
    )

    seu <- CreateSeuratObject(
      counts = counts,
      project = paste0(dataset_id, "_", sname),
      min.cells = 0,
      min.features = 0
    )

    seu$dataset <- dataset_id
    seu$sample_id <- sname
    seu$donor_id <- sname
    seu$diagnosis_std <- ifelse(grepl("^AD", sname), "AD", "Control")
    seu$diagnosis_std <- factor(seu$diagnosis_std, levels = c("Control", "AD"))

    # Make cell names globally unique before merging.
    seu <- RenameCells(seu, add.cell.id = sname)

    seu_list[[sname]] <- seu

    rm(counts, seu)
    gc(verbose = FALSE)
  }

  cat("\nMerging", length(seu_list), "samples for", dataset_id, "\n")
  seu_merged <- merge(seu_list[[1]], y = seu_list[-1])
  rm(seu_list)
  gc(verbose = FALSE)

  qc_params <- list(
    min_features = PARAMS$sn_min_features,
    max_features = Inf,
    max_counts = Inf,
    max_mito_pct = PARAMS$sn_max_mito_pct
  )

  seu_merged <- basic_qc_filter(seu_merged, dataset_id, qc_params)
  seu_merged <- run_seurat_clustering(seu_merged, dataset_id)
  seu_merged <- add_module_scores_by_panel(seu_merged)
  seu_merged <- make_cluster_review_outputs(seu_merged, dataset_id)

  saveRDS(seu_merged, file.path(RDS_DIR, paste0(dataset_id, "_clustered_seurat.rds")))
  cat("Saved:", file.path(RDS_DIR, paste0(dataset_id, "_clustered_seurat.rds")), "\n\n")

  return(invisible(seu_merged))
}

# ------------------------------------------------------------
# Function: load_GSE188545_raw
# ------------------------------------------------------------
load_GSE188545_raw <- function() {

  dataset_id <- "GSE188545"
  raw_dir <- DATASETS$GSE188545$raw_dir

  cat("============================================================\n")
  cat("Loading", dataset_id, "\n")
  cat("============================================================\n")

  files <- list.files(raw_dir, full.names = FALSE)

  sample_names <- unique(
    sub(".*_(AD[0-9A-Za-z]+|HC[0-9A-Za-z]+)_.*", "\\1", files)
  )
  sample_names <- sample_names[grepl("^(AD|HC)", sample_names)]
  sample_names <- sort(sample_names)

  cat("Samples parsed:", length(sample_names), "\n")
  cat("AD samples:", sum(grepl("^AD", sample_names)), "\n")
  cat("Control samples:", sum(grepl("^HC", sample_names)), "\n\n")

  seu_list <- list()

  for (sname in sample_names) {

    cat("  Reading", sname, "\n")

    counts <- read_flat_10x_sample(
      raw_dir = raw_dir,
      sample_name = sname,
      dataset_id = dataset_id
    )

    seu <- CreateSeuratObject(
      counts = counts,
      project = paste0(dataset_id, "_", sname),
      min.cells = 0,
      min.features = 0
    )

    seu$dataset <- dataset_id
    seu$sample_id <- sname
    seu$donor_id <- sname
    seu$diagnosis_std <- ifelse(grepl("^AD", sname), "AD", "Control")
    seu$diagnosis_std <- factor(seu$diagnosis_std, levels = c("Control", "AD"))

    seu <- RenameCells(seu, add.cell.id = sname)

    seu_list[[sname]] <- seu

    rm(counts, seu)
    gc(verbose = FALSE)
  }

  cat("\nMerging", length(seu_list), "samples for", dataset_id, "\n")
  seu_merged <- merge(seu_list[[1]], y = seu_list[-1])
  rm(seu_list)
  gc(verbose = FALSE)

  qc_params <- list(
    min_features = PARAMS$sn_min_features,
    max_features = Inf,
    max_counts = Inf,
    max_mito_pct = PARAMS$sn_max_mito_pct
  )

  seu_merged <- basic_qc_filter(seu_merged, dataset_id, qc_params)
  seu_merged <- run_seurat_clustering(seu_merged, dataset_id)
  seu_merged <- add_module_scores_by_panel(seu_merged)
  seu_merged <- make_cluster_review_outputs(seu_merged, dataset_id)

  saveRDS(seu_merged, file.path(RDS_DIR, paste0(dataset_id, "_clustered_seurat.rds")))
  cat("Saved:", file.path(RDS_DIR, paste0(dataset_id, "_clustered_seurat.rds")), "\n\n")

  return(invisible(seu_merged))
}

# ------------------------------------------------------------
# Function: load_GSE160936_raw
# ------------------------------------------------------------
load_GSE160936_raw <- function() {

  dataset_id <- "GSE160936"
  raw_dir <- DATASETS$GSE160936$raw_dir
  meta_path <- DATASETS$GSE160936$parsed_meta

  cat("============================================================\n")
  cat("Loading", dataset_id, "\n")
  cat("============================================================\n")

  if (!file.exists(meta_path)) {
    stop("Parsed metadata missing: ", meta_path,
         "\nRun 03B_parse_GSE160936_metadata_AD61026_fixed.R first.")
  }

  sample_meta <- read.csv(meta_path, stringsAsFactors = FALSE)

  cat("Metadata rows:", nrow(sample_meta), "\n")
  cat("Unique donors:", length(unique(sample_meta$donor_id)), "\n")
  cat("Sample-level diagnosis table:\n")
  print(table(sample_meta$diagnosis_std))
  cat("\n")

  seu_list <- list()

  for (i in seq_len(nrow(sample_meta))) {

    row <- sample_meta[i, ]
    sname <- row$local_sample_id

    cat("  Reading", sname, " donor=", row$donor_id,
        " region=", row$region, "\n", sep = "")

    tar_path <- file.path(raw_dir, row$file)

    counts <- read_tarred_10x_sample(
      tar_path = tar_path,
      dataset_id = dataset_id,
      sample_id = sname
    )

    seu <- CreateSeuratObject(
      counts = counts,
      project = paste0(dataset_id, "_", sname),
      min.cells = 0,
      min.features = 0
    )

    seu$dataset <- dataset_id
    seu$sample_id <- sname
    seu$gsm_id <- row$gsm_id
    seu$donor_id <- row$donor_id
    seu$region <- row$region
    seu$disease_state <- row$disease_state
    seu$diagnosis_std <- factor(row$diagnosis_std, levels = c("Control", "AD"))
    seu$braak_stage <- row$braak_stage
    seu$age <- row$age
    seu$sex <- row$sex
    seu$rin <- row$rin

    seu <- RenameCells(seu, add.cell.id = sname)

    seu_list[[sname]] <- seu

    rm(counts, seu)
    gc(verbose = FALSE)
  }

  cat("\nMerging", length(seu_list), "samples for", dataset_id, "\n")
  seu_merged <- merge(seu_list[[1]], y = seu_list[-1])
  rm(seu_list)
  gc(verbose = FALSE)

  # Dataset-specific QC reflecting the original Smith-style thresholds
  # described during planning:
  #   nFeature_RNA >= 200
  #   nFeature_RNA <= 6000
  #   nCount_RNA   <= 25000
  #   percent.mt   <= 5
  qc_params <- list(
    min_features = 200,
    max_features = 6000,
    max_counts = 25000,
    max_mito_pct = 5
  )

  seu_merged <- basic_qc_filter(seu_merged, dataset_id, qc_params)
  seu_merged <- run_seurat_clustering(seu_merged, dataset_id)
  seu_merged <- add_module_scores_by_panel(seu_merged)
  seu_merged <- make_cluster_review_outputs(seu_merged, dataset_id)

  saveRDS(seu_merged, file.path(RDS_DIR, paste0(dataset_id, "_clustered_seurat.rds")))
  cat("Saved:", file.path(RDS_DIR, paste0(dataset_id, "_clustered_seurat.rds")), "\n\n")

  return(invisible(seu_merged))
}

# ============================================================
# SECTION 6: Run datasets one at a time
# ============================================================
# Important:
#   These objects can be large. We process, save, remove, and run gc()
#   after each dataset.

cat("Starting raw dataset clustering.\n")
cat("This may take a long time. Do not run other R scripts in parallel.\n\n")

seu_157827 <- load_GSE157827_raw()
rm(seu_157827)
gc(verbose = FALSE)

seu_188545 <- load_GSE188545_raw()
rm(seu_188545)
gc(verbose = FALSE)

seu_160936 <- load_GSE160936_raw()
rm(seu_160936)
gc(verbose = FALSE)

cat("============================================================\n")
cat("Stage 1 complete.\n\n")
cat("Next action:\n")
cat("  Open the following files and inspect marker evidence:\n")
cat("    Output/Tables/GSE157827_microglia_cluster_review_template.csv\n")
cat("    Output/Tables/GSE160936_microglia_cluster_review_template.csv\n")
cat("    Output/Tables/GSE188545_microglia_cluster_review_template.csv\n\n")
cat("  Also inspect the dot plots and UMAPs in Output/QC/.\n\n")
cat("  After review, fill keep_as_microglia = TRUE/FALSE for each cluster.\n")
cat("  Then we will run Stage 2 to extract microglia and pseudobulk.\n")
cat("============================================================\n")
