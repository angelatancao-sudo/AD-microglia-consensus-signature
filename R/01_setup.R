# ============================================================
# 01_setup.R
# ============================================================
# Project: AD Microglia Multi-Dataset Transcriptomic Analysis
#
# Purpose:
#   1. Define the restarted AD61026 project directory.
#   2. Create the expected folder structure.
#   3. Load/install required R packages.
#   4. Save one central project setup object used by all later scripts.
#
# Important:
#   This script does NOT run any analysis.
#   It only prepares paths, package checks, dataset registry, and parameters.
#
# Run order:
#   source("R/01_Setup_AD61026.R")
#
# Restart note:
#   This version fixes the earlier GSE157827 tar_file/raw_dir mismatch.
#   GSE157827, GSE160936, and GSE188545 are raw-folder datasets.
# ============================================================

# ============================================================
# SECTION 1: Define project root directory
# ============================================================

# ------------------------------------------------------------
# Project folder. Set the AD_PROJECT_ROOT environment variable to the
# folder that contains Data/, RDS/, Output/ and Supplementary_Tables/,
# or start R with that folder as the working directory.
# ------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("AD_PROJECT_ROOT", unset = getwd())

ROOT_DIR <- PROJECT_ROOT  # all other folders are created inside this folder

DATA_DIR    <- file.path(ROOT_DIR, "Data")
SCRIPTS_DIR <- file.path(ROOT_DIR, "Scripts")
OUTPUT_DIR  <- file.path(ROOT_DIR, "Output")
RDS_DIR     <- file.path(ROOT_DIR, "RDS")
QC_DIR      <- file.path(OUTPUT_DIR, "QC")
DE_DIR      <- file.path(OUTPUT_DIR, "DE")
FIG_DIR     <- file.path(OUTPUT_DIR, "Figures")
TAB_DIR     <- file.path(OUTPUT_DIR, "Tables")

cat("============================================================\n")
cat("  01 Setup: AD61026 project\n")
cat("============================================================\n\n")
cat("ROOT_DIR:", ROOT_DIR, "\n\n")

# ============================================================
# SECTION 2: Create folder structure
# ============================================================

dirs_to_create <- c(
  DATA_DIR,
  SCRIPTS_DIR,
  OUTPUT_DIR,
  RDS_DIR,
  QC_DIR,
  DE_DIR,
  FIG_DIR,
  TAB_DIR
)

cat("Creating/checking project folders...\n")
for (d in dirs_to_create) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    cat("  Created:", d, "\n")
  } else {
    cat("  Exists: ", d, "\n")
  }
}
cat("\n")

# ============================================================
# SECTION 3: Required packages
# ============================================================
# We keep package installation here so the rest of the scripts
# can assume the required packages are available.
#
# Notes:
#   - Seurat is used for standard single-nucleus QC, normalization,
#     PCA, clustering, UMAP, and cluster marker inspection.
#   - edgeR + limma are used for donor-level pseudobulk DE.
#   - fgsea + msigdbr are used for pathway analysis.
#   - RankProd is used for rank-based cross-cohort meta-analysis.

required_cran <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "tibble",
  "readr",
  "readxl",
  "stringr",
  "patchwork",
  "pheatmap",
  "RColorBrewer",
  "ggrepel",
  "writexl",
  "Matrix",
  "R.utils",
  "hdf5r",
  "Seurat"
)

required_bioc <- c(
  "limma",
  "edgeR",
  "fgsea",
  "msigdbr",
  "RankProd"
)

cat("Checking CRAN packages...\n")
for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing:", pkg, "\n")
    install.packages(pkg, quiet = TRUE)
  } else {
    cat("  Already installed:", pkg, "\n")
  }
}

cat("\nChecking Bioconductor packages...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing:", pkg, "\n")
    BiocManager::install(pkg, ask = FALSE, quiet = TRUE)
  } else {
    cat("  Already installed:", pkg, "\n")
  }
}

cat("\nLoading packages...\n")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(readr)
  library(readxl)
  library(stringr)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(writexl)
  library(Matrix)
  library(R.utils)
  library(hdf5r)
  library(Seurat)
  library(limma)
  library(edgeR)
  library(fgsea)
  library(msigdbr)
  library(RankProd)
})
cat("All packages loaded successfully.\n\n")

# ============================================================
# SECTION 4: Dataset registry
# ============================================================
# This is the central file/folder registry for the restarted project.
# If a file name changes, update it here and rerun this setup script.
#
# Cell-type handling plan:
#   SEA_AD    : author-provided Microglia/Immune H5AD subset
#   GSE174367 : author cell labels, microglia = MG
#   GSE157827 : raw 10x files, needs cluster-based annotation
#   GSE160936 : raw 10x tar.gz files, needs cluster-based annotation
#   GSE188545 : raw 10x files, needs cluster-based annotation
#   GSE243292 : already microglia-specific/pre-annotated H5AD

DATASETS <- list(

  SEA_AD = list(
    name         = "SEA-AD",
    accession    = "Allen Brain Atlas / AWS",
    donors       = 84,
    brain_region = "MTG + DLPFC",
    role         = "Discovery",
    annotation   = "Author Microglia-PVM subset",
    h5ad_file    = file.path(DATA_DIR,
                             "SEA-AD_Microglia-and-Immune_multi-regional_final-nuclei_AAIC.h5ad"),
    meta_file    = file.path(DATA_DIR,
                             "Sea-ad_cohort_donor_metadata_072524.xlsx")
  ),

  GSE174367 = list(
    name         = "Morabito 2021",
    accession    = "GSE174367",
    donors       = 18,
    brain_region = "Prefrontal Cortex",
    role         = "Replication 1",
    annotation   = "Author cell-type labels; MG retained",
    h5_file      = file.path(DATA_DIR,
                             "GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5"),
    meta_file    = file.path(DATA_DIR,
                             "GSE174367_snRNA-seq_cell_meta.csv.gz")
  ),

  GSE157827 = list(
    name         = "Lau 2020",
    accession    = "GSE157827",
    donors       = 21,
    brain_region = "Prefrontal Cortex",
    role         = "Replication 2",
    annotation   = "Raw 10x; cluster-based microglia annotation required",
    raw_dir      = file.path(DATA_DIR, "GSE157827_RAW")
  ),

  GSE160936 = list(
    name         = "Smith 2022",
    accession    = "GSE160936",
    donors       = 12,
    samples      = 24,
    brain_region = "Entorhinal cortex + somatosensory cortex",
    role         = "Replication 3",
    annotation   = "Raw 10x tar.gz; cluster-based microglia annotation required",
    raw_dir      = file.path(DATA_DIR, "GSE160936_RAW"),
    series_meta  = file.path(DATA_DIR, "GSE160936_series_metadata_lines.txt"),
    parsed_meta  = file.path(DATA_DIR, "GSE160936_sample_metadata_parsed.csv")
  ),

  GSE188545 = list(
    name         = "Zhang 2023",
    accession    = "GSE188545",
    donors       = 12,
    brain_region = "Middle temporal gyrus",
    role         = "Replication 4",
    annotation   = "Raw 10x; cluster-based microglia annotation required",
    raw_dir      = file.path(DATA_DIR, "GSE188545_RAW")
  ),

  GSE243292 = list(
    name         = "Nguyen 2023",
    accession    = "GSE243292",
    donors       = 15,
    brain_region = "DLPFC",
    role         = "Replication 5",
    annotation   = "Microglia-specific/pre-annotated H5AD",
    h5ad_file    = file.path(DATA_DIR,
                             "GSE243292_Microglia_GEO.h5ad")
  )
)

# ============================================================
# SECTION 5: Global analysis parameters
# ============================================================
# These parameters should be kept explicit and stable.
# Do not tune thresholds to recover a desired gene list.
#
# For raw unannotated datasets, the new GEO script will use
# cluster-based annotation rather than per-cell marker thresholding.

PARAMS <- list(

  # ----------------------------
  # Differential expression
  # ----------------------------
  fdr_threshold    = 0.05,
  logfc_threshold  = 0.25,

  # ----------------------------
  # Pseudobulk donor QC
  # ----------------------------
  # This is a minimum donor inclusion floor after microglia extraction.
  # We will report donor-level cell counts. If many donors have very low
  # microglial nuclei, we may later add a sensitivity analysis with a
  # higher cutoff rather than silently changing the primary rule.
  min_cells_per_donor = 10,

  # ----------------------------
  # Single-nucleus QC defaults for raw GEO clustering
  # ----------------------------
  # The raw GEO script can override these with dataset-specific published
  # thresholds when available.
  #
  # These are not used for SEA-AD or author-annotated datasets.
  sn_min_features = 200,
  sn_max_features = Inf,
  sn_max_counts   = Inf,
  sn_max_mito_pct = 5,
  gene_min_cells  = 3,

  # ----------------------------
  # Seurat clustering defaults
  # ----------------------------
  seurat_nfeatures  = 3000,
  seurat_npcs       = 30,
  seurat_resolution = 0.5,

  # ----------------------------
  # Microglia marker panel for cluster-level annotation
  # ----------------------------
  # These markers are used to annotate clusters, not to select individual
  # cells by single-marker positivity.
  microglia_markers = c("P2RY12", "CX3CR1", "CSF1R", "TMEM119",
                        "AIF1", "CTSS", "C3", "APBB1IP", "SALL1"),

  # ----------------------------
  # Negative marker panels for excluding non-microglial clusters
  # ----------------------------
  neuron_markers = c("RBFOX3", "SNAP25", "SLC17A7", "GAD1", "GAD2"),
  astro_markers  = c("AQP4", "GFAP", "SLC1A2", "ALDH1L1"),
  olig_markers   = c("MBP", "MOG", "PLP1", "MOBP"),
  opc_markers    = c("PDGFRA", "VCAN", "CSPG4"),
  vascular_markers = c("CLDN5", "FLT1", "PECAM1", "COL1A1", "DCN"),
  macrophage_pvm_markers = c("CD163", "MRC1", "MSR1", "LYVE1"),

  # ----------------------------
  # Biological sanity-check genes
  # ----------------------------
  homeostatic_genes = c("P2RY12", "CX3CR1", "MEF2C",
                        "TMEM119", "CSF1R", "SALL1",
                        "SORL1", "LPAR6"),

  disease_genes = c("SPP1", "FKBP5", "SLC11A1",
                    "APOE", "TYROBP", "CD163",
                    "LPL", "GPNMB", "HSP90B1"),

  key_up_genes   = c("CD163", "FKBP5", "SLC11A1",
                     "SPP1", "LPL", "GPNMB"),

  key_down_genes = c("CX3CR1", "P2RY12", "MEF2C",
                     "SORL1", "SALL1", "IGF1"),

  # ----------------------------
  # Plot colors
  # ----------------------------
  ad_color      = "#E74C3C",
  control_color = "#3498DB",
  up_color      = "#E74C3C",
  down_color    = "#3498DB",
  ns_color      = "grey70",

  # ----------------------------
  # Figure dimensions
  # ----------------------------
  fig_width_single = 7,
  fig_width_double = 12,
  fig_height       = 6
)

# ============================================================
# SECTION 6: Verify all expected input files/folders
# ============================================================

cat("Verifying registered input files and folders...\n")

all_files_present <- TRUE

check_path <- function(label, path, required = TRUE) {
  # Safely handle undefined paths.
  if (is.null(path) || length(path) == 0 || is.na(path) || path == "") {
    msg <- paste0("  NOT DEFINED: ", label)
    if (required) {
      cat(msg, "\n")
      all_files_present <<- FALSE
    } else {
      cat("  Optional not defined:", label, "\n")
    }
    return(invisible(FALSE))
  }

  present <- file.exists(path) || dir.exists(path)

  if (present) {
    cat("  OK :", label, "\n")
  } else {
    cat("  MISSING:", label, "\n")
    cat("       ->", path, "\n")
    if (required) all_files_present <<- FALSE
  }

  return(invisible(present))
}

# SEA-AD
check_path("SEA-AD H5AD", DATASETS$SEA_AD$h5ad_file)
check_path("SEA-AD donor metadata", DATASETS$SEA_AD$meta_file)

# GSE174367
check_path("GSE174367 H5 count matrix", DATASETS$GSE174367$h5_file)
check_path("GSE174367 cell metadata", DATASETS$GSE174367$meta_file)

# Raw GEO datasets
check_path("GSE157827 RAW folder", DATASETS$GSE157827$raw_dir)
check_path("GSE160936 RAW folder", DATASETS$GSE160936$raw_dir)
check_path("GSE160936 series metadata", DATASETS$GSE160936$series_meta)
check_path("GSE160936 parsed metadata", DATASETS$GSE160936$parsed_meta, required = FALSE)
check_path("GSE188545 RAW folder", DATASETS$GSE188545$raw_dir)

# GSE243292
check_path("GSE243292 microglia H5AD", DATASETS$GSE243292$h5ad_file)

cat("\n")

if (all_files_present) {
  cat("All required input files/folders are present.\n\n")
} else {
  cat("Some required files/folders are missing. Fix those before analysis.\n\n")
}

# ============================================================
# SECTION 7: Save setup object
# ============================================================

project_setup <- list(
  ROOT_DIR    = ROOT_DIR,
  DATA_DIR    = DATA_DIR,
  SCRIPTS_DIR = SCRIPTS_DIR,
  OUTPUT_DIR  = OUTPUT_DIR,
  RDS_DIR     = RDS_DIR,
  QC_DIR      = QC_DIR,
  DE_DIR      = DE_DIR,
  FIG_DIR     = FIG_DIR,
  TAB_DIR     = TAB_DIR,
  DATASETS    = DATASETS,
  PARAMS      = PARAMS
)

setup_out <- file.path(RDS_DIR, "project_setup.rds")
saveRDS(project_setup, setup_out)

cat("Saved setup object to:\n")
cat("  ", setup_out, "\n\n")

cat("Quick check:\n")
cat("  GSE157827 raw_dir:", DATASETS$GSE157827$raw_dir, "\n")
cat("  GSE160936 raw_dir:", DATASETS$GSE160936$raw_dir, "\n")
cat("  GSE188545 raw_dir:", DATASETS$GSE188545$raw_dir, "\n\n")

cat("Setup complete.\n")
cat("Next scripts:\n")
cat("  02_load_SEAAD_AD61026.R has already run successfully.\n")
cat("  03B_parse_GSE160936_metadata_AD61026_fixed.R should be rerun because the first parser misread age.\n")
cat("  Then run the new 03 GEO processing script once written.\n")
cat("============================================================\n")
