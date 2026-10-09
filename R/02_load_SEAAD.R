# ============================================================
# BLOCK 2: Load and Process SEA-AD Discovery Dataset
# ============================================================
# Dataset : Seattle Alzheimer's Disease Brain Cell Atlas
# Role    : Discovery dataset (84 donors, 42 AD / 42 Control)
# File    : SEA-AD_Microglia-and-Immune_multi-regional_final-nuclei_AAIC.h5ad
#
# Key findings from initial data exploration:
#   - All 240,651 cells are already microglia (Subclass = Microglia-PVM)
#   - Diagnosis column: "Cognitive.Status" (Dementia / No dementia)
#   - Donor column: "Donor.ID" (84 unique donors)
#   - H5AD categorical columns use signed 8-bit integer codes
#     with spaces converted to dots by as.data.frame()
#   - All 84 donors pass QC (min 512 cells per donor)
#
# Output:
#   RDS/SEAAD_counts_raw.rds          raw sparse matrix
#   RDS/SEAAD_meta_raw.rds            raw cell metadata
#   RDS/SEAAD_pseudobulk_counts.rds   pseudobulk matrix
#   RDS/SEAAD_pseudobulk_meta.rds     donor metadata
#   QC/SEAAD_cells_per_donor.pdf      QC plot
#
# Runtime: ~5-10 minutes (3GB file)
# IMPORTANT: Run 01_setup.R before this script
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

setup <- readRDS(file.path(PROJECT_ROOT, "RDS/project_setup.rds"))
list2env(setup, envir = .GlobalEnv)

library(hdf5r)
library(Matrix)
library(dplyr)
library(ggplot2)
library(readxl)

cat("=================================================\n")
cat("  BLOCK 2: SEA-AD Discovery Dataset\n")
cat("=================================================\n")
cat("Expected: 84 donors, 42 AD / 42 Control\n\n")

# ============================================================
# SECTION 2: Read H5AD file
# ============================================================
# H5AD is the AnnData format from Python single-cell tools.
# We read it directly using hdf5r without needing Python.
#
# IMPORTANT: Categorical columns in this file use signed 8-bit
# integer codes stored with spaces in column names. When
# as.data.frame() is called, spaces become dots e.g.:
#   "Cognitive Status" -> "Cognitive.Status"
#   "Donor ID"         -> "Donor.ID"
#   "Brain Region"     -> "Brain.Region"

h5ad_path <- DATASETS$SEA_AD$h5ad_file

if (!file.exists(h5ad_path)) {
  stop("SEA-AD H5AD file not found at:\n", h5ad_path)
}

cat("File:", basename(h5ad_path), "\n")
cat("Size:", round(file.size(h5ad_path) / 1e9, 2), "GB\n")
cat("Opening file — please wait 5-10 minutes...\n\n")

h5 <- H5File$new(h5ad_path, mode = "r")

# ---- 2a: Read count matrix ----
# Stored as CSR sparse matrix with three arrays:
#   data    = non-zero count values
#   indices = gene index for each non-zero value (0-based)
#   indptr  = marks start of each cell's data

cat("Reading count matrix...\n")
X         <- h5[["X"]]
shape     <- h5[["X"]]$attr_open("shape")$read()
n_cells   <- shape[1]
n_genes   <- shape[2]

counts_raw <- sparseMatrix(
  i    = X[["indices"]][] + 1L,   # 0-based to 1-based
  p    = X[["indptr"]][],
  x    = as.numeric(X[["data"]][]),
  dims = c(n_genes, n_cells),      # genes x cells
  repr = "C"
)

# Assign gene symbols from var/index
rownames(counts_raw) <- h5[["var"]][["index"]][]
cat("  Matrix:", nrow(counts_raw), "genes x",
    ncol(counts_raw), "cells\n\n")

# ---- 2b: Read cell metadata ----
# Categorical columns are stored as groups with two subkeys:
#   categories = array of unique labels
#   codes      = integer index per cell pointing to categories
#
# CRITICAL: codes are signed 8-bit integers in this file.
# Values above 127 wrap to negative. We fix this by adding
# 256 to any negative code values before indexing.

cat("Reading cell metadata...\n")
obs      <- h5[["obs"]]
obs_list <- list()

for (col in names(obs)) {
  tryCatch({
    item <- obs[[col]]
    
    if (inherits(item, "H5Group")) {
      sub_names <- names(item)
      
      if (all(c("categories", "codes") %in% sub_names)) {
        # Standard categorical column
        categories <- item[["categories"]][]
        codes_raw  <- item[["codes"]][]
        
        # Fix signed 8-bit integer overflow
        codes_int <- as.integer(codes_raw)
        codes_int[codes_int < 0] <- codes_int[codes_int < 0] + 256L
        
        # Convert to 1-based indexing
        codes_1based <- codes_int + 1L
        codes_1based[codes_1based < 1 |
                       codes_1based > length(categories)] <- NA
        
        obs_list[[col]] <- categories[codes_1based]
        
      } else if ("Latino" %in% sub_names) {
        # Hispanic column has non-standard structure
        tryCatch({
          v <- item[["Latino"]][]
          if (length(v) == n_cells) obs_list[[col]] <- v
        }, error = function(e) {})
        
      } else {
        # Other group: try first subkey
        tryCatch({
          v <- item[[sub_names[1]]][]
          if (length(v) == n_cells) obs_list[[col]] <- v
        }, error = function(e) {})
      }
      
    } else {
      # Plain array column
      v <- item[]
      if (length(v) == n_cells) obs_list[[col]] <- v
    }
    
  }, error = function(e) {})
}

# Read cell barcodes
tryCatch({
  obs_list[["cell_id"]] <- obs[["_index"]][]
}, error = function(e) {})

# Close file — everything needed is now in memory
h5$close_all()

# Build metadata data frame
# NOTE: as.data.frame() converts spaces to dots in column names
# e.g. "Donor ID" becomes "Donor.ID"
meta_raw <- as.data.frame(obs_list, stringsAsFactors = FALSE)

cat("  Metadata:", nrow(meta_raw), "rows x",
    ncol(meta_raw), "columns\n\n")

# Save raw objects immediately
# These are saved so we never need to re-read the 3GB file
saveRDS(counts_raw, file.path(RDS_DIR, "SEAAD_counts_raw.rds"))
saveRDS(meta_raw,   file.path(RDS_DIR, "SEAAD_meta_raw.rds"))
cat("Raw objects saved to RDS/\n\n")

# ============================================================
# SECTION 3: Define column names
# ============================================================
# Column names use dots instead of spaces due to as.data.frame()
# Diagnosis labels confirmed from data inspection

DONOR_COL     <- "Donor.ID"
DX_COL        <- "Cognitive.Status"
CELLTYPE_COL  <- "Subclass"        # all cells = Microglia-PVM
REGION_COL    <- "Brain.Region"
BRAAK_COL     <- "Braak"
APOE_COL      <- "APOE.Genotype"

AD_LABEL      <- "Dementia"        # AD donors
CONTROL_LABEL <- "No dementia"     # control donors

# ============================================================
# SECTION 4: Verify key columns
# ============================================================

cat("=== KEY COLUMN VERIFICATION ===\n\n")

cat("Cognitive Status:\n")
print(table(meta_raw[[DX_COL]], useNA = "ifany"))

cat("\nSubclass (cell type):\n")
print(table(meta_raw[[CELLTYPE_COL]], useNA = "ifany"))

cat("\nBrain Region:\n")
print(table(meta_raw[[REGION_COL]], useNA = "ifany"))

cat("\nUnique donors:", length(unique(meta_raw[[DONOR_COL]])), "\n")

cat("\nBraak stage:\n")
print(table(meta_raw[[BRAAK_COL]], useNA = "ifany"))

cat("\nAPOE Genotype:\n")
print(table(meta_raw[[APOE_COL]], useNA = "ifany"))

# ============================================================
# SECTION 5: Filter to AD vs Control
# ============================================================
# All cells in this file are already microglia (Microglia-PVM)
# so no cell type filtering is needed.
# We only need to filter by diagnosis.

cat("\n=== FILTERING ===\n\n")

keep_dx   <- meta_raw[[DX_COL]] %in% c(AD_LABEL, CONTROL_LABEL)
counts_mg <- counts_raw[, keep_dx]
meta_mg   <- meta_raw[keep_dx, ]

cat("Cells after diagnosis filter:", ncol(counts_mg), "\n")
cat("Donors:", length(unique(meta_mg[[DONOR_COL]])), "\n")
print(table(meta_mg[[DX_COL]]))

# ============================================================
# SECTION 6: Quality control — cells per donor
# ============================================================
# Remove donors with fewer than min_cells_per_donor microglial
# nuclei to avoid noisy pseudobulk profiles.
# All 84 SEA-AD donors pass this threshold (min = 512 cells).

cat("\n=== QUALITY CONTROL ===\n\n")

cells_per_donor <- table(meta_mg[[DONOR_COL]])
cat("Cells per donor summary:\n")
print(summary(as.numeric(cells_per_donor)))

donors_keep <- names(cells_per_donor)[
  cells_per_donor >= PARAMS$min_cells_per_donor
]
donors_fail <- names(cells_per_donor)[
  cells_per_donor <  PARAMS$min_cells_per_donor
]

cat("Donors passing QC:", length(donors_keep), "\n")
cat("Donors removed   :", length(donors_fail), "\n")

keep_donors <- meta_mg[[DONOR_COL]] %in% donors_keep
counts_mg   <- counts_mg[, keep_donors]
meta_mg     <- meta_mg[keep_donors, ]

# QC plot: cells per donor
qc_df <- data.frame(
  donor   = names(cells_per_donor),
  n_cells = as.numeric(cells_per_donor),
  stringsAsFactors = FALSE
)
dx_lookup <- meta_raw[!duplicated(meta_raw[[DONOR_COL]]),
                      c(DONOR_COL, DX_COL)]
colnames(dx_lookup) <- c("donor", "diagnosis")
qc_df <- merge(qc_df, dx_lookup, by = "donor", all.x = TRUE)
qc_df$pass_qc <- qc_df$donor %in% donors_keep

ggplot(qc_df, aes(x     = reorder(donor, n_cells),
                  y     = n_cells,
                  fill  = diagnosis,
                  alpha = pass_qc)) +
  geom_col() +
  geom_hline(yintercept = PARAMS$min_cells_per_donor,
             linetype   = "dashed",
             color      = "black") +
  scale_fill_manual(
    values   = c("Dementia"    = PARAMS$ad_color,
                 "No dementia" = PARAMS$control_color),
    na.value = "grey70"
  ) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3),
                     guide  = "none") +
  coord_flip() +
  labs(title = "SEA-AD: Microglial nuclei per donor",
       x     = "Donor",
       y     = "Nuclei",
       fill  = "Diagnosis") +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_text(size = 7))

ggsave(file.path(QC_DIR, "SEAAD_cells_per_donor.pdf"),
       width = 8, height = 10)
cat("QC plot saved\n\n")

# ============================================================
# SECTION 7: Pseudobulk aggregation
# ============================================================
# Sum raw counts across all cells from the same donor.
# This produces one gene expression profile per donor,
# enabling proper donor-level statistical modeling with
# edgeR/limma-voom without pseudoreplication.

cat("=== PSEUDOBULK AGGREGATION ===\n\n")

donors  <- unique(meta_mg[[DONOR_COL]])
pb_list <- lapply(donors, function(d) {
  idx <- meta_mg[[DONOR_COL]] == d
  Matrix::rowSums(counts_mg[, idx, drop = FALSE])
})

counts_pb <- do.call(cbind, pb_list)
colnames(counts_pb) <- donors
counts_pb <- as.matrix(counts_pb)

cat("Pseudobulk matrix:", nrow(counts_pb), "genes x",
    ncol(counts_pb), "donors\n\n")

# ============================================================
# SECTION 8: Build donor metadata
# ============================================================

meta_pb <- meta_mg[!duplicated(meta_mg[[DONOR_COL]]), ]
meta_pb <- meta_pb[match(donors, meta_pb[[DONOR_COL]]), ]
rownames(meta_pb) <- donors

# Standardized diagnosis: "AD" and "Control"
# Control is set as reference level for DE modeling
meta_pb$diagnosis_std <- factor(
  ifelse(meta_pb[[DX_COL]] == AD_LABEL, "AD", "Control"),
  levels = c("Control", "AD")
)

cat("Donor breakdown:\n")
print(table(meta_pb$diagnosis_std))

# Merge richer Excel metadata (Braak, APOE, age, sex etc.)
donor_excel     <- read_excel(DATASETS$SEA_AD$meta_file)
donor_col_excel <- "Donor ID"   # confirmed column name in Excel

meta_pb <- merge(
  meta_pb, donor_excel,
  by.x  = DONOR_COL,
  by.y  = donor_col_excel,
  all.x = TRUE
)
rownames(meta_pb) <- meta_pb[[DONOR_COL]]
meta_pb <- meta_pb[match(colnames(counts_pb),
                         meta_pb[[DONOR_COL]]), ]

cat("Merged with Excel metadata:", nrow(meta_pb), "donors\n\n")

# ============================================================
# SECTION 9: Sanity check — key marker genes
# ============================================================
# Verify expected marker genes are present and expressed.
# Homeostatic markers should be detected (P2RY12, CX3CR1).
# Disease markers should also be present (SPP1, CD163).

cat("=== MARKER GENE CHECK ===\n\n")

all_markers <- c(PARAMS$key_down_genes, PARAMS$key_up_genes)
present     <- all_markers[all_markers %in% rownames(counts_pb)]
missing     <- all_markers[!all_markers %in% rownames(counts_pb)]

cat("Mean pseudobulk counts (all donors):\n")
if (length(present) > 0) {
  means <- rowMeans(counts_pb[present, ])
  print(round(sort(means, decreasing = TRUE), 1))
}
if (length(missing) > 0) {
  cat("Not in dataset:", paste(missing, collapse = ", "), "\n")
}

# ============================================================
# SECTION 10: Save outputs
# ============================================================

saveRDS(counts_pb, file.path(RDS_DIR, "SEAAD_pseudobulk_counts.rds"))
saveRDS(meta_pb,   file.path(RDS_DIR, "SEAAD_pseudobulk_meta.rds"))

cat("\n=================================================\n")
cat("  BLOCK 2 COMPLETE\n")
cat("=================================================\n")
cat("Donors :", ncol(counts_pb), "\n")
cat("AD     :", sum(meta_pb$diagnosis_std == "AD"), "\n")
cat("Control:", sum(meta_pb$diagnosis_std == "Control"), "\n")
cat("Genes  :", nrow(counts_pb), "\n")
cat("Saved  : RDS/SEAAD_pseudobulk_counts.rds\n")
cat("         RDS/SEAAD_pseudobulk_meta.rds\n")
cat("         QC/SEAAD_cells_per_donor.pdf\n")
cat("NEXT   : Run Scripts/03_load_GEO_datasets.R\n")
cat("=================================================\n")