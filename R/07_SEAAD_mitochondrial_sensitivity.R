# ============================================================
# 07_SEAAD_mitochondrial_sensitivity.R
# ============================================================
# Project: AD Microglia Multi-Dataset Transcriptomic Analysis
#
# Purpose:
#   Run a mitochondrial fraction sensitivity analysis for the rebuilt
#   AD61026 SEA-AD discovery model using a fair, directly comparable
#   gene universe.
#
# Why this v2 script is preferred:
#   The first MT sensitivity script correctly tested:
#
#      expression ~ mt_fraction + diagnosis
#
#   However, it allowed edgeR::filterByExpr() to re-filter genes using
#   the adjusted design. That produced a different adjusted gene universe
#   than the original SEA-AD DE table.
#
#   For sensitivity analysis, the cleanest comparison is:
#
#      original SEA-AD DE genes tested
#      versus
#      the same genes tested again after MT-fraction adjustment
#
#   Therefore this script restricts the adjusted model to the same gene
#   universe already present in:
#
#      Output/DE/SEAAD_DE_results.csv
#
# Main model:
#   Original model, already run in 04:
#      expression ~ diagnosis
#
#   Sensitivity model, run here:
#      expression ~ mt_fraction + diagnosis
#
# Main interpretation:
#   This is not meant to replace the primary diagnosis model.
#   It asks whether the direction and magnitude of the diagnosis signal
#   are stable after adjusting for donor-level mitochondrial fraction.
#
# Key outputs:
#   Output/DE/SEAAD_mito_sensitivity_summary_AD61026_v2.csv
#   Output/DE/SEAAD_DE_results_mito_adjusted_AD61026_v2.csv
#   Output/DE/SEAAD_mito_sensitivity_discovery_genes_AD61026_v2.csv
#   Output/DE/SEAAD_mito_sensitivity_strict48_AD61026_v2.csv
#   Output/DE/SEAAD_mito_sensitivity_expanded125_AD61026_v2.csv
#   Output/DE/SEAAD_mito_fraction_diagnosis_diagnostics_AD61026_v2.csv
#
# Recommended run command:
#   source("R/04B_SEAAD_mito_sensitivity_AD61026_v2_FAIR_UNIVERSE.R")
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

required_pkgs <- c(
  "edgeR", "limma", "Matrix",
  "dplyr", "tibble", "tidyr",
  "ggplot2", "ggrepel"
)

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall the missing package(s), restart R, and rerun this script."
  )
}

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})

# ------------------------------------------------------------
# Define output directories robustly.
# ------------------------------------------------------------
if (!exists("OUT_DIR")) {
  if (exists("TAB_DIR")) {
    OUT_DIR <- dirname(TAB_DIR)
  } else if (exists("DE_DIR")) {
    OUT_DIR <- dirname(DE_DIR)
  } else if (exists("FIG_DIR")) {
    OUT_DIR <- dirname(FIG_DIR)
  } else {
    OUT_DIR <- file.path(PROJECT_ROOT, "Output")
  }
}

if (!exists("DE_DIR")) {
  DE_DIR <- file.path(OUT_DIR, "DE")
}
if (!exists("FIG_DIR")) {
  FIG_DIR <- file.path(OUT_DIR, "Figures")
}
if (!exists("RDS_DIR")) {
  RDS_DIR <- file.path(PROJECT_ROOT, "RDS")
}

SOURCE_DIR <- file.path(OUT_DIR, "Manuscript_Source_Data")

dir.create(DE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("  04B v2: SEA-AD MT sensitivity, fair gene universe\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Use thresholds from setup if available.
# ------------------------------------------------------------
if (exists("PARAMS")) {
  FDR_THRESHOLD <- PARAMS$fdr_threshold
  LOGFC_THRESHOLD <- PARAMS$logfc_threshold
} else {
  FDR_THRESHOLD <- 0.05
  LOGFC_THRESHOLD <- 0.25
}

cat("Discovery threshold:\n")
cat("  FDR <", FDR_THRESHOLD, "\n")
cat("  |logFC| >", LOGFC_THRESHOLD, "\n\n")

# ============================================================
# SECTION 2: Helper functions
# ============================================================

read_required_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Required file missing: ", path)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

make_direction <- function(logFC) {
  dplyr::case_when(
    is.na(logFC) ~ NA_character_,
    logFC > 0 ~ "UP_in_AD",
    logFC < 0 ~ "DOWN_in_AD",
    TRUE ~ "NO_CHANGE"
  )
}

safe_percent_attenuation <- function(logFC_orig, logFC_adj) {
  ifelse(
    is.na(logFC_orig) | is.na(logFC_adj) | abs(logFC_orig) == 0,
    NA_real_,
    100 * (abs(logFC_orig) - abs(logFC_adj)) / abs(logFC_orig)
  )
}

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

safe_pearson <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "pearson"))
}

# ============================================================
# SECTION 3: Load SEA-AD pseudobulk and original DE results
# ============================================================

cat("Loading SEA-AD pseudobulk counts, metadata, and original DE...\n\n")

counts_path <- file.path(RDS_DIR, "SEAAD_pseudobulk_counts.rds")
meta_path   <- file.path(RDS_DIR, "SEAAD_pseudobulk_meta.rds")
de_path     <- file.path(DE_DIR, "SEAAD_DE_results.csv")

if (!file.exists(counts_path)) stop("Missing counts file: ", counts_path)
if (!file.exists(meta_path)) stop("Missing metadata file: ", meta_path)
if (!file.exists(de_path)) stop("Missing original DE file: ", de_path)

counts <- readRDS(counts_path)
meta <- readRDS(meta_path)
original_de <- read_required_csv(de_path)

if (inherits(counts, "sparseMatrix")) {
  counts <- as.matrix(counts)
}

# ------------------------------------------------------------
# Align metadata rows to count matrix columns.
# ------------------------------------------------------------
if (!identical(rownames(meta), colnames(counts))) {

  if ("pseudobulk_id" %in% colnames(meta) &&
      all(colnames(counts) %in% as.character(meta$pseudobulk_id))) {

    meta <- meta[match(colnames(counts), as.character(meta$pseudobulk_id)), , drop = FALSE]
    rownames(meta) <- meta$pseudobulk_id

  } else if ("donor_id" %in% colnames(meta) &&
             all(colnames(counts) %in% as.character(meta$donor_id))) {

    meta <- meta[match(colnames(counts), as.character(meta$donor_id)), , drop = FALSE]
    rownames(meta) <- colnames(counts)

  } else {
    stop("Could not align SEA-AD count columns and metadata rows.")
  }
}

if (!identical(rownames(meta), colnames(counts))) {
  stop("SEA-AD metadata and count matrix alignment failed.")
}

required_de_cols <- c("gene", "logFC", "P.Value", "adj.P.Val")
missing_de_cols <- setdiff(required_de_cols, colnames(original_de))

if (length(missing_de_cols) > 0) {
  stop("SEAAD_DE_results.csv is missing required columns: ",
       paste(missing_de_cols, collapse = ", "))
}

cat("SEA-AD pseudobulk dimensions:\n")
cat("  count genes:", nrow(counts), "\n")
cat("  donors:", ncol(counts), "\n")
cat("Original SEA-AD DE tested genes:", nrow(original_de), "\n\n")

# ============================================================
# SECTION 4: Calculate mitochondrial fraction
# ============================================================

cat("Calculating donor-level mitochondrial fraction...\n\n")

mt_genes <- rownames(counts)[grepl("^MT-", rownames(counts))]

if (length(mt_genes) == 0) {
  stop("No genes beginning with 'MT-' were found in the count matrix.")
}

total_counts <- colSums(counts)
mt_counts <- colSums(counts[mt_genes, , drop = FALSE])

meta$mt_counts <- as.numeric(mt_counts[rownames(meta)])
meta$total_counts <- as.numeric(total_counts[rownames(meta)])
meta$mt_fraction <- meta$mt_counts / meta$total_counts
meta$mt_percent <- 100 * meta$mt_fraction

cat("MT genes detected in pseudobulk count matrix:", length(mt_genes), "\n")
cat("MT fraction summary:\n")
print(summary(meta$mt_fraction))
cat("\n")

write.csv(
  meta %>%
    rownames_to_column("sample_id") %>%
    select(sample_id, everything()),
  file.path(DE_DIR, "SEAAD_mito_fraction_by_donor_AD61026_v2.csv"),
  row.names = FALSE
)

# ============================================================
# SECTION 5: Diagnosis and MT-fraction diagnostics
# ============================================================

cat("Checking whether MT fraction is associated with diagnosis...\n\n")

dx_col <- "diagnosis_std"

if (!dx_col %in% colnames(meta)) {
  stop("diagnosis_std column not found in SEA-AD metadata.")
}

dx <- factor(as.character(meta[[dx_col]]), levels = c("Control", "AD"))

if (any(is.na(dx))) {
  stop("Missing or unsupported diagnosis values in diagnosis_std.")
}

meta$dx_for_model <- dx
meta$dx_numeric <- ifelse(dx == "AD", 1, 0)

mt_group_summary <- meta %>%
  as.data.frame() %>%
  group_by(dx_for_model) %>%
  summarise(
    n = n(),
    mean_mt_fraction = mean(mt_fraction, na.rm = TRUE),
    median_mt_fraction = median(mt_fraction, na.rm = TRUE),
    sd_mt_fraction = sd(mt_fraction, na.rm = TRUE),
    min_mt_fraction = min(mt_fraction, na.rm = TRUE),
    max_mt_fraction = max(mt_fraction, na.rm = TRUE),
    .groups = "drop"
  )

mt_wilcox_p <- tryCatch(
  wilcox.test(mt_fraction ~ dx_for_model, data = meta)$p.value,
  error = function(e) NA_real_
)

mt_lm <- lm(dx_numeric ~ mt_fraction, data = meta)
mt_r2 <- summary(mt_lm)$r.squared
mt_vif_like <- 1 / (1 - mt_r2)

mt_diagnostics <- data.frame(
  metric = c(
    "MT genes detected",
    "Control donors",
    "AD donors",
    "Wilcoxon p for MT fraction by diagnosis",
    "R2 of diagnosis_numeric ~ mt_fraction",
    "VIF_like_for_diagnosis_due_to_mt_fraction"
  ),
  value = c(
    length(mt_genes),
    sum(dx == "Control"),
    sum(dx == "AD"),
    mt_wilcox_p,
    mt_r2,
    mt_vif_like
  ),
  stringsAsFactors = FALSE
)

write.csv(
  mt_group_summary,
  file.path(DE_DIR, "SEAAD_mito_fraction_by_diagnosis_summary_AD61026_v2.csv"),
  row.names = FALSE
)

write.csv(
  mt_diagnostics,
  file.path(DE_DIR, "SEAAD_mito_fraction_diagnosis_diagnostics_AD61026_v2.csv"),
  row.names = FALSE
)

cat("MT fraction diagnosis diagnostics:\n")
print(mt_diagnostics)
cat("\n")

# ============================================================
# SECTION 6: Load current strict and expanded signatures
# ============================================================

cat("Loading current strict and expanded signatures...\n\n")

strict_signature <- read_required_csv(
  file.path(SOURCE_DIR, "MainTable2_strict_consensus_signature_source.csv")
)

expanded_signature <- read_required_csv(
  file.path(SOURCE_DIR, "SuppTable6_expanded_drug_signature_source.csv")
)

if (!"gene" %in% colnames(strict_signature)) {
  stop("Strict signature source does not contain a gene column.")
}
if (!"gene" %in% colnames(expanded_signature)) {
  stop("Expanded signature source does not contain a gene column.")
}

strict_genes <- unique(strict_signature$gene)
expanded_genes <- unique(expanded_signature$gene)

strict_mt_genes <- strict_genes[grepl("^MT-", strict_genes)]
expanded_mt_genes <- expanded_genes[grepl("^MT-", expanded_genes)]

cat("Strict signature genes:", length(strict_genes), "\n")
cat("Expanded signature genes:", length(expanded_genes), "\n")
cat("Strict MT genes:", ifelse(length(strict_mt_genes) == 0, "none", paste(strict_mt_genes, collapse = ", ")), "\n")
cat("Expanded MT genes:", ifelse(length(expanded_mt_genes) == 0, "none", paste(expanded_mt_genes, collapse = ", ")), "\n\n")

# ============================================================
# SECTION 7: Run mito-adjusted model using the original gene universe
# ============================================================

cat("Running MT-adjusted model on the original SEA-AD tested gene universe...\n\n")

original_gene_universe <- unique(original_de$gene)
genes_for_adjusted_model <- intersect(original_gene_universe, rownames(counts))

missing_original_genes <- setdiff(original_gene_universe, rownames(counts))

if (length(missing_original_genes) > 0) {
  warning(
    length(missing_original_genes),
    " original DE genes were not found in the count matrix and will be omitted."
  )
}

counts_for_model <- counts[genes_for_adjusted_model, , drop = FALSE]

# DGEList from the SAME gene universe as original DE output.
dge <- edgeR::DGEList(counts = counts_for_model, samples = meta, group = dx)

# TMM normalization.
# edgeR now also uses normLibSizes(), but calcNormFactors() remains commonly used.
dge <- edgeR::calcNormFactors(dge, method = "TMM")

# Adjusted model.
# dxAD tests AD versus Control after adjustment for donor MT fraction.
design_adj <- model.matrix(~ mt_fraction + dx, data = meta)

cat("Adjusted design columns:\n")
print(colnames(design_adj))
cat("Design matrix condition number:", kappa(design_adj), "\n\n")

v <- limma::voom(dge, design_adj, plot = FALSE)

fit <- limma::lmFit(v, design_adj)
fit <- limma::eBayes(fit, robust = TRUE)

if (!"dxAD" %in% colnames(fit$coefficients)) {
  stop("Could not find dxAD coefficient in adjusted model.")
}

adjusted_de <- limma::topTable(
  fit,
  coef = "dxAD",
  number = Inf,
  sort.by = "none",
  adjust.method = "BH"
) %>%
  rownames_to_column("gene") %>%
  mutate(
    dataset = "SEAAD",
    model = "diagnosis_plus_mt_fraction",
    direction = make_direction(logFC),
    significant_discovery_threshold =
      adj.P.Val < FDR_THRESHOLD &
      abs(logFC) > LOGFC_THRESHOLD
  )

# Reorder adjusted DE to match original DE order.
adjusted_de <- adjusted_de %>%
  right_join(
    data.frame(gene = original_gene_universe, original_order = seq_along(original_gene_universe)),
    by = "gene"
  ) %>%
  arrange(original_order) %>%
  select(-original_order)

write.csv(
  adjusted_de,
  file.path(DE_DIR, "SEAAD_DE_results_mito_adjusted_AD61026_v2.csv"),
  row.names = FALSE
)

saveRDS(
  adjusted_de,
  file.path(RDS_DIR, "SEAAD_DE_results_mito_adjusted_AD61026_v2.rds")
)

cat("MT-adjusted DE complete.\n")
cat("  Adjusted genes tested:", nrow(adjusted_de), "\n")
cat("  Significant UP in AD:",
    sum(adjusted_de$adj.P.Val < FDR_THRESHOLD &
          adjusted_de$logFC > LOGFC_THRESHOLD, na.rm = TRUE), "\n")
cat("  Significant DOWN in AD:",
    sum(adjusted_de$adj.P.Val < FDR_THRESHOLD &
          adjusted_de$logFC < -LOGFC_THRESHOLD, na.rm = TRUE), "\n\n")

# ============================================================
# SECTION 8: Compare original and adjusted results
# ============================================================

cat("Comparing original and MT-adjusted results...\n\n")

comparison_all <- original_de %>%
  select(
    gene,
    logFC_orig = logFC,
    P.Value_orig = P.Value,
    FDR_orig = adj.P.Val
  ) %>%
  left_join(
    adjusted_de %>%
      select(
        gene,
        logFC_adj = logFC,
        P.Value_adj = P.Value,
        FDR_adj = adj.P.Val
      ),
    by = "gene"
  ) %>%
  mutate(
    direction_orig = make_direction(logFC_orig),
    direction_adj = make_direction(logFC_adj),
    same_direction = direction_orig == direction_adj,
    original_discovery_gene =
      FDR_orig < FDR_THRESHOLD &
      abs(logFC_orig) > LOGFC_THRESHOLD,
    retained_after_mito_adjustment_FDR =
      FDR_adj < FDR_THRESHOLD &
      abs(logFC_adj) > LOGFC_THRESHOLD,
    retained_after_mito_adjustment_nominal =
      P.Value_adj < 0.05 &
      abs(logFC_adj) > LOGFC_THRESHOLD,
    abs_adj_logFC_at_least_50pct_orig =
      abs(logFC_adj) >= 0.50 * abs(logFC_orig),
    abs_adj_logFC_at_least_75pct_orig =
      abs(logFC_adj) >= 0.75 * abs(logFC_orig),
    delta_logFC = logFC_adj - logFC_orig,
    percent_attenuation = safe_percent_attenuation(logFC_orig, logFC_adj),
    mt_gene = grepl("^MT-", gene),
    strict_primary_consensus = gene %in% strict_genes,
    expanded_signature = gene %in% expanded_genes
  )

discovery_sensitivity <- comparison_all %>%
  filter(original_discovery_gene) %>%
  arrange(FDR_orig)

strict_sensitivity <- comparison_all %>%
  filter(strict_primary_consensus) %>%
  arrange(FDR_orig)

expanded_sensitivity <- comparison_all %>%
  filter(expanded_signature) %>%
  arrange(FDR_orig)

mt_gene_sensitivity <- comparison_all %>%
  filter(mt_gene) %>%
  arrange(FDR_orig)

write.csv(
  comparison_all,
  file.path(DE_DIR, "SEAAD_mito_sensitivity_all_tested_genes_AD61026_v2.csv"),
  row.names = FALSE
)

write.csv(
  discovery_sensitivity,
  file.path(DE_DIR, "SEAAD_mito_sensitivity_discovery_genes_AD61026_v2.csv"),
  row.names = FALSE
)

write.csv(
  strict_sensitivity,
  file.path(DE_DIR, "SEAAD_mito_sensitivity_strict48_AD61026_v2.csv"),
  row.names = FALSE
)

write.csv(
  expanded_sensitivity,
  file.path(DE_DIR, "SEAAD_mito_sensitivity_expanded125_AD61026_v2.csv"),
  row.names = FALSE
)

write.csv(
  mt_gene_sensitivity,
  file.path(DE_DIR, "SEAAD_mito_sensitivity_all_MT_genes_AD61026_v2.csv"),
  row.names = FALSE
)

# ============================================================
# SECTION 9: Summary table
# ============================================================

summary_table <- data.frame(
  metric = c(
    "SEA-AD genes tested in original DE",
    "SEA-AD genes tested in MT-adjusted DE",
    "MT genes detected in pseudobulk counts",
    "Wilcoxon p for MT fraction by diagnosis",
    "R2 of diagnosis_numeric ~ mt_fraction",
    "VIF_like_for_diagnosis_due_to_mt_fraction",
    "Original SEA-AD discovery genes",
    "Original discovery genes same direction after MT adjustment",
    "Original discovery genes retained after MT adjustment, FDR<0.05 and |logFC|>0.25",
    "Original discovery genes retained after MT adjustment, nominal P<0.05 and |logFC|>0.25",
    "Original discovery genes with adjusted |logFC| >= 50 percent of original",
    "Original discovery genes with adjusted |logFC| >= 75 percent of original",
    "Spearman logFC correlation, all tested genes",
    "Pearson logFC correlation, all tested genes",
    "Spearman logFC correlation, original discovery genes",
    "Pearson logFC correlation, original discovery genes",
    "Strict 48-gene signature genes",
    "Strict signature MT genes",
    "Strict signature genes same direction after MT adjustment",
    "Strict signature genes retained after MT adjustment, FDR<0.05 and |logFC|>0.25",
    "Strict signature genes retained after MT adjustment, nominal P<0.05 and |logFC|>0.25",
    "Strict signature genes with adjusted |logFC| >= 50 percent of original",
    "Strict signature genes with adjusted |logFC| >= 75 percent of original",
    "Expanded 125-gene signature genes",
    "Expanded signature MT genes",
    "Expanded signature genes same direction after MT adjustment",
    "Expanded signature genes retained after MT adjustment, FDR<0.05 and |logFC|>0.25",
    "Expanded signature genes retained after MT adjustment, nominal P<0.05 and |logFC|>0.25",
    "Expanded signature genes with adjusted |logFC| >= 50 percent of original",
    "Expanded signature genes with adjusted |logFC| >= 75 percent of original",
    "All MT genes significant in original SEA-AD DE",
    "All MT genes retained after MT adjustment, FDR<0.05 and |logFC|>0.25",
    "All MT genes same direction after MT adjustment"
  ),
  value = c(
    nrow(original_de),
    nrow(adjusted_de),
    length(mt_genes),
    mt_wilcox_p,
    mt_r2,
    mt_vif_like,
    nrow(discovery_sensitivity),
    sum(discovery_sensitivity$same_direction, na.rm = TRUE),
    sum(discovery_sensitivity$retained_after_mito_adjustment_FDR, na.rm = TRUE),
    sum(discovery_sensitivity$retained_after_mito_adjustment_nominal, na.rm = TRUE),
    sum(discovery_sensitivity$abs_adj_logFC_at_least_50pct_orig, na.rm = TRUE),
    sum(discovery_sensitivity$abs_adj_logFC_at_least_75pct_orig, na.rm = TRUE),
    safe_spearman(comparison_all$logFC_orig, comparison_all$logFC_adj),
    safe_pearson(comparison_all$logFC_orig, comparison_all$logFC_adj),
    safe_spearman(discovery_sensitivity$logFC_orig, discovery_sensitivity$logFC_adj),
    safe_pearson(discovery_sensitivity$logFC_orig, discovery_sensitivity$logFC_adj),
    length(strict_genes),
    length(strict_mt_genes),
    sum(strict_sensitivity$same_direction, na.rm = TRUE),
    sum(strict_sensitivity$retained_after_mito_adjustment_FDR, na.rm = TRUE),
    sum(strict_sensitivity$retained_after_mito_adjustment_nominal, na.rm = TRUE),
    sum(strict_sensitivity$abs_adj_logFC_at_least_50pct_orig, na.rm = TRUE),
    sum(strict_sensitivity$abs_adj_logFC_at_least_75pct_orig, na.rm = TRUE),
    length(expanded_genes),
    length(expanded_mt_genes),
    sum(expanded_sensitivity$same_direction, na.rm = TRUE),
    sum(expanded_sensitivity$retained_after_mito_adjustment_FDR, na.rm = TRUE),
    sum(expanded_sensitivity$retained_after_mito_adjustment_nominal, na.rm = TRUE),
    sum(expanded_sensitivity$abs_adj_logFC_at_least_50pct_orig, na.rm = TRUE),
    sum(expanded_sensitivity$abs_adj_logFC_at_least_75pct_orig, na.rm = TRUE),
    sum(mt_gene_sensitivity$FDR_orig < FDR_THRESHOLD &
          abs(mt_gene_sensitivity$logFC_orig) > LOGFC_THRESHOLD,
        na.rm = TRUE),
    sum(mt_gene_sensitivity$retained_after_mito_adjustment_FDR, na.rm = TRUE),
    sum(mt_gene_sensitivity$same_direction, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  file.path(DE_DIR, "SEAAD_mito_sensitivity_summary_AD61026_v2.csv"),
  row.names = FALSE
)

# ============================================================
# SECTION 10: Save manuscript source data and draft plot
# ============================================================

mito_source <- comparison_all %>%
  mutate(
    analysis_group = case_when(
      strict_primary_consensus ~ "strict_primary_consensus_48",
      expanded_signature ~ "expanded_high_plus_moderate_125",
      original_discovery_gene ~ "SEAAD_discovery_gene",
      mt_gene ~ "MT_gene",
      TRUE ~ "other_tested_gene"
    )
  ) %>%
  arrange(
    desc(strict_primary_consensus),
    desc(expanded_signature),
    desc(original_discovery_gene),
    FDR_orig
  )

write.csv(
  mito_source,
  file.path(SOURCE_DIR, "SuppTable12_mito_sensitivity_source_v2.csv"),
  row.names = FALSE
)

signature_mito_source <- mito_source %>%
  filter(strict_primary_consensus | expanded_signature | mt_gene) %>%
  arrange(
    desc(strict_primary_consensus),
    desc(expanded_signature),
    desc(mt_gene),
    FDR_orig
  )

write.csv(
  signature_mito_source,
  file.path(SOURCE_DIR, "SuppTable12B_signature_and_MT_gene_mito_sensitivity_source_v2.csv"),
  row.names = FALSE
)

figure_source <- discovery_sensitivity %>%
  mutate(
    gene_group = case_when(
      strict_primary_consensus & mt_gene ~ "Strict_MT_gene",
      strict_primary_consensus ~ "Strict_consensus",
      expanded_signature & mt_gene ~ "Expanded_MT_gene",
      expanded_signature ~ "Expanded_only",
      mt_gene ~ "MT_gene",
      TRUE ~ "Other_discovery_gene"
    ),
    label_gene = strict_primary_consensus | mt_gene
  )

write.csv(
  figure_source,
  file.path(SOURCE_DIR, "Figure_MT_sensitivity_scatter_source_v2.csv"),
  row.names = FALSE
)

p_scatter <- ggplot(
  figure_source,
  aes(x = logFC_orig, y = logFC_adj)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.4) +
  geom_point(aes(shape = gene_group), alpha = 0.75, size = 2) +
  ggrepel::geom_text_repel(
    data = figure_source %>% filter(label_gene),
    aes(label = gene),
    size = 2.5,
    max.overlaps = 50,
    box.padding = 0.3
  ) +
  labs(
    title = "SEA-AD mitochondrial fraction sensitivity analysis",
    subtitle = "Original diagnosis model versus diagnosis + MT fraction model",
    x = "Original SEA-AD logFC",
    y = "MT-adjusted SEA-AD logFC",
    shape = "Gene group"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(FIG_DIR, "Fig_MT_sensitivity_logFC_scatter_draft_v2.pdf"),
  p_scatter,
  width = 8,
  height = 6
)

saveRDS(
  list(
    adjusted_de = adjusted_de,
    comparison_all = comparison_all,
    discovery_sensitivity = discovery_sensitivity,
    strict_sensitivity = strict_sensitivity,
    expanded_sensitivity = expanded_sensitivity,
    mt_gene_sensitivity = mt_gene_sensitivity,
    meta_with_mt_fraction = meta,
    mt_group_summary = mt_group_summary,
    mt_diagnostics = mt_diagnostics,
    summary_table = summary_table
  ),
  file.path(RDS_DIR, "SEAAD_mito_sensitivity_AD61026_v2.rds")
)

# ============================================================
# SECTION 11: Print completion summary
# ============================================================

cat("\n============================================================\n")
cat("04B v2 mitochondrial sensitivity analysis complete.\n\n")

cat("Summary:\n")
print(summary_table)

cat("\nKey interpretation checks:\n")
cat("  Strict signature MT genes:",
    ifelse(length(strict_mt_genes) == 0, "none", paste(strict_mt_genes, collapse = ", ")),
    "\n")
cat("  Expanded signature MT genes:",
    ifelse(length(expanded_mt_genes) == 0, "none", paste(expanded_mt_genes, collapse = ", ")),
    "\n")
cat("  Adjusted model used same gene universe as original DE:",
    nrow(adjusted_de) == nrow(original_de),
    "\n\n")

cat("Key output files:\n")
cat("  ", file.path(DE_DIR, "SEAAD_mito_sensitivity_summary_AD61026_v2.csv"), "\n")
cat("  ", file.path(DE_DIR, "SEAAD_mito_sensitivity_strict48_AD61026_v2.csv"), "\n")
cat("  ", file.path(DE_DIR, "SEAAD_mito_sensitivity_expanded125_AD61026_v2.csv"), "\n")
cat("  ", file.path(DE_DIR, "SEAAD_mito_fraction_diagnosis_diagnostics_AD61026_v2.csv"), "\n")
cat("  ", file.path(SOURCE_DIR, "SuppTable12_mito_sensitivity_source_v2.csv"), "\n\n")

cat("Next step:\n")
cat("  Upload SEAAD_mito_sensitivity_summary_AD61026_v2.csv.\n")
cat("============================================================\n")
