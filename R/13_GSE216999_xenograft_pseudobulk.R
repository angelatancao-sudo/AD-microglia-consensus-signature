# ============================================================
# GSE216999 mouse-level pseudobulk orthogonal support analysis
# Strict 48 AD microglia consensus signature vs APP-NLGF / APP-WT
# ============================================================
#
# Purpose:
# Reanalyze GSE216999 using mouse-level pseudobulk, not cell-level
# mean expression. This aligns the orthogonal support analysis with
# the main manuscript's donor/sample-level pseudobulk framework.
#
# Main comparison:
# APP-NLGF host brain vs APP-WT host brain
#
# Cell subset:
# mouse_genotype == APP-NLGF or APP-WT
# x_mg_genotype == WT
# x_mg_background == H9
# injection == None
# age == 6m
#
# Analysis:
# 1. Download raw count matrix files if needed
# 2. Read cell metadata
# 3. Subset cells to clean APP-NLGF vs APP-WT contrast
# 4. Aggregate raw counts to mouse-level pseudobulk
# 5. Run edgeR / limma-voom
# 6. Intersect with strict 48 consensus signature
# 7. Test direction concordance using a two-sided binomial sign test
# 8. Calculate Pearson and Spearman logFC correlations
# 9. Save clean outputs for manuscript and supplement
# ============================================================


# ============================================================
# 1. Load/install required packages
# ============================================================

# CRAN packages

# ------------------------------------------------------------
# Project folder. Set the AD_PROJECT_ROOT environment variable to the
# folder that contains Data/, RDS/, Output/ and Supplementary_Tables/,
# or start R with that folder as the working directory.
# ------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("AD_PROJECT_ROOT", unset = getwd())

cran_packages <- c(
  "Matrix",
  "dplyr",
  "readr",
  "stringr",
  "janitor",
  "tibble",
  "ggplot2",
  "openxlsx"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c("edgeR", "limma")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(Matrix)
library(dplyr)
library(readr)
library(stringr)
library(janitor)
library(tibble)
library(ggplot2)
library(openxlsx)
library(edgeR)
library(limma)


# ============================================================
# 2. Define folders, files, and analysis label
# ============================================================

outdir <- file.path(PROJECT_ROOT, "GSE216999_PSEUDOBULK_ORTHOGONAL_SUPPORT")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# This file was generated earlier when listing GEO supplementary files.
# It should contain columns including fname and url.
supp_file_list <- file.path(
  "GSE216999_FULL_MATRIX_ORTHOGONAL_SUPPORT",
  "GSE216999_GEO_supplementary_file_list_NO_DOWNLOAD.csv"
)

# If your file list is in the current working directory instead,
# change the path above accordingly.

# Your strict 48-gene covariate sensitivity file.
strict48_file <- file.path(PROJECT_ROOT, "RDS/SEAAD_Covariate_Sensitivity_TRUE_Strict48.csv")

# Output label used for all saved files.
analysis_label <- "PSEUDOBULK_Strict48_vs_GSE216999_APP_NLGF_vs_APP_WT"


# ============================================================
# 3. Define GSE216999 supplementary file list directly
# ============================================================
# This avoids needing the previously generated file:
# GSE216999_GEO_supplementary_file_list_NO_DOWNLOAD.csv
#
# The files are hosted on the NCBI GEO FTP server.
# We define them directly so the script can continue no matter
# where your previous file-list CSV was saved.
# ============================================================

geo_base_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE216nnn/GSE216999/suppl/"

supp_files <- tibble::tibble(
  fname = c(
    "GSE216999_Mancuso2022_NormExpr.mtx.gz",
    "GSE216999_Mancuso2022_NormExpr_cellbarcodes.txt.gz",
    "GSE216999_Mancuso2022_NormExpr_genes.txt.gz",
    "GSE216999_Mancuso2022_RNAcounts.mtx.gz",
    "GSE216999_Mancuso2022_RNAcounts_cellbarcodes.txt.gz",
    "GSE216999_Mancuso2022_RNAcounts_genes.txt.gz",
    "GSE216999_Mancuso2022_metadata.csv.gz"
  )
) %>%
  mutate(
    url = paste0(geo_base_url, fname)
  )

cat("\nGSE216999 supplementary file list defined directly:\n")
print(supp_files)

# Save this file list inside the pseudobulk output folder
# so the analysis folder is self-contained.

write_csv(
  supp_files,
  file.path(outdir, "GSE216999_GEO_supplementary_file_list_DIRECT.csv")
)



# ============================================================
# 4. Helper function to download one GEO supplementary file
# ============================================================

download_geo_file <- function(pattern, outdir, supp_files) {
  
  selected <- supp_files %>%
    filter(str_detect(fname, regex(pattern, ignore_case = TRUE)))
  
  if (nrow(selected) != 1) {
    stop(
      "Could not uniquely identify file for pattern: ",
      pattern,
      "\nMatched files:\n",
      paste(selected$fname, collapse = "\n")
    )
  }
  
  local_file <- file.path(outdir, selected$fname)
  
  if (!file.exists(local_file)) {
    cat("\nDownloading:\n")
    cat(selected$fname, "\n")
    
    download.file(
      url = selected$url,
      destfile = local_file,
      mode = "wb",
      method = "libcurl"
    )
    
  } else {
    cat("\nFile already exists locally, skipping download:\n")
    cat(selected$fname, "\n")
  }
  
  return(local_file)
}


# ============================================================
# 5. Download raw count matrix, gene, barcode, and metadata files
# ============================================================
# We use raw counts for pseudobulk aggregation.
# This is better aligned with the main donor-level pseudobulk framework.

counts_mtx_file <- download_geo_file(
  pattern = "RNAcounts\\.mtx\\.gz$",
  outdir = outdir,
  supp_files = supp_files
)

counts_barcodes_file <- download_geo_file(
  pattern = "RNAcounts_cellbarcodes\\.txt\\.gz$",
  outdir = outdir,
  supp_files = supp_files
)

counts_genes_file <- download_geo_file(
  pattern = "RNAcounts_genes\\.txt\\.gz$",
  outdir = outdir,
  supp_files = supp_files
)

metadata_file <- download_geo_file(
  pattern = "metadata\\.csv\\.gz$",
  outdir = outdir,
  supp_files = supp_files
)


# ============================================================
# 6. Helper function to parse GEO indexed text files
# ============================================================
# The GEO gene/barcode files look like:
#   "1" "AAACCCAAGTCAGGGT-1_lib3.4"
#
# The first field is just an index.
# The second field is the actual gene or barcode.
# This function extracts the second field.

parse_indexed_text_file <- function(x) {
  
  x <- as.character(x)
  x <- str_trim(x)
  x <- x[x != ""]
  
  second_quoted_field <- str_match(
    x,
    '^"?[0-9]+"?\\s+"?([^"]+)"?$'
  )[, 2]
  
  parsed <- ifelse(
    !is.na(second_quoted_field),
    second_quoted_field,
    x
  )
  
  parsed <- parsed %>%
    str_replace_all('"', "") %>%
    str_replace_all("'", "") %>%
    str_trim()
  
  return(parsed)
}


# ============================================================
# 7. Read metadata, genes, and barcodes
# ============================================================

cell_meta <- read_csv(metadata_file, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(
    cell_id = as.character(cell_id),
    cell_id = str_trim(cell_id),
    cell_id = str_replace_all(cell_id, '"', ""),
    cell_id = str_replace_all(cell_id, "'", "")
  )

cat("\nCell metadata dimensions:\n")
print(dim(cell_meta))

cat("\nCell metadata columns:\n")
print(colnames(cell_meta))

barcodes_raw <- read_lines(counts_barcodes_file)
genes_raw <- read_lines(counts_genes_file)

barcodes_clean <- parse_indexed_text_file(barcodes_raw)
genes_clean <- parse_indexed_text_file(genes_raw)

cat("\nParsed barcode count:", length(barcodes_clean), "\n")
cat("Parsed gene count:", length(genes_clean), "\n")

cat("\nFirst 5 parsed barcodes:\n")
print(head(barcodes_clean, 5))

cat("\nFirst 5 parsed genes:\n")
print(head(genes_clean, 5))


# ============================================================
# 8. Read raw count sparse matrix
# ============================================================
# This may take several minutes.

cat("\nReading raw count matrix. This may take several minutes...\n")

counts_mat <- readMM(gzfile(counts_mtx_file))

cat("\nRaw count matrix dimensions:\n")
print(dim(counts_mat))


# ============================================================
# 9. Fix one-line mismatch if present
# ============================================================
# As with the normalized matrix, the gene/barcode files may contain
# one extra header/index line. If so, remove the first entry.

if (length(genes_clean) == nrow(counts_mat) + 1) {
  cat("\nGene file has one extra line. Removing first gene entry:\n")
  print(genes_clean[1])
  genes_clean <- genes_clean[-1]
}

if (length(barcodes_clean) == ncol(counts_mat) + 1) {
  cat("\nBarcode file has one extra line. Removing first barcode entry:\n")
  print(barcodes_clean[1])
  barcodes_clean <- barcodes_clean[-1]
}

cat("\nCorrected barcode count:", length(barcodes_clean), "\n")
cat("Corrected gene count:", length(genes_clean), "\n")


# ============================================================
# 10. Assign row and column names to count matrix
# ============================================================

if (nrow(counts_mat) == length(genes_clean) && ncol(counts_mat) == length(barcodes_clean)) {
  
  cat("\nMatrix orientation detected: genes x cells\n")
  rownames(counts_mat) <- genes_clean
  colnames(counts_mat) <- barcodes_clean
  
} else if (nrow(counts_mat) == length(barcodes_clean) && ncol(counts_mat) == length(genes_clean)) {
  
  cat("\nMatrix orientation detected: cells x genes. Transposing to genes x cells.\n")
  counts_mat <- t(counts_mat)
  rownames(counts_mat) <- genes_clean
  colnames(counts_mat) <- barcodes_clean
  
} else {
  
  stop(
    "Matrix dimensions do not match after parsing.\n",
    "Matrix dimensions: ", paste(dim(counts_mat), collapse = " x "), "\n",
    "Genes: ", length(genes_clean), "\n",
    "Barcodes: ", length(barcodes_clean)
  )
}

cat("\nFinal raw count matrix dimensions:\n")
print(dim(counts_mat))


# ============================================================
# 11. Match metadata cells to matrix barcodes
# ============================================================

n_meta_in_matrix <- sum(cell_meta$cell_id %in% colnames(counts_mat))

cat("\nMetadata cells found in raw count matrix:", n_meta_in_matrix, "of", nrow(cell_meta), "\n")

if (n_meta_in_matrix == 0) {
  
  cat("\nFirst 10 metadata cell IDs:\n")
  print(head(cell_meta$cell_id, 10))
  
  cat("\nFirst 10 matrix barcodes:\n")
  print(head(colnames(counts_mat), 10))
  
  stop("No metadata cell IDs matched raw count matrix barcodes.")
}

cell_meta_aligned <- cell_meta %>%
  filter(cell_id %in% colnames(counts_mat)) %>%
  mutate(matrix_order = match(cell_id, colnames(counts_mat))) %>%
  arrange(matrix_order)

cat("\nAligned metadata dimensions:\n")
print(dim(cell_meta_aligned))


# ============================================================
# 12. Define the clean APP-NLGF vs APP-WT analysis subset
# ============================================================
# This removes other experimental perturbations such as:
# TREM2 altered xenografted microglia
# APOE genotype experiments
# injected oAb/Scr treatments
# APP-NLGF/ApoeKO hosts

main_meta <- cell_meta_aligned %>%
  filter(
    mouse_genotype %in% c("APP-NLGF", "APP-WT"),
    x_mg_genotype == "WT",
    x_mg_background == "H9",
    injection == "None",
    age == "6m"
  ) %>%
  mutate(
    contrast_group = case_when(
      mouse_genotype == "APP-NLGF" ~ "NLGF",
      mouse_genotype == "APP-WT" ~ "WT",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(contrast_group))

cat("\nSelected cell counts by contrast group:\n")
print(table(main_meta$contrast_group))

cat("\nSelected cell counts by mouse:\n")
print(table(main_meta$contrast_group, main_meta$mouse_id))

mouse_qc <- main_meta %>%
  count(mouse_id, contrast_group, name = "n_cells") %>%
  arrange(contrast_group, mouse_id)

write_csv(
  mouse_qc,
  file.path(outdir, paste0(analysis_label, "_mouse_cell_count_QC.csv"))
)

cat("\nMouse-level cell count QC:\n")
print(mouse_qc, n = Inf)

if (length(unique(main_meta$contrast_group)) < 2) {
  stop("The selected subset does not contain both NLGF and WT groups.")
}


# ============================================================
# 13. Optional minimum-cell filtering at mouse level
# ============================================================
# This avoids unstable pseudobulk profiles from mice with very few cells.
# Set min_cells_per_mouse to 0 to keep every mouse.
#
# A threshold of 50 is conservative but not too strict.
# Check the mouse QC output before changing this.

min_cells_per_mouse <- 50

mouse_keep <- mouse_qc %>%
  filter(n_cells >= min_cells_per_mouse) %>%
  pull(mouse_id)

main_meta <- main_meta %>%
  filter(mouse_id %in% mouse_keep)

mouse_qc_filtered <- main_meta %>%
  count(mouse_id, contrast_group, name = "n_cells") %>%
  arrange(contrast_group, mouse_id)

cat("\nMouse-level cell count QC after filtering:\n")
print(mouse_qc_filtered, n = Inf)

cat("\nNumber of mice retained by group:\n")
print(table(mouse_qc_filtered$contrast_group))

write_csv(
  mouse_qc_filtered,
  file.path(outdir, paste0(analysis_label, "_mouse_cell_count_QC_after_filtering.csv"))
)

if (length(unique(mouse_qc_filtered$contrast_group)) < 2) {
  stop("After mouse-level cell filtering, one group has no retained mice.")
}


# ============================================================
# 14. Aggregate raw counts to mouse-level pseudobulk
# ============================================================
# Matrix dimensions:
# counts_mat = genes x cells
#
# group_matrix = cells x mice
#
# pseudobulk_counts = counts_mat %*% group_matrix
# Result = genes x mice

main_cells <- main_meta$cell_id

# Subset count matrix to selected cells.
main_counts <- counts_mat[, main_cells, drop = FALSE]

# Create mouse factor in the same order as main_cells.
mouse_factor <- factor(main_meta$mouse_id)

# Sparse cell-to-mouse design matrix.
group_matrix <- sparseMatrix(
  i = seq_along(mouse_factor),
  j = as.integer(mouse_factor),
  x = 1,
  dims = c(length(mouse_factor), nlevels(mouse_factor)),
  dimnames = list(main_cells, levels(mouse_factor))
)

# Sum raw counts per mouse.
pseudobulk_counts <- main_counts %*% group_matrix

# Convert to matrix-like object with integer-ish counts.
pseudobulk_counts <- as.matrix(pseudobulk_counts)

cat("\nPseudobulk count matrix dimensions, genes x mice:\n")
print(dim(pseudobulk_counts))

# Create mouse-level sample metadata.
sample_meta <- mouse_qc_filtered %>%
  distinct(mouse_id, contrast_group, n_cells) %>%
  arrange(match(mouse_id, colnames(pseudobulk_counts)))

# Reorder pseudobulk matrix to match sample metadata.
pseudobulk_counts <- pseudobulk_counts[, sample_meta$mouse_id, drop = FALSE]

cat("\nFinal sample metadata:\n")
print(sample_meta, n = Inf)

write_csv(
  sample_meta,
  file.path(outdir, paste0(analysis_label, "_pseudobulk_sample_metadata.csv"))
)

# Save pseudobulk count matrix as CSV for record keeping.
# This file can be large but is useful for reproducibility.
write_csv(
  as.data.frame(pseudobulk_counts) %>%
    rownames_to_column("gene"),
  file.path(outdir, paste0(analysis_label, "_pseudobulk_counts_by_mouse.csv"))
)


# ============================================================
# 15. edgeR / limma-voom pseudobulk differential analysis
# Corrected version: avoids duplicated gene column error
# ============================================================

# Make group factor with WT as reference.
# This means the model coefficient groupNLGF represents APP-NLGF vs APP-WT.
group <- factor(sample_meta$contrast_group, levels = c("WT", "NLGF"))

# Create DGEList object from mouse-level pseudobulk raw counts.
dge <- DGEList(
  counts = pseudobulk_counts,
  group = group
)

# Store gene names in the DGEList.
# topTable() will carry this gene column into the output.
dge$genes <- data.frame(
  gene = rownames(pseudobulk_counts),
  stringsAsFactors = FALSE
)

# Filter lowly expressed genes using edgeR's expression filter.
# This keeps genes with enough counts across mouse-level pseudobulk samples.
keep <- filterByExpr(dge, group = group)

cat("\nGenes before filtering:", nrow(dge), "\n")
cat("Genes retained after filterByExpr:", sum(keep), "\n")

dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

# TMM normalization.
dge_filtered <- calcNormFactors(dge_filtered)

# Design matrix for NLGF vs WT.
design <- model.matrix(~ group)

cat("\nDesign matrix:\n")
print(design)

# voom transformation.
v <- voom(
  dge_filtered,
  design = design,
  plot = FALSE
)

# Fit limma model.
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Extract full results.
# Coefficient "groupNLGF" means APP-NLGF vs APP-WT.
#
# Important fix:
# We use row_id instead of gene in rownames_to_column()
# because topTable already contains a gene column from dge$genes.
pb_de_raw <- topTable(
  fit,
  coef = "groupNLGF",
  number = Inf,
  sort.by = "none"
) %>%
  as.data.frame() %>%
  rownames_to_column("row_id") %>%
  as_tibble() %>%
  clean_names()

cat("\nColumns in raw pseudobulk DE output:\n")
print(colnames(pb_de_raw))

# Decide which column contains gene names.
# Usually this will be "gene" from dge$genes.
# If not present, fall back to row_id.
gene_col_for_pb <- if ("gene" %in% colnames(pb_de_raw)) {
  "gene"
} else {
  "row_id"
}

cat("\nUsing this column as gene name:", gene_col_for_pb, "\n")

# Create clean pseudobulk DE table.
pb_de <- pb_de_raw %>%
  transmute(
    gene = as.character(.data[[gene_col_for_pb]]),
    gene_upper = str_to_upper(gene),
    pseudobulk_logFC = as.numeric(log_fc),
    pseudobulk_ave_expr = as.numeric(ave_expr),
    pseudobulk_t = as.numeric(t),
    pseudobulk_p_value = as.numeric(p_value),
    pseudobulk_adj_p_value = as.numeric(adj_p_val),
    pseudobulk_b = as.numeric(b)
  ) %>%
  mutate(
    pseudobulk_direction = case_when(
      pseudobulk_logFC > 0 ~ "NLGF_up",
      pseudobulk_logFC < 0 ~ "NLGF_down",
      TRUE ~ "neutral"
    ),
    pseudobulk_sign = case_when(
      pseudobulk_logFC > 0 ~ 1,
      pseudobulk_logFC < 0 ~ -1,
      TRUE ~ 0
    )
  )

write_csv(
  pb_de,
  file.path(outdir, paste0(analysis_label, "_FULL_pseudobulk_DE_all_genes.csv"))
)

cat("\nPseudobulk DE table saved.\n")
cat("Number of genes in pseudobulk DE table:", nrow(pb_de), "\n")


# ============================================================
# 16. Read strict 48 signature
# ============================================================

strict48_raw <- read_csv(strict48_file, show_col_types = FALSE) %>%
  clean_names()

cat("\nColumns in strict 48 file:\n")
print(colnames(strict48_raw))

find_first_matching_col <- function(df, patterns) {
  
  cols <- colnames(df)
  
  matched <- cols[
    str_detect(
      cols,
      regex(paste(patterns, collapse = "|"), ignore_case = TRUE)
    )
  ]
  
  if (length(matched) == 0) {
    stop(
      "Could not find matching column. Available columns are:\n",
      paste(cols, collapse = ", ")
    )
  }
  
  matched[1]
}

strict_gene_col <- find_first_matching_col(
  strict48_raw,
  patterns = c("^gene$", "gene_name", "symbol")
)

strict_logfc_col <- find_first_matching_col(
  strict48_raw,
  patterns = c(
    "diagnosis_only_log_fc",
    "seaad_log_fc",
    "log2fc",
    "log_fc",
    "logfc"
  )
)

strict48 <- strict48_raw %>%
  transmute(
    strict_gene = str_trim(as.character(.data[[strict_gene_col]])),
    gene_upper = str_to_upper(strict_gene),
    seaad_logFC = as.numeric(.data[[strict_logfc_col]])
  ) %>%
  filter(
    !is.na(strict_gene),
    strict_gene != "",
    !is.na(seaad_logFC)
  ) %>%
  mutate(
    seaad_direction = case_when(
      seaad_logFC > 0 ~ "AD_up",
      seaad_logFC < 0 ~ "AD_down",
      TRUE ~ "neutral"
    ),
    seaad_sign = case_when(
      seaad_logFC > 0 ~ 1,
      seaad_logFC < 0 ~ -1,
      TRUE ~ 0
    )
  )

cat("\nStrict 48 genes loaded:", nrow(strict48), "\n")


# ============================================================
# 17. Intersect strict 48 with pseudobulk DE table
# ============================================================

strict48_pb_overlap <- strict48 %>%
  left_join(
    pb_de,
    by = "gene_upper"
  ) %>%
  mutate(
    detected_in_pseudobulk = !is.na(pseudobulk_logFC),
    concordant_direction = case_when(
      detected_in_pseudobulk & seaad_sign == pseudobulk_sign ~ TRUE,
      detected_in_pseudobulk & seaad_sign != pseudobulk_sign ~ FALSE,
      TRUE ~ NA
    ),
    concordance_label = case_when(
      concordant_direction == TRUE ~ "Concordant",
      concordant_direction == FALSE ~ "Discordant",
      TRUE ~ "Not retained in pseudobulk DE table"
    )
  )

gene_level_output <- strict48_pb_overlap %>%
  arrange(
    desc(detected_in_pseudobulk),
    desc(concordant_direction),
    desc(abs(pseudobulk_logFC))
  ) %>%
  select(
    gene = strict_gene,
    seaad_logFC,
    seaad_direction,
    pseudobulk_logFC,
    pseudobulk_direction,
    pseudobulk_ave_expr,
    pseudobulk_p_value,
    pseudobulk_adj_p_value,
    detected_in_pseudobulk,
    concordant_direction,
    concordance_label
  )

write_csv(
  gene_level_output,
  file.path(outdir, paste0(analysis_label, "_strict48_gene_level_output.csv"))
)

cat("\nStrict 48 pseudobulk overlap table:\n")
print(gene_level_output, n = Inf)


# ============================================================
# 18. Concordance summary
# ============================================================

summary_tbl <- strict48_pb_overlap %>%
  summarize(
    strict48_total = n(),
    detected_in_pseudobulk = sum(detected_in_pseudobulk, na.rm = TRUE),
    not_detected_in_pseudobulk = strict48_total - detected_in_pseudobulk,
    concordant = sum(concordant_direction == TRUE, na.rm = TRUE),
    discordant = sum(concordant_direction == FALSE, na.rm = TRUE),
    concordance_rate = concordant / detected_in_pseudobulk
  )

write_csv(
  summary_tbl,
  file.path(outdir, paste0(analysis_label, "_summary.csv"))
)

cat("\nConcordance summary:\n")
print(summary_tbl)


# ============================================================
# 19. Direction-specific summary
# ============================================================

direction_summary <- strict48_pb_overlap %>%
  filter(detected_in_pseudobulk) %>%
  group_by(seaad_direction) %>%
  summarize(
    detected = n(),
    concordant = sum(concordant_direction == TRUE, na.rm = TRUE),
    discordant = sum(concordant_direction == FALSE, na.rm = TRUE),
    concordance_rate = concordant / detected,
    .groups = "drop"
  )

write_csv(
  direction_summary,
  file.path(outdir, paste0(analysis_label, "_direction_summary.csv"))
)

cat("\nDirection-specific summary:\n")
print(direction_summary)


# ============================================================
# 20. Two-sided binomial sign test
# ============================================================
# This tests whether directional concordance differs from 50%.
# Two-sided is used to avoid the reviewer question of why a one-sided
# test was selected.

n_concordant <- summary_tbl$concordant
n_detected <- summary_tbl$detected_in_pseudobulk

if (!is.na(n_detected) && n_detected > 0) {
  
  binom_result <- binom.test(
    x = n_concordant,
    n = n_detected,
    p = 0.5,
    alternative = "two.sided"
  )
  
  binom_tbl <- tibble(
    detected_in_pseudobulk = n_detected,
    concordant = n_concordant,
    discordant = summary_tbl$discordant,
    concordance_rate = n_concordant / n_detected,
    binomial_p_value_two_sided = binom_result$p.value,
    conf_low = binom_result$conf.int[1],
    conf_high = binom_result$conf.int[2]
  )
  
} else {
  
  binom_tbl <- tibble(
    detected_in_pseudobulk = NA_integer_,
    concordant = NA_integer_,
    discordant = NA_integer_,
    concordance_rate = NA_real_,
    binomial_p_value_two_sided = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_
  )
}

write_csv(
  binom_tbl,
  file.path(outdir, paste0(analysis_label, "_binomial_test_two_sided.csv"))
)

cat("\nTwo-sided binomial sign test:\n")
print(binom_tbl)


# ============================================================
# 21. Pearson and Spearman correlation
# ============================================================

cor_df <- strict48_pb_overlap %>%
  filter(
    detected_in_pseudobulk,
    !is.na(seaad_logFC),
    !is.na(pseudobulk_logFC)
  )

if (nrow(cor_df) >= 3) {
  
  pearson_test <- cor.test(
    cor_df$seaad_logFC,
    cor_df$pseudobulk_logFC,
    method = "pearson"
  )
  
  spearman_test <- cor.test(
    cor_df$seaad_logFC,
    cor_df$pseudobulk_logFC,
    method = "spearman",
    exact = FALSE
  )
  
  correlation_tbl <- tibble(
    n_genes = nrow(cor_df),
    pearson_r = unname(pearson_test$estimate),
    pearson_p_value = pearson_test$p.value,
    spearman_rho = unname(spearman_test$estimate),
    spearman_p_value = spearman_test$p.value
  )
  
} else {
  
  correlation_tbl <- tibble(
    n_genes = nrow(cor_df),
    pearson_r = NA_real_,
    pearson_p_value = NA_real_,
    spearman_rho = NA_real_,
    spearman_p_value = NA_real_
  )
}

write_csv(
  correlation_tbl,
  file.path(outdir, paste0(analysis_label, "_correlation_tests.csv"))
)

cat("\nCorrelation results:\n")
print(correlation_tbl)


# ============================================================
# 22. Scatter plot
# ============================================================

if (nrow(cor_df) > 0) {
  
  p <- ggplot(
    cor_df,
    aes(
      x = seaad_logFC,
      y = pseudobulk_logFC
    )
  ) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_point(size = 2.7, alpha = 0.85) +
    geom_text(
      aes(label = strict_gene),
      vjust = -0.7,
      size = 3
    ) +
    labs(
      title = "Mouse-level pseudobulk support of the strict 48-gene signature",
      subtitle = "SEA-AD AD-vs-control logFC compared with GSE216999 APP-NLGF-vs-APP-WT pseudobulk logFC",
      x = "SEA-AD AD vs control logFC",
      y = "GSE216999 APP-NLGF vs APP-WT pseudobulk logFC"
    ) +
    theme_classic(base_size = 12)
  
  print(p)
  
  ggsave(
    filename = file.path(outdir, paste0(analysis_label, "_logFC_correlation_scatter.pdf")),
    plot = p,
    width = 7,
    height = 5.5
  )
  
  ggsave(
    filename = file.path(outdir, paste0(analysis_label, "_logFC_correlation_scatter.png")),
    plot = p,
    width = 7,
    height = 5.5,
    dpi = 300
  )
}


# ============================================================
# 23. Create clean supplementary Excel workbook
# ============================================================

summary_final <- summary_tbl %>%
  bind_cols(
    binom_tbl %>%
      select(
        binomial_p_value_two_sided,
        conf_low,
        conf_high
      )
  ) %>%
  bind_cols(
    correlation_tbl %>%
      select(
        pearson_r,
        pearson_p_value,
        spearman_rho,
        spearman_p_value
      )
  ) %>%
  mutate(
    analysis = "GSE216999 mouse-level pseudobulk; APP-NLGF vs APP-WT host-brain contrast; WT H9 xenografted human microglia; age 6 months; no injection",
    interpretation = "Orthogonal support analysis only; not used to define the strict consensus signature",
    n_wt_mice = sum(sample_meta$contrast_group == "WT"),
    n_nlgf_mice = sum(sample_meta$contrast_group == "NLGF"),
    n_total_mice = nrow(sample_meta),
    n_wt_cells = sum(main_meta$contrast_group == "WT"),
    n_nlgf_cells = sum(main_meta$contrast_group == "NLGF"),
    n_total_cells = nrow(main_meta),
    min_cells_per_mouse = min_cells_per_mouse
  ) %>%
  select(
    analysis,
    interpretation,
    n_wt_mice,
    n_nlgf_mice,
    n_total_mice,
    n_wt_cells,
    n_nlgf_cells,
    n_total_cells,
    min_cells_per_mouse,
    strict48_total,
    detected_in_pseudobulk,
    not_detected_in_pseudobulk,
    concordant,
    discordant,
    concordance_rate,
    binomial_p_value_two_sided,
    conf_low,
    conf_high,
    pearson_r,
    pearson_p_value,
    spearman_rho,
    spearman_p_value
  )

direction_final <- direction_summary %>%
  arrange(seaad_direction)

gene_final <- gene_level_output %>%
  arrange(
    desc(detected_in_pseudobulk),
    desc(concordant_direction),
    seaad_direction,
    gene
  )

output_excel <- file.path(
  outdir,
  "Supplementary_Table_GSE216999_MouseLevel_Pseudobulk_Orthogonal_Support.xlsx"
)

wb <- createWorkbook()

addWorksheet(wb, "Summary")
addWorksheet(wb, "Mouse_QC")
addWorksheet(wb, "Direction_summary")
addWorksheet(wb, "Consensus48_gene_level")

writeData(wb, "Summary", summary_final)
writeData(wb, "Mouse_QC", mouse_qc_filtered)
writeData(wb, "Direction_summary", direction_final)
writeData(wb, "Consensus48_gene_level", gene_final)

setColWidths(wb, "Summary", cols = 1:ncol(summary_final), widths = "auto")
setColWidths(wb, "Mouse_QC", cols = 1:ncol(mouse_qc_filtered), widths = "auto")
setColWidths(wb, "Direction_summary", cols = 1:ncol(direction_final), widths = "auto")
setColWidths(wb, "Consensus48_gene_level", cols = 1:ncol(gene_final), widths = "auto")

freezePane(wb, "Summary", firstRow = TRUE)
freezePane(wb, "Mouse_QC", firstRow = TRUE)
freezePane(wb, "Direction_summary", firstRow = TRUE)
freezePane(wb, "Consensus48_gene_level", firstRow = TRUE)

saveWorkbook(
  wb,
  output_excel,
  overwrite = TRUE
)

cat("\nSupplementary Excel workbook saved to:\n")
cat(output_excel, "\n")


# ============================================================
# 24. Final message
# ============================================================

cat("\n============================================================\n")
cat("Mouse-level pseudobulk orthogonal support analysis complete.\n")
cat("Outputs saved in:\n")
cat(outdir, "\n")
cat("\nKey output files:\n")
cat(file.path(outdir, paste0(analysis_label, "_summary.csv")), "\n")
cat(file.path(outdir, paste0(analysis_label, "_binomial_test_two_sided.csv")), "\n")
cat(file.path(outdir, paste0(analysis_label, "_correlation_tests.csv")), "\n")
cat(file.path(outdir, paste0(analysis_label, "_direction_summary.csv")), "\n")
cat(file.path(outdir, paste0(analysis_label, "_strict48_gene_level_output.csv")), "\n")
cat(file.path(outdir, paste0(analysis_label, "_mouse_cell_count_QC_after_filtering.csv")), "\n")
cat(output_excel, "\n")
cat("============================================================\n")