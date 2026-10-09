# ============================================================
# SEA-AD covariate sensitivity analysis for TRUE strict 48 genes
# AD61026 project
#
# Purpose:
#   Test whether adjustment for age, sex, and PMI in SEA-AD changes
#   the direction of effect for the strict 48-gene AD microglial
#   consensus signature.
#
# Models compared:
#   1. Diagnosis-only model:
#        ~ diagnosis
#
#   2. Covariate-adjusted model:
#        ~ diagnosis + age + sex + PMI
#
# Outputs:
#   SEAAD_Covariate_Sensitivity_AllGenes.csv
#   SEAAD_Covariate_Sensitivity_TRUE_Strict48.csv
#
# Important:
#   This is a SEA-AD-only sensitivity check.
#   It does not redo the full cross-cohort replication pipeline.
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

# ------------------------------------------------------------
# Project folder. Set the AD_PROJECT_ROOT environment variable to the
# folder that contains Data/, RDS/, Output/ and Supplementary_Tables/,
# or start R with that folder as the working directory.
# ------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("AD_PROJECT_ROOT", unset = getwd())

library(edgeR)
library(limma)
library(dplyr)
library(readr)
library(readxl)
library(tibble)


# ------------------------------------------------------------
# 2. Set file paths
# ------------------------------------------------------------

# Folder containing SEA-AD RDS files
rds_dir <- file.path(PROJECT_ROOT, "RDS")

# Supplementary Table 2 file containing the TRUE Strict48 sheet
# If your file name is different, change this line only.
supp2_file <- file.path(PROJECT_ROOT, "Supplementary_Tables", "Supplementary_Table_2_DE_Consensus_Signatures.xlsx")

# Output folder
out_dir <- rds_dir


# ------------------------------------------------------------
# 3. Load SEA-AD pseudobulk counts and metadata
# ------------------------------------------------------------

seaad_counts <- readRDS(file.path(rds_dir, "SEAAD_pseudobulk_counts.rds"))
seaad_meta   <- readRDS(file.path(rds_dir, "SEAAD_pseudobulk_meta.rds"))

counts <- seaad_counts
meta   <- seaad_meta


# ------------------------------------------------------------
# 4. Load the TRUE strict 48-gene list from Supplementary Table 2
# ------------------------------------------------------------

strict48_table <- readxl::read_excel(
  supp2_file,
  sheet = "Consensus48"
)

# Check the sheet column names
cat("\nColumns in Strict48 sheet:\n")
print(names(strict48_table))

# Extract gene column
if ("gene" %in% names(strict48_table)) {
  strict48 <- strict48_table$gene
} else if ("Gene" %in% names(strict48_table)) {
  strict48 <- strict48_table$Gene
} else {
  stop("Could not find a column named gene or Gene in the Strict48 sheet.")
}

# Clean and de-duplicate gene names
strict48 <- unique(as.character(strict48))
strict48 <- strict48[!is.na(strict48) & strict48 != ""]

# Confirm that this is truly 48 genes
cat("\nNumber of strict consensus genes loaded:", length(strict48), "\n")

if (length(strict48) != 48) {
  stop(
    paste0(
      "The Strict48 gene list has ",
      length(strict48),
      " genes, not 48. Check the Strict48 sheet."
    )
  )
}


# ------------------------------------------------------------
# 5. Define SEA-AD metadata column names
# ------------------------------------------------------------
# These column names come from your seaad_meta object.
#
# Donor.ID matches the column names of seaad_counts.
# diagnosis_std contains Control / AD.
# Age.at.Death, Gender, and PMI.x are used as covariates.

sample_col    <- "Donor.ID"
diagnosis_col <- "diagnosis_std"
age_col       <- "Age.at.Death"
sex_col       <- "Gender"
pmi_col       <- "PMI.x"


# ------------------------------------------------------------
# 6. Basic input checks
# ------------------------------------------------------------

cat("\nCount matrix dimensions:\n")
print(dim(counts))

cat("\nMetadata dimensions:\n")
print(dim(meta))

required_cols <- c(sample_col, diagnosis_col, age_col, sex_col, pmi_col)

missing_cols <- setdiff(required_cols, names(meta))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "Missing required metadata columns: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

if (is.null(rownames(counts))) {
  stop("Count matrix must have gene symbols as row names.")
}

if (is.null(colnames(counts))) {
  stop("Count matrix must have donor/sample IDs as column names.")
}


# ------------------------------------------------------------
# 7. Confirm metadata matches count matrix
# ------------------------------------------------------------

n_matched <- sum(colnames(counts) %in% meta[[sample_col]])

cat("\nNumber of count columns found in metadata:", n_matched, "/", ncol(counts), "\n")

if (!all(colnames(counts) %in% meta[[sample_col]])) {
  missing_samples <- setdiff(colnames(counts), meta[[sample_col]])
  stop(
    paste0(
      "Some count matrix samples are missing from metadata:\n",
      paste(missing_samples, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 8. Align metadata to count matrix column order
# ------------------------------------------------------------

meta <- meta %>%
  dplyr::filter(.data[[sample_col]] %in% colnames(counts))

meta <- meta[match(colnames(counts), meta[[sample_col]]), ]

if (!all(meta[[sample_col]] == colnames(counts))) {
  stop("Metadata rows are not aligned with count matrix columns.")
}

cat("\nMetadata successfully aligned to count matrix.\n")


# ------------------------------------------------------------
# 9. Prepare diagnosis and covariates
# ------------------------------------------------------------

# Diagnosis factor
meta$diagnosis_model <- factor(meta[[diagnosis_col]])

cat("\nDiagnosis groups:\n")
print(table(meta$diagnosis_model, useNA = "ifany"))

# Set Control as reference so positive logFC means higher in AD
if ("Control" %in% levels(meta$diagnosis_model)) {
  meta$diagnosis_model <- relevel(meta$diagnosis_model, ref = "Control")
} else {
  stop("Control label not found in diagnosis_std.")
}

# Covariates
meta$age_model <- suppressWarnings(as.numeric(meta[[age_col]]))
meta$sex_model <- factor(meta[[sex_col]])
meta$pmi_model <- suppressWarnings(as.numeric(meta[[pmi_col]]))

cat("\nSex groups:\n")
print(table(meta$sex_model, useNA = "ifany"))

cat("\nMissing covariate counts:\n")
cat("Age missing:", sum(is.na(meta$age_model)), "\n")
cat("Sex missing:", sum(is.na(meta$sex_model)), "\n")
cat("PMI missing:", sum(is.na(meta$pmi_model)), "\n")


# ------------------------------------------------------------
# 10. Remove donors with missing covariates
# ------------------------------------------------------------

complete_covariates <- complete.cases(
  meta[, c("diagnosis_model", "age_model", "sex_model", "pmi_model")]
)

cat("\nSEA-AD donors before covariate filtering:", nrow(meta), "\n")
cat("SEA-AD donors after covariate filtering:", sum(complete_covariates), "\n")
cat("SEA-AD donors removed because of missing covariates:", sum(!complete_covariates), "\n")

meta_cov <- meta[complete_covariates, ]
counts_cov <- counts[, meta_cov[[sample_col]], drop = FALSE]


# ------------------------------------------------------------
# 11. Create DGEList and filter lowly expressed genes
# ------------------------------------------------------------

dge <- edgeR::DGEList(counts = counts_cov)

# Use diagnosis-only design for expression filtering
design_filter <- model.matrix(~ diagnosis_model, data = meta_cov)

keep <- edgeR::filterByExpr(dge, design = design_filter)

dge <- dge[keep, , keep.lib.sizes = FALSE]

cat("\nGenes retained after filterByExpr:", nrow(dge), "\n")


# ------------------------------------------------------------
# 12. Normalize counts using TMM
# ------------------------------------------------------------

dge <- edgeR::calcNormFactors(dge, method = "TMM")


# ------------------------------------------------------------
# 13. Model 1: diagnosis-only model
# ------------------------------------------------------------

design_unadjusted <- model.matrix(
  ~ diagnosis_model,
  data = meta_cov
)

cat("\nDiagnosis-only model design columns:\n")
print(colnames(design_unadjusted))

v_unadjusted <- limma::voom(
  dge,
  design_unadjusted,
  plot = FALSE
)

fit_unadjusted <- limma::lmFit(v_unadjusted, design_unadjusted)
fit_unadjusted <- limma::eBayes(fit_unadjusted)

# Identify diagnosis coefficient
coef_unadjusted <- grep("^diagnosis_model", colnames(design_unadjusted), value = TRUE)

if (length(coef_unadjusted) != 1) {
  stop("Could not identify a single diagnosis coefficient in the diagnosis-only model.")
}

cat("\nDiagnosis-only coefficient used:\n")
print(coef_unadjusted)

res_unadjusted <- limma::topTable(
  fit_unadjusted,
  coef = coef_unadjusted,
  number = Inf,
  sort.by = "none"
) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::rename(
    diagnosis_only_logFC = logFC,
    diagnosis_only_AveExpr = AveExpr,
    diagnosis_only_t = t,
    diagnosis_only_P.Value = P.Value,
    diagnosis_only_adj.P.Val = adj.P.Val,
    diagnosis_only_B = B
  )


# ------------------------------------------------------------
# 14. Model 2: covariate-adjusted model
# ------------------------------------------------------------

design_adjusted <- model.matrix(
  ~ diagnosis_model + age_model + sex_model + pmi_model,
  data = meta_cov
)

cat("\nCovariate-adjusted model design columns:\n")
print(colnames(design_adjusted))

v_adjusted <- limma::voom(
  dge,
  design_adjusted,
  plot = FALSE
)

fit_adjusted <- limma::lmFit(v_adjusted, design_adjusted)
fit_adjusted <- limma::eBayes(fit_adjusted)

# Identify diagnosis coefficient
coef_adjusted <- grep("^diagnosis_model", colnames(design_adjusted), value = TRUE)

if (length(coef_adjusted) != 1) {
  stop("Could not identify a single diagnosis coefficient in the adjusted model.")
}

cat("\nCovariate-adjusted diagnosis coefficient used:\n")
print(coef_adjusted)

res_adjusted <- limma::topTable(
  fit_adjusted,
  coef = coef_adjusted,
  number = Inf,
  sort.by = "none"
) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::rename(
    covariate_adjusted_logFC = logFC,
    covariate_adjusted_AveExpr = AveExpr,
    covariate_adjusted_t = t,
    covariate_adjusted_P.Value = P.Value,
    covariate_adjusted_adj.P.Val = adj.P.Val,
    covariate_adjusted_B = B
  )


# ------------------------------------------------------------
# 15. Combine diagnosis-only and adjusted results
# ------------------------------------------------------------

covariate_sensitivity_all <- res_unadjusted %>%
  dplyr::inner_join(
    res_adjusted,
    by = "gene"
  ) %>%
  dplyr::mutate(
    same_direction = sign(diagnosis_only_logFC) == sign(covariate_adjusted_logFC),
    abs_logFC_retention = abs(covariate_adjusted_logFC) / abs(diagnosis_only_logFC),
    direction_unadjusted = dplyr::case_when(
      diagnosis_only_logFC > 0 ~ "AD-upregulated",
      diagnosis_only_logFC < 0 ~ "AD-downregulated",
      TRUE ~ "No change"
    ),
    direction_adjusted = dplyr::case_when(
      covariate_adjusted_logFC > 0 ~ "AD-upregulated",
      covariate_adjusted_logFC < 0 ~ "AD-downregulated",
      TRUE ~ "No change"
    )
  )


# ------------------------------------------------------------
# 16. Extract TRUE strict 48-gene sensitivity results
# ------------------------------------------------------------

covariate_sensitivity_strict48 <- covariate_sensitivity_all %>%
  dplyr::filter(gene %in% strict48) %>%
  dplyr::mutate(
    strict48_gene = TRUE
  ) %>%
  dplyr::arrange(desc(abs(diagnosis_only_logFC)))

n_found <- nrow(covariate_sensitivity_strict48)

cat("\nStrict consensus genes found in model output:", n_found, "/ 48\n")

if (n_found != 48) {
  missing_strict_genes <- setdiff(strict48, covariate_sensitivity_all$gene)
  
  warning(
    paste0(
      "Only found ",
      n_found,
      " of 48 strict genes in model output. Missing genes: ",
      paste(missing_strict_genes, collapse = ", ")
    )
  )
}


# ------------------------------------------------------------
# 17. Summarize strict 48-gene sensitivity results
# ------------------------------------------------------------

n_same_direction <- sum(covariate_sensitivity_strict48$same_direction, na.rm = TRUE)

n_retained_50 <- sum(
  covariate_sensitivity_strict48$abs_logFC_retention >= 0.50,
  na.rm = TRUE
)

n_retained_75 <- sum(
  covariate_sensitivity_strict48$abs_logFC_retention >= 0.75,
  na.rm = TRUE
)

cat("\n============================================================\n")
cat("SEA-AD covariate sensitivity summary for TRUE strict 48\n")
cat("============================================================\n")
cat("Strict consensus genes found:", n_found, "/ 48\n")
cat("Same direction after age/sex/PMI adjustment:", n_same_direction, "/", n_found, "\n")
cat("At least 50% absolute logFC retained:", n_retained_50, "/", n_found, "\n")
cat("At least 75% absolute logFC retained:", n_retained_75, "/", n_found, "\n")
cat("============================================================\n")


# ------------------------------------------------------------
# 18. Show genes that changed direction, if any
# ------------------------------------------------------------

changed_direction <- covariate_sensitivity_strict48 %>%
  dplyr::filter(!same_direction) %>%
  dplyr::select(
    gene,
    diagnosis_only_logFC,
    covariate_adjusted_logFC,
    diagnosis_only_adj.P.Val,
    covariate_adjusted_adj.P.Val,
    abs_logFC_retention
  )

if (nrow(changed_direction) > 0) {
  cat("\nStrict consensus genes that changed direction after covariate adjustment:\n")
  print(changed_direction)
} else {
  cat("\nNo strict consensus genes changed direction after covariate adjustment.\n")
}


# ------------------------------------------------------------
# 19. Save output files
# ------------------------------------------------------------

out_all <- file.path(out_dir, "SEAAD_Covariate_Sensitivity_AllGenes.csv")
out_strict48 <- file.path(out_dir, "SEAAD_Covariate_Sensitivity_TRUE_Strict48.csv")

readr::write_csv(covariate_sensitivity_all, out_all)
readr::write_csv(covariate_sensitivity_strict48, out_strict48)

cat("\nSaved output files:\n")
cat(out_all, "\n")
cat(out_strict48, "\n")


# ------------------------------------------------------------
# 20. Save a one-row summary table for easy manuscript reporting
# ------------------------------------------------------------

summary_table <- tibble::tibble(
  analysis = "SEAAD_age_sex_PMI_covariate_sensitivity_TRUE_Strict48",
  strict_genes_expected = 48,
  strict_genes_found = n_found,
  same_direction_after_adjustment = n_same_direction,
  retained_50pct_abs_logFC = n_retained_50,
  retained_75pct_abs_logFC = n_retained_75,
  n_changed_direction = nrow(changed_direction)
)

out_summary <- file.path(out_dir, "SEAAD_Covariate_Sensitivity_TRUE_Strict48_Summary.csv")

readr::write_csv(summary_table, out_summary)

cat(out_summary, "\n")
