# ============================================================
# 15c_drug_prioritization_final.R
# ============================================================
# Project: AD Microglia Multi-Dataset Transcriptomic Analysis
#
# Purpose:
#   Build a sign-checked final single-drug prioritization table using
#   CLUE/CMap, iLINCS, DREIMT, and DGIdb.
#
# Why this version exists:
#   Each drug-prioritization platform uses a different scoring convention.
#   This script explicitly defines the correct reversal direction for each
#   tool before assigning support tiers.
#
# Correct sign conventions used here:
#
#   CLUE/CMap:
#     Reversal direction = negative normalized connectivity score.
#     More negative norm_cs = stronger reversal.
#     Strong support requires:
#       (median_norm_cs <= -0.70 OR mean_norm_cs <= -0.70)
#       AND n_negative >= 2
#     Rationale:
#       This avoids calling a drug strong based on only one negative profile.
#
#   iLINCS:
#     Reversal direction = Correlation == "-"
#     In the iLINCS connected perturbation output, zScore is a positive
#     strength statistic for the negative/discordant relationship.
#     Strong support requires:
#       Correlation == "-"
#       zScore > 3
#       pValue < 0.05
#
#   DREIMT:
#     Reversal/inhibitory direction = tau < 0.
#     Strong reversal direction used for screening:
#       tau < -75
#     FDR-supported DREIMT reversal:
#       tau < -75
#       AND either up_dr < 0.05 OR down_fdr < 0.05
#     Positive tau is NOT reversal in this analysis.
#
#   DGIdb:
#     DGIdb is NOT a transcriptomic reversal tool.
#     It is used only for target-level annotation and mechanistic support.
#
# Main output files:
#   final_single_drug_prioritization_SIGN_CHECKED_AD61026.csv
#   final_single_drug_curated_candidates_SIGN_CHECKED_AD61026.csv
#   final_single_drug_high_caution_SIGN_CHECKED_AD61026.csv
#   final_single_drug_tier_summary_SIGN_CHECKED_AD61026.csv
#
# Recommended run:
#   source("R/06D_final_single_drug_prioritization_AD61026_SIGN_CHECKED.R")
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

required_pkgs <- c("dplyr", "tibble", "tidyr", "stringr", "ggplot2")

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

if (!exists("OUT_DIR")) {
  OUT_DIR <- file.path(PROJECT_ROOT, "Output")
}

if (!exists("FIG_DIR")) {
  FIG_DIR <- file.path(OUT_DIR, "Figures")
}

DRUG_DIR <- file.path(OUT_DIR, "Drug_Prioritization")
PROCESSED_DIR <- file.path(DRUG_DIR, "processed_drug_results")
DREIMT_DIR <- file.path(DRUG_DIR, "DREIMT")

dir.create(DRUG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DREIMT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("  06D SIGN CHECKED: Final single-drug prioritization\n")
cat("============================================================\n\n")


# ============================================================
# SECTION 2: Locate required input files
# ============================================================

clue_path <- file.path(PROCESSED_DIR, "CLUE_CMap_negative_perturbagen_summary.csv")
ilincs_path <- file.path(PROCESSED_DIR, "iLINCS_negative_connected_perturbations.csv")
dgidb_summary_path <- file.path(DRUG_DIR, "dgidb_combined_drug_summary_AD61026.csv")

dreimt_candidate_paths <- c(
  file.path(DRUG_DIR, "microglia_project.csv"),
  file.path(DREIMT_DIR, "microglia_project.csv"),
  file.path(DRUG_DIR, "DREIMT_microglia_project.csv")
)

dreimt_path <- dreimt_candidate_paths[file.exists(dreimt_candidate_paths)][1]

required_files <- c(clue_path, ilincs_path, dgidb_summary_path)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "\nMissing required file(s):\n  ",
    paste(missing_files, collapse = "\n  "),
    "\n"
  )
}

if (is.na(dreimt_path) || length(dreimt_path) == 0) {
  stop(
    "\nDREIMT file not found. Copy microglia_project.csv into:\n  ",
    DRUG_DIR,
    "\n"
  )
}

cat("Input files found:\n")
cat("  CLUE:   ", clue_path, "\n")
cat("  iLINCS: ", ilincs_path, "\n")
cat("  DGIdb:  ", dgidb_summary_path, "\n")
cat("  DREIMT: ", dreimt_path, "\n\n")


# ============================================================
# SECTION 3: Helper functions
# ============================================================

clean_drug_name <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace_all(x, "\\s+", " ")
  x_lower <- str_to_lower(x)

  case_when(
    x_lower %in% c("pg 490", "pg-490", "triptolide") ~ "triptolide / pg 490",
    x_lower %in% c("trichostatin a, streptomyces sp.", "trichostatin a") ~ "trichostatin a",
    x_lower %in% c("cucurbitacin-i", "cucurbitacin i") ~ "cucurbitacin i",
    x_lower %in% c("alvocidib", "flavopiridol") ~ "alvocidib",
    x_lower %in% c("bms-387032", "sns-032") ~ "bms-387032",
    TRUE ~ x_lower
  )
}

safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

collapse_unique <- function(x, max_items = 12) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "NA"]
  x <- unique(x)
  if (length(x) == 0) return(NA_character_)
  paste(head(x, max_items), collapse = "; ")
}

standardize_score <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    return(ifelse(is.na(x), NA_real_, 1))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

classify_drug_class <- function(drug, targets, moa) {
  combined <- str_to_lower(paste(drug, targets, moa, sep = " "))

  case_when(
    str_detect(combined, "cdk|cyclin|cell cycle") ~ "CDK/cell-cycle inhibitor",
    str_detect(combined, "hdac|histone deacetylase") ~ "HDAC/epigenetic inhibitor",
    str_detect(combined, "hsp90") ~ "HSP90 inhibitor",
    str_detect(combined, "top2|top1|topoisomerase|epirubicin|doxorubicin|teniposide|dna binding|dna damage|chromomycin") ~
      "DNA/topoisomerase/cytotoxic",
    str_detect(combined, "rna polymerase|triptolide") ~ "RNA polymerase/toxicity-prone",
    str_detect(combined, "proteasome|psmb|psma|mg-132|bortezomib") ~ "Proteasome inhibitor",
    str_detect(combined, "map2k|mek|trametinib|selumetinib|pd-0325901|mapk|p38") ~ "MAPK pathway inhibitor",
    str_detect(combined, "pi3k|akt|mtor|pik3") ~ "PI3K/AKT/mTOR pathway",
    str_detect(combined, "flt1|kdr|vegf|vegfr|sunitinib|axitinib|cediranib|sorafenib|rtk") ~
      "VEGF/RTK inhibitor",
    str_detect(combined, "stat3|jak") ~ "JAK/STAT pathway",
    str_detect(combined, "hmgcr|statin") ~ "Statin/lipid metabolism",
    str_detect(combined, "fnta|fntb|tipifarnib|farnesyltransferase") ~ "Farnesyltransferase/RAS pathway",
    str_detect(combined, "aldh2|alda") ~ "Mitochondrial/ALDH2 related",
    str_detect(combined, "glucocorticoid|nr3c1|dexamethasone|beclomethasone|ciclesonide|prednisone") ~
      "Glucocorticoid receptor related",
    str_detect(combined, "serotonin|htr1|htr2|htr3") ~ "Serotonergic/neurotransmitter related",
    str_detect(combined, "glutamate|nmda|grin|memantine") ~ "Glutamatergic/neurotransmitter related",
    str_detect(combined, "phosphodiesterase|pde|cilostazol") ~ "Phosphodiesterase/cAMP related",
    str_detect(combined, "metformin|ampk") ~ "Metabolic/AMPK related",
    TRUE ~ "Other/unclear"
  )
}

flag_translation_caution <- function(drug_class, drug, targets, moa) {
  combined <- str_to_lower(paste(drug_class, drug, targets, moa, sep = " "))

  case_when(
    str_detect(combined, "cdk|cell-cycle|topoisomerase|cytotoxic|proteasome|hsp90|dna/topoisomerase|dna binding|rna polymerase|chromomycin|triptolide") ~
      "High caution: cytotoxic/toxicity-prone mechanism",
    str_detect(combined, "hdac|pi3k|akt|mtor|vegf|rtk") ~
      "Moderate caution: broad signaling/oncology-associated mechanism",
    str_detect(combined, "statin|farnesyltransferase|aldh2|glucocorticoid|mapk|phosphodiesterase|serotonergic|glutamatergic|metabolic") ~
      "More interpretable but still requires AD/microglia validation",
    TRUE ~ "Needs manual review"
  )
}


# ============================================================
# SECTION 4: Load raw tables
# ============================================================

cat("Loading raw result tables...\n\n")

clue <- read.csv(clue_path, stringsAsFactors = FALSE, check.names = FALSE)
ilincs <- read.csv(ilincs_path, stringsAsFactors = FALSE, check.names = FALSE)
dgidb_summary <- read.csv(dgidb_summary_path, stringsAsFactors = FALSE, check.names = FALSE)
dreimt <- read.csv(dreimt_path, stringsAsFactors = FALSE, check.names = FALSE)

cat("Rows loaded:\n")
cat("  CLUE negative perturbagen summary:", nrow(clue), "\n")
cat("  iLINCS negative perturbations:", nrow(ilincs), "\n")
cat("  DGIdb drug summary:", nrow(dgidb_summary), "\n")
cat("  DREIMT rows:", nrow(dreimt), "\n\n")


# ============================================================
# SECTION 5: CLUE/CMap sign-checked support
# ============================================================

cat("Standardizing CLUE/CMap evidence with sign-checked criteria...\n")

required_clue_cols <- c("pert_iname", "min_norm_cs", "median_norm_cs", "mean_norm_cs", "n_negative")

if (!all(required_clue_cols %in% colnames(clue))) {
  stop(
    "CLUE file is missing expected columns: ",
    paste(setdiff(required_clue_cols, colnames(clue)), collapse = ", ")
  )
}

clue_std <- clue %>%
  mutate(
    drug_clean = clean_drug_name(pert_iname),
    clue_present = TRUE,
    clue_min_norm_cs = suppressWarnings(as.numeric(min_norm_cs)),
    clue_median_norm_cs = suppressWarnings(as.numeric(median_norm_cs)),
    clue_mean_norm_cs = suppressWarnings(as.numeric(mean_norm_cs)),
    clue_n_negative = suppressWarnings(as.numeric(n_negative)),
    clue_max_fdr_q_nlog10 =
      if ("max_fdr_q_nlog10" %in% colnames(.)) suppressWarnings(as.numeric(max_fdr_q_nlog10)) else NA_real_,
    clue_moa = if ("moa" %in% colnames(.)) moa else NA_character_,
    clue_targets = if ("targets" %in% colnames(.)) targets else NA_character_
  ) %>%
  group_by(drug_clean) %>%
  summarise(
    clue_present = TRUE,
    clue_min_norm_cs = safe_min(clue_min_norm_cs),
    clue_median_norm_cs = safe_min(clue_median_norm_cs),
    clue_mean_norm_cs = safe_mean(clue_mean_norm_cs),
    clue_n_negative = safe_max(clue_n_negative),
    clue_max_fdr_q_nlog10 = safe_max(clue_max_fdr_q_nlog10),
    clue_moa = collapse_unique(clue_moa),
    clue_targets = collapse_unique(clue_targets),
    .groups = "drop"
  ) %>%
  mutate(
    # Any reversal is negative norm_cs, but not all any-reversal hits should be prioritized.
    clue_any_reversal = !is.na(clue_min_norm_cs) & clue_min_norm_cs < 0,

    # Strong support requires repeated negative profiles or consistent negative summary.
    clue_strong_support =
      !is.na(clue_n_negative) &
      clue_n_negative >= 2 &
      (
        (!is.na(clue_median_norm_cs) & clue_median_norm_cs <= -0.70) |
          (!is.na(clue_mean_norm_cs) & clue_mean_norm_cs <= -0.70)
      )
  )

cat("  CLUE drugs with any negative score:", sum(clue_std$clue_any_reversal, na.rm = TRUE), "\n")
cat("  CLUE strong support drugs:", sum(clue_std$clue_strong_support, na.rm = TRUE), "\n\n")


# ============================================================
# SECTION 6: iLINCS sign-checked support
# ============================================================

cat("Standardizing iLINCS evidence with sign-checked criteria...\n")

required_ilincs_cols <- c("Perturbagen", "Correlation", "NoOfSignatures", "pValue", "zScore")

if (!all(required_ilincs_cols %in% colnames(ilincs))) {
  stop(
    "iLINCS file is missing expected columns: ",
    paste(setdiff(required_ilincs_cols, colnames(ilincs)), collapse = ", ")
  )
}

ilincs_std <- ilincs %>%
  filter(Correlation == "-") %>%
  mutate(
    drug_clean = clean_drug_name(Perturbagen),
    ilincs_negative_correlation = TRUE,
    ilincs_pvalue = suppressWarnings(as.numeric(pValue)),
    ilincs_zscore = suppressWarnings(as.numeric(zScore)),
    ilincs_n_signatures = suppressWarnings(as.numeric(NoOfSignatures)),
    ilincs_targets = if ("GeneTargets" %in% colnames(.)) GeneTargets else NA_character_
  ) %>%
  group_by(drug_clean) %>%
  summarise(
    ilincs_present = TRUE,
    ilincs_negative_correlation = TRUE,
    ilincs_best_pvalue = safe_min(ilincs_pvalue),
    ilincs_best_zscore = safe_max(ilincs_zscore),
    ilincs_n_signatures = safe_max(ilincs_n_signatures),
    ilincs_targets = collapse_unique(ilincs_targets),
    .groups = "drop"
  ) %>%
  mutate(
    ilincs_strong_support =
      !is.na(ilincs_best_zscore) &
      !is.na(ilincs_best_pvalue) &
      ilincs_best_zscore > 3 &
      ilincs_best_pvalue < 0.05
  )

cat("  iLINCS negative-correlation drugs:", nrow(ilincs_std), "\n")
cat("  iLINCS strong support drugs:", sum(ilincs_std$ilincs_strong_support, na.rm = TRUE), "\n\n")


# ============================================================
# SECTION 7: DREIMT sign-checked support
# ============================================================

cat("Standardizing DREIMT evidence with sign-checked criteria...\n")

required_dreimt_cols <- c("drug_name", "up_dr", "down_fdr", "tau", "drug_status")

if (!all(required_dreimt_cols %in% colnames(dreimt))) {
  stop(
    "DREIMT file is missing expected columns: ",
    paste(setdiff(required_dreimt_cols, colnames(dreimt)), collapse = ", ")
  )
}

dreimt_clean <- dreimt %>%
  mutate(
    drug_clean = clean_drug_name(drug_name),
    dreimt_tau = suppressWarnings(as.numeric(tau)),
    dreimt_up_dr = suppressWarnings(as.numeric(up_dr)),
    dreimt_down_fdr = suppressWarnings(as.numeric(down_fdr)),
    dreimt_status = drug_status,
    dreimt_moa = if ("drug_moa" %in% colnames(.)) drug_moa else NA_character_,
    dreimt_targets = if ("drug_target_gene_names" %in% colnames(.)) drug_target_gene_names else NA_character_,
    dreimt_source_db = if ("drug_source_db" %in% colnames(.)) drug_source_db else NA_character_,

    # Correct sign:
    #   tau < -75 is reversal/inhibitory support.
    #   tau > +75 would be same-direction/mimic and is not used.
    dreimt_tau_reversal = !is.na(dreimt_tau) & dreimt_tau < -75,
    dreimt_tau_mimic = !is.na(dreimt_tau) & dreimt_tau > 75,

    dreimt_fdr_supported_reversal =
      dreimt_tau_reversal &
      (
        (!is.na(dreimt_up_dr) & dreimt_up_dr < 0.05) |
          (!is.na(dreimt_down_fdr) & dreimt_down_fdr < 0.05)
      ),

    dreimt_approved = str_to_upper(dreimt_status) == "APPROVED"
  )

dreimt_std <- dreimt_clean %>%
  group_by(drug_clean) %>%
  summarise(
    dreimt_present = TRUE,
    dreimt_n_rows = n(),
    dreimt_best_tau = safe_min(dreimt_tau),
    dreimt_best_up_dr = safe_min(dreimt_up_dr),
    dreimt_best_down_fdr = safe_min(dreimt_down_fdr),
    dreimt_tau_support = any(dreimt_tau_reversal, na.rm = TRUE),
    dreimt_strong_support = any(dreimt_fdr_supported_reversal, na.rm = TRUE),
    dreimt_mimic_flag = any(dreimt_tau_mimic, na.rm = TRUE),
    dreimt_any_approved = any(dreimt_approved, na.rm = TRUE),
    dreimt_status = collapse_unique(dreimt_status),
    dreimt_moa = collapse_unique(dreimt_moa),
    dreimt_targets = collapse_unique(dreimt_targets),
    dreimt_source_db = collapse_unique(dreimt_source_db),
    .groups = "drop"
  )

dreimt_summary_path <- file.path(DRUG_DIR, "dreimt_drug_summary_SIGN_CHECKED_AD61026.csv")
write.csv(dreimt_std, dreimt_summary_path, row.names = FALSE)

cat("  DREIMT tau-reversal drugs, tau < -75:", sum(dreimt_std$dreimt_tau_support, na.rm = TRUE), "\n")
cat("  DREIMT FDR-supported reversal drugs:", sum(dreimt_std$dreimt_strong_support, na.rm = TRUE), "\n")
cat("  DREIMT mimic-flag drugs, tau > 75:", sum(dreimt_std$dreimt_mimic_flag, na.rm = TRUE), "\n\n")


# ============================================================
# SECTION 8: DGIdb target-level annotation
# ============================================================

cat("Standardizing DGIdb evidence as target-level annotation only...\n")

if (!"drug" %in% colnames(dgidb_summary)) {
  stop("DGIdb summary file does not contain a column named 'drug'.")
}

dgidb_std <- dgidb_summary %>%
  mutate(
    drug_clean = clean_drug_name(drug),
    dgidb_present = TRUE,
    dgidb_n_genes =
      if ("n_genes" %in% colnames(.)) suppressWarnings(as.numeric(n_genes)) else NA_real_,
    dgidb_relevant_interactions =
      if ("n_relevant_direction_interactions" %in% colnames(.)) {
        suppressWarnings(as.numeric(n_relevant_direction_interactions))
      } else {
        NA_real_
      },
    dgidb_target_genes =
      if ("target_genes" %in% colnames(.)) target_genes else NA_character_,
    dgidb_interaction_types =
      if ("interaction_types" %in% colnames(.)) interaction_types else NA_character_
  ) %>%
  group_by(drug_clean) %>%
  summarise(
    dgidb_present = TRUE,
    dgidb_n_genes = safe_max(dgidb_n_genes),
    dgidb_relevant_interactions = safe_max(dgidb_relevant_interactions),
    dgidb_target_genes = collapse_unique(dgidb_target_genes),
    dgidb_interaction_types = collapse_unique(dgidb_interaction_types),
    .groups = "drop"
  )

cat("  DGIdb annotated drugs:", nrow(dgidb_std), "\n\n")


# ============================================================
# SECTION 9: Merge sign-checked evidence
# ============================================================

cat("Merging sign-checked evidence layers...\n\n")

all_drugs <- unique(c(
  clue_std$drug_clean,
  ilincs_std$drug_clean,
  dreimt_std$drug_clean,
  dgidb_std$drug_clean
))

final_drugs <- data.frame(drug_clean = all_drugs, stringsAsFactors = FALSE) %>%
  left_join(clue_std, by = "drug_clean") %>%
  left_join(ilincs_std, by = "drug_clean") %>%
  left_join(dreimt_std, by = "drug_clean") %>%
  left_join(dgidb_std, by = "drug_clean") %>%
  mutate(
    clue_strong_support = ifelse(is.na(clue_strong_support), FALSE, clue_strong_support),
    clue_any_reversal = ifelse(is.na(clue_any_reversal), FALSE, clue_any_reversal),

    ilincs_strong_support = ifelse(is.na(ilincs_strong_support), FALSE, ilincs_strong_support),
    ilincs_negative_correlation = ifelse(is.na(ilincs_negative_correlation), FALSE, ilincs_negative_correlation),

    dreimt_tau_support = ifelse(is.na(dreimt_tau_support), FALSE, dreimt_tau_support),
    dreimt_strong_support = ifelse(is.na(dreimt_strong_support), FALSE, dreimt_strong_support),
    dreimt_mimic_flag = ifelse(is.na(dreimt_mimic_flag), FALSE, dreimt_mimic_flag),
    dreimt_any_approved = ifelse(is.na(dreimt_any_approved), FALSE, dreimt_any_approved),

    dgidb_present = ifelse(is.na(dgidb_present), FALSE, dgidb_present),

    n_transcriptomic_tools =
      as.integer(clue_strong_support) +
      as.integer(ilincs_strong_support) +
      as.integer(dreimt_tau_support),

    n_transcriptomic_tools_fdr_strict =
      as.integer(clue_strong_support) +
      as.integer(ilincs_strong_support) +
      as.integer(dreimt_strong_support),

    n_total_support_sources =
      n_transcriptomic_tools +
      as.integer(dgidb_present),

    combined_targets = paste(
      ifelse(is.na(clue_targets), "", clue_targets),
      ifelse(is.na(ilincs_targets), "", ilincs_targets),
      ifelse(is.na(dreimt_targets), "", dreimt_targets),
      ifelse(is.na(dgidb_target_genes), "", dgidb_target_genes),
      sep = "; "
    ),

    combined_moa = paste(
      ifelse(is.na(clue_moa), "", clue_moa),
      ifelse(is.na(dreimt_moa), "", dreimt_moa),
      ifelse(is.na(dgidb_interaction_types), "", dgidb_interaction_types),
      sep = "; "
    ),

    drug_class = classify_drug_class(drug_clean, combined_targets, combined_moa),

    translation_caution =
      flag_translation_caution(drug_class, drug_clean, combined_targets, combined_moa)
  )


# ============================================================
# SECTION 10: Score and assign tiers
# ============================================================

final_drugs <- final_drugs %>%
  mutate(
    # Strength components use correct sign conventions:
    #   CLUE: multiply by -1 so more negative becomes larger strength.
    #   iLINCS: positive zScore among negative-correlation rows is stronger.
    #   DREIMT: multiply tau by -1 so more negative tau becomes larger strength.
    clue_strength_component = standardize_score(-1 * clue_median_norm_cs),
    ilincs_strength_component = standardize_score(ilincs_best_zscore),
    dreimt_strength_component = standardize_score(-1 * dreimt_best_tau),
    dgidb_strength_component = standardize_score(dgidb_relevant_interactions),

    high_caution =
      translation_caution == "High caution: cytotoxic/toxicity-prone mechanism",

    moderate_caution =
      translation_caution == "Moderate caution: broad signaling/oncology-associated mechanism",

    final_priority_score =
      2.0 * n_transcriptomic_tools +
      1.0 * n_transcriptomic_tools_fdr_strict +
      0.4 * as.integer(dgidb_present) +
      0.5 * as.integer(dreimt_any_approved) +
      0.5 * ifelse(is.na(clue_strength_component), 0, clue_strength_component) +
      0.5 * ifelse(is.na(ilincs_strength_component), 0, ilincs_strength_component) +
      0.5 * ifelse(is.na(dreimt_strength_component), 0, dreimt_strength_component) +
      0.3 * ifelse(is.na(dgidb_strength_component), 0, dgidb_strength_component) -
      1.2 * as.integer(high_caution) -
      0.5 * as.integer(moderate_caution),

    final_tier = case_when(
      n_transcriptomic_tools == 3 &
        dreimt_strong_support &
        !high_caution ~
        "Tier 1A: all 3 transcriptomic tools, DREIMT FDR-supported, not high caution",

      n_transcriptomic_tools == 3 &
        !high_caution ~
        "Tier 1B: all 3 transcriptomic tools, not high caution",

      n_transcriptomic_tools == 3 &
        high_caution ~
        "High-caution: all 3 transcriptomic tools but toxicity-prone",

      n_transcriptomic_tools == 2 &
        !high_caution ~
        "Tier 2: two transcriptomic tools, not high caution",

      n_transcriptomic_tools == 2 &
        high_caution ~
        "High-caution: two transcriptomic tools but toxicity-prone",

      n_transcriptomic_tools == 1 &
        dgidb_present &
        !high_caution ~
        "Tier 3: one transcriptomic tool plus DGIdb target support",

      n_transcriptomic_tools == 1 &
        !high_caution ~
        "Tier 4: one transcriptomic tool only",

      dgidb_present ~
        "Target-level only: DGIdb annotation",

      TRUE ~ "Background/low priority"
    )
  ) %>%
  arrange(
    desc(n_transcriptomic_tools),
    desc(n_transcriptomic_tools_fdr_strict),
    desc(dreimt_any_approved),
    high_caution,
    desc(final_priority_score),
    drug_clean
  )


# ============================================================
# SECTION 11: Save output tables
# ============================================================

final_path <- file.path(DRUG_DIR, "final_single_drug_prioritization_SIGN_CHECKED_AD61026.csv")
curated_path <- file.path(DRUG_DIR, "final_single_drug_curated_candidates_SIGN_CHECKED_AD61026.csv")
high_caution_path <- file.path(DRUG_DIR, "final_single_drug_high_caution_SIGN_CHECKED_AD61026.csv")
tier_summary_path <- file.path(DRUG_DIR, "final_single_drug_tier_summary_SIGN_CHECKED_AD61026.csv")
class_summary_path <- file.path(DRUG_DIR, "final_single_drug_class_summary_SIGN_CHECKED_AD61026.csv")

write.csv(final_drugs, final_path, row.names = FALSE)

curated_candidates <- final_drugs %>%
  filter(
    final_tier %in% c(
      "Tier 1A: all 3 transcriptomic tools, DREIMT FDR-supported, not high caution",
      "Tier 1B: all 3 transcriptomic tools, not high caution",
      "Tier 2: two transcriptomic tools, not high caution",
      "Tier 3: one transcriptomic tool plus DGIdb target support"
    )
  ) %>%
  filter(!str_detect(drug_clean, "^brd-|^mls|^chemb|^ac1|^spectrum|^prestwick|^san")) %>%
  arrange(
    final_tier,
    desc(n_transcriptomic_tools),
    desc(n_transcriptomic_tools_fdr_strict),
    desc(dreimt_any_approved),
    desc(final_priority_score)
  )

write.csv(curated_candidates, curated_path, row.names = FALSE)

high_caution_candidates <- final_drugs %>%
  filter(str_detect(final_tier, "High-caution")) %>%
  arrange(desc(n_transcriptomic_tools), desc(final_priority_score))

write.csv(high_caution_candidates, high_caution_path, row.names = FALSE)

tier_summary <- final_drugs %>%
  group_by(final_tier) %>%
  summarise(
    n_drugs = n(),
    n_clue = sum(clue_strong_support, na.rm = TRUE),
    n_ilincs = sum(ilincs_strong_support, na.rm = TRUE),
    n_dreimt_tau = sum(dreimt_tau_support, na.rm = TRUE),
    n_dreimt_fdr = sum(dreimt_strong_support, na.rm = TRUE),
    n_dgidb = sum(dgidb_present, na.rm = TRUE),
    examples = paste(head(drug_clean, 10), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_drugs))

write.csv(tier_summary, tier_summary_path, row.names = FALSE)

class_summary <- final_drugs %>%
  group_by(drug_class, translation_caution) %>%
  summarise(
    n_drugs = n(),
    n_transcriptomic_3 = sum(n_transcriptomic_tools == 3, na.rm = TRUE),
    n_transcriptomic_2 = sum(n_transcriptomic_tools == 2, na.rm = TRUE),
    n_transcriptomic_1 = sum(n_transcriptomic_tools == 1, na.rm = TRUE),
    examples = paste(head(drug_clean, 10), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_transcriptomic_3), desc(n_transcriptomic_2), desc(n_drugs))

write.csv(class_summary, class_summary_path, row.names = FALSE)


# ============================================================
# SECTION 12: Save simple draft figure
# ============================================================

plot_df <- curated_candidates %>%
  slice_head(n = 30) %>%
  mutate(drug_clean = factor(drug_clean, levels = rev(drug_clean)))

if (nrow(plot_df) > 0) {
  p <- ggplot(plot_df, aes(x = drug_clean, y = final_priority_score)) +
    geom_col() +
    coord_flip() +
    labs(
      title = "Sign-checked final single-drug candidates",
      subtitle = "Correct reversal direction used for CLUE/CMap, iLINCS, and DREIMT",
      x = "Drug candidate",
      y = "Priority score"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(
    file.path(FIG_DIR, "Fig_final_drug_evidence_SIGN_CHECKED_AD61026_draft.pdf"),
    p,
    width = 8,
    height = 7
  )
}


# ============================================================
# SECTION 13: Print summary
# ============================================================

cat("============================================================\n")
cat("06D SIGN CHECKED complete.\n\n")

cat("Final merged table rows:", nrow(final_drugs), "\n")
cat("Curated candidate rows:", nrow(curated_candidates), "\n")
cat("High-caution candidate rows:", nrow(high_caution_candidates), "\n\n")

cat("Support counts using sign-checked criteria:\n")
cat("  CLUE strong support:", sum(final_drugs$clue_strong_support, na.rm = TRUE), "\n")
cat("  iLINCS strong support:", sum(final_drugs$ilincs_strong_support, na.rm = TRUE), "\n")
cat("  DREIMT tau support:", sum(final_drugs$dreimt_tau_support, na.rm = TRUE), "\n")
cat("  DREIMT FDR-supported reversal:", sum(final_drugs$dreimt_strong_support, na.rm = TRUE), "\n")
cat("  DGIdb annotation:", sum(final_drugs$dgidb_present, na.rm = TRUE), "\n\n")

cat("Top sign-checked curated candidates:\n")
print(
  curated_candidates %>%
    select(
      drug_clean,
      final_tier,
      n_transcriptomic_tools,
      n_transcriptomic_tools_fdr_strict,
      dreimt_any_approved,
      drug_class,
      translation_caution,
      clue_median_norm_cs,
      ilincs_best_zscore,
      dreimt_best_tau,
      dreimt_best_up_dr,
      dreimt_best_down_fdr,
      final_priority_score
    ) %>%
    slice_head(n = 25)
)

cat("\nHigh-caution transcriptomic hits are saved separately and should not be treated as primary translational candidates.\n\n")

cat("Output files:\n")
cat("  ", final_path, "\n")
cat("  ", curated_path, "\n")
cat("  ", high_caution_path, "\n")
cat("  ", tier_summary_path, "\n")
cat("  ", class_summary_path, "\n")
cat("  ", dreimt_summary_path, "\n")
cat("============================================================\n")
