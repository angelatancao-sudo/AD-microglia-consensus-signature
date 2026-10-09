# ============================================================
# 15d_drug_name_QC.R
# ============================================================
# Purpose:
#   Perform a final QC pass on the sign-checked single-drug table and
#   collapse likely duplicate drug names caused by punctuation differences.
#
# Why this is needed:
#   Different tools may write the same drug with different punctuation:
#      ALDA-1 vs ALDA 1
#      CGP-60474 vs CGP 60474
#      GSK-1059615 vs GSK 1059615
#
#   These should be treated as the same perturbagen/drug when counting
#   cross-tool support.
#
# Input:
#   Output/Drug_Prioritization/final_single_drug_prioritization_SIGN_CHECKED_AD61026.csv
#
# Outputs:
#   Output/Drug_Prioritization/final_single_drug_prioritization_QC_COLLAPSED_AD61026.csv
#   Output/Drug_Prioritization/final_single_drug_curated_candidates_QC_COLLAPSED_AD61026.csv
#   Output/Drug_Prioritization/final_single_drug_high_caution_QC_COLLAPSED_AD61026.csv
#   Output/Drug_Prioritization/final_single_drug_QC_summary_AD61026.csv
#
# Recommended run:
#   source("R/06D_QC_and_collapse_drug_names_AD61026.R")
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

required_pkgs <- c("dplyr", "stringr", "tibble", "tidyr")

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

if (!exists("OUT_DIR")) {
  OUT_DIR <- file.path(PROJECT_ROOT, "Output")
}

DRUG_DIR <- file.path(OUT_DIR, "Drug_Prioritization")

input_path <- file.path(
  DRUG_DIR,
  "final_single_drug_prioritization_SIGN_CHECKED_AD61026.csv"
)

if (!file.exists(input_path)) {
  stop("Missing input file: ", input_path)
}

cat("============================================================\n")
cat("  06D QC: collapse duplicate drug names and verify signs\n")
cat("============================================================\n\n")

cat("Input file:\n")
cat("  ", input_path, "\n\n")


# ============================================================
# SECTION 2: Helper functions
# ============================================================

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  x <- as.character(x)
  x <- tolower(x)
  x %in% c("true", "t", "1", "yes")
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

collapse_unique <- function(x, max_items = 20) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "NA"]
  x <- unique(x)
  if (length(x) == 0) return(NA_character_)
  paste(head(x, max_items), collapse = "; ")
}

make_drug_key <- function(x) {
  # Remove punctuation, spaces, and capitalization differences.
  # This treats ALDA-1 and ALDA 1 as the same key.
  x <- tolower(as.character(x))
  gsub("[^a-z0-9]", "", x)
}

choose_display_name <- function(x) {
  # Choose a readable display name from all variants in a duplicate group.
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)

  # Prefer names that are not BRD/CHEMBL-like screening IDs.
  x_preferred <- x[!str_detect(x, "^brd-|^brd |^chemb|^mls|^ac1|^spectrum|^prestwick|^san")]
  if (length(x_preferred) == 0) x_preferred <- x

  # Prefer the shortest readable version.
  x_preferred[which.min(nchar(x_preferred))]
}

classify_tier <- function(n_tools, n_tools_strict, dgidb_present, high_caution, dreimt_strong_support) {
  dplyr::case_when(
    n_tools == 3 & dreimt_strong_support & !high_caution ~
      "Tier 1A: all 3 transcriptomic tools, DREIMT FDR-supported, not high caution",

    n_tools == 3 & !high_caution ~
      "Tier 1B: all 3 transcriptomic tools, not high caution",

    n_tools == 3 & high_caution ~
      "High-caution: all 3 transcriptomic tools but toxicity-prone",

    n_tools == 2 & !high_caution ~
      "Tier 2: two transcriptomic tools, not high caution",

    n_tools == 2 & high_caution ~
      "High-caution: two transcriptomic tools but toxicity-prone",

    n_tools == 1 & high_caution ~
      "High-caution: one transcriptomic tool but toxicity-prone",

    n_tools == 1 & dgidb_present & !high_caution ~
      "Tier 3: one transcriptomic tool plus DGIdb target support",

    n_tools == 1 & !high_caution ~
      "Tier 4: one transcriptomic tool only",

    n_tools == 0 & dgidb_present ~
      "Target-level only: DGIdb annotation",

    TRUE ~ "Background/low priority"
  )
}


# ============================================================
# SECTION 3: Load sign-checked table and verify sign criteria
# ============================================================

df <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)

# Convert important logical columns safely because CSVs sometimes store
# TRUE/FALSE as text.
logical_cols <- c(
  "clue_strong_support",
  "clue_any_reversal",
  "ilincs_strong_support",
  "ilincs_negative_correlation",
  "dreimt_tau_support",
  "dreimt_strong_support",
  "dreimt_mimic_flag",
  "dreimt_any_approved",
  "dgidb_present",
  "high_caution",
  "moderate_caution"
)

for (cc in logical_cols) {
  if (cc %in% colnames(df)) {
    df[[cc]] <- as_logical_safe(df[[cc]])
  }
}

# Sign and criterion checks.
qc_clue_violations <- df %>%
  filter(
    clue_strong_support,
    !(
      clue_n_negative >= 2 &
        (
          clue_median_norm_cs <= -0.70 |
            clue_mean_norm_cs <= -0.70
        )
    )
  )

qc_ilincs_violations <- df %>%
  filter(
    ilincs_strong_support,
    !(ilincs_best_zscore > 3 & ilincs_best_pvalue < 0.05)
  )

qc_dreimt_tau_violations <- df %>%
  filter(
    dreimt_tau_support,
    !(dreimt_best_tau < -75)
  )

qc_dreimt_strong_violations <- df %>%
  filter(
    dreimt_strong_support,
    !(
      dreimt_best_tau < -75 &
        (
          dreimt_best_up_dr < 0.05 |
            dreimt_best_down_fdr < 0.05
        )
    )
  )

recomputed_n_tools <-
  as.integer(df$clue_strong_support) +
  as.integer(df$ilincs_strong_support) +
  as.integer(df$dreimt_tau_support)

n_tool_mismatch <- sum(recomputed_n_tools != df$n_transcriptomic_tools, na.rm = TRUE)

cat("Sign/criterion QC checks:\n")
cat("  CLUE strong-support violations:", nrow(qc_clue_violations), "\n")
cat("  iLINCS strong-support violations:", nrow(qc_ilincs_violations), "\n")
cat("  DREIMT tau-support violations:", nrow(qc_dreimt_tau_violations), "\n")
cat("  DREIMT FDR-support violations:", nrow(qc_dreimt_strong_violations), "\n")
cat("  Transcriptomic tool count mismatches:", n_tool_mismatch, "\n\n")


# ============================================================
# SECTION 4: Collapse duplicate drug names
# ============================================================

df2 <- df %>%
  mutate(
    drug_key = make_drug_key(drug_clean)
  )

duplicate_groups <- df2 %>%
  group_by(drug_key) %>%
  summarise(
    n_name_variants = n_distinct(drug_clean),
    name_variants = collapse_unique(drug_clean),
    any_support = any(n_total_support_sources > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_name_variants > 1, any_support)

duplicate_path <- file.path(
  DRUG_DIR,
  "final_single_drug_duplicate_name_groups_QC_AD61026.csv"
)

write.csv(duplicate_groups, duplicate_path, row.names = FALSE)

collapsed <- df2 %>%
  group_by(drug_key) %>%
  summarise(
    drug_clean = choose_display_name(drug_clean),
    name_variants = collapse_unique(drug_clean),

    clue_present = any(clue_present, na.rm = TRUE),
    clue_any_reversal = any(clue_any_reversal, na.rm = TRUE),
    clue_strong_support = any(clue_strong_support, na.rm = TRUE),
    clue_min_norm_cs = safe_min(clue_min_norm_cs),
    clue_median_norm_cs = safe_min(clue_median_norm_cs),
    clue_mean_norm_cs = safe_mean(clue_mean_norm_cs),
    clue_n_negative = safe_max(clue_n_negative),
    clue_max_fdr_q_nlog10 = safe_max(clue_max_fdr_q_nlog10),
    clue_moa = collapse_unique(clue_moa),
    clue_targets = collapse_unique(clue_targets),

    ilincs_present = any(ilincs_present, na.rm = TRUE),
    ilincs_negative_correlation = any(ilincs_negative_correlation, na.rm = TRUE),
    ilincs_strong_support = any(ilincs_strong_support, na.rm = TRUE),
    ilincs_best_pvalue = safe_min(ilincs_best_pvalue),
    ilincs_best_zscore = safe_max(ilincs_best_zscore),
    ilincs_n_signatures = safe_max(ilincs_n_signatures),
    ilincs_targets = collapse_unique(ilincs_targets),

    dreimt_present = any(dreimt_present, na.rm = TRUE),
    dreimt_n_rows = safe_max(dreimt_n_rows),
    dreimt_best_tau = safe_min(dreimt_best_tau),
    dreimt_best_up_dr = safe_min(dreimt_best_up_dr),
    dreimt_best_down_fdr = safe_min(dreimt_best_down_fdr),
    dreimt_tau_support = any(dreimt_tau_support, na.rm = TRUE),
    dreimt_strong_support = any(dreimt_strong_support, na.rm = TRUE),
    dreimt_mimic_flag = any(dreimt_mimic_flag, na.rm = TRUE),
    dreimt_any_approved = any(dreimt_any_approved, na.rm = TRUE),
    dreimt_status = collapse_unique(dreimt_status),
    dreimt_moa = collapse_unique(dreimt_moa),
    dreimt_targets = collapse_unique(dreimt_targets),
    dreimt_source_db = collapse_unique(dreimt_source_db),

    dgidb_present = any(dgidb_present, na.rm = TRUE),
    dgidb_n_genes = safe_max(dgidb_n_genes),
    dgidb_relevant_interactions = safe_max(dgidb_relevant_interactions),
    dgidb_target_genes = collapse_unique(dgidb_target_genes),
    dgidb_interaction_types = collapse_unique(dgidb_interaction_types),

    combined_targets = collapse_unique(combined_targets),
    combined_moa = collapse_unique(combined_moa),
    drug_class = collapse_unique(drug_class),
    translation_caution = collapse_unique(translation_caution),

    .groups = "drop"
  ) %>%
  mutate(
    high_caution = str_detect(translation_caution, "High caution"),
    moderate_caution = str_detect(translation_caution, "Moderate caution"),

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

    final_tier = classify_tier(
      n_tools = n_transcriptomic_tools,
      n_tools_strict = n_transcriptomic_tools_fdr_strict,
      dgidb_present = dgidb_present,
      high_caution = high_caution,
      dreimt_strong_support = dreimt_strong_support
    ),

    # Recalculate a simple priority score after name collapsing.
    clue_strength_component = ifelse(is.na(clue_median_norm_cs), 0, -1 * clue_median_norm_cs),
    ilincs_strength_component = ifelse(is.na(ilincs_best_zscore), 0, ilincs_best_zscore / 10),
    dreimt_strength_component = ifelse(is.na(dreimt_best_tau), 0, -1 * dreimt_best_tau / 100),
    dgidb_strength_component = ifelse(is.na(dgidb_relevant_interactions), 0, dgidb_relevant_interactions / 10),

    final_priority_score =
      2.0 * n_transcriptomic_tools +
      1.0 * n_transcriptomic_tools_fdr_strict +
      0.4 * as.integer(dgidb_present) +
      0.5 * as.integer(dreimt_any_approved) +
      0.4 * clue_strength_component +
      0.4 * ilincs_strength_component +
      0.4 * dreimt_strength_component +
      0.2 * dgidb_strength_component -
      1.2 * as.integer(high_caution) -
      0.5 * as.integer(moderate_caution)
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
# SECTION 5: Build final review tables
# ============================================================

collapsed_path <- file.path(
  DRUG_DIR,
  "final_single_drug_prioritization_QC_COLLAPSED_AD61026.csv"
)

curated_path <- file.path(
  DRUG_DIR,
  "final_single_drug_curated_candidates_QC_COLLAPSED_AD61026.csv"
)

high_caution_path <- file.path(
  DRUG_DIR,
  "final_single_drug_high_caution_QC_COLLAPSED_AD61026.csv"
)

tier_summary_path <- file.path(
  DRUG_DIR,
  "final_single_drug_tier_summary_QC_COLLAPSED_AD61026.csv"
)

qc_summary_path <- file.path(
  DRUG_DIR,
  "final_single_drug_QC_summary_AD61026.csv"
)

write.csv(collapsed, collapsed_path, row.names = FALSE)

curated <- collapsed %>%
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
    desc(n_transcriptomic_tools),
    desc(n_transcriptomic_tools_fdr_strict),
    desc(dreimt_any_approved),
    desc(final_priority_score)
  )

write.csv(curated, curated_path, row.names = FALSE)

high_caution <- collapsed %>%
  filter(str_detect(final_tier, "High-caution")) %>%
  arrange(
    desc(n_transcriptomic_tools),
    desc(n_transcriptomic_tools_fdr_strict),
    desc(final_priority_score)
  )

write.csv(high_caution, high_caution_path, row.names = FALSE)

tier_summary <- collapsed %>%
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

qc_summary <- tibble::tibble(
  qc_item = c(
    "Input rows before name collapsing",
    "Rows after name collapsing",
    "Duplicate name groups with support",
    "CLUE strong-support criterion violations",
    "iLINCS strong-support criterion violations",
    "DREIMT tau-support criterion violations",
    "DREIMT FDR-support criterion violations",
    "Transcriptomic tool count mismatches"
  ),
  value = c(
    nrow(df),
    nrow(collapsed),
    nrow(duplicate_groups),
    nrow(qc_clue_violations),
    nrow(qc_ilincs_violations),
    nrow(qc_dreimt_tau_violations),
    nrow(qc_dreimt_strong_violations),
    n_tool_mismatch
  )
)

write.csv(qc_summary, qc_summary_path, row.names = FALSE)


# ============================================================
# SECTION 6: Print completion summary
# ============================================================

cat("QC complete.\n\n")

print(qc_summary)

cat("\nTier summary after name collapsing:\n")
print(tier_summary)

cat("\nTop curated candidates after name collapsing:\n")
print(
  curated %>%
    select(
      drug_clean,
      name_variants,
      final_tier,
      n_transcriptomic_tools,
      n_transcriptomic_tools_fdr_strict,
      n_total_support_sources,
      dreimt_any_approved,
      drug_class,
      translation_caution,
      clue_median_norm_cs,
      ilincs_best_zscore,
      dreimt_best_tau,
      dgidb_present,
      final_priority_score
    ) %>%
    slice_head(n = 30)
)

cat("\nOutput files:\n")
cat("  ", collapsed_path, "\n")
cat("  ", curated_path, "\n")
cat("  ", high_caution_path, "\n")
cat("  ", tier_summary_path, "\n")
cat("  ", duplicate_path, "\n")
cat("  ", qc_summary_path, "\n")
cat("============================================================\n")
