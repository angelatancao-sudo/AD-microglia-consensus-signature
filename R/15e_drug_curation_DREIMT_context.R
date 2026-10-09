# ============================================================
# 15e_drug_curation_DREIMT_context.R
# ============================================================
# Project: AD Microglia Multi-Dataset Transcriptomic Analysis
#
# Purpose:
#   Add a stricter, literature-aligned DREIMT context classification
#   before finalizing the single-drug list for pair prioritization.
#
# Why this script is needed:
#   DREIMT can return multiple rows for the same drug. A drug may have:
#     - reversal rows: tau < -80
#     - mimic/same-direction rows: tau > +80
#     - both reversal and mimic rows across different perturbation contexts
#
#   Therefore, we should not simply count any negative tau row and ignore
#   positive tau rows. This script classifies DREIMT support as:
#
#     1. clean_reversal
#        tau < -80 present, no tau > +80 mimic row
#
#     2. mixed_reversal_and_mimic
#        tau < -80 present and tau > +80 present
#
#     3. discordant_mimic_only
#        tau > +80 present, no tau < -80 reversal row
#
#     4. no_strong_dreimt_signal
#        no tau < -80 and no tau > +80
#
#   It also separates FDR-supported DREIMT reversal from tau-only reversal.
#
# Correct sign conventions:
#   CLUE/CMap:
#     negative normalized connectivity = reversal
#
#   iLINCS:
#     Correlation = "-" = reversal
#     zScore > 3 and pValue < 0.05 = strong support
#
#   DREIMT:
#     tau < -80 = reversal / inhibitory direction
#     tau > +80 = mimic / same-direction
#     tau near 0 = no strong direction
#
#   DGIdb:
#     target-level annotation only
#     not transcriptomic reversal
#
# Inputs:
#   Output/Drug_Prioritization/final_single_drug_prioritization_QC_COLLAPSED_AD61026.csv
#   Output/Drug_Prioritization/microglia_project.csv
#
# Outputs:
#   Output/Drug_Prioritization/final_single_drug_DREIMT_CONTEXT_CLASSIFIED_AD61026.csv
#   Output/Drug_Prioritization/final_single_drug_pair_pool_PRIMARY_AD61026.csv
#   Output/Drug_Prioritization/final_single_drug_pair_pool_CONTEXT_DEPENDENT_AD61026.csv
#   Output/Drug_Prioritization/final_single_drug_pair_pool_EXPLORATORY_AD61026.csv
#   Output/Drug_Prioritization/final_single_drug_DREIMT_context_summary_AD61026.csv
#   Output/Drug_Prioritization/final_single_drug_FINAL_CURATION_SUMMARY_AD61026.csv
#
# Recommended run:
#   source("R/06F_DREIMT_context_classification_and_final_drug_curation_AD61026.R")
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

# Note: the "pair pool" files written by this script are the three candidate
# categories (primary, context-dependent, exploratory) used for Table 4. The
# names come from an earlier drug-pair analysis that is not in the published
# paper; no drug pairs are generated here.

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

final_table_path <- file.path(
  DRUG_DIR,
  "final_single_drug_prioritization_QC_COLLAPSED_AD61026.csv"
)

dreimt_path_candidates <- c(
  file.path(DRUG_DIR, "microglia_project.csv"),
  file.path(DRUG_DIR, "DREIMT", "microglia_project.csv")
)

dreimt_path <- dreimt_path_candidates[file.exists(dreimt_path_candidates)][1]

if (!file.exists(final_table_path)) {
  stop("Missing final QC-collapsed drug table: ", final_table_path)
}

if (is.na(dreimt_path) || length(dreimt_path) == 0) {
  stop("Missing DREIMT raw file. Expected microglia_project.csv in Drug_Prioritization folder.")
}

cat("============================================================\n")
cat("  06F: DREIMT context classification and final curation\n")
cat("============================================================\n\n")

cat("Final QC-collapsed table:\n")
cat("  ", final_table_path, "\n")
cat("DREIMT raw file:\n")
cat("  ", dreimt_path, "\n\n")


# ============================================================
# SECTION 2: Helper functions
# ============================================================

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  x <- tolower(as.character(x))
  x %in% c("true", "t", "1", "yes")
}

clean_drug_name <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace_all(x, "\\s+", " ")
  x_lower <- str_to_lower(x)

  dplyr::case_when(
    x_lower %in% c("pg 490", "pg-490", "triptolide") ~ "triptolide / pg 490",
    x_lower %in% c("trichostatin a, streptomyces sp.", "trichostatin a") ~ "trichostatin a",
    x_lower %in% c("cucurbitacin-i", "cucurbitacin i") ~ "cucurbitacin i",
    x_lower %in% c("alvocidib", "flavopiridol") ~ "alvocidib",
    x_lower %in% c("bms-387032", "sns-032") ~ "bms-387032",
    x_lower %in% c("alda 1", "alda-1") ~ "alda-1",
    x_lower %in% c("valproic-acid", "valproic acid", "valproate", "sodium valproate", "divalproex") ~ "valproic acid",
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

collapse_unique <- function(x, max_items = 20) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "NA"]
  x <- unique(x)
  if (length(x) == 0) return(NA_character_)
  paste(head(x, max_items), collapse = "; ")
}

make_drug_key <- function(x) {
  x <- tolower(as.character(x))
  gsub("[^a-z0-9]", "", x)
}


# ============================================================
# SECTION 3: Load and normalize tables
# ============================================================

final_drugs <- read.csv(final_table_path, stringsAsFactors = FALSE, check.names = FALSE)
dreimt_raw <- read.csv(dreimt_path, stringsAsFactors = FALSE, check.names = FALSE)

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
  if (cc %in% colnames(final_drugs)) {
    final_drugs[[cc]] <- as_logical_safe(final_drugs[[cc]])
  }
}

final_drugs <- final_drugs %>%
  mutate(
    drug_clean_original = drug_clean,
    drug_clean = clean_drug_name(drug_clean),
    drug_key = make_drug_key(drug_clean)
  )

dreimt_raw <- dreimt_raw %>%
  mutate(
    drug_clean = clean_drug_name(drug_name),
    drug_key = make_drug_key(drug_clean),
    tau = suppressWarnings(as.numeric(tau)),
    up_dr = suppressWarnings(as.numeric(up_dr)),
    down_fdr = suppressWarnings(as.numeric(down_fdr)),
    dreimt_row_reversal_80 = !is.na(tau) & tau < -80,
    dreimt_row_mimic_80 = !is.na(tau) & tau > 80,
    dreimt_row_fdr_reversal_80 =
      dreimt_row_reversal_80 &
      (
        (!is.na(up_dr) & up_dr < 0.05) |
          (!is.na(down_fdr) & down_fdr < 0.05)
      )
  )


# ============================================================
# SECTION 4: Build DREIMT context classification
# ============================================================

dreimt_context <- dreimt_raw %>%
  group_by(drug_key) %>%
  summarise(
    drug_clean_dreimt = first(drug_clean),
    dreimt_n_rows_raw = n(),

    dreimt_has_reversal_80 = any(dreimt_row_reversal_80, na.rm = TRUE),
    dreimt_has_mimic_80 = any(dreimt_row_mimic_80, na.rm = TRUE),
    dreimt_has_fdr_reversal_80 = any(dreimt_row_fdr_reversal_80, na.rm = TRUE),

    dreimt_best_reversal_tau_80 =
      ifelse(any(dreimt_row_reversal_80, na.rm = TRUE),
             min(tau[dreimt_row_reversal_80], na.rm = TRUE),
             NA_real_),

    dreimt_best_reversal_up_dr_80 =
      ifelse(any(dreimt_row_reversal_80, na.rm = TRUE),
             up_dr[which.min(ifelse(dreimt_row_reversal_80, tau, Inf))],
             NA_real_),

    dreimt_best_reversal_down_fdr_80 =
      ifelse(any(dreimt_row_reversal_80, na.rm = TRUE),
             down_fdr[which.min(ifelse(dreimt_row_reversal_80, tau, Inf))],
             NA_real_),

    dreimt_strongest_mimic_tau_80 =
      ifelse(any(dreimt_row_mimic_80, na.rm = TRUE),
             max(tau[dreimt_row_mimic_80], na.rm = TRUE),
             NA_real_),

    dreimt_status_raw = collapse_unique(drug_status),
    dreimt_moa_raw = collapse_unique(drug_moa),
    dreimt_targets_raw = collapse_unique(drug_target_gene_names),

    .groups = "drop"
  ) %>%
  mutate(
    dreimt_context_class = case_when(
      dreimt_has_reversal_80 & !dreimt_has_mimic_80 ~ "clean_reversal",
      dreimt_has_reversal_80 & dreimt_has_mimic_80 ~ "mixed_reversal_and_mimic",
      !dreimt_has_reversal_80 & dreimt_has_mimic_80 ~ "discordant_mimic_only",
      TRUE ~ "no_strong_dreimt_signal"
    ),

    dreimt_support_level = case_when(
      dreimt_context_class == "clean_reversal" & dreimt_has_fdr_reversal_80 ~
        "clean_FDR_supported_reversal",
      dreimt_context_class == "clean_reversal" & !dreimt_has_fdr_reversal_80 ~
        "clean_tau_only_reversal",
      dreimt_context_class == "mixed_reversal_and_mimic" & dreimt_has_fdr_reversal_80 ~
        "mixed_FDR_supported_reversal",
      dreimt_context_class == "mixed_reversal_and_mimic" & !dreimt_has_fdr_reversal_80 ~
        "mixed_tau_only_reversal",
      dreimt_context_class == "discordant_mimic_only" ~
        "mimic_only_not_reversal",
      TRUE ~
        "no_strong_dreimt_signal"
    )
  )


# ============================================================
# SECTION 5: Merge context into final drug table
# ============================================================

curated_all <- final_drugs %>%
  left_join(dreimt_context, by = "drug_key") %>%
  mutate(
    dreimt_context_class = ifelse(
      is.na(dreimt_context_class),
      "no_dreimt_row",
      dreimt_context_class
    ),

    dreimt_support_level = ifelse(
      is.na(dreimt_support_level),
      "no_dreimt_row",
      dreimt_support_level
    ),

    # Recompute transcriptomic support using literature-aligned DREIMT tau < -80.
    n_transcriptomic_tools_80 =
      as.integer(clue_strong_support) +
      as.integer(ilincs_strong_support) +
      as.integer(dreimt_has_reversal_80 %in% TRUE),

    n_transcriptomic_tools_80_clean =
      as.integer(clue_strong_support) +
      as.integer(ilincs_strong_support) +
      as.integer(dreimt_context_class == "clean_reversal"),

    dreimt_penalty_note = case_when(
      dreimt_context_class == "mixed_reversal_and_mimic" ~
        "DREIMT context-dependent: both reversal and mimic rows",
      dreimt_context_class == "discordant_mimic_only" ~
        "DREIMT-discordant: mimic only, no reversal row",
      dreimt_context_class == "clean_reversal" ~
        "DREIMT clean reversal",
      dreimt_context_class == "no_strong_dreimt_signal" ~
        "No strong DREIMT reversal or mimic signal",
      dreimt_context_class == "no_dreimt_row" ~
        "No DREIMT row found",
      TRUE ~ "Manual review"
    ),

    final_curation_category = case_when(
      high_caution ~
        "Exclude/high-caution despite computational support",

      n_transcriptomic_tools_80 >= 3 &
        dreimt_context_class == "clean_reversal" ~
        "Primary: all 3 tools with clean DREIMT reversal",

      n_transcriptomic_tools_80 >= 3 &
        dreimt_context_class == "mixed_reversal_and_mimic" ~
        "Context-dependent: all 3 tools but DREIMT mixed",

      n_transcriptomic_tools_80 == 2 &
        dreimt_context_class == "clean_reversal" &
        dreimt_has_fdr_reversal_80 ~
        "Primary: two tools with clean FDR-supported DREIMT reversal",

      n_transcriptomic_tools_80 == 2 &
        dreimt_context_class == "clean_reversal" &
        !dreimt_has_fdr_reversal_80 ~
        "Primary/secondary: two tools with clean tau-only DREIMT reversal",

      n_transcriptomic_tools_80 == 2 &
        dreimt_context_class %in% c("no_dreimt_row", "no_strong_dreimt_signal") ~
        "Secondary: two tools without DREIMT reversal",

      n_transcriptomic_tools_80 == 2 &
        dreimt_context_class == "discordant_mimic_only" ~
        "Secondary/DREIMT-discordant: two tools but mimic-only DREIMT",

      n_transcriptomic_tools_80 == 1 &
        dgidb_present ~
        "Exploratory: one transcriptomic tool plus DGIdb",

      n_transcriptomic_tools_80 == 1 ~
        "Low/secondary: one transcriptomic tool only",

      dgidb_present ~
        "Target-only annotation",

      TRUE ~
        "Background/low priority"
    ),

    pair_pool_category = case_when(
      final_curation_category %in% c(
        "Primary: all 3 tools with clean DREIMT reversal",
        "Primary: two tools with clean FDR-supported DREIMT reversal",
        "Primary/secondary: two tools with clean tau-only DREIMT reversal",
        "Secondary: two tools without DREIMT reversal"
      ) ~ "primary_pair_pool",

      final_curation_category %in% c(
        "Context-dependent: all 3 tools but DREIMT mixed",
        "Secondary/DREIMT-discordant: two tools but mimic-only DREIMT"
      ) ~ "context_dependent_or_secondary_pair_pool",

      final_curation_category == "Exploratory: one transcriptomic tool plus DGIdb" ~
        "exploratory_pair_pool",

      TRUE ~ "not_for_primary_pair_design"
    )
  ) %>%
  arrange(
    factor(
      pair_pool_category,
      levels = c(
        "primary_pair_pool",
        "context_dependent_or_secondary_pair_pool",
        "exploratory_pair_pool",
        "not_for_primary_pair_design"
      )
    ),
    desc(n_transcriptomic_tools_80),
    desc(dreimt_has_fdr_reversal_80),
    desc(dreimt_any_approved),
    desc(final_priority_score),
    drug_clean
  )


# ============================================================
# SECTION 6: Create final candidate pools
# ============================================================

primary_pool <- curated_all %>%
  filter(pair_pool_category == "primary_pair_pool") %>%
  filter(!str_detect(drug_clean, "^brd-|^mls|^chemb|^ac1|^spectrum|^prestwick|^san"))

context_dependent_pool <- curated_all %>%
  filter(pair_pool_category == "context_dependent_or_secondary_pair_pool") %>%
  filter(!str_detect(drug_clean, "^brd-|^mls|^chemb|^ac1|^spectrum|^prestwick|^san"))

exploratory_pool <- curated_all %>%
  filter(pair_pool_category == "exploratory_pair_pool") %>%
  filter(!str_detect(drug_clean, "^brd-|^mls|^chemb|^ac1|^spectrum|^prestwick|^san"))

# Current working drugs that we have been discussing.
working_drugs <- clean_drug_name(c(
  "atorvastatin",
  "cilostazol",
  "metformin",
  "memantine",
  "prednisone",
  "vx-745",
  "tipifarnib",
  "alda-1",
  "colforsin",
  "valproic acid",
  "cediranib",
  "triptolide / pg 490"
))

working_review <- curated_all %>%
  filter(drug_clean %in% working_drugs) %>%
  select(
    drug_clean,
    drug_clean_original,
    n_transcriptomic_tools_80,
    n_transcriptomic_tools_80_clean,
    clue_strong_support,
    clue_median_norm_cs,
    clue_n_negative,
    ilincs_strong_support,
    ilincs_best_zscore,
    ilincs_best_pvalue,
    dreimt_context_class,
    dreimt_support_level,
    dreimt_best_reversal_tau_80,
    dreimt_best_reversal_up_dr_80,
    dreimt_best_reversal_down_fdr_80,
    dreimt_strongest_mimic_tau_80,
    dgidb_present,
    dreimt_any_approved,
    high_caution,
    drug_class,
    translation_caution,
    dreimt_penalty_note,
    final_curation_category,
    pair_pool_category,
    final_priority_score
  ) %>%
  arrange(
    factor(
      pair_pool_category,
      levels = c(
        "primary_pair_pool",
        "context_dependent_or_secondary_pair_pool",
        "exploratory_pair_pool",
        "not_for_primary_pair_design"
      )
    ),
    desc(n_transcriptomic_tools_80),
    drug_clean
  )


# ============================================================
# SECTION 7: Summaries and output files
# ============================================================

context_summary <- curated_all %>%
  group_by(dreimt_context_class, dreimt_support_level) %>%
  summarise(
    n_drugs = n(),
    n_primary_pair_pool = sum(pair_pool_category == "primary_pair_pool", na.rm = TRUE),
    n_context_dependent = sum(pair_pool_category == "context_dependent_or_secondary_pair_pool", na.rm = TRUE),
    n_exploratory = sum(pair_pool_category == "exploratory_pair_pool", na.rm = TRUE),
    examples = paste(head(drug_clean, 12), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_drugs))

curation_summary <- curated_all %>%
  group_by(final_curation_category, pair_pool_category) %>%
  summarise(
    n_drugs = n(),
    n_clue = sum(clue_strong_support, na.rm = TRUE),
    n_ilincs = sum(ilincs_strong_support, na.rm = TRUE),
    n_dreimt_reversal_80 = sum(dreimt_has_reversal_80 %in% TRUE, na.rm = TRUE),
    n_dreimt_clean = sum(dreimt_context_class == "clean_reversal", na.rm = TRUE),
    n_dreimt_mixed = sum(dreimt_context_class == "mixed_reversal_and_mimic", na.rm = TRUE),
    n_dreimt_mimic_only = sum(dreimt_context_class == "discordant_mimic_only", na.rm = TRUE),
    n_dgidb = sum(dgidb_present, na.rm = TRUE),
    examples = paste(head(drug_clean, 12), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(pair_pool_category, desc(n_drugs))

full_path <- file.path(
  DRUG_DIR,
  "final_single_drug_DREIMT_CONTEXT_CLASSIFIED_AD61026.csv"
)

primary_path <- file.path(
  DRUG_DIR,
  "final_single_drug_pair_pool_PRIMARY_AD61026.csv"
)

context_path <- file.path(
  DRUG_DIR,
  "final_single_drug_pair_pool_CONTEXT_DEPENDENT_AD61026.csv"
)

exploratory_path <- file.path(
  DRUG_DIR,
  "final_single_drug_pair_pool_EXPLORATORY_AD61026.csv"
)

working_review_path <- file.path(
  DRUG_DIR,
  "final_single_drug_WORKING_REVIEW_AD61026.csv"
)

context_summary_path <- file.path(
  DRUG_DIR,
  "final_single_drug_DREIMT_context_summary_AD61026.csv"
)

curation_summary_path <- file.path(
  DRUG_DIR,
  "final_single_drug_FINAL_CURATION_SUMMARY_AD61026.csv"
)

write.csv(curated_all, full_path, row.names = FALSE)
write.csv(primary_pool, primary_path, row.names = FALSE)
write.csv(context_dependent_pool, context_path, row.names = FALSE)
write.csv(exploratory_pool, exploratory_path, row.names = FALSE)
write.csv(working_review, working_review_path, row.names = FALSE)
write.csv(context_summary, context_summary_path, row.names = FALSE)
write.csv(curation_summary, curation_summary_path, row.names = FALSE)


# ============================================================
# SECTION 8: Print summary
# ============================================================

cat("DREIMT context classification complete.\n\n")

cat("DREIMT context summary:\n")
print(context_summary)

cat("\nFinal curation summary:\n")
print(curation_summary)

cat("\nWorking drug review:\n")
print(working_review)

cat("\nOutput files:\n")
cat("  ", full_path, "\n")
cat("  ", primary_path, "\n")
cat("  ", context_path, "\n")
cat("  ", exploratory_path, "\n")
cat("  ", working_review_path, "\n")
cat("  ", context_summary_path, "\n")
cat("  ", curation_summary_path, "\n")

cat("\nScientific interpretation:\n")
cat("  Use the primary pair pool for main drug-pair prioritization.\n")
cat("  Use context-dependent drugs only if mechanistically justified.\n")
cat("  Use exploratory drugs only in a separate exploratory section or supplement.\n")
cat("============================================================\n")
