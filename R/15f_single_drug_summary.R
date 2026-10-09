# ============================================================
# 15f_single_drug_summary.R
# ============================================================
# Purpose:
#   Build the final curated single-drug candidate table used for Table 4
#   and Fig 6: candidate category (primary, context-dependent, exploratory),
#   transcriptomic tool support, DREIMT quality, DGIdb annotation, caution
#   class, and biological program annotations.
#
# This code was taken from an earlier drug-pair analysis script that is
# not part of the published study. Only the single-drug steps are kept;
# no drug pairs are generated. Some column names (manual_pool,
# pair_use_note) keep their original names because scripts 16 and 17a
# read them.
#
# Input:  Output/Drug_Prioritization/final_single_drug_DREIMT_CONTEXT_CLASSIFIED_AD61026.csv (from 15e)
# Output: Output/Drug_Prioritization/final_single_drug_STRICT_PAIR_INPUT_SUMMARY_AD61026.csv
#         (file name kept for compatibility with scripts 16 and 17a)
# ============================================================

# ------------------------------------------------------------
# Project folder. Set the AD_PROJECT_ROOT environment variable to the
# folder that contains Data/, RDS/, Output/ and Supplementary_Tables/,
# or start R with that folder as the working directory.
# ------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("AD_PROJECT_ROOT", unset = getwd())

# ============================================================
# SECTION 1: Load setup and packages
# ============================================================

setup_path <- file.path(PROJECT_ROOT, "RDS/project_setup.rds")

if (!file.exists(setup_path)) {
  stop("project_setup.rds not found. Run 01_setup.R first.")
}

setup <- readRDS(setup_path)
list2env(setup, envir = .GlobalEnv)

required_pkgs <- c("dplyr", "tidyr", "stringr", "tibble", "ggplot2")

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})

if (!exists("OUT_DIR")) {
  OUT_DIR <- file.path(PROJECT_ROOT, "Output")
}

if (!exists("FIG_DIR")) {
  FIG_DIR <- file.path(OUT_DIR, "Figures")
}

DRUG_DIR <- file.path(OUT_DIR, "Drug_Prioritization")
dir.create(DRUG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(
  DRUG_DIR,
  "final_single_drug_DREIMT_CONTEXT_CLASSIFIED_AD61026.csv"
)

if (!file.exists(input_path)) {
  stop("Missing input file: ", input_path)
}

cat("============================================================\n")
cat("  15f: Final single-drug candidate summary\n")
cat("============================================================\n\n")
cat("Input file:\n")
cat("  ", input_path, "\n\n")


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



# ============================================================
# SECTION 3: Load finalized single-drug table
# ============================================================

single_drugs <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)

logical_cols <- c(
  "clue_strong_support",
  "ilincs_strong_support",
  "dgidb_present",
  "high_caution",
  "moderate_caution",
  "dreimt_any_approved",
  "dreimt_has_reversal_80",
  "dreimt_has_mimic_80",
  "dreimt_has_fdr_reversal_80"
)

for (cc in logical_cols) {
  if (cc %in% colnames(single_drugs)) {
    single_drugs[[cc]] <- as_logical_safe(single_drugs[[cc]])
  }
}

single_drugs <- single_drugs %>%
  mutate(drug_clean = clean_drug_name(drug_clean))


# ============================================================
# SECTION 4: Freeze manually curated single-drug pools
# ============================================================

# Primary candidates.
primary_pair_pool <- clean_drug_name(c(
  "cilostazol",
  "vx-745",
  "metformin",
  "memantine",
  "prednisone",
  "valproic acid",
  "alda-1",
  "tipifarnib"
))

# Context-dependent drug: strong multi-tool signal but mixed DREIMT reversal/mimic.
context_dependent_pool <- clean_drug_name(c(
  "atorvastatin"
))

# Exploratory drug: target-supported or one-tool-plus-target evidence.
exploratory_pool <- clean_drug_name(c(
  "cediranib"
))

# High-caution example retained for tracking but excluded from the candidate list.
excluded_from_primary_pair_design <- clean_drug_name(c(
  "triptolide / pg 490"
))

allowed_drugs <- c(
  primary_pair_pool,
  context_dependent_pool,
  exploratory_pool
)

manual_pool_table <- tibble(
  drug_clean = allowed_drugs,
  manual_pool = case_when(
    drug_clean %in% primary_pair_pool ~ "primary",
    drug_clean %in% context_dependent_pool ~ "context_dependent",
    drug_clean %in% exploratory_pool ~ "exploratory",
    TRUE ~ "not_selected"
  )
)


# ============================================================
# SECTION 5: Manual mechanism and program annotations
# ============================================================

# Program labels are categorical biology annotations. They are not numerical weights.

drug_annotation <- tribble(
  ~drug_clean,       ~manual_mechanism_label,                           ~programs_targeted,                                           ~manual_caution_class,        ~pair_use_note,
  "cilostazol",      "PDE3/cAMP signaling; vascular-inflammatory modulation", "inflammatory_signaling; receptor_cAMP_vascular",         "moderate",                  "Clean DREIMT FDR-supported reversal; approved drug but vascular/platelet biology requires safety review.",
  "vx-745",          "p38 MAPK pathway inhibition",                    "inflammatory_signaling; stress_MAPK_signaling",               "moderate_experimental",     "Clean DREIMT FDR-supported reversal; experimental p38 inhibitor.",
  "metformin",       "AMPK/metabolic and mitochondrial stress modulation", "metabolic_mitochondrial; inflammatory_signaling",            "low_to_moderate",           "Clean DREIMT tau-only reversal; approved drug; metabolic mechanism is plausible.",
  "memantine",       "NMDA/glutamatergic signaling modulation",        "receptor_neuroimmune; stress_excitotoxicity",                  "low_to_moderate",           "Clean DREIMT tau-only reversal; neurologically relevant approved AD drug.",
  "prednisone",      "Glucocorticoid receptor-mediated immunosuppression", "inflammatory_signaling; stress_response",                    "moderate_broad",            "Clean DREIMT tau-only reversal; broad steroid effect may dominate interpretation.",
  "valproic acid",   "HDAC/GABA/sodium-channel related broad modulation", "stress_proteostasis_epigenetic; receptor_neuroimmune",        "moderate_broad",            "Clean DREIMT tau-only reversal; broad HDAC/GABA/sodium-channel biology.",
  "alda-1",          "ALDH2 activation / oxidative aldehyde detoxification", "metabolic_mitochondrial; oxidative_stress",                  "moderate_experimental",     "CLUE+iLINCS support; experimental mitochondrial/oxidative-stress candidate.",
  "tipifarnib",      "Farnesyltransferase/RAS pathway modulation",      "receptor_RAS_signaling; inflammatory_signaling",               "moderate_experimental",     "CLUE+iLINCS support; no strong DREIMT reversal.",
  "atorvastatin",    "HMGCR/statin lipid-inflammatory modulation",      "lipid_metabolism; inflammatory_signaling",                     "context_dependent",         "CLUE+iLINCS plus mixed DREIMT reversal/mimic; keep as context-dependent.",
  "cediranib",       "VEGF/RTK/FLT1 pathway inhibition",               "receptor_RTK_vascular; inflammatory_signaling",                "exploratory",               "CLUE+DGIdb only; target-supported exploratory candidate."
) %>%
  mutate(drug_clean = clean_drug_name(drug_clean))


# ============================================================
# SECTION 6: Build final single-drug table
# ============================================================

pair_drugs <- single_drugs %>%
  filter(drug_clean %in% allowed_drugs) %>%
  left_join(manual_pool_table, by = "drug_clean") %>%
  left_join(drug_annotation, by = "drug_clean") %>%
  mutate(
    manual_pool = ifelse(is.na(manual_pool), "not_selected", manual_pool),

    # Transcriptomic tool count:
    #   CLUE + iLINCS + DREIMT reversal.
    #   DGIdb is NOT included because it is target annotation only.
    transcriptomic_tool_count = n_transcriptomic_tools_80,

    has_clean_dreimt_reversal =
      dreimt_support_level %in% c(
        "clean_FDR_supported_reversal",
        "clean_tau_only_reversal"
      ),

    has_clean_fdr_dreimt_reversal =
      dreimt_support_level == "clean_FDR_supported_reversal",

    has_clean_tau_only_dreimt_reversal =
      dreimt_support_level == "clean_tau_only_reversal",

    has_mixed_dreimt =
      dreimt_support_level %in% c(
        "mixed_FDR_supported_reversal",
        "mixed_tau_only_reversal"
      ),

    has_no_strong_dreimt =
      dreimt_support_level %in% c(
        "no_dreimt_row",
        "no_strong_dreimt_signal"
      ),

    dreimt_quality_label = case_when(
      has_clean_fdr_dreimt_reversal ~ "clean FDR-supported DREIMT reversal",
      has_clean_tau_only_dreimt_reversal ~ "clean tau-only DREIMT reversal",
      has_mixed_dreimt ~ "mixed DREIMT reversal/mimic",
      dreimt_support_level == "mimic_only_not_reversal" ~ "DREIMT mimic only",
      TRUE ~ "no strong DREIMT reversal"
    )
  )

missing_allowed <- setdiff(allowed_drugs, pair_drugs$drug_clean)

if (length(missing_allowed) > 0) {
  warning("These allowed drugs were not found in the finalized single-drug table: ",
          paste(missing_allowed, collapse = ", "))
}



# ============================================================
# SECTION 7: Save the single-drug summary
# ============================================================

single_drug_summary <- pair_drugs %>%
  select(
    drug_clean,
    manual_pool,
    transcriptomic_tool_count,
    clue_strong_support,
    ilincs_strong_support,
    dreimt_support_level,
    dreimt_quality_label,
    dgidb_present,
    manual_caution_class,
    programs_targeted,
    pair_use_note
  ) %>%
  arrange(manual_pool, desc(transcriptomic_tool_count), drug_clean)


single_drug_summary_path <- file.path(
  DRUG_DIR,
  "final_single_drug_STRICT_PAIR_INPUT_SUMMARY_AD61026.csv"
)
write.csv(single_drug_summary, single_drug_summary_path, row.names = FALSE)
cat("Saved:", single_drug_summary_path, "\n")
