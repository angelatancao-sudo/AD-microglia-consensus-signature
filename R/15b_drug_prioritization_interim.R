# ============================================================
# 15b_drug_prioritization_interim.R
# ============================================================
# Project: AD Microglia Multi-Dataset Transcriptomic Analysis
#
# Purpose:
#   Build an interim drug-prioritization table before DREIMT arrives,
#   using already available evidence:
#      1. CLUE/CMap negative perturbagen summary
#      2. iLINCS negative connected perturbations
#      3. DGIdb target-level support
#
# Why this fixed version is needed:
#   The earlier 06C script continued even when the CLUE/iLINCS processed
#   files were missing. This fixed version stops early with a clear message
#   if those files are not in the expected folder.
#
# Required local files:
#   Put these two processed files in:
#     <project folder>/Output/Drug_Prioritization/processed_drug_results/
#
#   Required:
#     CLUE_CMap_negative_perturbagen_summary.csv
#     iLINCS_negative_connected_perturbations.csv
#
#   DGIdb file should already be in:
#     <project folder>/Output/Drug_Prioritization/
#
#   Required:
#     dgidb_combined_drug_summary_AD61026.csv
#
# Recommended run:
#   source("R/06C_interim_drug_prioritization_before_DREIMT_AD61026_FIXED.R")
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

dir.create(DRUG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("  06C FIXED: Interim drug prioritization before DREIMT\n")
cat("============================================================\n\n")

cat("Drug directory:\n")
cat("  ", DRUG_DIR, "\n")
cat("Processed result directory:\n")
cat("  ", PROCESSED_DIR, "\n\n")

# ============================================================
# SECTION 2: Confirm required files exist
# ============================================================

clue_path <- file.path(PROCESSED_DIR, "CLUE_CMap_negative_perturbagen_summary.csv")
ilincs_path <- file.path(PROCESSED_DIR, "iLINCS_negative_connected_perturbations.csv")
dgidb_summary_path <- file.path(DRUG_DIR, "dgidb_combined_drug_summary_AD61026.csv")

required_files <- c(clue_path, ilincs_path, dgidb_summary_path)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "\nMissing required file(s):\n  ",
    paste(missing_files, collapse = "\n  "),
    "\n\nPlease copy the processed CLUE and iLINCS CSV files into:\n  ",
    PROCESSED_DIR,
    "\n\nThen rerun this script.\n"
  )
}

# ============================================================
# SECTION 3: Helper functions
# ============================================================

clean_drug_name <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace_all(x, "\\s+", " ")
  x_lower <- str_to_lower(x)

  case_when(
    x_lower %in% c("pg 490", "triptolide") ~ "triptolide / pg 490",
    x_lower %in% c("trichostatin a, streptomyces sp.", "trichostatin a") ~ "trichostatin a",
    x_lower %in% c("cucurbitacin-i", "cucurbitacin i") ~ "cucurbitacin i",
    x_lower %in% c("alvocidib", "flavopiridol") ~ "alvocidib",
    TRUE ~ x_lower
  )
}

classify_drug_class <- function(drug, targets) {
  combined <- str_to_lower(paste(drug, targets, sep = " "))

  case_when(
    str_detect(combined, "cdk|cyclin") ~ "CDK/cell-cycle inhibitor",
    str_detect(combined, "hdac") ~ "HDAC/epigenetic inhibitor",
    str_detect(combined, "hsp90") ~ "HSP90 inhibitor",
    str_detect(combined, "top2|top1|epirubicin|doxorubicin|teniposide") ~ "DNA/topoisomerase/cytotoxic",
    str_detect(combined, "proteasome|psmb|psma|mg-132|bortezomib") ~ "Proteasome inhibitor",
    str_detect(combined, "map2k|mek|trametinib|selumetinib|pd-0325901") ~ "MEK/MAPK pathway inhibitor",
    str_detect(combined, "pi3k|akt|mtor|pik3") ~ "PI3K/AKT/mTOR pathway",
    str_detect(combined, "flt1|kdr|vegf|sunitinib|axitinib|cediranib|sorafenib") ~ "VEGF/RTK inhibitor",
    str_detect(combined, "stat3|jak") ~ "JAK/STAT pathway",
    str_detect(combined, "hmgcr|statin") ~ "Statin/lipid metabolism",
    str_detect(combined, "fnta|fntb|tipifarnib") ~ "Farnesyltransferase/RAS pathway",
    str_detect(combined, "aldh2|alda") ~ "Mitochondrial/ALDH2 related",
    str_detect(combined, "glucocorticoid|nr3c1|dexamethasone|beclomethasone|ciclesonide") ~ "Glucocorticoid receptor related",
    TRUE ~ "Other/unclear"
  )
}

flag_translation_caution <- function(drug_class, drug, targets) {
  combined <- str_to_lower(paste(drug_class, drug, targets, sep = " "))

  case_when(
    str_detect(combined, "cdk|cell-cycle|topoisomerase|cytotoxic|proteasome|hsp90|dna") ~
      "High caution: cytotoxic/oncology-like mechanism",
    str_detect(combined, "hdac|pi3k|akt|mtor|vegf|rtk") ~
      "Moderate caution: broad signaling/oncology mechanism",
    str_detect(combined, "statin|farnesyltransferase|aldh2|glucocorticoid") ~
      "More interpretable but still requires AD/microglia validation",
    TRUE ~ "Needs manual review"
  )
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

# ============================================================
# SECTION 4: Load input files
# ============================================================

cat("Loading input files...\n\n")

clue <- read.csv(clue_path, stringsAsFactors = FALSE, check.names = FALSE)
ilincs <- read.csv(ilincs_path, stringsAsFactors = FALSE, check.names = FALSE)
dgidb_summary <- read.csv(dgidb_summary_path, stringsAsFactors = FALSE, check.names = FALSE)

cat("Rows loaded:\n")
cat("  CLUE negative perturbagen summary:", nrow(clue), "\n")
cat("  iLINCS negative connected perturbations:", nrow(ilincs), "\n")
cat("  DGIdb drug summary:", nrow(dgidb_summary), "\n\n")

# ============================================================
# SECTION 5: Standardize CLUE results
# ============================================================

cat("Standardizing CLUE results...\n")

clue_cols <- colnames(clue)

drug_col <- clue_cols[str_to_lower(clue_cols) %in%
                        c("perturbagen", "pert_iname", "drug", "compound", "name", "drug_clean")]

if (length(drug_col) == 0) {
  drug_col <- clue_cols[str_detect(str_to_lower(clue_cols), "pert|drug|compound|name")][1]
} else {
  drug_col <- drug_col[1]
}

score_col <- clue_cols[str_to_lower(clue_cols) %in%
                         c("mean_tau", "best_tau", "tau", "score", "connectivity",
                           "normalized_connectivity_score", "clue_best_negative_score")]

if (length(score_col) == 0) {
  score_col <- clue_cols[str_detect(str_to_lower(clue_cols), "tau|score|connect")][1]
} else {
  score_col <- score_col[1]
}

count_col <- clue_cols[str_to_lower(clue_cols) %in%
                         c("n_signatures", "n_sig", "n", "n_rows", "signature_count",
                           "clue_n_signatures")]

if (length(count_col) == 0) {
  count_col <- clue_cols[str_detect(str_to_lower(clue_cols), "count|signature|n_")][1]
} else {
  count_col <- count_col[1]
}

if (is.na(drug_col)) {
  stop("Could not identify drug/perturbagen column in CLUE file.")
}

clue_std <- clue %>%
  mutate(
    drug_raw = .data[[drug_col]],
    drug_clean = clean_drug_name(drug_raw),
    clue_present = TRUE,
    clue_score_raw = if (!is.na(score_col)) suppressWarnings(as.numeric(.data[[score_col]])) else NA_real_,
    clue_n_signatures_raw = if (!is.na(count_col)) suppressWarnings(as.numeric(.data[[count_col]])) else NA_real_
  ) %>%
  group_by(drug_clean) %>%
  summarise(
    clue_present = TRUE,
    clue_best_negative_score = suppressWarnings(min(clue_score_raw, na.rm = TRUE)),
    clue_mean_negative_score = suppressWarnings(mean(clue_score_raw, na.rm = TRUE)),
    clue_n_signatures = suppressWarnings(max(clue_n_signatures_raw, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    clue_best_negative_score = ifelse(is.infinite(clue_best_negative_score), NA_real_, clue_best_negative_score),
    clue_mean_negative_score = ifelse(is.nan(clue_mean_negative_score), NA_real_, clue_mean_negative_score),
    clue_n_signatures = ifelse(is.infinite(clue_n_signatures), NA_real_, clue_n_signatures)
  )

cat("  Standardized CLUE drugs:", nrow(clue_std), "\n\n")

# ============================================================
# SECTION 6: Standardize iLINCS results
# ============================================================

cat("Standardizing iLINCS results...\n")

required_ilincs_cols <- c("Perturbagen", "Correlation", "NoOfSignatures", "pValue", "zScore")

if (!all(required_ilincs_cols %in% colnames(ilincs))) {
  stop("iLINCS file is missing expected columns: ",
       paste(setdiff(required_ilincs_cols, colnames(ilincs)), collapse = ", "))
}

ilincs_std <- ilincs %>%
  filter(Correlation == "-") %>%
  mutate(
    drug_raw = Perturbagen,
    drug_clean = clean_drug_name(drug_raw),
    targets = if ("GeneTargets" %in% colnames(.)) GeneTargets else NA_character_,
    ilincs_present = TRUE,
    ilincs_pvalue = suppressWarnings(as.numeric(pValue)),
    ilincs_zscore = suppressWarnings(as.numeric(zScore)),
    ilincs_n_signatures = suppressWarnings(as.numeric(NoOfSignatures))
  ) %>%
  group_by(drug_clean) %>%
  summarise(
    ilincs_present = TRUE,
    ilincs_best_pvalue = min(ilincs_pvalue, na.rm = TRUE),
    ilincs_best_zscore = max(ilincs_zscore, na.rm = TRUE),
    ilincs_n_signatures = max(ilincs_n_signatures, na.rm = TRUE),
    ilincs_targets = paste(unique(na.omit(targets)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    ilincs_best_pvalue = ifelse(is.infinite(ilincs_best_pvalue), NA_real_, ilincs_best_pvalue),
    ilincs_best_zscore = ifelse(is.infinite(ilincs_best_zscore), NA_real_, ilincs_best_zscore),
    ilincs_n_signatures = ifelse(is.infinite(ilincs_n_signatures), NA_real_, ilincs_n_signatures)
  )

cat("  Standardized iLINCS drugs:", nrow(ilincs_std), "\n\n")

# ============================================================
# SECTION 7: Standardize DGIdb results
# ============================================================

cat("Standardizing DGIdb results...\n")

if (!"drug" %in% colnames(dgidb_summary)) {
  stop("DGIdb summary file does not contain a column named 'drug'.")
}

dgidb_std <- dgidb_summary %>%
  mutate(
    drug_raw = drug,
    drug_clean = clean_drug_name(drug_raw),
    dgidb_present = TRUE,
    dgidb_n_genes = if ("n_genes" %in% colnames(.)) suppressWarnings(as.numeric(n_genes)) else NA_real_,
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
    dgidb_n_genes = max(dgidb_n_genes, na.rm = TRUE),
    dgidb_relevant_interactions = max(dgidb_relevant_interactions, na.rm = TRUE),
    dgidb_target_genes = paste(unique(na.omit(dgidb_target_genes)), collapse = "; "),
    dgidb_interaction_types = paste(unique(na.omit(dgidb_interaction_types)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    dgidb_n_genes = ifelse(is.infinite(dgidb_n_genes), NA_real_, dgidb_n_genes),
    dgidb_relevant_interactions =
      ifelse(is.infinite(dgidb_relevant_interactions), NA_real_, dgidb_relevant_interactions)
  )

cat("  Standardized DGIdb drugs:", nrow(dgidb_std), "\n\n")

# ============================================================
# SECTION 8: Merge evidence and score
# ============================================================

cat("Merging evidence sources...\n\n")

all_drugs <- unique(c(clue_std$drug_clean, ilincs_std$drug_clean, dgidb_std$drug_clean))

interim <- data.frame(drug_clean = all_drugs, stringsAsFactors = FALSE) %>%
  left_join(clue_std, by = "drug_clean") %>%
  left_join(ilincs_std, by = "drug_clean") %>%
  left_join(dgidb_std, by = "drug_clean") %>%
  mutate(
    clue_present = ifelse(is.na(clue_present), FALSE, clue_present),
    ilincs_present = ifelse(is.na(ilincs_present), FALSE, ilincs_present),
    dgidb_present = ifelse(is.na(dgidb_present), FALSE, dgidb_present),
    dreimt_present = NA,
    dreimt_score = NA_real_,
    n_transcriptomic_tools_pre_DREIMT =
      as.integer(clue_present) + as.integer(ilincs_present),
    n_evidence_sources_pre_DREIMT =
      as.integer(clue_present) + as.integer(ilincs_present) + as.integer(dgidb_present),
    display_name = drug_clean,
    combined_targets = paste(
      ifelse(is.na(ilincs_targets), "", ilincs_targets),
      ifelse(is.na(dgidb_target_genes), "", dgidb_target_genes),
      sep = "; "
    ),
    drug_class = classify_drug_class(display_name, combined_targets),
    translation_caution = flag_translation_caution(drug_class, display_name, combined_targets),
    clue_component = ifelse(clue_present, 1, 0),
    ilincs_z_component = standardize_score(ilincs_best_zscore),
    ilincs_nsig_component = standardize_score(log10(ilincs_n_signatures + 1)),
    dgidb_component = standardize_score(dgidb_relevant_interactions),
    cytotoxic_penalty = case_when(
      translation_caution == "High caution: cytotoxic/oncology-like mechanism" ~ 1,
      translation_caution == "Moderate caution: broad signaling/oncology mechanism" ~ 0.5,
      TRUE ~ 0
    ),
    interim_score_pre_DREIMT =
      2.0 * n_transcriptomic_tools_pre_DREIMT +
      0.8 * clue_component +
      1.0 * ifelse(is.na(ilincs_z_component), 0, ilincs_z_component) +
      0.5 * ifelse(is.na(ilincs_nsig_component), 0, ilincs_nsig_component) +
      0.8 * ifelse(is.na(dgidb_component), 0, dgidb_component) -
      0.8 * cytotoxic_penalty
  ) %>%
  arrange(
    desc(n_transcriptomic_tools_pre_DREIMT),
    desc(n_evidence_sources_pre_DREIMT),
    desc(interim_score_pre_DREIMT),
    ilincs_best_pvalue
  )

# ============================================================
# SECTION 9: Save output tables and draft plot
# ============================================================

interim_path <- file.path(DRUG_DIR, "interim_drug_prioritization_pre_DREIMT_AD61026.csv")
top_review_path <- file.path(DRUG_DIR, "interim_drug_prioritization_top_review_pre_DREIMT_AD61026.csv")
class_path <- file.path(DRUG_DIR, "drug_class_review_pre_DREIMT_AD61026.csv")

write.csv(interim, interim_path, row.names = FALSE)

top_review <- interim %>%
  filter(
    n_transcriptomic_tools_pre_DREIMT >= 1,
    !str_detect(display_name, "^brd-|^mls|^chemb|^ac1|^spectrum|^prestwick|^san")
  ) %>%
  slice_head(n = 100)

write.csv(top_review, top_review_path, row.names = FALSE)

drug_class_review <- interim %>%
  group_by(drug_class, translation_caution) %>%
  summarise(
    n_drugs = n(),
    n_clue = sum(clue_present, na.rm = TRUE),
    n_ilincs = sum(ilincs_present, na.rm = TRUE),
    n_dgidb = sum(dgidb_present, na.rm = TRUE),
    examples = paste(head(display_name, 8), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_drugs))

write.csv(drug_class_review, class_path, row.names = FALSE)

plot_df <- top_review %>%
  slice_head(n = 30) %>%
  mutate(display_name = factor(display_name, levels = rev(display_name)))

if (nrow(plot_df) > 0) {
  p <- ggplot(plot_df, aes(x = display_name, y = interim_score_pre_DREIMT)) +
    geom_col() +
    coord_flip() +
    labs(
      title = "Interim drug prioritization before DREIMT",
      subtitle = "Preliminary score from CLUE/CMap, iLINCS, and DGIdb only",
      x = "Candidate",
      y = "Interim score"
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(
    file.path(FIG_DIR, "Fig_interim_drug_evidence_pre_DREIMT_draft.pdf"),
    p,
    width = 8,
    height = 7
  )
}

# ============================================================
# SECTION 10: Completion summary
# ============================================================

cat("============================================================\n")
cat("06C FIXED complete.\n\n")

cat("Interim table rows:", nrow(interim), "\n")
cat("Top-review table rows:", nrow(top_review), "\n\n")

cat("Top 15 interim candidates before DREIMT:\n")
print(
  top_review %>%
    select(
      display_name,
      n_transcriptomic_tools_pre_DREIMT,
      n_evidence_sources_pre_DREIMT,
      drug_class,
      translation_caution,
      ilincs_best_pvalue,
      ilincs_best_zscore,
      clue_best_negative_score,
      interim_score_pre_DREIMT
    ) %>%
    slice_head(n = 15)
)

cat("\nOutput files:\n")
cat("  ", interim_path, "\n")
cat("  ", top_review_path, "\n")
cat("  ", class_path, "\n")

cat("\nNext step after DREIMT arrives:\n")
cat("  Add DREIMT results and rerun final scoring in 06D.\n")
cat("============================================================\n")
