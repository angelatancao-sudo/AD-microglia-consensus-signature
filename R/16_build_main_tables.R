# ============================================================
# 16_build_main_tables.R
# ============================================================
# Project: AD61026 AD microglia multi-dataset transcriptomic analysis
#
# Purpose:
#   Build final manuscript tables from the cleaned display index.
#
# Main principle:
#   Tables should be generated from frozen source-data files in:
#     Output/Manuscript_Display_Index_AD61026
#
# Why this script is separate from the source-record script:
#   Script 09/09B organized and renamed source files.
#   This script creates the actual manuscript-ready table outputs.
#
# Output:
#   Output/Manuscript_Final_Tables_AD61026/
#     Main_Table_1__cohort_characteristics.csv
#     Main_Table_2__strict_48_consensus_signature.csv
#     Main_Table_3__pathway_summary.csv
#     Main_Table_4__final_single_drug_candidates.csv
#     Main_Tables_AD61026.xlsx
#     Supplementary_Table_S1_*.xlsx
#     ...
#     Supplementary_Table_S10_*.xlsx
#     Table_QC_summary_AD61026.csv
#
# Important scientific notes:
#   - Strict 48 genes = primary biological consensus signature.
#   - Expanded 117 non-mitochondrial genes = perturbational drug-query signature.
#   - DGIdb is target-level annotation only, not transcriptomic reversal.
#   - Pair-level supporting analyses are not validation.
#
# Recommended run:
#   source("R/10_build_manuscript_tables_AD61026.R")
# ============================================================


# ============================================================
# SECTION 1: Setup
# ============================================================

# ------------------------------------------------------------
# Project folder. Set the AD_PROJECT_ROOT environment variable to the
# folder that contains Data/, RDS/, Output/ and Supplementary_Tables/,
# or start R with that folder as the working directory.
# ------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("AD_PROJECT_ROOT", unset = getwd())

# Sections 8 and 10 build an earlier ten-table supplement (S1-S10) that
# included drug-pair analyses not in the published paper. The published
# S1-S5 Tables were assembled from the analysis outputs. Skipped by default.
RUN_LEGACY_SUPPLEMENT <- FALSE

setup_path <- file.path(PROJECT_ROOT, "RDS/project_setup.rds")

if (file.exists(setup_path)) {
  setup <- readRDS(setup_path)
  list2env(setup, envir = .GlobalEnv)
}

if (!exists("BASE_DIR")) {
  BASE_DIR <- PROJECT_ROOT
}
if (!exists("OUT_DIR")) {
  OUT_DIR <- file.path(BASE_DIR, "Output")
}

DISPLAY_DIR <- file.path(OUT_DIR, "Manuscript_Display_Index_AD61026")
DISPLAY_MANIFEST_DIR <- file.path(DISPLAY_DIR, "00_MANIFESTS")
DISPLAY_MAIN_DIR <- file.path(DISPLAY_DIR, "Main_Tables")
DISPLAY_FIG_DIR <- file.path(DISPLAY_DIR, "Figures")
DISPLAY_SUPP_DIR <- file.path(DISPLAY_DIR, "Supplementary_Source_Archive")

FINAL_TABLE_DIR <- file.path(OUT_DIR, "Manuscript_Final_Tables_AD61026")
FINAL_MAIN_DIR <- file.path(FINAL_TABLE_DIR, "Main_Tables")
FINAL_SUPP_DIR <- file.path(FINAL_TABLE_DIR, "Supplementary_Tables")
FINAL_QC_DIR <- file.path(FINAL_TABLE_DIR, "QC")

dir.create(FINAL_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FINAL_MAIN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FINAL_SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FINAL_QC_DIR, recursive = TRUE, showWarnings = FALSE)

required_pkgs <- c("dplyr", "stringr", "tibble", "readr", "tidyr", "purrr")

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "),
       "\nInstall them before running this script.")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(readr)
  library(tidyr)
  library(purrr)
})

has_writexl <- requireNamespace("writexl", quietly = TRUE)

cat("============================================================\n")
cat("  10: Build manuscript tables\n")
cat("============================================================\n\n")

cat("Using display index folder:\n")
cat("  ", DISPLAY_DIR, "\n\n")

cat("Final table output folder:\n")
cat("  ", FINAL_TABLE_DIR, "\n\n")


# ============================================================
# SECTION 2: Helper functions
# ============================================================

find_file_anywhere <- function(file_name) {
  # Search clean display folder first, then broader project output.
  roots <- c(DISPLAY_DIR, OUT_DIR, BASE_DIR)
  roots <- roots[file.exists(roots)]

  hits <- unlist(lapply(roots, function(root) {
    list.files(
      path = root,
      pattern = paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", file_name), "$"),
      recursive = TRUE,
      full.names = TRUE
    )
  }))

  hits <- unique(hits)

  if (length(hits) == 0) {
    return(NA_character_)
  }

  hits[order(file.info(hits)$mtime, decreasing = TRUE)][1]
}

read_csv_safely <- function(file_name, required = TRUE) {
  path <- find_file_anywhere(file_name)

  if (is.na(path) || !file.exists(path)) {
    msg <- paste0("Could not find source file: ", file_name)
    if (required) stop(msg) else warning(msg)
    return(tibble())
  }

  readr::read_csv(path, show_col_types = FALSE, guess_max = 100000) %>%
    mutate(.source_file = basename(path))
}

standardize_bool <- function(x) {
  # Converts TRUE/FALSE or true/false to Yes/No text for manuscript tables.
  case_when(
    is.na(x) ~ "",
    x %in% TRUE ~ "Yes",
    x %in% FALSE ~ "No",
    str_to_lower(as.character(x)) %in% c("true", "yes", "1") ~ "Yes",
    str_to_lower(as.character(x)) %in% c("false", "no", "0") ~ "No",
    TRUE ~ as.character(x)
  )
}

format_num <- function(x, digits = 3) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "", as.character(round(x, digits)))
}

format_p <- function(x, digits = 2) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(
    is.na(x),
    "",
    ifelse(
      x < 0.001,
      formatC(x, format = "e", digits = digits),
      as.character(signif(x, digits + 1))
    )
  )
}

safe_sheet_name <- function(x) {
  # Excel sheet names must be <= 31 characters and avoid special characters.
  x <- str_replace_all(x, "[:\\\\/?*\\[\\]]", "_")
  x <- str_replace_all(x, "\\s+", "_")
  x <- substr(x, 1, 31)
  x
}

make_excel_safe <- function(df, max_chars = 32000) {
  # Excel cannot store >32,767 characters in one cell.
  # CSV files retain full text, but xlsx sheets are truncated if needed.
  df %>%
    mutate(
      across(
        where(is.character),
        ~ ifelse(
          !is.na(.x) & nchar(.x) > max_chars,
          paste0(substr(.x, 1, max_chars), " ... [truncated for Excel export]"),
          .x
        )
      )
    )
}

write_table_outputs <- function(df, csv_path) {
  readr::write_csv(df, csv_path)
  cat("  Saved:", csv_path, "\n")
}

qc_rows <- list()

add_qc <- function(table_id, check_name, value, expected = NA_character_, status = "INFO", note = "") {
  qc_rows[[length(qc_rows) + 1]] <<- tibble(
    table_id = table_id,
    check_name = check_name,
    value = as.character(value),
    expected = as.character(expected),
    status = status,
    note = note
  )
}


# ============================================================
# SECTION 3: Main Table 1 — Cohort characteristics
# ============================================================
# Purpose:
#   Provide dataset-level sample and pseudobulk summary.
#
# Source:
#   Main_Table_1__cohort_characteristics_source.csv
#   or MainTable1_cohort_summary_source.csv
#
# QC:
#   Make sure this table is dataset-level and has no duplicated cohort rows.

cat("Building Main Table 1: Cohort characteristics...\n")

table1_source <- read_csv_safely(
  "Main_Table_1__cohort_characteristics_source.csv",
  required = FALSE
)

if (nrow(table1_source) == 0) {
  table1_source <- read_csv_safely("MainTable1_cohort_summary_source.csv", required = TRUE)
}

table1 <- table1_source %>%
  select(-any_of(".source_file")) %>%
  distinct()

# Prefer a clean ordering if a dataset/cohort column exists.
dataset_col <- intersect(
  c("Dataset", "dataset", "Cohort", "cohort", "GEO", "study"),
  colnames(table1)
)[1]

if (!is.na(dataset_col)) {
  table1 <- table1 %>%
    arrange(.data[[dataset_col]])
}

write_table_outputs(
  table1,
  file.path(FINAL_MAIN_DIR, "Main_Table_1__cohort_characteristics.csv")
)

add_qc(
  "Main_Table_1",
  "Number of rows",
  nrow(table1),
  expected = "One row per cohort/dataset",
  status = ifelse(nrow(table1) > 0, "PASS", "FAIL")
)

if (!is.na(dataset_col)) {
  add_qc(
    "Main_Table_1",
    "Dataset/cohort column used for ordering",
    dataset_col,
    expected = "Dataset-level table",
    status = "INFO"
  )
}


# ============================================================
# SECTION 4: Main Table 2 — Strict 48-gene consensus signature
# ============================================================
# Purpose:
#   Present the primary biological consensus signature.
#
# Source:
#   final_strict_consensus_high_or_rankprod_AD61026.csv
#
# Important expected checks:
#   - 48 genes total
#   - 24 AD-up and 24 AD-down
#   - 0 mitochondrial genes beginning with MT-
#
# Display columns:
#   Gene, Direction, SEA-AD logFC, SEA-AD FDR,
#   confidence, same-direction datasets, nominal replication count,
#   mean replication logFC, RankProd pfp, RankProd AveFC,
#   evidence basis.

cat("Building Main Table 2: Strict 48-gene consensus signature...\n")

strict48 <- read_csv_safely(
  "final_strict_consensus_high_or_rankprod_AD61026.csv",
  required = TRUE
)

table2 <- strict48 %>%
  mutate(
    Direction = case_when(
      discovery_direction %in% c("UP", "UP_in_AD", "AD-upregulated") ~ "AD-upregulated",
      discovery_direction %in% c("DOWN", "DOWN_in_AD", "AD-downregulated") ~ "AD-downregulated",
      discovery_logFC > 0 ~ "AD-upregulated",
      discovery_logFC < 0 ~ "AD-downregulated",
      TRUE ~ as.character(discovery_direction)
    ),
    Evidence_basis = case_when(
      high_confidence_vote_counting %in% TRUE &
        rankprod_pfp_lt_0.05_same_direction %in% TRUE ~
        "High-confidence vote-counting + RankProd",
      high_confidence_vote_counting %in% TRUE ~
        "High-confidence vote-counting",
      rankprod_pfp_lt_0.05_same_direction %in% TRUE ~
        "RankProd same-direction support",
      TRUE ~ "Strict consensus"
    )
  ) %>%
  transmute(
    Gene = gene,
    Direction,
    `SEA-AD logFC` = round(discovery_logFC, 3),
    `SEA-AD FDR` = format_p(discovery_adj.P.Val),
    Confidence = confidence,
    `Replication datasets same direction` = n_same_direction,
    `Nominally significant replication datasets` = n_nominal_p_lt_0.05,
    `Mean replication logFC` = round(mean_replication_logFC, 3),
    `RankProd pfp` = format_p(RankProd_pfp),
    `RankProd AveFC` = format_num(RankProd_AveFC, 3),
    Evidence_basis,
    Recommended_use = recommended_use
  ) %>%
  mutate(
    Direction_sort = ifelse(Direction == "AD-upregulated", 1, 0),
    Abs_logFC_sort = abs(suppressWarnings(as.numeric(`SEA-AD logFC`)))
  ) %>%
  arrange(Direction_sort, desc(Abs_logFC_sort), Gene) %>%
  select(-Direction_sort, -Abs_logFC_sort)

write_table_outputs(
  table2,
  file.path(FINAL_MAIN_DIR, "Main_Table_2__strict_48_consensus_signature.csv")
)

n_strict <- nrow(table2)
n_up <- sum(table2$Direction == "AD-upregulated", na.rm = TRUE)
n_down <- sum(table2$Direction == "AD-downregulated", na.rm = TRUE)
n_mt <- sum(str_detect(table2$Gene, "^MT-"), na.rm = TRUE)

add_qc("Main_Table_2", "Strict consensus gene count", n_strict, "48", ifelse(n_strict == 48, "PASS", "FAIL"))
add_qc("Main_Table_2", "AD-upregulated gene count", n_up, "24", ifelse(n_up == 24, "PASS", "FAIL"))
add_qc("Main_Table_2", "AD-downregulated gene count", n_down, "24", ifelse(n_down == 24, "PASS", "FAIL"))
add_qc("Main_Table_2", "Mitochondrial MT- genes", n_mt, "0", ifelse(n_mt == 0, "PASS", "FAIL"))


# ============================================================
# SECTION 5: Main Table 3 — Pathway enrichment summary
# ============================================================
# Purpose:
#   Provide selected representative pathway results.
#
# Primary pathway evidence:
#   SEA-AD genome-wide GSEA by moderated t-statistic.
#
# Supportive pathway evidence:
#   ORA of strict and expanded signatures.
#
# This table uses current GSEA files directly, not manually typed values.

cat("Building Main Table 3: Pathway enrichment summary...\n")

gsea_hallmark <- read_csv_safely("GSEA_SEAAD_Hallmark.csv", required = TRUE)
gsea_kegg <- read_csv_safely("GSEA_SEAAD_KEGG.csv", required = TRUE)
gsea_gobp <- read_csv_safely("GSEA_SEAAD_GO_BP.csv", required = TRUE)

make_gsea_table <- function(df, collection_label, n_keep) {
  df %>%
    filter(!is.na(padj), padj < 0.05) %>%
    arrange(padj, desc(abs(NES))) %>%
    slice_head(n = n_keep) %>%
    transmute(
      Collection = collection_label,
      Pathway = pathway_label,
      Direction = direction,
      NES = round(NES, 3),
      `P value` = format_p(pval),
      FDR = format_p(padj),
      `Gene set size` = size,
      `Leading-edge genes` = leadingEdge
    )
}

table3 <- bind_rows(
  make_gsea_table(gsea_hallmark, "Hallmark", 10),
  make_gsea_table(gsea_kegg, "KEGG Legacy", 6),
  make_gsea_table(gsea_gobp, "GO Biological Process", 8)
)

write_table_outputs(
  table3,
  file.path(FINAL_MAIN_DIR, "Main_Table_3__pathway_summary.csv")
)

add_qc("Main_Table_3", "Number of selected pathway rows", nrow(table3), "24", ifelse(nrow(table3) == 24, "PASS", "CHECK"))
add_qc("Main_Table_3", "Hallmark rows", sum(table3$Collection == "Hallmark"), "10", ifelse(sum(table3$Collection == "Hallmark") == 10, "PASS", "CHECK"))
add_qc("Main_Table_3", "KEGG rows", sum(table3$Collection == "KEGG Legacy"), "6", ifelse(sum(table3$Collection == "KEGG Legacy") == 6, "PASS", "CHECK"))


# ============================================================
# SECTION 6: Main Table 4 — Final prioritized single-drug candidates
# ============================================================
# Purpose:
#   Present the final curated single-drug candidates.
#
# Important correction:
#   The earlier Main Table 4 source index pointed to the curation summary
#   because drug prioritization was still evolving. For manuscript Table 4,
#   the detailed candidate source is:
#     final_single_drug_STRICT_PAIR_INPUT_SUMMARY_AD61026.csv
#
# This table separates transcriptomic reversal tools from DGIdb target annotation.

cat("Building Main Table 4: Final single-drug candidates...\n")

single_drug <- read_csv_safely(
  "final_single_drug_STRICT_PAIR_INPUT_SUMMARY_AD61026.csv",
  required = TRUE
)

table4 <- single_drug %>%
  mutate(
    Pair_pool = case_when(
      manual_pool == "primary" ~ "Primary",
      manual_pool == "context_dependent" ~ "Context-dependent",
      manual_pool == "exploratory" ~ "Exploratory",
      TRUE ~ as.character(manual_pool)
    ),
    CLUE_CMap = standardize_bool(clue_strong_support),
    iLINCS = standardize_bool(ilincs_strong_support),
    DGIdb_target_annotation = standardize_bool(dgidb_present),
    Sort_pool = case_when(
      manual_pool == "primary" ~ 1,
      manual_pool == "context_dependent" ~ 2,
      manual_pool == "exploratory" ~ 3,
      TRUE ~ 4
    )
  ) %>%
  transmute(
    Drug = drug_clean,
    `Category` = Pair_pool,
    `Transcriptomic reversal tools (n)` = transcriptomic_tool_count,
    `CLUE/CMap reversal` = CLUE_CMap,
    `iLINCS reversal` = iLINCS,
    `DREIMT support` = dreimt_quality_label,
    `DGIdb target annotation` = DGIdb_target_annotation,
    `Caution class` = manual_caution_class,
    `Biological programs` = programs_targeted,
    `Use note` = pair_use_note,
    Sort_pool
  ) %>%
  arrange(Sort_pool, desc(`Transcriptomic reversal tools (n)`), Drug) %>%
  select(-Sort_pool)

write_table_outputs(
  table4,
  file.path(FINAL_MAIN_DIR, "Main_Table_4__final_single_drug_candidates.csv")
)

add_qc("Main_Table_4", "Candidate drug count", nrow(table4), "10", ifelse(nrow(table4) == 10, "PASS", "CHECK"))
add_qc("Main_Table_4", "Primary candidate drugs", sum(table4$`Category` == "Primary"), "7 or more", ifelse(sum(table4$`Category` == "Primary") >= 7, "PASS", "CHECK"))
add_qc("Main_Table_4", "DGIdb interpretation", "DGIdb kept as target annotation only", "Not transcriptomic reversal", "PASS")


# ============================================================
# SECTION 7: Write combined main-table workbook
# ============================================================

if (has_writexl) {
  main_xlsx_path <- file.path(FINAL_MAIN_DIR, "Main_Tables_AD61026.xlsx")

  writexl::write_xlsx(
    list(
      Table1_Cohorts = make_excel_safe(table1),
      Table2_Consensus48 = make_excel_safe(table2),
      Table3_Pathways = make_excel_safe(table3),
      Table4_Drugs = make_excel_safe(table4)
    ),
    path = main_xlsx_path
  )

  cat("  Saved main-table workbook:", main_xlsx_path, "\n")
} else {
  warning("Package writexl is not installed. Main-table xlsx workbook was not created.")
}


# ============================================================
# SECTION 8: Supplementary Table builder
# ============================================================
if (RUN_LEGACY_SUPPLEMENT) {
# Purpose:
#   Build one Excel workbook per recommended supplementary table.
#
# Why not one source file per supplementary table?
#   Each supplementary table may contain several sheets. This keeps the
#   manuscript supplement list clean while preserving detailed evidence.
#
# Source list:
#   Output/Manuscript_Display_Index_AD61026/00_MANIFESTS/
#     AD61026_supplementary_recommended_index.csv

cat("Building Supplementary Tables S1-S10...\n")

supp_index_path <- file.path(DISPLAY_MANIFEST_DIR, "AD61026_supplementary_recommended_index.csv")

if (!file.exists(supp_index_path)) {
  stop("Missing supplementary recommended index: ", supp_index_path)
}

supp_index <- readr::read_csv(supp_index_path, show_col_types = FALSE)

read_source_for_supp <- function(file_name) {
  path <- find_file_anywhere(file_name)

  if (is.na(path) || !file.exists(path)) {
    return(tibble(
      missing_source_file = file_name,
      note = "Source file was not found when building supplementary workbook."
    ))
  }

  out <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE, guess_max = 100000),
    error = function(e) {
      tibble(
        source_file = file_name,
        read_error = as.character(e$message)
      )
    }
  )

  out %>%
    mutate(.source_file = basename(path))
}

supp_qc <- list()

for (i in seq_len(nrow(supp_index))) {

  if (supp_index$include_in_manuscript[i] != "Yes") {
    next
  }

  supp_id <- supp_index$supplement_id[i]
  supp_title <- supp_index$recommended_title[i]

  source_files <- unlist(strsplit(supp_index$source_files_to_use[i], ";"))
  source_files <- str_trim(source_files)
  source_files <- source_files[source_files != ""]

  cat("  Building ", supp_id, ": ", supp_title, "\n", sep = "")

  sheets <- list()

  # Add index/readme sheet.
  sheets[["README"]] <- tibble(
    Supplementary_Table = supp_id,
    Title = supp_title,
    Note = supp_index$note[i],
    Source_files = paste(source_files, collapse = "; "),
    Important_wording = case_when(
      supp_id == "Supplementary_Table_S7" ~ "DGIdb is target-level annotation only.",
      supp_id == "Supplementary_Table_S10" ~ "Pair-level analyses are supporting analyses, not validation.",
      TRUE ~ "Generated from AD61026 frozen source-data record."
    )
  )

  for (j in seq_along(source_files)) {
    df <- read_source_for_supp(source_files[j])
    sheet_base <- tools::file_path_sans_ext(basename(source_files[j]))
    sheet_name <- safe_sheet_name(paste0("S", j, "_", sheet_base))

    # Excel cannot handle empty names or duplicates.
    if (sheet_name %in% names(sheets)) {
      sheet_name <- safe_sheet_name(paste0(sheet_name, "_", j))
    }

    sheets[[sheet_name]] <- make_excel_safe(df)

    supp_qc[[length(supp_qc) + 1]] <- tibble(
      supplementary_table = supp_id,
      source_file = source_files[j],
      sheet_name = sheet_name,
      n_rows = nrow(df),
      n_cols = ncol(df),
      status = ifelse("missing_source_file" %in% colnames(df) | "read_error" %in% colnames(df), "CHECK", "PASS")
    )
  }

  if (has_writexl) {
    supp_xlsx_path <- file.path(
      FINAL_SUPP_DIR,
      paste0(supp_id, "__", str_replace_all(supp_title, "[^A-Za-z0-9]+", "_"), ".xlsx")
    )

    writexl::write_xlsx(sheets, path = supp_xlsx_path)
    cat("    Saved:", supp_xlsx_path, "\n")
  } else {
    # Fallback: save each sheet as CSV if writexl is not available.
    supp_csv_dir <- file.path(
      FINAL_SUPP_DIR,
      paste0(supp_id, "__CSV_sheets")
    )
    dir.create(supp_csv_dir, recursive = TRUE, showWarnings = FALSE)

    for (nm in names(sheets)) {
      readr::write_csv(
        sheets[[nm]],
        file.path(supp_csv_dir, paste0(nm, ".csv"))
      )
    }

    warning("writexl not installed. Saved CSV sheets for ", supp_id)
  }
}

supp_qc <- bind_rows(supp_qc)

readr::write_csv(
  supp_qc,
  file.path(FINAL_QC_DIR, "Supplementary_Table_Build_QC_AD61026.csv")
)


} else {
  supp_qc <- tibble(supplementary_table = character(), sheet_name = character(),
                    n_rows = integer(), n_cols = integer(), status = character(),
                    source_file = character())
}

# ============================================================
# SECTION 9: Final QC summary
# ============================================================

table_qc <- bind_rows(qc_rows)

readr::write_csv(
  table_qc,
  file.path(FINAL_QC_DIR, "Main_Table_QC_summary_AD61026.csv")
)

combined_qc <- bind_rows(
  table_qc %>%
    transmute(
      table_or_supplement = table_id,
      check_name,
      value,
      expected,
      status,
      note
    ),
  supp_qc %>%
    transmute(
      table_or_supplement = supplementary_table,
      check_name = paste0("Source sheet: ", sheet_name),
      value = paste0(n_rows, " rows x ", n_cols, " columns"),
      expected = "Source file readable",
      status,
      note = source_file
    )
)

readr::write_csv(
  combined_qc,
  file.path(FINAL_QC_DIR, "Table_QC_summary_AD61026.csv")
)


# ============================================================
# SECTION 10: Table legend starter text
# ============================================================
if (RUN_LEGACY_SUPPLEMENT) {
# This gives a plain text legend draft to help later manuscript writing.

legend_text <- c(
  "AD61026 manuscript table legend notes",
  "=====================================",
  "",
  "Table 1. Cohort characteristics and donor-level pseudobulk summary.",
  "This table summarizes the discovery and replication cohorts used for AD microglial pseudobulk differential expression analysis.",
  "",
  "Table 2. Strict 48-gene AD microglial consensus signature.",
  "Genes were included if they met high-confidence vote-counting criteria or had same-direction RankProd support among SEA-AD discovery genes. Direction refers to AD versus control in SEA-AD.",
  "",
  "Table 3. Pathway enrichment summary.",
  "Pathways were selected from SEA-AD genome-wide GSEA using moderated t-statistics. ORA results are provided in the supplementary tables as supportive analyses.",
  "",
  "Table 4. Final prioritized single-drug candidates.",
  "CLUE/CMap, iLINCS, and DREIMT are transcriptomic reversal support layers. DGIdb is target-level annotation only and should not be interpreted as transcriptomic reversal.",
  "",
  "Supplementary Table S10 note.",
  "Pair-level analyses provide supporting context for feasibility and mechanistic plausibility. They do not validate synergy, efficacy, or clinical safety."
)

writeLines(
  legend_text,
  con = file.path(FINAL_TABLE_DIR, "Table_legend_starter_text_AD61026.txt")
)


}

# ============================================================
# SECTION 11: Print summary
# ============================================================

cat("\n============================================================\n")
cat("  Manuscript table build complete\n")
cat("============================================================\n\n")

cat("Main table files:\n")
cat("  ", file.path(FINAL_MAIN_DIR, "Main_Table_1__cohort_characteristics.csv"), "\n")
cat("  ", file.path(FINAL_MAIN_DIR, "Main_Table_2__strict_48_consensus_signature.csv"), "\n")
cat("  ", file.path(FINAL_MAIN_DIR, "Main_Table_3__pathway_summary.csv"), "\n")
cat("  ", file.path(FINAL_MAIN_DIR, "Main_Table_4__final_single_drug_candidates.csv"), "\n\n")

if (has_writexl) {
  cat("Main table workbook:\n")
  cat("  ", file.path(FINAL_MAIN_DIR, "Main_Tables_AD61026.xlsx"), "\n\n")
}

cat("Supplementary table folder:\n")
cat("  ", FINAL_SUPP_DIR, "\n\n")

cat("QC files:\n")
cat("  ", file.path(FINAL_QC_DIR, "Main_Table_QC_summary_AD61026.csv"), "\n")
cat("  ", file.path(FINAL_QC_DIR, "Supplementary_Table_Build_QC_AD61026.csv"), "\n")
cat("  ", file.path(FINAL_QC_DIR, "Table_QC_summary_AD61026.csv"), "\n\n")

cat("Main table QC summary:\n")
print(table_qc)

if (any(table_qc$status == "FAIL")) {
  cat("\nWARNING: At least one main table QC check failed. Review before using tables.\n")
} else {
  cat("\nNo main table QC failures detected.\n")
}

cat("============================================================\n")
