# ============================================================
# Supplementary Table 7
# GSE216999 mouse-level pseudobulk orthogonal support analysis
# ============================================================
#
# Purpose:
# Create a clean supplementary Excel workbook for the manuscript.
#
# This table summarizes the orthogonal support analysis comparing
# the strict 48-gene AD microglia consensus signature with GSE216999
# mouse-level pseudobulk APP-NLGF vs APP-WT effects.
#
# Input files:
#   PSEUDOBULK_Strict48_vs_GSE216999_APP_NLGF_vs_APP_WT_summary.csv
#   PSEUDOBULK_Strict48_vs_GSE216999_APP_NLGF_vs_APP_WT_binomial_test_two_sided.csv
#   PSEUDOBULK_Strict48_vs_GSE216999_APP_NLGF_vs_APP_WT_correlation_tests.csv
#   PSEUDOBULK_Strict48_vs_GSE216999_APP_NLGF_vs_APP_WT_direction_summary.csv
#   PSEUDOBULK_Strict48_vs_GSE216999_APP_NLGF_vs_APP_WT_strict48_gene_level_output.csv
#   PSEUDOBULK_Strict48_vs_GSE216999_APP_NLGF_vs_APP_WT_mouse_cell_count_QC_after_filtering.csv
#
# Output:
#   Supplementary_Table_7_GSE216999_MouseLevel_Pseudobulk_Orthogonal_Support.xlsx
# ============================================================


# ============================================================
# 1. Load required packages
# ============================================================

# ------------------------------------------------------------
# Project folder. Set the AD_PROJECT_ROOT environment variable to the
# folder that contains Data/, RDS/, Output/ and Supplementary_Tables/,
# or start R with that folder as the working directory.
# ------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("AD_PROJECT_ROOT", unset = getwd())

packages_needed <- c(
  "dplyr",
  "readr",
  "janitor",
  "openxlsx",
  "stringr"
)

for (pkg in packages_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(dplyr)
library(readr)
library(janitor)
library(openxlsx)
library(stringr)


# ============================================================
# 2. Define input and output folders
# ============================================================
# Change input_dir if your files are saved somewhere else.
#

input_dir <- file.path(PROJECT_ROOT, "GSE216999_PSEUDOBULK_ORTHOGONAL_SUPPORT")

# Output folder.
# This keeps the final supplementary table in the same folder.
output_dir <- input_dir

# Analysis file prefix.
analysis_label <- "PSEUDOBULK_Strict48_vs_GSE216999_APP_NLGF_vs_APP_WT"

# Final Excel workbook name.
output_excel <- file.path(
  output_dir,
  "Supplementary_Table_7_GSE216999_MouseLevel_Pseudobulk_Orthogonal_Support.xlsx"
)


# ============================================================
# 3. Define input file paths
# ============================================================

summary_file <- file.path(
  input_dir,
  paste0(analysis_label, "_summary.csv")
)

binomial_file <- file.path(
  input_dir,
  paste0(analysis_label, "_binomial_test_two_sided.csv")
)

correlation_file <- file.path(
  input_dir,
  paste0(analysis_label, "_correlation_tests.csv")
)

direction_file <- file.path(
  input_dir,
  paste0(analysis_label, "_direction_summary.csv")
)

gene_file <- file.path(
  input_dir,
  paste0(analysis_label, "_strict48_gene_level_output.csv")
)

mouse_qc_file <- file.path(
  input_dir,
  paste0(analysis_label, "_mouse_cell_count_QC_after_filtering.csv")
)


# ============================================================
# 4. Check that all required files exist
# ============================================================

required_files <- c(
  summary_file,
  binomial_file,
  correlation_file,
  direction_file,
  gene_file,
  mouse_qc_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  cat("\nThe following required files were not found:\n")
  print(missing_files)
  stop("Please check input_dir or move the files into the expected folder.")
}

cat("\nAll required files found.\n")


# ============================================================
# 5. Read input tables
# ============================================================

summary_tbl <- read_csv(summary_file, show_col_types = FALSE) %>%
  clean_names()

binomial_tbl <- read_csv(binomial_file, show_col_types = FALSE) %>%
  clean_names()

correlation_tbl <- read_csv(correlation_file, show_col_types = FALSE) %>%
  clean_names()

direction_tbl <- read_csv(direction_file, show_col_types = FALSE) %>%
  clean_names()

gene_tbl <- read_csv(gene_file, show_col_types = FALSE) %>%
  clean_names()

mouse_qc_tbl <- read_csv(mouse_qc_file, show_col_types = FALSE) %>%
  clean_names()


# ============================================================
# 6. Build clean Summary sheet
# ============================================================
# This combines:
#   overall concordance summary
#   two-sided binomial sign test
#   Pearson/Spearman correlation
#   mouse and cell counts
#
# The key manuscript numbers should be:
#   7 total mice
#   5 APP-NLGF mice
#   2 APP-WT mice
#   8,555 cells
#   28 strict genes detected
#   23/28 concordant
#   two-sided binomial P = 9.12e-4
#   Pearson r = 0.523, P = 0.00431
#   Spearman rho = 0.643, P = 2.27e-4

n_wt_mice <- mouse_qc_tbl %>%
  filter(contrast_group == "WT") %>%
  distinct(mouse_id) %>%
  nrow()

n_nlgf_mice <- mouse_qc_tbl %>%
  filter(contrast_group == "NLGF") %>%
  distinct(mouse_id) %>%
  nrow()

n_total_mice <- mouse_qc_tbl %>%
  distinct(mouse_id) %>%
  nrow()

n_wt_cells <- mouse_qc_tbl %>%
  filter(contrast_group == "WT") %>%
  summarize(n = sum(n_cells, na.rm = TRUE)) %>%
  pull(n)

n_nlgf_cells <- mouse_qc_tbl %>%
  filter(contrast_group == "NLGF") %>%
  summarize(n = sum(n_cells, na.rm = TRUE)) %>%
  pull(n)

n_total_cells <- sum(mouse_qc_tbl$n_cells, na.rm = TRUE)

summary_final <- summary_tbl %>%
  bind_cols(
    binomial_tbl %>%
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
    n_wt_mice = n_wt_mice,
    n_nlgf_mice = n_nlgf_mice,
    n_total_mice = n_total_mice,
    n_wt_cells = n_wt_cells,
    n_nlgf_cells = n_nlgf_cells,
    n_total_cells = n_total_cells,
    min_cells_per_mouse = 50
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


# ============================================================
# 7. Build clean Mouse_QC sheet
# ============================================================

mouse_qc_final <- mouse_qc_tbl %>%
  arrange(contrast_group, mouse_id) %>%
  select(
    mouse_id,
    contrast_group,
    n_cells
  )


# ============================================================
# 8. Build clean Direction_summary sheet
# ============================================================

direction_final <- direction_tbl %>%
  arrange(seaad_direction) %>%
  select(
    seaad_direction,
    detected,
    concordant,
    discordant,
    concordance_rate
  )


# ============================================================
# 9. Build clean Strict48_gene_level sheet
# ============================================================
# This provides the gene-level evidence used in the manuscript.
#
# It includes:
#   SEA-AD logFC and direction
#   GSE216999 pseudobulk logFC and direction
#   pseudobulk P value and adjusted P value
#   concordance status

gene_final <- gene_tbl %>%
  select(
    gene,
    seaad_log_fc,
    seaad_direction,
    pseudobulk_log_fc,
    pseudobulk_direction,
    pseudobulk_ave_expr,
    pseudobulk_p_value,
    pseudobulk_adj_p_value,
    detected_in_pseudobulk,
    concordant_direction,
    concordance_label
  ) %>%
  arrange(
    desc(detected_in_pseudobulk),
    desc(concordant_direction),
    seaad_direction,
    gene
  )


# ============================================================
# 10. Build Legend sheet
# ============================================================

legend_tbl <- tibble::tibble(
  item = c(
    "Supplementary Table 7 title",
    "Dataset",
    "Analysis subset",
    "Pseudobulk unit",
    "Differential-expression method",
    "Concordance definition",
    "Statistical tests",
    "Important interpretation",
    "Sheet: Summary",
    "Sheet: Mouse_QC",
    "Sheet: Direction_summary",
    "Sheet: Strict48_gene_level"
  ),
  description = c(
    "Orthogonal support analysis using GSE216999 mouse-level pseudobulk human microglia xenograft data.",
    "GSE216999 from Mancuso et al.; human microglia xenografted into amyloid model and control mouse brains.",
    "WT H9 xenografted human microglia from 6-month APP-NLGF and APP-WT host brains with no injection treatment.",
    "Cells were aggregated to mouse-level pseudobulk profiles by summing raw counts for each gene within each mouse.",
    "edgeR normalization followed by limma-voom; APP-WT was used as the reference group.",
    "A strict consensus gene was considered concordant if the SEA-AD AD-versus-control logFC and GSE216999 APP-NLGF-versus-APP-WT pseudobulk logFC had the same sign.",
    "Directional concordance was assessed using a two-sided binomial sign test. Pearson and Spearman correlations were calculated across overlapping detected genes.",
    "This analysis was used as orthogonal support only and was not used to define the strict 48-gene consensus signature.",
    "Overall concordance, mouse/cell counts, binomial sign-test result, and Pearson/Spearman correlation statistics.",
    "Cell counts per mouse after filtering.",
    "Concordance summarized separately for SEA-AD AD-upregulated and AD-downregulated genes.",
    "Gene-level SEA-AD logFC, GSE216999 pseudobulk logFC, differential-expression statistics, and concordance status."
  )
)


# ============================================================
# 11. Create Excel workbook
# ============================================================

wb <- createWorkbook()

addWorksheet(wb, "Legend")
addWorksheet(wb, "Summary")
addWorksheet(wb, "Mouse_QC")
addWorksheet(wb, "Direction_summary")
addWorksheet(wb, "Consensus48_gene_level")

writeData(wb, "Legend", legend_tbl)
writeData(wb, "Summary", summary_final)
writeData(wb, "Mouse_QC", mouse_qc_final)
writeData(wb, "Direction_summary", direction_final)
writeData(wb, "Consensus48_gene_level", gene_final)


# ============================================================
# 12. Format workbook
# ============================================================

# Header style.
header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

# Wrap text style.
wrap_style <- createStyle(
  wrapText = TRUE,
  valign = "top"
)

# Apply header style to all sheets.
sheet_names <- names(wb)

for (sheet in sheet_names) {
  addStyle(
    wb,
    sheet = sheet,
    style = header_style,
    rows = 1,
    cols = 1:100,
    gridExpand = TRUE,
    stack = TRUE
  )
  
  freezePane(wb, sheet = sheet, firstRow = TRUE)
}

# Set readable widths.
setColWidths(wb, "Legend", cols = 1, widths = 28)
setColWidths(wb, "Legend", cols = 2, widths = 95)
addStyle(
  wb,
  sheet = "Legend",
  style = wrap_style,
  rows = 1:(nrow(legend_tbl) + 1),
  cols = 1:2,
  gridExpand = TRUE,
  stack = TRUE
)

setColWidths(wb, "Summary", cols = 1:ncol(summary_final), widths = "auto")
setColWidths(wb, "Mouse_QC", cols = 1:ncol(mouse_qc_final), widths = "auto")
setColWidths(wb, "Direction_summary", cols = 1:ncol(direction_final), widths = "auto")
setColWidths(wb, "Consensus48_gene_level", cols = 1:ncol(gene_final), widths = "auto")


# ============================================================
# 13. Save workbook
# ============================================================

saveWorkbook(
  wb,
  output_excel,
  overwrite = TRUE
)

cat("\nSupplementary Table 7 saved to:\n")
cat(output_excel, "\n")


# ============================================================
# 14. Print key manuscript numbers for quick check
# ============================================================

cat("\nKey Supplementary Table 7 values:\n")
cat("Total mice:", n_total_mice, "\n")
cat("APP-NLGF mice:", n_nlgf_mice, "\n")
cat("APP-WT mice:", n_wt_mice, "\n")
cat("Total cells:", n_total_cells, "\n")
cat("Strict 48 genes detected:", summary_final$detected_in_pseudobulk, "\n")
cat("Concordant:", summary_final$concordant, "\n")
cat("Discordant:", summary_final$discordant, "\n")
cat("Concordance rate:", summary_final$concordance_rate, "\n")
cat("Two-sided binomial P:", summary_final$binomial_p_value_two_sided, "\n")
cat("Pearson r:", summary_final$pearson_r, "\n")
cat("Pearson P:", summary_final$pearson_p_value, "\n")
cat("Spearman rho:", summary_final$spearman_rho, "\n")
cat("Spearman P:", summary_final$spearman_p_value, "\n")