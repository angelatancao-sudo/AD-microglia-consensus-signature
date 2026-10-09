# ============================================================
# 06_consensus_signature_and_source_data.R
# ============================================================
# Project: AD Microglia Multi-Dataset Transcriptomic Analysis
#
# Purpose:
#   Prepare clean source-data files for future manuscript tables,
#   supplementary tables, and figures.
#
# Important:
#   This script does NOT polish final manuscript figures yet.
#   It organizes the source data first, so later figures/tables can
#   be generated reproducibly from saved CSV files.
#
# Why this is different from the old 07 script:
#   The old script generated final figures and some hard-coded tables.
#   The new AD61026 workflow is still evolving through pathway analysis
#   and drug prioritization, so this script focuses on source-data
#   organization rather than final publication formatting.
#
# Run after:
#   04_DE_analysis_AD61026_FIELD_ALIGNED.R
#   and after RankProd has successfully produced:
#     rankprod_primary_results.csv
#
# Recommended run command:
#   source("R/07_prepare_manuscript_source_data_AD61026_FIXED.R")
#
# Note:
#   Run with source(). Do not paste the whole script line-by-line into
#   the console, because errors can cause later lines to continue running
#   without earlier objects being created.
#
# Main outputs:
#   Output/Manuscript_Source_Data/
#      README_manifest.csv
#      MainTable1_cohort_summary_source.csv
#      MainTable2_strict_consensus_signature_source.csv
#      SuppTable1_pseudobulk_input_check_source.csv
#      SuppTable2_microglia_cluster_decisions_source.csv
#      SuppTable3_SEAAD_DE_full_source.csv
#      SuppTable4_replication_assessment_all_discovery_genes_source.csv
#      SuppTable5_RankProd_primary_results_source.csv
#      SuppTable6_expanded_drug_signature_source.csv
#      SuppTable7_sensitivity_comparison_source.csv
#      SuppTable8_GSE243292_supplementary_support_source.csv
#      Figure1_workflow_source.csv
#      Figure3_SEAAD_volcano_source.csv
#      Figure4_replication_heatmap_strict_source.csv
#      Figure4_replication_heatmap_expanded_source.csv
#      Figure5_signature_summary_source.csv
#
# Key definitions used here:
#   Strict primary consensus signature:
#       high-confidence vote-counting genes
#       OR
#       RankProd pfp < 0.05 genes that are also SEA-AD discovery genes
#       and match the SEA-AD discovery direction.
#
#   Expanded drug sensitivity signature:
#       high-confidence genes
#       OR
#       moderate-confidence supportive genes.
#
#   GSE243292:
#       kept as supplementary only, not part of primary replication
#       confidence definitions.
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

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

cat("============================================================\n")
cat("  07: Prepare manuscript source data\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Define the main Output directory robustly.
# Some earlier setup scripts define TAB_DIR, DE_DIR, FIG_DIR, and RDS_DIR
# but may not define OUT_DIR explicitly. This fallback prevents the
# "object 'OUT_DIR' not found" error.
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

SOURCE_DIR <- file.path(OUT_DIR, "Manuscript_Source_Data")
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Output directory:\n")
cat("  ", OUT_DIR, "\n")
cat("Source-data directory:\n")
cat("  ", SOURCE_DIR, "\n\n")

# ============================================================
# SECTION 2: Helper functions
# ============================================================

# ------------------------------------------------------------
# Read a CSV file if it exists.
# If required = TRUE, the script stops if the file is missing.
# If required = FALSE, the script continues and returns NULL.
# ------------------------------------------------------------
read_csv_if_exists <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) {
      stop("Required file missing: ", path)
    } else {
      warning("Optional file missing: ", path)
      return(NULL)
    }
  }

  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

# ------------------------------------------------------------
# Create a manifest so every source-data file is traceable.
# ------------------------------------------------------------
manifest_rows <- list()

add_manifest <- function(file_name,
                         intended_use,
                         source_inputs,
                         notes = "") {
  manifest_rows[[length(manifest_rows) + 1]] <<- data.frame(
    file_name = file_name,
    intended_use = intended_use,
    source_inputs = source_inputs,
    notes = notes,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Write one source-data file and record it in the manifest.
# ------------------------------------------------------------
write_source <- function(df,
                         file_name,
                         intended_use,
                         source_inputs,
                         notes = "") {

  out_path <- file.path(SOURCE_DIR, file_name)
  write.csv(df, out_path, row.names = FALSE)

  add_manifest(
    file_name = file_name,
    intended_use = intended_use,
    source_inputs = source_inputs,
    notes = notes
  )

  cat("Saved:", out_path, "\n")
}

# ------------------------------------------------------------
# Find the moderate-confidence file, since the name changed slightly
# between 04 script versions.
# ------------------------------------------------------------
get_moderate_file <- function() {
  f_supportive <- file.path(DE_DIR, "moderate_confidence_genes_SUPPORTIVE.csv")
  f_primary <- file.path(DE_DIR, "moderate_confidence_genes_PRIMARY.csv")

  if (file.exists(f_supportive)) return(f_supportive)
  if (file.exists(f_primary)) return(f_primary)

  stop("Could not find a moderate-confidence gene file.")
}

# ------------------------------------------------------------
# Convert RankProd direction labels to simple UP/DOWN labels.
# ------------------------------------------------------------
standardize_rankprod_direction <- function(x) {
  case_when(
    x %in% c("UP", "Up", "up", "UP_in_AD", "AD_upregulated", "Upregulated") ~ "UP",
    x %in% c("DOWN", "Down", "down", "DOWN_in_AD", "AD_downregulated", "Downregulated") ~ "DOWN",
    str_detect(as.character(x), regex("^UP", ignore_case = TRUE)) ~ "UP",
    str_detect(as.character(x), regex("^DOWN", ignore_case = TRUE)) ~ "DOWN",
    TRUE ~ NA_character_
  )
}

# ============================================================
# SECTION 3: Load current AD61026 outputs
# ============================================================

cat("Loading current AD61026 output files...\n\n")

pseudobulk_check <- read_csv_if_exists(
  file.path(TAB_DIR, "03F_pseudobulk_input_check_summary.csv")
)

seaad_de <- read_csv_if_exists(
  file.path(DE_DIR, "SEAAD_DE_results.csv")
)

seaad_discovery_genes <- read_csv_if_exists(
  file.path(DE_DIR, "SEAAD_discovery_genes.csv")
)

replication_primary <- read_csv_if_exists(
  file.path(DE_DIR, "replication_assessment_PRIMARY.csv")
)

high_conf <- read_csv_if_exists(
  file.path(DE_DIR, "high_confidence_genes_PRIMARY.csv")
)

moderate_conf <- read_csv_if_exists(
  get_moderate_file()
)

rankprod_results <- read_csv_if_exists(
  file.path(DE_DIR, "rankprod_primary_results.csv"),
  required = FALSE
)

rankprod_up <- read_csv_if_exists(
  file.path(DE_DIR, "rankprod_primary_AD_upregulated.csv"),
  required = FALSE
)

rankprod_down <- read_csv_if_exists(
  file.path(DE_DIR, "rankprod_primary_AD_downregulated.csv"),
  required = FALSE
)

supp_gse243292 <- read_csv_if_exists(
  file.path(DE_DIR, "supplementary_support_GSE243292.csv"),
  required = FALSE
)

sensitivity_comparison <- read_csv_if_exists(
  file.path(DE_DIR, "sensitivity_comparison_summary.csv"),
  required = FALSE
)

# ============================================================
# SECTION 4: Reconstruct strict and expanded signatures
# ============================================================
# We recompute these here so the source data are reproducible from
# the 04 outputs and do not depend on manually edited files.

cat("Reconstructing strict and expanded gene signatures...\n\n")

# High-confidence genes from vote-counting.
high_set <- unique(high_conf$gene)

# Moderate-confidence genes from vote-counting.
moderate_set <- unique(moderate_conf$gene)

# RankProd support is used only if the gene:
#   1. has pfp < 0.05 in RankProd,
#   2. is one of the SEA-AD discovery genes assessed in the replication table,
#   3. has the same direction as the SEA-AD discovery result.
if (!is.null(rankprod_results)) {

  rankprod_sig <- rankprod_results %>%
    mutate(
      rankprod_direction_short =
        standardize_rankprod_direction(rankprod_direction)
    ) %>%
    filter(!is.na(rankprod_direction_short),
           pfp < 0.05)

  # Join to SEA-AD discovery-anchored replication table.
  rankprod_same_direction <- replication_primary %>%
    select(
      gene,
      discovery_logFC,
      discovery_adj.P.Val,
      discovery_direction,
      confidence,
      n_same_direction,
      n_nominal_p_lt_0.05,
      mean_replication_logFC
    ) %>%
    inner_join(
      rankprod_sig %>%
        select(
          gene,
          RankProd_direction = rankprod_direction,
          RankProd_direction_short = rankprod_direction_short,
          RankProd_pfp = pfp,
          RankProd_pval = pval,
          RankProd_AveFC = AveFC
        ),
      by = "gene"
    ) %>%
    filter(RankProd_direction_short == discovery_direction) %>%
    arrange(RankProd_pfp) %>%
    group_by(gene) %>%
    slice(1) %>%
    ungroup()

} else {

  rankprod_sig <- NULL

  rankprod_same_direction <- data.frame(
    gene = character(),
    discovery_logFC = numeric(),
    discovery_adj.P.Val = numeric(),
    discovery_direction = character(),
    confidence = character(),
    n_same_direction = integer(),
    n_nominal_p_lt_0.05 = integer(),
    mean_replication_logFC = numeric(),
    RankProd_direction = character(),
    RankProd_direction_short = character(),
    RankProd_pfp = numeric(),
    RankProd_pval = numeric(),
    RankProd_AveFC = numeric(),
    stringsAsFactors = FALSE
  )
}

rankprod_same_set <- unique(rankprod_same_direction$gene)

# Strict signature:
# high-confidence vote-counting OR same-direction RankProd support.
strict_set <- union(high_set, rankprod_same_set)

# Expanded drug sensitivity signature:
# high-confidence OR moderate-confidence supportive genes.
expanded_set <- union(high_set, moderate_set)

# Build strict table.
strict_signature <- replication_primary %>%
  filter(gene %in% strict_set) %>%
  select(
    gene,
    discovery_logFC,
    discovery_adj.P.Val,
    discovery_direction,
    confidence,
    n_same_direction,
    n_nominal_p_lt_0.05,
    mean_replication_logFC,
    everything()
  ) %>%
  mutate(
    high_confidence_vote_counting = gene %in% high_set,
    rankprod_pfp_lt_0.05_same_direction = gene %in% rankprod_same_set,
    strict_consensus_source = case_when(
      high_confidence_vote_counting & rankprod_pfp_lt_0.05_same_direction ~
        "High_confidence_and_RankProd",
      high_confidence_vote_counting ~
        "High_confidence_only",
      rankprod_pfp_lt_0.05_same_direction ~
        "RankProd_same_direction_only",
      TRUE ~ "Other"
    )
  ) %>%
  left_join(
    rankprod_same_direction %>%
      select(gene, RankProd_direction, RankProd_pfp,
             RankProd_pval, RankProd_AveFC),
    by = "gene"
  ) %>%
  arrange(discovery_direction, discovery_adj.P.Val, gene)

# Build expanded table.
expanded_signature <- replication_primary %>%
  filter(gene %in% expanded_set) %>%
  select(
    gene,
    discovery_logFC,
    discovery_adj.P.Val,
    discovery_direction,
    confidence,
    n_same_direction,
    n_nominal_p_lt_0.05,
    mean_replication_logFC,
    everything()
  ) %>%
  mutate(
    high_confidence_vote_counting = gene %in% high_set,
    moderate_confidence_supportive = gene %in% moderate_set,
    rankprod_pfp_lt_0.05_same_direction = gene %in% rankprod_same_set,
    strict_primary_consensus = gene %in% strict_set,
    recommended_use = "Expanded_drug_sensitivity_signature"
  ) %>%
  left_join(
    rankprod_same_direction %>%
      select(gene, RankProd_direction, RankProd_pfp,
             RankProd_pval, RankProd_AveFC),
    by = "gene"
  ) %>%
  arrange(discovery_direction, discovery_adj.P.Val, gene)

# Summary of signature construction.
rankprod_sig_n <- if (!is.null(rankprod_sig)) {
  length(unique(rankprod_sig$gene))
} else {
  0
}

signature_summary <- data.frame(
  category = c(
    "SEA-AD discovery genes assessed",
    "High-confidence vote-counting genes",
    "Moderate-confidence supportive genes",
    "RankProd pfp < 0.05 genes across primary replication datasets",
    "RankProd pfp < 0.05 genes also SEA-AD discovery and same direction",
    "Overlap: high-confidence and same-direction RankProd",
    "Strict primary consensus: high-confidence OR same-direction RankProd",
    "Expanded drug sensitivity signature: high + moderate"
  ),
  n_genes = c(
    nrow(replication_primary),
    length(high_set),
    length(moderate_set),
    rankprod_sig_n,
    length(rankprod_same_set),
    length(intersect(high_set, rankprod_same_set)),
    length(strict_set),
    length(expanded_set)
  ),
  stringsAsFactors = FALSE
)

# Save the reconstructed signatures.
write_source(
  strict_signature,
  "MainTable2_strict_consensus_signature_source.csv",
  intended_use = "Main Table 2: strict primary consensus signature.",
  source_inputs = "replication_assessment_PRIMARY.csv; high_confidence_genes_PRIMARY.csv; rankprod_primary_results.csv",
  notes = "Strict consensus is high-confidence vote-counting OR same-direction RankProd pfp < 0.05 among SEA-AD discovery genes."
)

write_source(
  expanded_signature,
  "SuppTable6_expanded_drug_signature_source.csv",
  intended_use = "Supplementary table / drug input source: expanded high + moderate signature.",
  source_inputs = "replication_assessment_PRIMARY.csv; high_confidence_genes_PRIMARY.csv; moderate_confidence_genes_SUPPORTIVE.csv",
  notes = "Expanded signature is for drug-prioritization sensitivity, not the strict biological consensus."
)

write_source(
  signature_summary,
  "Figure5_signature_summary_source.csv",
  intended_use = "Figure/source summary of strict and expanded signature construction.",
  source_inputs = "Current 04 DE/replication/RankProd outputs.",
  notes = "Useful for manuscript text and workflow figure labels."
)

write_source(
  rankprod_same_direction,
  "SuppTable5B_RankProd_same_direction_discovery_source.csv",
  intended_use = "Supplementary RankProd support table restricted to SEA-AD discovery genes with same-direction support.",
  source_inputs = "rankprod_primary_results.csv; replication_assessment_PRIMARY.csv",
  notes = "These genes contribute to the strict signature if not already high-confidence."
)

# ============================================================
# SECTION 5: Cohort and pseudobulk source tables
# ============================================================

cat("Preparing cohort and pseudobulk source tables...\n\n")

write_source(
  pseudobulk_check,
  "SuppTable1_pseudobulk_input_check_source.csv",
  intended_use = "Supplementary Table 1: pseudobulk sample counts, AD/control donor counts, metadata alignment, and microglia cell counts.",
  source_inputs = "Output/Tables/03F_pseudobulk_input_check_summary.csv",
  notes = "Audit table showing all pseudobulk inputs were present and aligned before DE."
)

main_table1 <- pseudobulk_check %>%
  filter(analysis %in% c("DISCOVERY", "PRIMARY", "SUPPLEMENTARY")) %>%
  mutate(
    role = case_when(
      dataset == "SEAAD" ~ "Discovery",
      dataset == "GSE243292" ~ "Supplementary support",
      TRUE ~ "Primary replication"
    ),
    microglia_source = case_when(
      dataset == "SEAAD" ~ "Author SEA-AD Microglia/Immune subset",
      dataset == "GSE174367" ~ "Author cell-type annotation, MG",
      dataset %in% c("GSE157827", "GSE160936", "GSE188545") ~
        "Cluster-validated microglia from raw snRNA-seq",
      dataset == "GSE243292" ~ "Microglia-specific/preannotated H5AD",
      TRUE ~ NA_character_
    ),
    analysis_note = case_when(
      dataset == "GSE160936" ~ "EC and SSC combined by donor before DE",
      dataset == "GSE243292" ~ "A+T+ as AD, A-T- as Control, A+T- excluded; supplementary due to 2 controls",
      dataset == "GSE188545" ~ "One donor has low extracted microglia count; retained with QC flag",
      dataset == "GSE157827" ~ "One original AD donor not represented after validated microglia extraction",
      TRUE ~ ""
    )
  ) %>%
  select(
    dataset,
    role,
    analysis,
    n_pseudobulk_samples,
    n_AD,
    n_Control,
    n_genes,
    min_microglia_cells,
    median_microglia_cells,
    max_microglia_cells,
    microglia_source,
    analysis_note,
    status,
    metadata_aligned
  )

write_source(
  main_table1,
  "MainTable1_cohort_summary_source.csv",
  intended_use = "Main Table 1: cohort and pseudobulk summary.",
  source_inputs = "Output/Tables/03F_pseudobulk_input_check_summary.csv",
  notes = "Use this as the main cohort table source."
)

figure1_source <- main_table1 %>%
  select(
    dataset,
    role,
    n_pseudobulk_samples,
    n_AD,
    n_Control,
    microglia_source,
    analysis_note
  )

write_source(
  figure1_source,
  "Figure1_workflow_source.csv",
  intended_use = "Figure 1 workflow/cohort count source data.",
  source_inputs = "MainTable1_cohort_summary_source.csv",
  notes = "Use this to label the study-design workflow figure."
)

# ============================================================
# SECTION 6: Microglia cluster decision source table
# ============================================================
# This table documents cluster inclusion/exclusion decisions for
# the three raw unannotated GEO datasets.

cat("Preparing microglia cluster decision source table...\n\n")

cluster_decisions <- data.frame(
  dataset = c(
    "GSE157827", "GSE157827", "GSE157827",
    "GSE160936", "GSE160936", "GSE160936", "GSE160936", "GSE160936",
    "GSE188545", "GSE188545", "GSE188545", "GSE188545", "GSE188545"
  ),
  cluster_id = c(
    "7", "28", "29",
    "1", "2", "6", "9", "27",
    "7", "18", "23", "29", "32"
  ),
  final_decision = c(
    "Primary", "Exclude", "Exclude",
    "Primary", "Primary", "Primary", "Primary", "Sensitivity_only",
    "Primary", "Primary", "Exclude", "Sensitivity_only", "Exclude"
  ),
  include_primary = c(
    TRUE, FALSE, FALSE,
    TRUE, TRUE, TRUE, TRUE, FALSE,
    TRUE, TRUE, FALSE, FALSE, FALSE
  ),
  include_sensitivity = c(
    TRUE, FALSE, FALSE,
    TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, FALSE, TRUE, FALSE
  ),
  rationale = c(
    "Clean microglia-enriched cluster with broad donor representation and markers including LRMDA, PLXDC2, DOCK8, APBB1IP, C3, and CSF1R.",
    "Excluded because of strong astrocyte-like signal.",
    "Excluded because of OPC-like or mixed signal including VCAN.",
    "Homeostatic microglia-like cluster with P2RY12, APBB1IP, CSF1R, and DOCK8 support.",
    "Microglia-enriched cluster with PLXDC2, LRMDA, RUNX1, APBB1IP, MEF2C, and DOCK8 support.",
    "Microglia/immune-enriched AD-associated cluster with MS4A6A, FKBP5, DOCK8, PLXDC2, and ARHGAP15 support.",
    "Disease-associated microglia-like cluster with SPP1, SLC11A1, LRMDA, TNFRSF1B, and PTPRC support.",
    "Small cell-cycle/proliferating microglia-like cluster; retained only for sensitivity.",
    "Clear microglia-enriched cluster with LRMDA, ARHGAP24, DOCK8, PLXDC2, APBB1IP, and CSF1R support.",
    "AD-enriched disease-associated microglia cluster with SPP1, APBB1IP, CD74, DOCK8, PTPRC, and SLC11A1 support.",
    "Excluded because of astrocyte-like signal including SLC1A3 and AQP4.",
    "Small cluster with microglial markers but high oligodendrocyte/myelin signal; sensitivity only.",
    "Excluded because top markers indicate lymphoid/T-cell-like identity including SKAP1, ITK, IL7R, and CCL5."
  ),
  stringsAsFactors = FALSE
)

write_source(
  cluster_decisions,
  "SuppTable2_microglia_cluster_decisions_source.csv",
  intended_use = "Supplementary Table 2: validated microglia cluster annotation decisions.",
  source_inputs = "03D validation outputs: candidate validation summaries, top candidate markers, dot plots.",
  notes = "Primary clusters are used for main DE; sensitivity-only clusters are used only in cluster-definition sensitivity analysis."
)

# ============================================================
# SECTION 7: DE, replication, RankProd, supplementary, sensitivity
# ============================================================

cat("Preparing DE and replication source tables...\n\n")

write_source(
  seaad_de,
  "SuppTable3_SEAAD_DE_full_source.csv",
  intended_use = "Supplementary Table 3: full SEA-AD discovery DE results.",
  source_inputs = "Output/DE/SEAAD_DE_results.csv",
  notes = "Positive logFC means higher in AD."
)

write_source(
  seaad_discovery_genes,
  "SuppTable3B_SEAAD_discovery_genes_source.csv",
  intended_use = "Supplementary source: SEA-AD discovery genes at FDR < 0.05 and |logFC| > 0.25.",
  source_inputs = "Output/DE/SEAAD_discovery_genes.csv",
  notes = "These genes were tested for replication."
)

write_source(
  replication_primary,
  "SuppTable4_replication_assessment_all_discovery_genes_source.csv",
  intended_use = "Supplementary Table 4: replication assessment for all SEA-AD discovery genes.",
  source_inputs = "Output/DE/replication_assessment_PRIMARY.csv",
  notes = "High and moderate are operational vote-counting categories across the 4 primary replication datasets."
)

if (!is.null(rankprod_results)) {
  write_source(
    rankprod_results,
    "SuppTable5_RankProd_primary_results_source.csv",
    intended_use = "Supplementary Table 5: full RankProd primary replication results.",
    source_inputs = "Output/DE/rankprod_primary_results.csv",
    notes = "RankProd was run across the 4 primary replication datasets only."
  )
}

if (!is.null(rankprod_up)) {
  write_source(
    rankprod_up,
    "SuppTable5C_RankProd_primary_AD_upregulated_source.csv",
    intended_use = "Supplementary RankProd AD-upregulated results.",
    source_inputs = "Output/DE/rankprod_primary_AD_upregulated.csv",
    notes = "Optional directional RankProd source table."
  )
}

if (!is.null(rankprod_down)) {
  write_source(
    rankprod_down,
    "SuppTable5D_RankProd_primary_AD_downregulated_source.csv",
    intended_use = "Supplementary RankProd AD-downregulated results.",
    source_inputs = "Output/DE/rankprod_primary_AD_downregulated.csv",
    notes = "Optional directional RankProd source table."
  )
}

if (!is.null(sensitivity_comparison)) {
  write_source(
    sensitivity_comparison,
    "SuppTable7_sensitivity_comparison_source.csv",
    intended_use = "Supplementary Table 7: primary versus sensitivity cluster-definition comparison.",
    source_inputs = "Output/DE/sensitivity_comparison_summary.csv",
    notes = "Assesses whether including borderline microglia-like clusters changes confidence categories."
  )
}

if (!is.null(supp_gse243292)) {
  write_source(
    supp_gse243292,
    "SuppTable8_GSE243292_supplementary_support_source.csv",
    intended_use = "Supplementary Table 8: GSE243292 supplementary support.",
    source_inputs = "Output/DE/supplementary_support_GSE243292.csv",
    notes = "GSE243292 is supplementary because only 2 control donors remain after A/T filtering."
  )
}

# ============================================================
# SECTION 8: Figure source-data files
# ============================================================

cat("Preparing figure source-data files...\n\n")

# ------------------------------------------------------------
# Figure 3 source: SEA-AD volcano plot.
# ------------------------------------------------------------
volcano_source <- seaad_de %>%
  mutate(
    significance_group = case_when(
      adj.P.Val < 0.05 & logFC > 0.25 ~ "AD_up_FDR_lt_0.05_logFC_gt_0.25",
      adj.P.Val < 0.05 & logFC < -0.25 ~ "AD_down_FDR_lt_0.05_logFC_lt_minus_0.25",
      TRUE ~ "Not_significant"
    ),
    in_strict_consensus = gene %in% strict_set,
    in_expanded_signature = gene %in% expanded_set,
    minus_log10_p = -log10(P.Value)
  )

write_source(
  volcano_source,
  "Figure3_SEAAD_volcano_source.csv",
  intended_use = "Figure source: SEA-AD discovery volcano plot.",
  source_inputs = "Output/DE/SEAAD_DE_results.csv; strict/expanded signatures reconstructed in 07.",
  notes = "Use in_strict_consensus or in_expanded_signature to choose labels."
)

# ------------------------------------------------------------
# Figure 4 source: replication heatmap for strict signature.
# ------------------------------------------------------------
make_heatmap_source <- function(replication_df, genes, signature_label) {

  dataset_cols <- grep("_logFC$", colnames(replication_df), value = TRUE)

  heatmap_df <- replication_df %>%
    filter(gene %in% genes) %>%
    select(gene, discovery_logFC, discovery_direction, confidence,
           all_of(dataset_cols)) %>%
    pivot_longer(
      cols = all_of(dataset_cols),
      names_to = "dataset",
      values_to = "logFC"
    ) %>%
    mutate(
      dataset = gsub("_logFC$", "", dataset),
      signature = signature_label
    )

  # Add SEA-AD discovery as its own dataset column for plotting.
  seaad_df <- replication_df %>%
    filter(gene %in% genes) %>%
    transmute(
      gene = gene,
      discovery_direction = discovery_direction,
      confidence = confidence,
      dataset = "SEAAD",
      logFC = discovery_logFC,
      signature = signature_label
    )

  bind_rows(seaad_df, heatmap_df) %>%
    arrange(gene, dataset)
}

heatmap_strict <- make_heatmap_source(
  replication_df = replication_primary,
  genes = strict_set,
  signature_label = "strict_primary_consensus"
)

heatmap_expanded <- make_heatmap_source(
  replication_df = replication_primary,
  genes = expanded_set,
  signature_label = "expanded_drug_sensitivity_signature"
)

write_source(
  heatmap_strict,
  "Figure4_replication_heatmap_strict_source.csv",
  intended_use = "Figure source: replication logFC heatmap for strict consensus genes.",
  source_inputs = "replication_assessment_PRIMARY.csv; strict signature reconstructed in 07.",
  notes = "Rows are genes; columns are SEAAD and the 4 primary replication datasets."
)

write_source(
  heatmap_expanded,
  "Figure4_replication_heatmap_expanded_source.csv",
  intended_use = "Figure source: replication logFC heatmap for expanded high+moderate genes.",
  source_inputs = "replication_assessment_PRIMARY.csv; expanded signature reconstructed in 07.",
  notes = "Useful for supplementary heatmap or drug-signature visualization."
)

# ------------------------------------------------------------
# Figure 5 source: signature category counts.
# ------------------------------------------------------------
strict_category_counts <- strict_signature %>%
  count(discovery_direction, strict_consensus_source, name = "n_genes") %>%
  mutate(signature = "strict_primary_consensus")

expanded_category_counts <- expanded_signature %>%
  count(discovery_direction, confidence, name = "n_genes") %>%
  mutate(signature = "expanded_drug_sensitivity_signature") %>%
  rename(category = confidence)

strict_category_counts2 <- strict_category_counts %>%
  rename(category = strict_consensus_source)

signature_category_source <- bind_rows(
  strict_category_counts2,
  expanded_category_counts
)

write_source(
  signature_category_source,
  "Figure5_signature_category_counts_source.csv",
  intended_use = "Figure source: bar plot or stacked summary of signature categories.",
  source_inputs = "Strict and expanded signatures reconstructed in 07.",
  notes = "Can be used to show UP/DOWN and evidence-source composition."
)

# ============================================================
# SECTION 9: Optional source placeholders for future modules
# ============================================================
# These are not written as data files yet, because pathway and drug
# prioritization should be rerun after the final gene signatures are
# confirmed. The manifest records what should be added later.

add_manifest(
  file_name = "Future_MainTable3_pathway_results_source.csv",
  intended_use = "Main/supplementary pathway tables after 05 pathway analysis is rerun.",
  source_inputs = "Future 05_pathway_analysis_AD61026 outputs.",
  notes = "Do not reuse old hard-coded pathway values; regenerate from the new 48-gene strict signature and SEA-AD ranked list."
)

add_manifest(
  file_name = "Future_MainTable4_drug_prioritization_source.csv",
  intended_use = "Main drug-prioritization table after 06 drug prioritization is rerun.",
  source_inputs = "Future 06_drug_prioritization_AD61026 outputs.",
  notes = "Use strict 48-gene signature as primary query and expanded 125-gene signature as sensitivity query."
)

add_manifest(
  file_name = "Future_Figure6_drug_prioritization_source.csv",
  intended_use = "Drug bubble/heatmap source data after drug prioritization is rerun.",
  source_inputs = "Future 06 drug outputs.",
  notes = "Only compounds supported by current AD61026 signatures should be included."
)

# ============================================================
# SECTION 10: Write manifest and final summary
# ============================================================

manifest <- bind_rows(manifest_rows)

write.csv(
  manifest,
  file.path(SOURCE_DIR, "README_manifest.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("07 manuscript source-data preparation complete.\n\n")
cat("Source-data folder:\n")
cat("  ", SOURCE_DIR, "\n\n")
cat("Key counts:\n")
cat("  SEA-AD discovery genes assessed:", nrow(replication_primary), "\n")
cat("  High-confidence genes:", length(high_set), "\n")
cat("  Moderate-confidence genes:", length(moderate_set), "\n")
cat("  Same-direction RankProd-supported SEA-AD discovery genes:",
    length(rankprod_same_set), "\n")
cat("  Strict primary consensus genes:", length(strict_set), "\n")
cat("  Expanded drug sensitivity genes:", length(expanded_set), "\n\n")
cat("Manifest saved:\n")
cat("  ", file.path(SOURCE_DIR, "README_manifest.csv"), "\n")
cat("============================================================\n")
