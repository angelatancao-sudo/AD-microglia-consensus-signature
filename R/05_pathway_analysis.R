# ============================================================
# 05_pathway_analysis.R
# ============================================================
# Project: AD Microglia Multi-Dataset Transcriptomic Analysis
#
# Purpose:
#   Rerun pathway analysis using the rebuilt AD61026 framework.
#
# Why this v2 script exists:
#   fgsea is now installed, but the KEGG subcollection name used in the
#   earlier script is not available in your current msigdbr version.
#   This script handles that gracefully:
#     - Hallmark is always attempted.
#     - GO BP is always attempted.
#     - KEGG is attempted using several possible modern/legacy names.
#     - If KEGG is unavailable, the script skips KEGG instead of failing.
#
# Important:
#   Run this script with source(). Do NOT paste it line by line.
#
# Recommended run command:
#   source("R/05_pathway_analysis_AD61026_FIELD_ALIGNED_v2.R")
#
# Main inputs:
#   Output/DE/SEAAD_DE_results.csv
#   Output/Manuscript_Source_Data/MainTable2_strict_consensus_signature_source.csv
#   Output/Manuscript_Source_Data/SuppTable6_expanded_drug_signature_source.csv
#
# Main outputs:
#   Output/Pathway/
#   Output/Manuscript_Source_Data/
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

required_pkgs <- c(
  "dplyr", "tidyr", "tibble", "stringr",
  "ggplot2", "fgsea", "msigdbr"
)

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall the missing package(s), restart R, and rerun this script with source()."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(fgsea)
  library(msigdbr)
})

cat("============================================================\n")
cat("  05: Pathway analysis for AD61026, field-aligned v2\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# Define output directories robustly.
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

if (!exists("FIG_DIR")) {
  FIG_DIR <- file.path(OUT_DIR, "Figures")
}
if (!exists("RDS_DIR")) {
  RDS_DIR <- file.path(PROJECT_ROOT, "RDS")
}

PATHWAY_DIR <- file.path(OUT_DIR, "Pathway")
SOURCE_DIR  <- file.path(OUT_DIR, "Manuscript_Source_Data")

dir.create(PATHWAY_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Pathway output directory:\n")
cat("  ", PATHWAY_DIR, "\n")
cat("Manuscript source-data directory:\n")
cat("  ", SOURCE_DIR, "\n\n")

# ============================================================
# SECTION 2: Helper functions
# ============================================================

read_required_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Required file missing: ", path)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

clean_pathway_label <- function(x) {
  x %>%
    gsub("^HALLMARK_", "", .) %>%
    gsub("^KEGG_", "", .) %>%
    gsub("^GOBP_", "", .) %>%
    gsub("^REACTOME_", "", .) %>%
    gsub("_", " ", .) %>%
    stringr::str_to_title()
}

# ------------------------------------------------------------
# Empty result tables. These prevent "object not found" errors if
# a collection such as KEGG is unavailable in the installed msigdbr.
# ------------------------------------------------------------
empty_gsea_table <- function(collection_label) {
  data.frame(
    collection = character(),
    pathway = character(),
    pathway_label = character(),
    direction = character(),
    NES = numeric(),
    pval = numeric(),
    padj = numeric(),
    size = integer(),
    leadingEdge = character(),
    stringsAsFactors = FALSE
  )
}

empty_ora_table <- function() {
  data.frame(
    signature = character(),
    query_direction = character(),
    collection = character(),
    pathway = character(),
    pathway_label = character(),
    query_size = integer(),
    pathway_size_in_background = integer(),
    overlap_n = integer(),
    overlap_genes = character(),
    odds_ratio = numeric(),
    p_value = numeric(),
    FDR = numeric(),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Try msigdbr using both modern and older argument names.
# This avoids failure across msigdbr versions.
# ------------------------------------------------------------
try_msigdbr <- function(collection, subcollection = NULL) {

  # First try the modern argument names.
  modern_try <- tryCatch(
    {
      if (is.null(subcollection)) {
        msigdbr::msigdbr(
          species = "Homo sapiens",
          collection = collection
        )
      } else {
        msigdbr::msigdbr(
          species = "Homo sapiens",
          collection = collection,
          subcollection = subcollection
        )
      }
    },
    error = function(e) NULL
  )

  if (!is.null(modern_try) && nrow(modern_try) > 0) {
    return(modern_try)
  }

  # Then try the older argument names.
  old_try <- tryCatch(
    {
      if (is.null(subcollection)) {
        msigdbr::msigdbr(
          species = "Homo sapiens",
          category = collection
        )
      } else {
        msigdbr::msigdbr(
          species = "Homo sapiens",
          category = collection,
          subcategory = subcollection
        )
      }
    },
    error = function(e) NULL
  )

  if (!is.null(old_try) && nrow(old_try) > 0) {
    return(old_try)
  }

  return(NULL)
}

# ------------------------------------------------------------
# Convert msigdbr data frame to named list for fgsea/ORA.
# ------------------------------------------------------------
msig_to_gene_sets <- function(msig_df) {

  if (is.null(msig_df) || nrow(msig_df) == 0) {
    return(list())
  }

  if (!all(c("gs_name", "gene_symbol") %in% colnames(msig_df))) {
    stop("msigdbr output does not contain expected columns gs_name and gene_symbol.")
  }

  gene_sets <- split(msig_df$gene_symbol, msig_df$gs_name)
  gene_sets <- lapply(gene_sets, unique)

  return(gene_sets)
}

# ------------------------------------------------------------
# Fetch a gene-set collection.
# ------------------------------------------------------------
fetch_gene_sets <- function(collection,
                            subcollection = NULL,
                            collection_label = collection) {

  msig_df <- try_msigdbr(collection, subcollection)

  gene_sets <- msig_to_gene_sets(msig_df)

  cat(collection_label, "gene sets loaded:", length(gene_sets), "\n")

  return(gene_sets)
}

# ------------------------------------------------------------
# Fetch KEGG robustly.
#
# In recent MSigDB/msigdbr versions, KEGG may be absent or renamed due
# to licensing/version changes. We try several possible names and skip
# KEGG if none are available.
# ------------------------------------------------------------
fetch_kegg_gene_sets <- function() {

  kegg_candidates <- c(
    "CP:KEGG",
    "CP:KEGG_LEGACY",
    "CP:KEGG_MEDICUS",
    "KEGG",
    "KEGG_LEGACY",
    "KEGG_MEDICUS"
  )

  for (candidate in kegg_candidates) {
    cat("Trying KEGG subcollection:", candidate, "\n")
    msig_df <- try_msigdbr("C2", candidate)

    gene_sets <- msig_to_gene_sets(msig_df)

    if (length(gene_sets) > 0) {
      cat("KEGG gene sets loaded using:", candidate, "\n")
      cat("KEGG gene sets loaded:", length(gene_sets), "\n")
      return(gene_sets)
    }
  }

  warning(
    "No KEGG subcollection was available in the installed msigdbr version. ",
    "The script will continue with Hallmark and GO BP only."
  )

  return(list())
}

# ------------------------------------------------------------
# Save available msigdbr collections for documentation/debugging.
# ------------------------------------------------------------
save_msig_collection_table <- function() {

  collections_tbl <- tryCatch(
    msigdbr::msigdbr_collections(),
    error = function(e) NULL
  )

  if (!is.null(collections_tbl)) {
    write.csv(
      collections_tbl,
      file.path(PATHWAY_DIR, "msigdbr_available_collections.csv"),
      row.names = FALSE
    )
    cat("Saved available msigdbr collections table.\n")
  }
}

# ------------------------------------------------------------
# Run preranked GSEA with fgsea.
# ------------------------------------------------------------
run_gsea <- function(ranks, gene_sets, collection_label) {

  if (length(gene_sets) == 0) {
    cat("Skipping GSEA for", collection_label, "because no gene sets were loaded.\n")
    return(empty_gsea_table(collection_label))
  }

  ranks <- ranks[is.finite(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)

  set.seed(42)

  fgsea_res <- fgsea::fgsea(
    pathways = gene_sets,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 0
  ) %>%
    as.data.frame()

  if (nrow(fgsea_res) == 0) {
    return(empty_gsea_table(collection_label))
  }

  fgsea_res <- fgsea_res %>%
    arrange(padj, desc(abs(NES))) %>%
    mutate(
      collection = collection_label,
      pathway_label = clean_pathway_label(pathway),
      direction = ifelse(NES > 0, "Enriched_in_AD", "Enriched_in_Control"),
      leadingEdge = sapply(leadingEdge, paste, collapse = ";")
    ) %>%
    select(
      collection,
      pathway,
      pathway_label,
      direction,
      NES,
      pval,
      padj,
      size,
      leadingEdge
    )

  return(fgsea_res)
}

# ------------------------------------------------------------
# Run over-representation analysis using Fisher exact test.
# ------------------------------------------------------------
run_ora <- function(query_genes,
                    background_genes,
                    gene_sets,
                    collection_label,
                    signature_label,
                    direction_label) {

  if (length(gene_sets) == 0) {
    cat("Skipping ORA for", collection_label, "because no gene sets were loaded.\n")
    return(empty_ora_table())
  }

  query_genes <- unique(query_genes)
  background_genes <- unique(background_genes)
  query_genes <- intersect(query_genes, background_genes)

  if (length(query_genes) < 3) {
    warning(
      "Very small ORA query for ",
      signature_label, " / ", direction_label,
      ": fewer than 3 genes after background filtering."
    )
  }

  ora_list <- lapply(names(gene_sets), function(pathway_name) {

    pathway_genes <- intersect(unique(gene_sets[[pathway_name]]), background_genes)

    if (length(pathway_genes) < 5 || length(pathway_genes) > 1000) {
      return(NULL)
    }

    overlap_genes <- intersect(query_genes, pathway_genes)

    a <- length(overlap_genes)
    b <- length(query_genes) - a
    c <- length(pathway_genes) - a
    d <- length(background_genes) - a - b - c

    mat <- matrix(c(a, b, c, d), nrow = 2)
    fisher <- fisher.test(mat, alternative = "greater")

    data.frame(
      signature = signature_label,
      query_direction = direction_label,
      collection = collection_label,
      pathway = pathway_name,
      pathway_label = clean_pathway_label(pathway_name),
      query_size = length(query_genes),
      pathway_size_in_background = length(pathway_genes),
      overlap_n = a,
      overlap_genes = paste(sort(overlap_genes), collapse = ";"),
      odds_ratio = unname(fisher$estimate),
      p_value = fisher$p.value,
      stringsAsFactors = FALSE
    )
  })

  ora <- bind_rows(ora_list)

  if (nrow(ora) == 0) {
    return(empty_ora_table())
  }

  ora <- ora %>%
    mutate(FDR = p.adjust(p_value, method = "BH")) %>%
    arrange(FDR, p_value, desc(overlap_n))

  return(ora)
}

# ------------------------------------------------------------
# Count significant rows safely.
# ------------------------------------------------------------
count_fdr_sig <- function(df, fdr_col) {
  if (is.null(df) || nrow(df) == 0 || !fdr_col %in% colnames(df)) {
    return(0)
  }
  sum(df[[fdr_col]] < 0.05, na.rm = TRUE)
}

# ------------------------------------------------------------
# Write source data and track it in a manifest.
# ------------------------------------------------------------
source_manifest <- list()

write_source <- function(df,
                         file_name,
                         intended_use,
                         source_inputs,
                         notes = "") {

  out_path <- file.path(SOURCE_DIR, file_name)
  write.csv(df, out_path, row.names = FALSE)

  source_manifest[[length(source_manifest) + 1]] <<- data.frame(
    file_name = file_name,
    intended_use = intended_use,
    source_inputs = source_inputs,
    notes = notes,
    stringsAsFactors = FALSE
  )

  cat("Saved source:", out_path, "\n")
}

# ============================================================
# SECTION 3: Load current signature and DE files
# ============================================================

cat("Loading DE and signature source files...\n\n")

seaad_de <- read_required_csv(
  file.path(DE_DIR, "SEAAD_DE_results.csv")
)

strict_signature <- read_required_csv(
  file.path(SOURCE_DIR, "MainTable2_strict_consensus_signature_source.csv")
)

expanded_signature <- read_required_csv(
  file.path(SOURCE_DIR, "SuppTable6_expanded_drug_signature_source.csv")
)

required_de_cols <- c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")
missing_de_cols <- setdiff(required_de_cols, colnames(seaad_de))

if (length(missing_de_cols) > 0) {
  stop(
    "SEAAD_DE_results.csv is missing required columns: ",
    paste(missing_de_cols, collapse = ", ")
  )
}

required_sig_cols <- c("gene", "discovery_direction")

missing_strict_cols <- setdiff(required_sig_cols, colnames(strict_signature))
missing_expanded_cols <- setdiff(required_sig_cols, colnames(expanded_signature))

if (length(missing_strict_cols) > 0) {
  stop(
    "Strict signature file is missing required columns: ",
    paste(missing_strict_cols, collapse = ", ")
  )
}

if (length(missing_expanded_cols) > 0) {
  stop(
    "Expanded signature file is missing required columns: ",
    paste(missing_expanded_cols, collapse = ", ")
  )
}

cat("SEA-AD DE genes loaded:", nrow(seaad_de), "\n")
cat("Strict consensus genes loaded:", nrow(strict_signature), "\n")
cat("Expanded signature genes loaded:", nrow(expanded_signature), "\n\n")

# ============================================================
# SECTION 4: Prepare ranked SEA-AD gene list for GSEA
# ============================================================

cat("Preparing SEA-AD ranked gene list using limma moderated t-statistic...\n\n")

rank_df <- seaad_de %>%
  filter(!is.na(gene),
         !is.na(t),
         is.finite(t)) %>%
  group_by(gene) %>%
  arrange(desc(abs(t)), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

ranks <- rank_df$t
names(ranks) <- rank_df$gene
ranks <- sort(ranks, decreasing = TRUE)

background_genes <- unique(rank_df$gene)

cat("Ranked genes for GSEA:", length(ranks), "\n")
cat("Background genes for ORA:", length(background_genes), "\n\n")

ranked_source <- rank_df %>%
  select(gene, logFC, AveExpr, t, P.Value, adj.P.Val) %>%
  arrange(desc(t))

write.csv(
  ranked_source,
  file.path(PATHWAY_DIR, "SEAAD_ranked_gene_list_by_moderated_t.csv"),
  row.names = FALSE
)

# ============================================================
# SECTION 5: Load MSigDB gene sets
# ============================================================

cat("Loading MSigDB gene sets through msigdbr...\n\n")

save_msig_collection_table()

hallmark_sets <- fetch_gene_sets(
  collection = "H",
  subcollection = NULL,
  collection_label = "Hallmark"
)

kegg_sets <- fetch_kegg_gene_sets()

gobp_sets <- fetch_gene_sets(
  collection = "C5",
  subcollection = "GO:BP",
  collection_label = "GO_BP"
)

# ============================================================
# SECTION 6: Run SEA-AD preranked GSEA
# ============================================================

cat("\nRunning SEA-AD preranked GSEA...\n\n")

gsea_hallmark <- run_gsea(ranks, hallmark_sets, "Hallmark")
gsea_kegg     <- run_gsea(ranks, kegg_sets, "KEGG")
gsea_gobp     <- run_gsea(ranks, gobp_sets, "GO_BP")

gsea_all <- bind_rows(gsea_hallmark, gsea_kegg, gsea_gobp) %>%
  arrange(padj, desc(abs(NES)))

write.csv(gsea_hallmark, file.path(PATHWAY_DIR, "GSEA_SEAAD_Hallmark.csv"), row.names = FALSE)
write.csv(gsea_kegg,     file.path(PATHWAY_DIR, "GSEA_SEAAD_KEGG.csv"), row.names = FALSE)
write.csv(gsea_gobp,     file.path(PATHWAY_DIR, "GSEA_SEAAD_GO_BP.csv"), row.names = FALSE)
write.csv(gsea_all,      file.path(PATHWAY_DIR, "GSEA_SEAAD_all_collections.csv"), row.names = FALSE)

saveRDS(
  list(
    hallmark = gsea_hallmark,
    kegg = gsea_kegg,
    gobp = gsea_gobp,
    all = gsea_all
  ),
  file.path(RDS_DIR, "pathway_GSEA_SEAAD_AD61026.rds")
)

cat("GSEA complete.\n")
cat("  Hallmark FDR<0.05:", count_fdr_sig(gsea_hallmark, "padj"), "\n")
cat("  KEGG FDR<0.05:", count_fdr_sig(gsea_kegg, "padj"), "\n")
cat("  GO BP FDR<0.05:", count_fdr_sig(gsea_gobp, "padj"), "\n\n")

# ============================================================
# SECTION 7: Run ORA for strict and expanded signatures
# ============================================================

cat("Running ORA for strict and expanded signatures...\n\n")

strict_all <- unique(strict_signature$gene)

strict_up <- strict_signature %>%
  filter(discovery_direction == "UP") %>%
  pull(gene) %>%
  unique()

strict_down <- strict_signature %>%
  filter(discovery_direction == "DOWN") %>%
  pull(gene) %>%
  unique()

expanded_all <- unique(expanded_signature$gene)

expanded_up <- expanded_signature %>%
  filter(discovery_direction == "UP") %>%
  pull(gene) %>%
  unique()

expanded_down <- expanded_signature %>%
  filter(discovery_direction == "DOWN") %>%
  pull(gene) %>%
  unique()

cat("Strict signature genes:\n")
cat("  all:", length(strict_all), " up:", length(strict_up), " down:", length(strict_down), "\n")
cat("Expanded signature genes:\n")
cat("  all:", length(expanded_all), " up:", length(expanded_up), " down:", length(expanded_down), "\n\n")

gene_set_registry <- list(
  Hallmark = hallmark_sets,
  KEGG = kegg_sets,
  GO_BP = gobp_sets
)

run_ora_across_collections <- function(query_genes,
                                       signature_label,
                                       direction_label) {

  bind_rows(lapply(names(gene_set_registry), function(coll) {
    run_ora(
      query_genes = query_genes,
      background_genes = background_genes,
      gene_sets = gene_set_registry[[coll]],
      collection_label = coll,
      signature_label = signature_label,
      direction_label = direction_label
    )
  }))
}

ora_strict_all <- run_ora_across_collections(
  strict_all,
  "strict_primary_consensus_48",
  "all"
)

ora_strict_by_direction <- bind_rows(
  run_ora_across_collections(strict_up, "strict_primary_consensus_48", "UP_in_AD"),
  run_ora_across_collections(strict_down, "strict_primary_consensus_48", "DOWN_in_AD")
)

ora_expanded_all <- run_ora_across_collections(
  expanded_all,
  "expanded_high_plus_moderate_125",
  "all"
)

ora_expanded_by_direction <- bind_rows(
  run_ora_across_collections(expanded_up, "expanded_high_plus_moderate_125", "UP_in_AD"),
  run_ora_across_collections(expanded_down, "expanded_high_plus_moderate_125", "DOWN_in_AD")
)

write.csv(ora_strict_all, file.path(PATHWAY_DIR, "ORA_strict_consensus_all.csv"), row.names = FALSE)
write.csv(ora_strict_by_direction, file.path(PATHWAY_DIR, "ORA_strict_consensus_by_direction.csv"), row.names = FALSE)
write.csv(ora_expanded_all, file.path(PATHWAY_DIR, "ORA_expanded_signature_all.csv"), row.names = FALSE)
write.csv(ora_expanded_by_direction, file.path(PATHWAY_DIR, "ORA_expanded_signature_by_direction.csv"), row.names = FALSE)

saveRDS(
  list(
    strict_all = ora_strict_all,
    strict_by_direction = ora_strict_by_direction,
    expanded_all = ora_expanded_all,
    expanded_by_direction = ora_expanded_by_direction
  ),
  file.path(RDS_DIR, "pathway_ORA_AD61026.rds")
)

cat("ORA complete.\n")
cat("  Strict all FDR<0.05:", count_fdr_sig(ora_strict_all, "FDR"), "\n")
cat("  Strict by direction FDR<0.05:", count_fdr_sig(ora_strict_by_direction, "FDR"), "\n")
cat("  Expanded all FDR<0.05:", count_fdr_sig(ora_expanded_all, "FDR"), "\n")
cat("  Expanded by direction FDR<0.05:", count_fdr_sig(ora_expanded_by_direction, "FDR"), "\n\n")

# ============================================================
# SECTION 8: Build concise pathway summary table
# ============================================================

cat("Building concise pathway summary table...\n\n")

top_gsea <- gsea_all %>%
  filter(padj < 0.05) %>%
  arrange(padj, desc(abs(NES))) %>%
  group_by(collection) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  transmute(
    analysis = "SEAAD_GSEA_moderated_t",
    signature = "SEAAD_full_ranked_DE",
    query_direction = direction,
    collection,
    pathway,
    pathway_label,
    statistic = NES,
    p_value = pval,
    FDR = padj,
    gene_count = size,
    genes = leadingEdge
  )

top_ora_strict <- bind_rows(ora_strict_all, ora_strict_by_direction) %>%
  filter(FDR < 0.25) %>%
  arrange(FDR, p_value, desc(overlap_n)) %>%
  group_by(signature, query_direction, collection) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  transmute(
    analysis = "ORA_strict_consensus",
    signature,
    query_direction,
    collection,
    pathway,
    pathway_label,
    statistic = odds_ratio,
    p_value,
    FDR,
    gene_count = overlap_n,
    genes = overlap_genes
  )

top_ora_expanded <- bind_rows(ora_expanded_all, ora_expanded_by_direction) %>%
  filter(FDR < 0.25) %>%
  arrange(FDR, p_value, desc(overlap_n)) %>%
  group_by(signature, query_direction, collection) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  transmute(
    analysis = "ORA_expanded_signature",
    signature,
    query_direction,
    collection,
    pathway,
    pathway_label,
    statistic = odds_ratio,
    p_value,
    FDR,
    gene_count = overlap_n,
    genes = overlap_genes
  )

pathway_summary <- bind_rows(top_gsea, top_ora_strict, top_ora_expanded) %>%
  arrange(analysis, FDR, desc(gene_count))

write.csv(
  pathway_summary,
  file.path(PATHWAY_DIR, "pathway_summary_top_terms.csv"),
  row.names = FALSE
)

# ============================================================
# SECTION 9: Create draft pathway review figures
# ============================================================

cat("Creating draft pathway review figures...\n\n")

figure_gsea_source <- gsea_all %>%
  filter(padj < 0.05) %>%
  arrange(padj, desc(abs(NES))) %>%
  group_by(collection) %>%
  slice_head(n = 12) %>%
  ungroup() %>%
  mutate(
    pathway_label = factor(pathway_label, levels = rev(unique(pathway_label))),
    minus_log10_FDR = -log10(padj)
  )

if (nrow(figure_gsea_source) > 0) {

  p_gsea <- ggplot(
    figure_gsea_source,
    aes(x = NES, y = pathway_label, size = minus_log10_FDR)
  ) +
    geom_point(alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    facet_grid(collection ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = "SEA-AD microglia pathway remodeling",
      subtitle = "Preranked GSEA using limma moderated t-statistic",
      x = "Normalized enrichment score",
      y = NULL,
      size = expression(-log[10]~FDR)
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    )

  ggsave(
    file.path(FIG_DIR, "Fig_GSEA_SEAAD_moderated_t_draft.pdf"),
    p_gsea,
    width = 9,
    height = 8
  )
}

figure_ora_strict_source <- bind_rows(ora_strict_all, ora_strict_by_direction) %>%
  filter(FDR < 0.25) %>%
  arrange(FDR, p_value, desc(overlap_n)) %>%
  slice_head(n = 25) %>%
  mutate(
    pathway_label = factor(pathway_label, levels = rev(unique(pathway_label))),
    minus_log10_FDR = -log10(FDR)
  )

if (nrow(figure_ora_strict_source) > 0) {

  p_ora_strict <- ggplot(
    figure_ora_strict_source,
    aes(x = minus_log10_FDR, y = pathway_label)
  ) +
    geom_col() +
    facet_grid(collection ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = "Pathway over-representation in strict consensus genes",
      subtitle = "48-gene strict primary consensus signature",
      x = expression(-log[10]~FDR),
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    )

  ggsave(
    file.path(FIG_DIR, "Fig_ORA_strict_consensus_draft.pdf"),
    p_ora_strict,
    width = 8,
    height = 7
  )
}

# ============================================================
# SECTION 10: Save manuscript source-data files
# ============================================================

cat("Saving pathway source-data files for manuscript tables/figures...\n\n")

write_source(
  pathway_summary,
  "MainTable3_pathway_summary_source.csv",
  intended_use = "Main Table 3 or pathway results table.",
  source_inputs = "SEAAD_DE_results.csv; strict and expanded signature source files; msigdbr gene sets.",
  notes = "Concise top pathway table. Full GSEA and ORA results are saved separately."
)

write_source(
  gsea_all,
  "SuppTable9_GSEA_full_source.csv",
  intended_use = "Supplementary Table 9: full preranked GSEA results.",
  source_inputs = "SEAAD_DE_results.csv ranked by limma moderated t-statistic.",
  notes = "Positive NES indicates enrichment among genes higher in AD. KEGG may be empty if unavailable in msigdbr."
)

write_source(
  bind_rows(ora_strict_all, ora_strict_by_direction),
  "SuppTable10_ORA_strict_consensus_source.csv",
  intended_use = "Supplementary Table 10: ORA for 48-gene strict primary consensus signature.",
  source_inputs = "MainTable2_strict_consensus_signature_source.csv.",
  notes = "Uses SEA-AD tested genes as background."
)

write_source(
  bind_rows(ora_expanded_all, ora_expanded_by_direction),
  "SuppTable11_ORA_expanded_signature_source.csv",
  intended_use = "Supplementary Table 11: ORA for 125-gene expanded high+moderate signature.",
  source_inputs = "SuppTable6_expanded_drug_signature_source.csv.",
  notes = "Supportive/sensitivity analysis; useful before drug prioritization."
)

write_source(
  figure_gsea_source,
  "Figure6_GSEA_dotplot_source.csv",
  intended_use = "Figure source: SEA-AD GSEA dot plot.",
  source_inputs = "GSEA_SEAAD_all_collections.csv.",
  notes = "Draft plot saved as Fig_GSEA_SEAAD_moderated_t_draft.pdf."
)

write_source(
  figure_ora_strict_source,
  "Figure7_ORA_strict_source.csv",
  intended_use = "Figure source: strict consensus ORA bar plot.",
  source_inputs = "ORA_strict_consensus_all.csv; ORA_strict_consensus_by_direction.csv.",
  notes = "Draft plot saved as Fig_ORA_strict_consensus_draft.pdf."
)

# ============================================================
# SECTION 11: Run summary and manifest
# ============================================================

pathway_run_summary <- data.frame(
  item = c(
    "SEA-AD ranked genes",
    "Strict consensus genes",
    "Strict consensus UP genes",
    "Strict consensus DOWN genes",
    "Expanded signature genes",
    "Expanded signature UP genes",
    "Expanded signature DOWN genes",
    "Hallmark gene sets loaded",
    "KEGG gene sets loaded",
    "GO BP gene sets loaded",
    "Hallmark GSEA FDR<0.05",
    "KEGG GSEA FDR<0.05",
    "GO BP GSEA FDR<0.05",
    "Strict ORA all FDR<0.05",
    "Strict ORA by direction FDR<0.05",
    "Expanded ORA all FDR<0.05",
    "Expanded ORA by direction FDR<0.05"
  ),
  value = c(
    length(ranks),
    length(strict_all),
    length(strict_up),
    length(strict_down),
    length(expanded_all),
    length(expanded_up),
    length(expanded_down),
    length(hallmark_sets),
    length(kegg_sets),
    length(gobp_sets),
    count_fdr_sig(gsea_hallmark, "padj"),
    count_fdr_sig(gsea_kegg, "padj"),
    count_fdr_sig(gsea_gobp, "padj"),
    count_fdr_sig(ora_strict_all, "FDR"),
    count_fdr_sig(ora_strict_by_direction, "FDR"),
    count_fdr_sig(ora_expanded_all, "FDR"),
    count_fdr_sig(ora_expanded_by_direction, "FDR")
  ),
  stringsAsFactors = FALSE
)

write.csv(
  pathway_run_summary,
  file.path(PATHWAY_DIR, "pathway_run_summary.csv"),
  row.names = FALSE
)

source_manifest_df <- bind_rows(source_manifest)

write.csv(
  source_manifest_df,
  file.path(SOURCE_DIR, "README_pathway_source_manifest.csv"),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("05 pathway analysis complete.\n\n")
cat("Key outputs:\n")
cat("  ", file.path(PATHWAY_DIR, "pathway_run_summary.csv"), "\n")
cat("  ", file.path(PATHWAY_DIR, "pathway_summary_top_terms.csv"), "\n")
cat("  ", file.path(SOURCE_DIR, "MainTable3_pathway_summary_source.csv"), "\n")
cat("  ", file.path(SOURCE_DIR, "SuppTable9_GSEA_full_source.csv"), "\n")
cat("  ", file.path(SOURCE_DIR, "SuppTable10_ORA_strict_consensus_source.csv"), "\n")
cat("  ", file.path(SOURCE_DIR, "SuppTable11_ORA_expanded_signature_source.csv"), "\n\n")
cat("Next step:\n")
cat("  Upload pathway_run_summary.csv and pathway_summary_top_terms.csv.\n")
cat("============================================================\n")
