# ============================================================
# 15a_drug_query_inputs.R
# ============================================================
# Project: AD Microglia Multi-Dataset Transcriptomic Analysis
#
# Purpose:
#   Prepare drug-prioritization query files for the rebuilt AD61026
#   analysis and run an optional DGIdb target-level query.
#
# Why this script comes before final drug ranking:
#   DREIMT, CMap, and iLINCS usually require web/manual submission.
#   This script creates clean, versioned input files for those tools.
#   It also queries DGIdb programmatically when internet access works.
#
# Current signature hierarchy:
#   Primary drug-query signature:
#      strict 48-gene consensus signature
#
#   Sensitivity drug-query signature:
#      expanded 125-gene high + moderate signature
#      plus a non-mitochondrial version because the expanded set
#      includes MT genes.
#
# Important interpretation:
#   - The strict 48-gene signature is the primary drug query.
#   - The expanded 125-gene signature is supportive/sensitivity only.
#   - The expanded non-MT signature is useful because the MT sensitivity
#     analysis showed mitochondrial genes are strongly attenuated after
#     MT-fraction adjustment.
#   - Do not hard-code old DREIMT/CMap/iLINCS results from AD6626 or
#     older 75-gene frameworks.
#
# Expected inputs:
#   Output/Manuscript_Source_Data/MainTable2_strict_consensus_signature_source.csv
#   Output/Manuscript_Source_Data/SuppTable6_expanded_drug_signature_source.csv
#   Output/Manuscript_Source_Data/SuppTable12_mito_sensitivity_source_v2.csv
#      optional but recommended
#
# Main outputs:
#   Output/Drug_Prioritization/
#      drug_query_gene_list_manifest_AD61026.csv
#      drug_query_genes_long_AD61026.csv
#      strict48_up_genes.txt
#      strict48_down_genes.txt
#      expanded125_up_genes.txt
#      expanded125_down_genes.txt
#      expanded125_nonMT_up_genes.txt
#      expanded125_nonMT_down_genes.txt
#      dgidb_strict48_interactions.csv
#      dgidb_expanded125_nonMT_interactions.csv
#
# Recommended run command:
#   source("R/06A_prepare_drug_query_inputs_AD61026.R")
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
  "dplyr", "tibble", "tidyr", "stringr",
  "ggplot2", "httr", "jsonlite"
)

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall the missing package(s), restart R, and rerun this script."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(httr)
  library(jsonlite)
})

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

if (!exists("RDS_DIR")) {
  RDS_DIR <- file.path(PROJECT_ROOT, "RDS")
}
if (!exists("FIG_DIR")) {
  FIG_DIR <- file.path(OUT_DIR, "Figures")
}

SOURCE_DIR <- file.path(OUT_DIR, "Manuscript_Source_Data")
DRUG_DIR   <- file.path(OUT_DIR, "Drug_Prioritization")

dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DRUG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("  06A: Prepare AD61026 drug-prioritization inputs\n")
cat("============================================================\n\n")

cat("Drug-prioritization output directory:\n")
cat("  ", DRUG_DIR, "\n\n")

# ============================================================
# SECTION 2: Helper functions
# ============================================================

read_required_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Required file missing: ", path)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_gene_list <- function(genes, path) {
  genes <- genes[!is.na(genes) & genes != ""]
  genes <- unique(genes)
  writeLines(genes, path)
  cat("Saved", length(genes), "genes to", path, "\n")
}

standardize_direction <- function(x) {
  x <- as.character(x)
  case_when(
    x %in% c("UP", "Up", "up", "UP_in_AD", "AD_upregulated", "AD-upregulated") ~ "UP",
    x %in% c("DOWN", "Down", "down", "DOWN_in_AD", "AD_downregulated", "AD-downregulated") ~ "DOWN",
    TRUE ~ NA_character_
  )
}

is_mito_gene <- function(gene) {
  grepl("^MT-", gene)
}

# ------------------------------------------------------------
# Very simple symbol-category flag.
# This does not replace platform recognition results from CMap/iLINCS/
# DREIMT, but helps us anticipate likely database coverage issues.
# ------------------------------------------------------------
classify_gene_symbol <- function(gene) {
  case_when(
    grepl("^MT-", gene) ~ "mitochondrial_encoded",
    grepl("^LINC", gene) ~ "lncRNA_like",
    grepl("-AS[0-9]*$", gene) ~ "antisense_lncRNA_like",
    grepl("^AC[0-9]", gene) ~ "uncharacterized_locus_like",
    grepl("^AL[0-9]", gene) ~ "uncharacterized_locus_like",
    grepl("\\.", gene) ~ "uncharacterized_locus_like",
    TRUE ~ "protein_coding_or_standard_symbol_like"
  )
}

# ============================================================
# SECTION 3: Load current strict and expanded signatures
# ============================================================

cat("Loading strict and expanded signatures...\n\n")

strict_path <- file.path(
  SOURCE_DIR,
  "MainTable2_strict_consensus_signature_source.csv"
)

expanded_path <- file.path(
  SOURCE_DIR,
  "SuppTable6_expanded_drug_signature_source.csv"
)

strict_signature <- read_required_csv(strict_path)
expanded_signature <- read_required_csv(expanded_path)

required_cols <- c("gene", "discovery_direction")

missing_strict <- setdiff(required_cols, colnames(strict_signature))
missing_expanded <- setdiff(required_cols, colnames(expanded_signature))

if (length(missing_strict) > 0) {
  stop("Strict signature is missing columns: ",
       paste(missing_strict, collapse = ", "))
}

if (length(missing_expanded) > 0) {
  stop("Expanded signature is missing columns: ",
       paste(missing_expanded, collapse = ", "))
}

# ------------------------------------------------------------
# Standardize and annotate signatures.
# ------------------------------------------------------------
strict_signature <- strict_signature %>%
  mutate(
    signature = "strict_primary_consensus_48",
    direction = standardize_direction(discovery_direction),
    mt_gene = is_mito_gene(gene),
    symbol_class = classify_gene_symbol(gene)
  )

expanded_signature <- expanded_signature %>%
  mutate(
    signature = "expanded_high_plus_moderate_125",
    direction = standardize_direction(discovery_direction),
    mt_gene = is_mito_gene(gene),
    symbol_class = classify_gene_symbol(gene)
  )

if (any(is.na(strict_signature$direction))) {
  stop("Some strict genes have unrecognized discovery_direction values.")
}

if (any(is.na(expanded_signature$direction))) {
  stop("Some expanded genes have unrecognized discovery_direction values.")
}

cat("Strict signature genes:", nrow(strict_signature), "\n")
cat("  UP:", sum(strict_signature$direction == "UP"), "\n")
cat("  DOWN:", sum(strict_signature$direction == "DOWN"), "\n")
cat("  MT genes:", sum(strict_signature$mt_gene), "\n\n")

cat("Expanded signature genes:", nrow(expanded_signature), "\n")
cat("  UP:", sum(expanded_signature$direction == "UP"), "\n")
cat("  DOWN:", sum(expanded_signature$direction == "DOWN"), "\n")
cat("  MT genes:", sum(expanded_signature$mt_gene), "\n\n")

# ============================================================
# SECTION 4: Add MT sensitivity information if available
# ============================================================

cat("Checking for MT sensitivity source table...\n\n")

mito_source_path <- file.path(
  SOURCE_DIR,
  "SuppTable12_mito_sensitivity_source_v2.csv"
)

if (file.exists(mito_source_path)) {

  mito_source <- read.csv(
    mito_source_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  mito_cols_to_keep <- intersect(
    c(
      "gene",
      "logFC_orig",
      "FDR_orig",
      "logFC_adj",
      "FDR_adj",
      "same_direction",
      "retained_after_mito_adjustment_FDR",
      "retained_after_mito_adjustment_nominal",
      "percent_attenuation"
    ),
    colnames(mito_source)
  )

  mito_source_small <- mito_source %>%
    select(all_of(mito_cols_to_keep))

  strict_signature <- strict_signature %>%
    left_join(mito_source_small, by = "gene")

  expanded_signature <- expanded_signature %>%
    left_join(mito_source_small, by = "gene")

  cat("MT sensitivity annotations added.\n\n")

} else {
  cat("MT sensitivity source table not found. Continuing without MT annotations.\n\n")
}

# ============================================================
# SECTION 5: Create query lists
# ============================================================

cat("Creating query gene lists...\n\n")

strict_up <- strict_signature %>%
  filter(direction == "UP") %>%
  pull(gene) %>%
  unique()

strict_down <- strict_signature %>%
  filter(direction == "DOWN") %>%
  pull(gene) %>%
  unique()

expanded_up <- expanded_signature %>%
  filter(direction == "UP") %>%
  pull(gene) %>%
  unique()

expanded_down <- expanded_signature %>%
  filter(direction == "DOWN") %>%
  pull(gene) %>%
  unique()

expanded_nonMT_up <- expanded_signature %>%
  filter(direction == "UP", !mt_gene) %>%
  pull(gene) %>%
  unique()

expanded_nonMT_down <- expanded_signature %>%
  filter(direction == "DOWN", !mt_gene) %>%
  pull(gene) %>%
  unique()

# Standard-symbol-like versions are useful for tools that reject lncRNAs
# or uncharacterized loci.
strict_standard_up <- strict_signature %>%
  filter(direction == "UP",
         symbol_class == "protein_coding_or_standard_symbol_like") %>%
  pull(gene) %>%
  unique()

strict_standard_down <- strict_signature %>%
  filter(direction == "DOWN",
         symbol_class == "protein_coding_or_standard_symbol_like") %>%
  pull(gene) %>%
  unique()

expanded_nonMT_standard_up <- expanded_signature %>%
  filter(direction == "UP",
         !mt_gene,
         symbol_class == "protein_coding_or_standard_symbol_like") %>%
  pull(gene) %>%
  unique()

expanded_nonMT_standard_down <- expanded_signature %>%
  filter(direction == "DOWN",
         !mt_gene,
         symbol_class == "protein_coding_or_standard_symbol_like") %>%
  pull(gene) %>%
  unique()

# ============================================================
# SECTION 6: Export text files for DREIMT, CMap, and iLINCS
# ============================================================

cat("Exporting text files for external drug-query tools...\n\n")

write_gene_list(strict_up, file.path(DRUG_DIR, "strict48_up_genes.txt"))
write_gene_list(strict_down, file.path(DRUG_DIR, "strict48_down_genes.txt"))

write_gene_list(expanded_up, file.path(DRUG_DIR, "expanded125_up_genes.txt"))
write_gene_list(expanded_down, file.path(DRUG_DIR, "expanded125_down_genes.txt"))

write_gene_list(expanded_nonMT_up, file.path(DRUG_DIR, "expanded125_nonMT_up_genes.txt"))
write_gene_list(expanded_nonMT_down, file.path(DRUG_DIR, "expanded125_nonMT_down_genes.txt"))

write_gene_list(strict_standard_up, file.path(DRUG_DIR, "strict48_standard_symbol_up_genes.txt"))
write_gene_list(strict_standard_down, file.path(DRUG_DIR, "strict48_standard_symbol_down_genes.txt"))

write_gene_list(expanded_nonMT_standard_up, file.path(DRUG_DIR, "expanded125_nonMT_standard_symbol_up_genes.txt"))
write_gene_list(expanded_nonMT_standard_down, file.path(DRUG_DIR, "expanded125_nonMT_standard_symbol_down_genes.txt"))

# ------------------------------------------------------------
# Tool-specific copies. These are intentionally duplicated so it is
# obvious which files were submitted to each tool.
# ------------------------------------------------------------
write_gene_list(strict_up, file.path(DRUG_DIR, "DREIMT_PRIMARY_strict48_up.txt"))
write_gene_list(strict_down, file.path(DRUG_DIR, "DREIMT_PRIMARY_strict48_down.txt"))

write_gene_list(expanded_nonMT_up, file.path(DRUG_DIR, "DREIMT_SENSITIVITY_expanded125_nonMT_up.txt"))
write_gene_list(expanded_nonMT_down, file.path(DRUG_DIR, "DREIMT_SENSITIVITY_expanded125_nonMT_down.txt"))

write_gene_list(strict_up, file.path(DRUG_DIR, "CMap_PRIMARY_strict48_up.txt"))
write_gene_list(strict_down, file.path(DRUG_DIR, "CMap_PRIMARY_strict48_down.txt"))

write_gene_list(expanded_nonMT_up, file.path(DRUG_DIR, "CMap_SENSITIVITY_expanded125_nonMT_up.txt"))
write_gene_list(expanded_nonMT_down, file.path(DRUG_DIR, "CMap_SENSITIVITY_expanded125_nonMT_down.txt"))

write_gene_list(strict_up, file.path(DRUG_DIR, "iLINCS_PRIMARY_strict48_up.txt"))
write_gene_list(strict_down, file.path(DRUG_DIR, "iLINCS_PRIMARY_strict48_down.txt"))

write_gene_list(expanded_nonMT_up, file.path(DRUG_DIR, "iLINCS_SENSITIVITY_expanded125_nonMT_up.txt"))
write_gene_list(expanded_nonMT_down, file.path(DRUG_DIR, "iLINCS_SENSITIVITY_expanded125_nonMT_down.txt"))

# ============================================================
# SECTION 7: Save query manifest and long gene table
# ============================================================

query_manifest <- data.frame(
  query_name = c(
    "strict48_primary",
    "expanded125_sensitivity",
    "expanded125_nonMT_sensitivity",
    "strict48_standard_symbol_primary",
    "expanded125_nonMT_standard_symbol_sensitivity"
  ),
  up_n = c(
    length(strict_up),
    length(expanded_up),
    length(expanded_nonMT_up),
    length(strict_standard_up),
    length(expanded_nonMT_standard_up)
  ),
  down_n = c(
    length(strict_down),
    length(expanded_down),
    length(expanded_nonMT_down),
    length(strict_standard_down),
    length(expanded_nonMT_standard_down)
  ),
  total_n = c(
    length(c(strict_up, strict_down)),
    length(c(expanded_up, expanded_down)),
    length(c(expanded_nonMT_up, expanded_nonMT_down)),
    length(c(strict_standard_up, strict_standard_down)),
    length(c(expanded_nonMT_standard_up, expanded_nonMT_standard_down))
  ),
  intended_use = c(
    "Primary drug-prioritization query.",
    "Sensitivity query, includes MT genes.",
    "Sensitivity query excluding MT genes.",
    "Primary query subset likely to have better platform recognition.",
    "Sensitivity non-MT subset likely to have better platform recognition."
  ),
  notes = c(
    "Use as primary CMap/iLINCS/DREIMT query when accepted.",
    "Useful to assess stability but interpret cautiously due MT genes.",
    "Preferred sensitivity query if MT genes dominate output.",
    "Do not replace primary unless platform recognition is poor; report recognized gene count.",
    "Use if DREIMT/CMap/iLINCS rejects many nonstandard symbols."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  query_manifest,
  file.path(DRUG_DIR, "drug_query_gene_list_manifest_AD61026.csv"),
  row.names = FALSE
)

drug_query_genes_long <- bind_rows(
  strict_signature,
  expanded_signature
) %>%
  distinct(signature, gene, .keep_all = TRUE) %>%
  mutate(
    recommended_primary_query = signature == "strict_primary_consensus_48",
    recommended_sensitivity_query = signature == "expanded_high_plus_moderate_125",
    recommended_sensitivity_nonMT_query =
      signature == "expanded_high_plus_moderate_125" & !mt_gene,
    dgidb_reversal_rule = case_when(
      direction == "UP" ~ "Prioritize inhibitors, antagonists, blockers, suppressors",
      direction == "DOWN" ~ "Prioritize activators, agonists, stimulators",
      TRUE ~ NA_character_
    )
  )

write.csv(
  drug_query_genes_long,
  file.path(DRUG_DIR, "drug_query_genes_long_AD61026.csv"),
  row.names = FALSE
)

cat("\nQuery manifest:\n")
print(query_manifest)
cat("\n")

# ============================================================
# SECTION 8: DGIdb query helper
# ============================================================
# DGIdb v5 uses GraphQL at:
#   https://dgidb.org/api/graphql
#
# If the local network blocks this API, the script saves all query
# files and continues. You can rerun this section later or use saved
# files for manual querying.

query_dgidb <- function(gene_table,
                        query_label,
                        output_prefix) {

  genes <- unique(gene_table$gene)
  genes <- genes[!is.na(genes) & genes != ""]

  if (length(genes) == 0) {
    warning("No genes supplied for DGIdb query: ", query_label)
    return(list(
      interactions = data.frame(),
      gene_summary = data.frame(),
      drug_summary = data.frame(),
      relevant = data.frame()
    ))
  }

  cat("------------------------------------------------------------\n")
  cat("DGIdb query:", query_label, "\n")
  cat("Genes:", length(genes), "\n")
  cat("------------------------------------------------------------\n")

  conn_test <- tryCatch(
    httr::GET("https://dgidb.org", httr::timeout(10)),
    error = function(e) NULL
  )

  if (is.null(conn_test) || httr::status_code(conn_test) != 200) {
    warning("DGIdb not reachable. Skipping API query for ", query_label)
    return(list(
      interactions = data.frame(),
      gene_summary = data.frame(),
      drug_summary = data.frame(),
      relevant = data.frame()
    ))
  }

  query_genes <- paste0('"', genes, '"', collapse = ", ")

  graphql_query <- paste0(
    '{ genes(names: [', query_genes, ']) {
        nodes {
          name
          interactions {
            drug { name conceptId }
            interactionTypes { type directionality }
            publications { pmid }
            sources { fullName }
          }
        }
      }
    }'
  )

  response <- tryCatch(
    {
      httr::POST(
        url = "https://dgidb.org/api/graphql",
        body = list(query = graphql_query),
        encode = "json",
        httr::timeout(120)
      )
    },
    error = function(e) {
      warning("DGIdb request error for ", query_label, ": ", e$message)
      NULL
    }
  )

  if (is.null(response)) {
    return(list(
      interactions = data.frame(),
      gene_summary = data.frame(),
      drug_summary = data.frame(),
      relevant = data.frame()
    ))
  }

  cat("DGIdb response status:", httr::status_code(response), "\n")

  if (httr::status_code(response) != 200) {
    warning("DGIdb returned non-200 status for ", query_label)
    return(list(
      interactions = data.frame(),
      gene_summary = data.frame(),
      drug_summary = data.frame(),
      relevant = data.frame()
    ))
  }

  raw_text <- httr::content(response, "text", encoding = "UTF-8")
  result <- jsonlite::fromJSON(raw_text, flatten = TRUE)

  if ("errors" %in% names(result)) {
    print(result$errors)
    warning("DGIdb GraphQL error for ", query_label)
    return(list(
      interactions = data.frame(),
      gene_summary = data.frame(),
      drug_summary = data.frame(),
      relevant = data.frame()
    ))
  }

  genes_data <- result$data$genes$nodes

  if (is.null(genes_data) || nrow(genes_data) == 0) {
    warning("No genes returned by DGIdb for ", query_label)
    return(list(
      interactions = data.frame(),
      gene_summary = data.frame(),
      drug_summary = data.frame(),
      relevant = data.frame()
    ))
  }

  interactions_list <- lapply(seq_len(nrow(genes_data)), function(i) {

    gene_name <- genes_data$name[i]
    ints <- genes_data$interactions[[i]]

    if (is.null(ints) || !is.data.frame(ints) || nrow(ints) == 0) {
      return(NULL)
    }

    if (!"drug.name" %in% colnames(ints)) {
      return(NULL)
    }

    valid <- !is.na(ints$drug.name) & ints$drug.name != ""
    ints <- ints[valid, , drop = FALSE]

    if (nrow(ints) == 0) {
      return(NULL)
    }

    drug_name <- ints$drug.name

    drug_id <- if ("drug.conceptId" %in% colnames(ints)) {
      ints$drug.conceptId
    } else {
      rep(NA_character_, nrow(ints))
    }

    int_types <- sapply(seq_len(nrow(ints)), function(j) {
      x <- ints$interactionTypes[[j]]
      if (is.null(x) || !is.data.frame(x) || nrow(x) == 0) return(NA_character_)
      paste(unique(x$type), collapse = ";")
    })

    int_dirs <- sapply(seq_len(nrow(ints)), function(j) {
      x <- ints$interactionTypes[[j]]
      if (is.null(x) || !is.data.frame(x) || nrow(x) == 0) return(NA_character_)
      paste(unique(x$directionality), collapse = ";")
    })

    n_pubs <- sapply(seq_len(nrow(ints)), function(j) {
      x <- ints$publications[[j]]
      if (is.null(x) || !is.data.frame(x)) 0L else nrow(x)
    })

    sources <- sapply(seq_len(nrow(ints)), function(j) {
      x <- ints$sources[[j]]
      if (is.null(x) || !is.data.frame(x) || nrow(x) == 0) return(NA_character_)
      paste(unique(x$fullName), collapse = ";")
    })

    data.frame(
      query_label = query_label,
      gene = gene_name,
      drug = drug_name,
      drug_id = drug_id,
      interaction_type = int_types,
      directionality = int_dirs,
      n_publications = n_pubs,
      sources = sources,
      stringsAsFactors = FALSE
    )
  })

  interactions_df <- bind_rows(interactions_list)

  if (nrow(interactions_df) == 0) {
    warning("No usable interactions found for ", query_label)
    return(list(
      interactions = data.frame(),
      gene_summary = data.frame(),
      drug_summary = data.frame(),
      relevant = data.frame()
    ))
  }

  direction_map <- gene_table %>%
    select(gene, direction, mt_gene, symbol_class) %>%
    distinct()

  interactions_df <- interactions_df %>%
    left_join(direction_map, by = "gene") %>%
    mutate(
      gene_direction = case_when(
        direction == "UP" ~ "AD_upregulated",
        direction == "DOWN" ~ "AD_downregulated",
        TRUE ~ NA_character_
      ),
      reversal_rule = case_when(
        direction == "UP" ~ "drug should inhibit or antagonize target",
        direction == "DOWN" ~ "drug should activate or agonize target",
        TRUE ~ NA_character_
      ),
      therapeutically_relevant_direction = case_when(
        direction == "UP" &
          grepl("inhibitor|antagonist|blocker|suppressor|negative modulator",
                interaction_type, ignore.case = TRUE) ~ TRUE,
        direction == "DOWN" &
          grepl("activator|agonist|stimulator|positive modulator",
                interaction_type, ignore.case = TRUE) ~ TRUE,
        TRUE ~ FALSE
      )
    )

  gene_summary <- interactions_df %>%
    group_by(query_label, gene, gene_direction, direction, mt_gene, symbol_class) %>%
    summarise(
      n_drugs = n_distinct(drug),
      n_interactions = n(),
      n_relevant_direction_interactions =
        sum(therapeutically_relevant_direction, na.rm = TRUE),
      interaction_types = paste(unique(na.omit(interaction_type)), collapse = "; "),
      .groups = "drop"
    ) %>%
    arrange(desc(n_relevant_direction_interactions), desc(n_drugs))

  drug_summary <- interactions_df %>%
    group_by(query_label, drug) %>%
    summarise(
      n_genes = n_distinct(gene),
      target_genes = paste(unique(gene), collapse = ", "),
      target_directions = paste(unique(na.omit(gene_direction)), collapse = ", "),
      n_relevant_direction_interactions =
        sum(therapeutically_relevant_direction, na.rm = TRUE),
      interaction_types = paste(unique(na.omit(interaction_type)), collapse = "; "),
      total_publications = sum(n_publications, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_relevant_direction_interactions),
            desc(n_genes),
            desc(total_publications))

  relevant <- interactions_df %>%
    filter(therapeutically_relevant_direction) %>%
    arrange(gene_direction, gene, drug)

  write.csv(
    interactions_df,
    file.path(DRUG_DIR, paste0(output_prefix, "_interactions.csv")),
    row.names = FALSE
  )

  write.csv(
    gene_summary,
    file.path(DRUG_DIR, paste0(output_prefix, "_gene_summary.csv")),
    row.names = FALSE
  )

  write.csv(
    drug_summary,
    file.path(DRUG_DIR, paste0(output_prefix, "_drug_summary.csv")),
    row.names = FALSE
  )

  write.csv(
    relevant,
    file.path(DRUG_DIR, paste0(output_prefix, "_relevant_direction_interactions.csv")),
    row.names = FALSE
  )

  saveRDS(
    list(
      interactions = interactions_df,
      gene_summary = gene_summary,
      drug_summary = drug_summary,
      relevant = relevant
    ),
    file.path(RDS_DIR, paste0(output_prefix, ".rds"))
  )

  cat("DGIdb summary for", query_label, ":\n")
  cat("  Interactions:", nrow(interactions_df), "\n")
  cat("  Genes with interactions:", n_distinct(interactions_df$gene), "of", length(genes), "\n")
  cat("  Unique drugs:", n_distinct(interactions_df$drug), "\n")
  cat("  Directionally relevant interactions:", nrow(relevant), "\n\n")

  return(list(
    interactions = interactions_df,
    gene_summary = gene_summary,
    drug_summary = drug_summary,
    relevant = relevant
  ))
}

# ============================================================
# SECTION 9: Run optional DGIdb queries
# ============================================================

cat("Running optional DGIdb queries if API is reachable...\n\n")

strict_dgidb <- query_dgidb(
  gene_table = strict_signature,
  query_label = "strict48_primary",
  output_prefix = "dgidb_strict48"
)

expanded_nonMT_table <- expanded_signature %>%
  filter(!mt_gene)

expanded_nonMT_dgidb <- query_dgidb(
  gene_table = expanded_nonMT_table,
  query_label = "expanded125_nonMT_sensitivity",
  output_prefix = "dgidb_expanded125_nonMT"
)

# ============================================================
# SECTION 10: Combined DGIdb summary if available
# ============================================================

combined_interactions <- bind_rows(
  strict_dgidb$interactions,
  expanded_nonMT_dgidb$interactions
)

combined_drug_summary <- bind_rows(
  strict_dgidb$drug_summary,
  expanded_nonMT_dgidb$drug_summary
)

combined_relevant <- bind_rows(
  strict_dgidb$relevant,
  expanded_nonMT_dgidb$relevant
)

if (nrow(combined_interactions) > 0) {

  write.csv(
    combined_interactions,
    file.path(DRUG_DIR, "dgidb_combined_interactions_AD61026.csv"),
    row.names = FALSE
  )

  write.csv(
    combined_drug_summary,
    file.path(DRUG_DIR, "dgidb_combined_drug_summary_AD61026.csv"),
    row.names = FALSE
  )

  write.csv(
    combined_relevant,
    file.path(DRUG_DIR, "dgidb_combined_relevant_direction_interactions_AD61026.csv"),
    row.names = FALSE
  )

  # Simple review plot: top strict-primary target genes by relevant interactions.
  strict_gene_plot_df <- strict_dgidb$gene_summary

  if (nrow(strict_gene_plot_df) > 0) {

    p_gene <- strict_gene_plot_df %>%
      slice_max(order_by = n_relevant_direction_interactions, n = 20, with_ties = FALSE) %>%
      ggplot(
        aes(
          x = reorder(gene, n_relevant_direction_interactions),
          y = n_relevant_direction_interactions
        )
      ) +
      geom_col() +
      coord_flip() +
      labs(
        title = "DGIdb directionally relevant interactions",
        subtitle = "Strict 48-gene primary consensus signature",
        x = "Gene",
        y = "Number of directionally relevant interactions"
      ) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold")
      )

    ggsave(
      file.path(FIG_DIR, "Fig_DGIdb_strict48_relevant_targets_draft.pdf"),
      p_gene,
      width = 8,
      height = 6
    )
  }
}

# ============================================================
# SECTION 11: Manual upload instructions table
# ============================================================

manual_tool_instructions <- data.frame(
  tool = c(
    "DREIMT primary",
    "DREIMT sensitivity",
    "CMap primary",
    "CMap sensitivity",
    "iLINCS primary",
    "iLINCS sensitivity"
  ),
  up_file = c(
    "DREIMT_PRIMARY_strict48_up.txt",
    "DREIMT_SENSITIVITY_expanded125_nonMT_up.txt",
    "CMap_PRIMARY_strict48_up.txt",
    "CMap_SENSITIVITY_expanded125_nonMT_up.txt",
    "iLINCS_PRIMARY_strict48_up.txt",
    "iLINCS_SENSITIVITY_expanded125_nonMT_up.txt"
  ),
  down_file = c(
    "DREIMT_PRIMARY_strict48_down.txt",
    "DREIMT_SENSITIVITY_expanded125_nonMT_down.txt",
    "CMap_PRIMARY_strict48_down.txt",
    "CMap_SENSITIVITY_expanded125_nonMT_down.txt",
    "iLINCS_PRIMARY_strict48_down.txt",
    "iLINCS_SENSITIVITY_expanded125_nonMT_down.txt"
  ),
  interpretation = c(
    "Primary immune-enriched perturbational query.",
    "Sensitivity query excluding MT genes from expanded signature.",
    "Primary L1000 Touchstone reversal query.",
    "Sensitivity L1000 query excluding MT genes.",
    "Primary LINCS/iLINCS reversal query.",
    "Sensitivity iLINCS query excluding MT genes."
  ),
  result_file_to_save_later = c(
    "dreimt_primary_strict48_results.csv",
    "dreimt_sensitivity_expanded125_nonMT_results.csv",
    "cmap_primary_strict48_results.csv",
    "cmap_sensitivity_expanded125_nonMT_results.csv",
    "ilincs_primary_strict48_results.csv",
    "ilincs_sensitivity_expanded125_nonMT_results.csv"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  manual_tool_instructions,
  file.path(DRUG_DIR, "manual_drug_tool_submission_instructions_AD61026.csv"),
  row.names = FALSE
)

# ============================================================
# SECTION 12: Completion summary
# ============================================================

cat("\n============================================================\n")
cat("06A drug query preparation complete.\n\n")

cat("Main query files are in:\n")
cat("  ", DRUG_DIR, "\n\n")

cat("Primary query:\n")
cat("  strict48_up_genes.txt:", length(strict_up), "genes\n")
cat("  strict48_down_genes.txt:", length(strict_down), "genes\n\n")

cat("Sensitivity query excluding MT genes:\n")
cat("  expanded125_nonMT_up_genes.txt:", length(expanded_nonMT_up), "genes\n")
cat("  expanded125_nonMT_down_genes.txt:", length(expanded_nonMT_down), "genes\n\n")

cat("DGIdb outputs:\n")
if (nrow(combined_interactions) > 0) {
  cat("  DGIdb query completed and results were saved.\n")
  cat("  Combined interactions:", nrow(combined_interactions), "\n")
  cat("  Directionally relevant interactions:", nrow(combined_relevant), "\n\n")
} else {
  cat("  DGIdb results were not generated, likely due to API/network limits.\n")
  cat("  Query files were still generated for manual use.\n\n")
}

cat("Next step:\n")
cat("  Upload drug_query_gene_list_manifest_AD61026.csv and any DGIdb output files.\n")
cat("============================================================\n")
