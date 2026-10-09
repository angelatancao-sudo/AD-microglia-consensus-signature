# ============================================================
# 04_differential_expression_and_replication.R
# ============================================================
# Project: AD Microglia Multi-Dataset Transcriptomic Analysis
#
# Purpose:
#   Run donor-level pseudobulk differential expression after the
#   field-aligned microglia extraction rebuild.
#
# Key design:
#
#   1. Discovery dataset
#        SEA-AD Microglia/Immune donor-level pseudobulk
#
#   2. Primary replication datasets
#        GSE174367  author-annotated microglia
#        GSE157827  cluster-validated microglia, PRIMARY clusters
#        GSE160936  cluster-validated microglia, PRIMARY clusters
#        GSE188545  cluster-validated microglia, PRIMARY clusters
#
#   3. Supplementary dataset
#        GSE243292  microglia-specific/preannotated H5AD
#        Important: this has only 2 control donors after A/T filtering,
#        so it is reported as supplementary support, not used to define
#        primary replication confidence.
#
#   4. Sensitivity analysis
#        Re-run the full primary replication assessment while holding
#        GSE174367 and GSE157827 constant and using expanded cluster
#        definitions for the datasets with validated borderline
#        microglia-like clusters: GSE160936 and GSE188545.
#
# Field-practice principles used here:
#   - donor-level pseudobulk counts, not cell-level DE
#   - raw counts are summed by donor
#   - edgeR TMM normalization
#   - limma-voom mean-variance modeling
#   - empirical Bayes moderation
#   - no covariates are added when cohort metadata are incomplete or
#     sample size is too small to support stable covariate adjustment
#
# Interpretation:
#   Positive logFC = higher expression in AD than Control.
#
# Run after:
#   03F_check_pseudobulk_inputs_before_DE_AD61026.R shows OK
#   for SEAAD, GSE174367, GSE157827, GSE160936, GSE188545,
#   and GSE243292.
#
# Outputs:
#   Output/DE/SEAAD_DE_results.csv
#   Output/DE/[dataset]_DE_results_PRIMARY.csv
#   Output/DE/GSE243292_DE_results_SUPPLEMENTARY.csv
#   Output/DE/replication_assessment_PRIMARY.csv
#   Output/DE/high_confidence_genes_PRIMARY.csv
#   Output/DE/moderate_confidence_genes_PRIMARY.csv
#   Output/DE/supplementary_support_GSE243292.csv
#   Output/DE/final_consensus_gene_list_PRIMARY.csv
#   Output/DE/sensitivity_replication_assessment.csv
#   Output/DE/sensitivity_comparison_summary.csv
#
#   RDS/all_DE_results_PRIMARY.rds
#   RDS/final_consensus_gene_list_PRIMARY.rds
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
  library(edgeR)
  library(limma)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})

# RankProd is used for the rank-based replication meta-analysis.
# If it is not installed, the script will still run all DE and replication
# assessment steps, but RankProd outputs will be skipped.
HAS_RANKPROD <- requireNamespace("RankProd", quietly = TRUE)
if (HAS_RANKPROD) {
  suppressPackageStartupMessages(library(RankProd))
} else {
  warning("RankProd package not installed. RankProd meta-analysis will be skipped.")
}

dir.create(DE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("  04: Donor-level pseudobulk DE analysis\n")
cat("============================================================\n\n")

# ============================================================
# SECTION 2: Dataset registry
# ============================================================

DISCOVERY <- list(
  id = "SEAAD",
  label = "SEA-AD",
  role = "Discovery",
  counts_file = file.path(RDS_DIR, "SEAAD_pseudobulk_counts.rds"),
  meta_file = file.path(RDS_DIR, "SEAAD_pseudobulk_meta.rds"),
  dx_col = "diagnosis_std"
)

PRIMARY_REPLICATION <- list(
  GSE174367 = list(
    id = "GSE174367",
    label = "GSE174367",
    role = "Primary replication",
    counts_file = file.path(RDS_DIR, "GSE174367_pseudobulk_counts.rds"),
    meta_file = file.path(RDS_DIR, "GSE174367_pseudobulk_meta.rds"),
    dx_col = "diagnosis_std"
  ),
  GSE157827 = list(
    id = "GSE157827",
    label = "GSE157827 validated microglia primary",
    role = "Primary replication",
    counts_file = file.path(RDS_DIR, "GSE157827_microglia_PRIMARY_pseudobulk_counts.rds"),
    meta_file = file.path(RDS_DIR, "GSE157827_microglia_PRIMARY_pseudobulk_meta.rds"),
    dx_col = "diagnosis_std"
  ),
  GSE160936 = list(
    id = "GSE160936",
    label = "GSE160936 validated microglia primary",
    role = "Primary replication",
    counts_file = file.path(RDS_DIR, "GSE160936_microglia_PRIMARY_pseudobulk_counts.rds"),
    meta_file = file.path(RDS_DIR, "GSE160936_microglia_PRIMARY_pseudobulk_meta.rds"),
    dx_col = "diagnosis_std"
  ),
  GSE188545 = list(
    id = "GSE188545",
    label = "GSE188545 validated microglia primary",
    role = "Primary replication",
    counts_file = file.path(RDS_DIR, "GSE188545_microglia_PRIMARY_pseudobulk_counts.rds"),
    meta_file = file.path(RDS_DIR, "GSE188545_microglia_PRIMARY_pseudobulk_meta.rds"),
    dx_col = "diagnosis_std"
  )
)

SUPPLEMENTARY <- list(
  GSE243292 = list(
    id = "GSE243292",
    label = "GSE243292 supplementary A/T-defined microglia",
    role = "Supplementary support",
    counts_file = file.path(RDS_DIR, "GSE243292_pseudobulk_counts.rds"),
    meta_file = file.path(RDS_DIR, "GSE243292_pseudobulk_meta.rds"),
    dx_col = "diagnosis_std",
    note = "A+T+ treated as AD, A-T- as Control, A+T- excluded; only 2 controls after filtering."
  )
)

SENSITIVITY_REPLICATION <- list(
  GSE174367 = PRIMARY_REPLICATION$GSE174367,
  GSE157827 = PRIMARY_REPLICATION$GSE157827,
  GSE160936 = list(
    id = "GSE160936",
    label = "GSE160936 validated microglia sensitivity",
    role = "Sensitivity replication",
    counts_file = file.path(RDS_DIR, "GSE160936_microglia_SENSITIVITY_pseudobulk_counts.rds"),
    meta_file = file.path(RDS_DIR, "GSE160936_microglia_SENSITIVITY_pseudobulk_meta.rds"),
    dx_col = "diagnosis_std"
  ),
  GSE188545 = list(
    id = "GSE188545",
    label = "GSE188545 validated microglia sensitivity",
    role = "Sensitivity replication",
    counts_file = file.path(RDS_DIR, "GSE188545_microglia_SENSITIVITY_pseudobulk_counts.rds"),
    meta_file = file.path(RDS_DIR, "GSE188545_microglia_SENSITIVITY_pseudobulk_meta.rds"),
    dx_col = "diagnosis_std"
  )
)

KEY_GENES_DOWN <- PARAMS$key_down_genes
KEY_GENES_UP <- PARAMS$key_up_genes
KEY_GENES_ALL <- unique(c(KEY_GENES_DOWN, KEY_GENES_UP))

# ============================================================
# SECTION 3: Utility functions
# ============================================================

# ------------------------------------------------------------
# Function: load and align pseudobulk count matrix and metadata
# ------------------------------------------------------------
load_pseudobulk <- function(ds) {

  if (!file.exists(ds$counts_file)) {
    stop("Missing counts file for ", ds$id, ": ", ds$counts_file)
  }
  if (!file.exists(ds$meta_file)) {
    stop("Missing metadata file for ", ds$id, ": ", ds$meta_file)
  }

  counts <- readRDS(ds$counts_file)
  meta <- readRDS(ds$meta_file)

  # Convert to ordinary matrix for edgeR/limma stability.
  # These pseudobulk matrices are small enough for this to be safe.
  if (inherits(counts, "sparseMatrix")) {
    counts <- as.matrix(counts)
  }

  # Check alignment. Newer files use rownames(meta) = colnames(counts).
  if (!identical(rownames(meta), colnames(counts))) {

    if ("pseudobulk_id" %in% colnames(meta) &&
        all(colnames(counts) %in% as.character(meta$pseudobulk_id))) {
      meta <- meta[match(colnames(counts), as.character(meta$pseudobulk_id)), , drop = FALSE]
      rownames(meta) <- meta$pseudobulk_id

    } else if ("donor_id" %in% colnames(meta) &&
               all(colnames(counts) %in% as.character(meta$donor_id))) {
      meta <- meta[match(colnames(counts), as.character(meta$donor_id)), , drop = FALSE]
      rownames(meta) <- colnames(counts)

    } else {
      stop("Could not align count columns and metadata rows for ", ds$id)
    }
  }

  if (!identical(rownames(meta), colnames(counts))) {
    stop("Alignment failed for ", ds$id)
  }

  return(list(counts = counts, meta = meta))
}

# ------------------------------------------------------------
# Function: run edgeR-voom-limma DE
# ------------------------------------------------------------
run_pseudobulk_DE <- function(counts,
                              meta,
                              dx_col,
                              dataset_id,
                              dataset_label,
                              plot_voom = FALSE) {

  cat("------------------------------------------------------------\n")
  cat("DE: ", dataset_id, " — ", dataset_label, "\n", sep = "")
  cat("------------------------------------------------------------\n")

  if (!dx_col %in% colnames(meta)) {
    stop("Diagnosis column ", dx_col, " not found for ", dataset_id)
  }

  dx <- factor(as.character(meta[[dx_col]]), levels = c("Control", "AD"))

  if (any(is.na(dx))) {
    stop("Missing or unsupported diagnosis values in ", dataset_id)
  }

  n_ad <- sum(dx == "AD")
  n_control <- sum(dx == "Control")

  cat("  Donors:", ncol(counts), "\n")
  cat("  AD:", n_ad, " Control:", n_control, "\n")

  if (n_ad < 2 || n_control < 2) {
    warning(dataset_id, " has fewer than 2 donors in one group. DE may be unstable.")
  }

  # Optional low-cell warning for transparency.
  if ("n_microglia_cells" %in% colnames(meta)) {
    min_cells <- min(meta$n_microglia_cells, na.rm = TRUE)
    med_cells <- median(meta$n_microglia_cells, na.rm = TRUE)
    cat("  Microglia cells per donor: min=", min_cells,
        " median=", med_cells, "\n", sep = "")

    if (min_cells < 50) {
      warning(dataset_id, " has at least one donor with <50 microglia nuclei. ",
              "Keeping it for primary analysis but flagging for interpretation.")
    }
  }

  # Create edgeR object from raw pseudobulk counts.
  dge <- DGEList(counts = counts, samples = meta, group = dx)

  # filterByExpr uses the experimental design to keep genes with
  # enough expression to be statistically testable.
  design <- model.matrix(~ 0 + dx)
  colnames(design) <- c("Control", "AD")

  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  cat("  Genes before filtering:", nrow(counts), "\n")
  cat("  Genes after filtering :", nrow(dge), "\n")

  # TMM normalization for compositional differences in library size.
  dge <- calcNormFactors(dge, method = "TMM")

  # Voom estimates the mean-variance trend and observation weights.
  v <- voom(dge, design, plot = plot_voom)

  # Fit linear model and contrast AD minus Control.
  contrast <- makeContrasts(AD_vs_Control = AD - Control, levels = design)
  fit <- lmFit(v, design)
  fit <- contrasts.fit(fit, contrast)
  fit <- eBayes(fit, robust = TRUE)

  res <- topTable(
    fit,
    coef = "AD_vs_Control",
    number = Inf,
    sort.by = "P",
    adjust.method = "BH"
  ) %>%
    rownames_to_column("gene") %>%
    mutate(
      dataset = dataset_id,
      dataset_label = dataset_label,
      direction = case_when(
        logFC > 0 ~ "UP_in_AD",
        logFC < 0 ~ "DOWN_in_AD",
        TRUE ~ "NO_CHANGE"
      ),
      significant_discovery_threshold =
        adj.P.Val < PARAMS$fdr_threshold &
        abs(logFC) > PARAMS$logfc_threshold
    )

  n_up <- sum(res$adj.P.Val < PARAMS$fdr_threshold &
                res$logFC > PARAMS$logfc_threshold, na.rm = TRUE)
  n_down <- sum(res$adj.P.Val < PARAMS$fdr_threshold &
                  res$logFC < -PARAMS$logfc_threshold, na.rm = TRUE)

  cat("  Significant genes at FDR<", PARAMS$fdr_threshold,
      " and |logFC|>", PARAMS$logfc_threshold, ":\n", sep = "")
  cat("    UP in AD:", n_up, "\n")
  cat("    DOWN in AD:", n_down, "\n\n")

  return(res)
}

# ------------------------------------------------------------
# Function: run DE for a registry list
# ------------------------------------------------------------
run_DE_for_registry <- function(registry, suffix) {

  results_list <- list()

  for (ds_id in names(registry)) {

    ds <- registry[[ds_id]]
    pb <- load_pseudobulk(ds)

    res <- run_pseudobulk_DE(
      counts = pb$counts,
      meta = pb$meta,
      dx_col = ds$dx_col,
      dataset_id = ds$id,
      dataset_label = ds$label
    )

    results_list[[ds_id]] <- res

    write.csv(
      res,
      file.path(DE_DIR, paste0(ds_id, "_DE_results_", suffix, ".csv")),
      row.names = FALSE
    )

    rm(pb, res)
    gc(verbose = FALSE)
  }

  return(results_list)
}

# ------------------------------------------------------------
# Function: volcano plot for discovery dataset
# ------------------------------------------------------------
make_discovery_volcano <- function(res) {

  res <- res %>%
    mutate(
      volcano_group = case_when(
        adj.P.Val < PARAMS$fdr_threshold &
          logFC > PARAMS$logfc_threshold ~ "UP",
        adj.P.Val < PARAMS$fdr_threshold &
          logFC < -PARAMS$logfc_threshold ~ "DOWN",
        TRUE ~ "NS"
      ),
      neg_log10_p = pmin(-log10(P.Value), 20),
      label = ifelse(gene %in% KEY_GENES_ALL, gene, "")
    )

  p <- ggplot(res, aes(x = logFC, y = neg_log10_p, color = volcano_group)) +
    geom_point(size = 0.6, alpha = 0.65) +
    geom_text_repel(
      aes(label = label),
      size = 2.6,
      max.overlaps = 25,
      color = "black",
      fontface = "italic"
    ) +
    scale_color_manual(
      values = c(
        "UP" = PARAMS$up_color,
        "DOWN" = PARAMS$down_color,
        "NS" = PARAMS$ns_color
      )
    ) +
    geom_vline(
      xintercept = c(-PARAMS$logfc_threshold, PARAMS$logfc_threshold),
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.4
    ) +
    geom_hline(
      yintercept = -log10(PARAMS$fdr_threshold),
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.4
    ) +
    labs(
      title = "SEA-AD discovery: AD vs Control microglia",
      subtitle = "Donor-level pseudobulk edgeR-voom-limma",
      x = expression(log[2]~"fold change (AD / Control)"),
      y = expression(-log[10]~"p-value"),
      color = ""
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggsave(
    file.path(FIG_DIR, "SEAAD_discovery_volcano.pdf"),
    p,
    width = 7,
    height = 6
  )
}

# ============================================================
# SECTION 4: Discovery DE
# ============================================================

cat("============================================================\n")
cat("STAGE 1: SEA-AD discovery DE\n")
cat("============================================================\n\n")

pb_disc <- load_pseudobulk(DISCOVERY)

disc_res <- run_pseudobulk_DE(
  counts = pb_disc$counts,
  meta = pb_disc$meta,
  dx_col = DISCOVERY$dx_col,
  dataset_id = DISCOVERY$id,
  dataset_label = DISCOVERY$label
)

write.csv(
  disc_res,
  file.path(DE_DIR, "SEAAD_DE_results.csv"),
  row.names = FALSE
)
saveRDS(disc_res, file.path(RDS_DIR, "SEAAD_DE_results.rds"))

make_discovery_volcano(disc_res)

discovery_genes <- disc_res %>%
  filter(adj.P.Val < PARAMS$fdr_threshold,
         abs(logFC) > PARAMS$logfc_threshold) %>%
  arrange(adj.P.Val)

write.csv(
  discovery_genes,
  file.path(DE_DIR, "SEAAD_discovery_genes.csv"),
  row.names = FALSE
)
saveRDS(discovery_genes, file.path(RDS_DIR, "discovery_genes.rds"))

cat("Discovery genes:", nrow(discovery_genes), "\n")
cat("  UP in AD:", sum(discovery_genes$logFC > 0), "\n")
cat("  DOWN in AD:", sum(discovery_genes$logFC < 0), "\n\n")

rm(pb_disc)
gc(verbose = FALSE)

# ============================================================
# SECTION 5: Primary replication DE
# ============================================================

cat("============================================================\n")
cat("STAGE 2: Primary replication DE\n")
cat("============================================================\n\n")

primary_rep_results <- run_DE_for_registry(
  registry = PRIMARY_REPLICATION,
  suffix = "PRIMARY"
)

all_DE_results_PRIMARY <- c(list(SEAAD = disc_res), primary_rep_results)
saveRDS(
  all_DE_results_PRIMARY,
  file.path(RDS_DIR, "all_DE_results_PRIMARY.rds")
)

# ============================================================
# SECTION 6: Supplementary DE
# ============================================================

cat("============================================================\n")
cat("STAGE 3: Supplementary DE\n")
cat("============================================================\n\n")

supp_results <- run_DE_for_registry(
  registry = SUPPLEMENTARY,
  suffix = "SUPPLEMENTARY"
)

saveRDS(
  supp_results,
  file.path(RDS_DIR, "supplementary_DE_results.rds")
)

# ============================================================
# SECTION 7: Primary replication assessment
# ============================================================
# Because the primary replication set contains 4 independent cohorts,
# the conservative replication rule is:
#
#   High confidence:
#      same direction in at least 3 of 4 primary replication datasets
#      AND nominal p < 0.05 in at least 2 of 4 datasets
#
#   Moderate confidence:
#      same direction in at least 3 of 4 datasets
#      OR nominal p < 0.05 in at least 2 of 4 datasets
#
# GSE243292 is not used to define this confidence category because it
# has only 2 control donors after A/T filtering. It is assessed later
# as supplementary support.

assess_replication <- function(discovery_genes,
                               rep_results,
                               label = "PRIMARY") {

  rep_dataset_ids <- names(rep_results)

  out <- lapply(discovery_genes$gene, function(g) {

    disc_row <- discovery_genes[discovery_genes$gene == g, ][1, ]
    disc_dir <- ifelse(disc_row$logFC > 0, "UP", "DOWN")

    same_direction <- c()
    nominal_sig <- c()
    rep_logfc <- c()
    rep_p <- c()

    per_dataset <- list()

    for (ds_id in rep_dataset_ids) {

      res <- rep_results[[ds_id]]
      row <- res[res$gene == g, ]

      if (nrow(row) > 0) {
        obs_logfc <- row$logFC[1]
        obs_dir <- ifelse(obs_logfc > 0, "UP", "DOWN")
        obs_p <- row$P.Value[1]

        same_direction <- c(same_direction, obs_dir == disc_dir)
        nominal_sig <- c(nominal_sig, obs_p < 0.05)
        rep_logfc <- c(rep_logfc, obs_logfc)
        rep_p <- c(rep_p, obs_p)

        per_dataset[[paste0(ds_id, "_logFC")]] <- obs_logfc
        per_dataset[[paste0(ds_id, "_P.Value")]] <- obs_p
        per_dataset[[paste0(ds_id, "_same_direction")]] <- obs_dir == disc_dir

      } else {
        same_direction <- c(same_direction, NA)
        nominal_sig <- c(nominal_sig, NA)
        rep_logfc <- c(rep_logfc, NA)
        rep_p <- c(rep_p, NA)

        per_dataset[[paste0(ds_id, "_logFC")]] <- NA_real_
        per_dataset[[paste0(ds_id, "_P.Value")]] <- NA_real_
        per_dataset[[paste0(ds_id, "_same_direction")]] <- NA
      }
    }

    n_available <- sum(!is.na(same_direction))
    n_same <- sum(same_direction, na.rm = TRUE)
    n_nominal <- sum(nominal_sig, na.rm = TRUE)

    high_conf <- n_same >= 3 && n_nominal >= 2
    moderate_conf <- !high_conf && (n_same >= 3 || n_nominal >= 2)

    confidence <- ifelse(
      high_conf,
      "High",
      ifelse(moderate_conf, "Moderate", "Not_replicated")
    )

    base <- data.frame(
      analysis = label,
      gene = g,
      discovery_logFC = disc_row$logFC,
      discovery_adj.P.Val = disc_row$adj.P.Val,
      discovery_direction = disc_dir,
      n_primary_datasets_available = n_available,
      n_same_direction = n_same,
      n_nominal_p_lt_0.05 = n_nominal,
      mean_replication_logFC = mean(rep_logfc, na.rm = TRUE),
      confidence = confidence,
      stringsAsFactors = FALSE
    )

    bind_cols(base, as.data.frame(per_dataset, check.names = FALSE))
  })

  bind_rows(out)
}

cat("============================================================\n")
cat("STAGE 4: Primary replication assessment\n")
cat("============================================================\n\n")

replication_primary <- assess_replication(
  discovery_genes = discovery_genes,
  rep_results = primary_rep_results,
  label = "PRIMARY"
)

write.csv(
  replication_primary,
  file.path(DE_DIR, "replication_assessment_PRIMARY.csv"),
  row.names = FALSE
)
saveRDS(
  replication_primary,
  file.path(RDS_DIR, "replication_assessment_PRIMARY.rds")
)

high_conf <- replication_primary %>%
  filter(confidence == "High") %>%
  arrange(discovery_adj.P.Val)

moderate_conf <- replication_primary %>%
  filter(confidence == "Moderate") %>%
  arrange(discovery_adj.P.Val)

write.csv(
  high_conf,
  file.path(DE_DIR, "high_confidence_genes_PRIMARY.csv"),
  row.names = FALSE
)
write.csv(
  moderate_conf,
  file.path(DE_DIR, "moderate_confidence_genes_PRIMARY.csv"),
  row.names = FALSE
)

cat("Primary replication confidence:\n")
print(table(replication_primary$confidence))
cat("\n")

# ============================================================
# SECTION 8: Supplementary support from GSE243292
# ============================================================

cat("============================================================\n")
cat("STAGE 5: Supplementary GSE243292 support\n")
cat("============================================================\n\n")

if ("GSE243292" %in% names(supp_results)) {

  supp_res <- supp_results$GSE243292

  supplementary_support <- replication_primary %>%
    select(gene, discovery_logFC, discovery_direction, confidence) %>%
    left_join(
      supp_res %>%
        select(gene, GSE243292_logFC = logFC,
               GSE243292_P.Value = P.Value,
               GSE243292_adj.P.Val = adj.P.Val),
      by = "gene"
    ) %>%
    mutate(
      GSE243292_direction = case_when(
        GSE243292_logFC > 0 ~ "UP",
        GSE243292_logFC < 0 ~ "DOWN",
        TRUE ~ NA_character_
      ),
      GSE243292_same_direction =
        GSE243292_direction == discovery_direction,
      GSE243292_nominal_p_lt_0.05 =
        GSE243292_P.Value < 0.05
    )

  write.csv(
    supplementary_support,
    file.path(DE_DIR, "supplementary_support_GSE243292.csv"),
    row.names = FALSE
  )

  cat("Supplementary support table saved.\n\n")
}

# ============================================================
# SECTION 9: RankProd primary replication meta-analysis
# ============================================================
# This is a rank-based consensus analysis using only the 4 primary
# replication datasets. GSE243292 is excluded because of its highly
# imbalanced 8 AD / 2 Control design.

run_rankprod_primary <- function(rep_results) {

  if (!HAS_RANKPROD) {
    cat("RankProd skipped because package is not installed.\n\n")
    return(NULL)
  }

  cat("============================================================\n")
  cat("STAGE 6: RankProd primary replication meta-analysis\n")
  cat("============================================================\n\n")

  # Use genes tested in all primary replication datasets.
  common_genes <- Reduce(
    intersect,
    lapply(rep_results, function(res) res$gene)
  )

  common_genes <- sort(common_genes)

  cat("Genes tested in all primary replication datasets:",
      length(common_genes), "\n")

  fc_matrix <- sapply(names(rep_results), function(ds_id) {
    res <- rep_results[[ds_id]]
    res$logFC[match(common_genes, res$gene)]
  })

  rownames(fc_matrix) <- common_genes
  colnames(fc_matrix) <- names(rep_results)

  # Save the input matrix for transparency.
  write.csv(
    data.frame(gene = rownames(fc_matrix), fc_matrix, check.names = FALSE),
    file.path(DE_DIR, "rankprod_primary_logFC_matrix.csv"),
    row.names = FALSE
  )

  cat("Running RankProducts on logFC matrix...\n")
  cat("This may take several minutes.\n\n")

  set.seed(42)
  rp_out <- RankProducts(
    data = fc_matrix,
    cl = rep(1, ncol(fc_matrix)),
    logged = TRUE,
    na.rm = TRUE,
    plot = FALSE,
    rand = 42,
    gene.names = common_genes
  )

  # RankProd returns two directions. We retain both and label by AveFC sign.
  rp_down <- data.frame(
    gene = common_genes,
    RP = rp_out$RPs[, 1],
    pfp = rp_out$pfp[, 1],
    pval = rp_out$pval[, 1],
    AveFC = rp_out$AveFC[, 1],
    rankprod_direction = "DOWN_in_AD",
    stringsAsFactors = FALSE
  ) %>%
    arrange(pfp)

  rp_up <- data.frame(
    gene = common_genes,
    RP = rp_out$RPs[, 2],
    pfp = rp_out$pfp[, 2],
    pval = rp_out$pval[, 2],
    AveFC = rp_out$AveFC[, 1],
    rankprod_direction = "UP_in_AD",
    stringsAsFactors = FALSE
  ) %>%
    arrange(pfp)

  write.csv(
    rp_up,
    file.path(DE_DIR, "rankprod_primary_AD_upregulated.csv"),
    row.names = FALSE
  )
  write.csv(
    rp_down,
    file.path(DE_DIR, "rankprod_primary_AD_downregulated.csv"),
    row.names = FALSE
  )

  rp_combined <- bind_rows(rp_up, rp_down) %>%
    arrange(pfp)

  write.csv(
    rp_combined,
    file.path(DE_DIR, "rankprod_primary_results.csv"),
    row.names = FALSE
  )
  saveRDS(
    rp_combined,
    file.path(RDS_DIR, "rankprod_primary_results.rds")
  )

  cat("RankProd complete.\n")
  cat("  UP pfp<0.05:",
      sum(rp_up$pfp < 0.05, na.rm = TRUE), "\n")
  cat("  DOWN pfp<0.05:",
      sum(rp_down$pfp < 0.05, na.rm = TRUE), "\n\n")

  return(list(up = rp_up, down = rp_down, combined = rp_combined))
}

rankprod_primary <- run_rankprod_primary(primary_rep_results)

# ============================================================
# SECTION 10: Final primary consensus gene list
# ============================================================
# The primary consensus list is intentionally conservative:
#   - genes with High replication confidence
#   - plus RankProd pfp<0.05 genes, if RankProd ran successfully
#
# Moderate-confidence genes are saved separately and should be described
# as supportive but not part of the strict final consensus list unless
# you decide otherwise later.

cat("============================================================\n")
cat("STAGE 7: Final primary consensus gene list\n")
cat("============================================================\n\n")

final_sources <- list()

if (nrow(high_conf) > 0) {
  final_sources[["High_confidence_replication"]] <- high_conf %>%
    transmute(
      gene,
      direction = discovery_direction,
      source = "High_confidence_replication"
    )
}

if (!is.null(rankprod_primary)) {
  rp_sig_up <- rankprod_primary$up %>%
    filter(pfp < 0.05) %>%
    transmute(gene, direction = "UP", source = "RankProd_primary_pfp_lt_0.05")

  rp_sig_down <- rankprod_primary$down %>%
    filter(pfp < 0.05) %>%
    transmute(gene, direction = "DOWN", source = "RankProd_primary_pfp_lt_0.05")

  final_sources[["RankProd"]] <- bind_rows(rp_sig_up, rp_sig_down)
}

if (length(final_sources) > 0) {
  final_long <- bind_rows(final_sources)

  final_consensus <- final_long %>%
    group_by(gene, direction) %>%
    summarise(
      sources = paste(sort(unique(source)), collapse = ";"),
      n_sources = n_distinct(source),
      .groups = "drop"
    ) %>%
    arrange(direction, gene)

} else {
  final_consensus <- data.frame(
    gene = character(),
    direction = character(),
    sources = character(),
    n_sources = integer()
  )
}

write.csv(
  final_consensus,
  file.path(DE_DIR, "final_consensus_gene_list_PRIMARY.csv"),
  row.names = FALSE
)
saveRDS(
  final_consensus,
  file.path(RDS_DIR, "final_consensus_gene_list_PRIMARY.rds")
)

write.csv(
  moderate_conf,
  file.path(DE_DIR, "moderate_confidence_genes_SUPPORTIVE.csv"),
  row.names = FALSE
)

cat("Final primary consensus genes:", nrow(final_consensus), "\n")
if (nrow(final_consensus) > 0) {
  print(table(final_consensus$direction))
}
cat("\n")

# ============================================================
# SECTION 11: Sensitivity analysis for cluster definitions
# ============================================================

cat("============================================================\n")
cat("STAGE 8: Sensitivity analysis using expanded cluster calls\n")
cat("============================================================\n\n")

sensitivity_rep_results <- run_DE_for_registry(
  registry = SENSITIVITY_REPLICATION,
  suffix = "SENSITIVITY"
)

sensitivity_assessment <- assess_replication(
  discovery_genes = discovery_genes,
  rep_results = sensitivity_rep_results,
  label = "SENSITIVITY"
)

write.csv(
  sensitivity_assessment,
  file.path(DE_DIR, "sensitivity_replication_assessment.csv"),
  row.names = FALSE
)

saveRDS(
  sensitivity_assessment,
  file.path(RDS_DIR, "sensitivity_replication_assessment.rds")
)

sensitivity_comparison <- replication_primary %>%
  select(gene, primary_confidence = confidence,
         primary_n_same_direction = n_same_direction,
         primary_n_nominal = n_nominal_p_lt_0.05) %>%
  left_join(
    sensitivity_assessment %>%
      select(gene, sensitivity_confidence = confidence,
             sensitivity_n_same_direction = n_same_direction,
             sensitivity_n_nominal = n_nominal_p_lt_0.05),
    by = "gene"
  ) %>%
  mutate(
    confidence_changed = primary_confidence != sensitivity_confidence
  )

write.csv(
  sensitivity_comparison,
  file.path(DE_DIR, "sensitivity_comparison_summary.csv"),
  row.names = FALSE
)

cat("Sensitivity confidence table:\n")
print(table(sensitivity_assessment$confidence))
cat("\n")
cat("Genes with changed confidence:",
    sum(sensitivity_comparison$confidence_changed, na.rm = TRUE), "\n\n")

# ============================================================
# SECTION 12: Save run summary
# ============================================================

run_summary <- data.frame(
  item = c(
    "SEAAD discovery donors",
    "Primary replication datasets",
    "Supplementary datasets",
    "Discovery genes",
    "High confidence primary genes",
    "Moderate confidence supportive genes",
    "Final primary consensus genes",
    "Sensitivity confidence changes"
  ),
  value = c(
    ncol(load_pseudobulk(DISCOVERY)$counts),
    length(PRIMARY_REPLICATION),
    length(SUPPLEMENTARY),
    nrow(discovery_genes),
    nrow(high_conf),
    nrow(moderate_conf),
    nrow(final_consensus),
    sum(sensitivity_comparison$confidence_changed, na.rm = TRUE)
  )
)

write.csv(
  run_summary,
  file.path(DE_DIR, "04_DE_run_summary.csv"),
  row.names = FALSE
)

cat("============================================================\n")
cat("04 DE analysis complete.\n\n")
cat("Main outputs:\n")
cat("  Output/DE/SEAAD_DE_results.csv\n")
cat("  Output/DE/replication_assessment_PRIMARY.csv\n")
cat("  Output/DE/high_confidence_genes_PRIMARY.csv\n")
cat("  Output/DE/moderate_confidence_genes_SUPPORTIVE.csv\n")
cat("  Output/DE/supplementary_support_GSE243292.csv\n")
cat("  Output/DE/final_consensus_gene_list_PRIMARY.csv\n")
cat("  Output/DE/sensitivity_comparison_summary.csv\n\n")
cat("Next step:\n")
cat("  Review 04_DE_run_summary.csv and final_consensus_gene_list_PRIMARY.csv.\n")
cat("============================================================\n")
