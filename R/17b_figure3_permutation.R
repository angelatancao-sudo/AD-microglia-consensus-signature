# ============================================================
# Figure 7: Empirical permutation analysis of cross-cohort
# replication
#
# Purpose:
#   1. Shuffle gene labels independently within each replication
#      cohort.
#   2. Preserve the observed logFC/P-value pairing within each
#      cohort and the number of genes available in each cohort.
#   3. Reapply the exact primary replication criteria.
#   4. Generate null distributions for:
#        A. High-confidence replication
#        B. High + Moderate replication
#   5. Compare the observed counts with the permutation null.
#
# Input:
#   replication_assessment_PRIMARY.rds
#
# Output:
#   Figure7_empirical_permutation.jpg
#   permutation_results_PRIMARY.rds
#   permutation_summary_PRIMARY.txt
# ============================================================

# ------------------------------------------------------------
# Project folder. Set the AD_PROJECT_ROOT environment variable to the
# folder that contains Data/, RDS/, Output/ and Supplementary_Tables/,
# or start R with that folder as the working directory.
# ------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("AD_PROJECT_ROOT", unset = getwd())

rm(list = ls())

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(patchwork)

# ------------------------------------------------------------
# 2. Settings
# ------------------------------------------------------------

set.seed(61026)

n_perm <- 10000

input_file <- file.path(PROJECT_ROOT, "RDS", "replication_assessment_PRIMARY.rds")

dir.create(file.path(PROJECT_ROOT, "Output", "Figures"), recursive = TRUE, showWarnings = FALSE)
output_figure <- file.path(PROJECT_ROOT, "Output", "Figures", "Fig3_empirical_permutation.jpg")
output_rds <- file.path(PROJECT_ROOT, "RDS", "permutation_results_PRIMARY.rds")
output_summary <- "permutation_summary_PRIMARY.txt"

# ------------------------------------------------------------
# 3. Load replication assessment
# ------------------------------------------------------------

rep <- readRDS(input_file)

stopifnot(
  is.data.frame(rep),
  nrow(rep) == 504
)

# ------------------------------------------------------------
# 4. Define the four replication cohorts
# ------------------------------------------------------------

cohorts <- c(
  "GSE174367",
  "GSE157827",
  "GSE160936",
  "GSE188545"
)

# ------------------------------------------------------------
# 5. Confirm required columns
# ------------------------------------------------------------

required_cols <- c(
  "gene",
  "discovery_direction",
  paste0(cohorts, "_logFC"),
  paste0(cohorts, "_P.Value")
)

missing_cols <- setdiff(required_cols, colnames(rep))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 6. Observed replication counts
#
# Primary criteria:
#
# High-confidence:
#   >= 3 same-direction cohorts
#   AND
#   >= 2 nominal P < 0.05
#
# Moderate:
#   >= 3 same-direction cohorts
#   OR
#   >= 2 nominal P < 0.05
# ------------------------------------------------------------

same_direction_matrix <- as.matrix(
  rep[, paste0(cohorts, "_same_direction")]
)

pvalue_matrix <- as.matrix(
  rep[, paste0(cohorts, "_P.Value")]
)

n_same_direction <- rowSums(
  same_direction_matrix == TRUE,
  na.rm = TRUE
)

n_nominal_p <- rowSums(
  pvalue_matrix < 0.05,
  na.rm = TRUE
)

observed_high <- sum(
  n_same_direction >= 3 &
    n_nominal_p >= 2,
  na.rm = TRUE
)

observed_moderate <- sum(
  (n_same_direction >= 3 |
     n_nominal_p >= 2) &
    !(n_same_direction >= 3 &
        n_nominal_p >= 2),
  na.rm = TRUE
)

observed_combined <- observed_high + observed_moderate

cat("\nObserved results\n")
cat("----------------------------\n")
cat("High-confidence:", observed_high, "\n")
cat("Moderate:", observed_moderate, "\n")
cat("Combined:", observed_combined, "\n\n")

# These should be:
# High = 36
# Moderate = 89
# Combined = 125

stopifnot(observed_high == 36)
stopifnot(observed_moderate == 89)
stopifnot(observed_combined == 125)

# ------------------------------------------------------------
# 7. Prepare cohort-specific data
#
# For each cohort:
#   - retain only genes with observed data
#   - retain the paired logFC and P value
#
# Gene labels are then permuted independently within each
# cohort. This preserves:
#   - the number of available genes
#   - the observed pairing between logFC and P value
# while disrupting the correspondence of gene identities
# across cohorts.
# ------------------------------------------------------------

cohort_data <- list()

for (cohort in cohorts) {
  
  logfc_col <- paste0(cohort, "_logFC")
  p_col <- paste0(cohort, "_P.Value")
  
  tmp <- data.frame(
    gene = rep$gene,
    logFC = rep[[logfc_col]],
    P.Value = rep[[p_col]]
  )
  
  tmp <- tmp[
    !is.na(tmp$logFC) &
      !is.na(tmp$P.Value),
    ,
    drop = FALSE
  ]
  
  cohort_data[[cohort]] <- tmp
  
  cat(
    cohort,
    ":",
    nrow(tmp),
    "genes available\n"
  )
}

# ------------------------------------------------------------
# 8. Permutation function
# ------------------------------------------------------------

run_one_permutation <- function() {
  
  # Start with fixed SEA-AD discovery genes/directions
  discovery_direction <- rep$discovery_direction
  
  # Matrices to hold permuted results
  perm_same <- matrix(
    FALSE,
    nrow = nrow(rep),
    ncol = length(cohorts)
  )
  
  perm_p <- matrix(
    NA_real_,
    nrow = nrow(rep),
    ncol = length(cohorts)
  )
  
  colnames(perm_same) <- cohorts
  colnames(perm_p) <- cohorts
  
  # --------------------------------------------------------
  # Independently permute gene labels within each cohort
  # --------------------------------------------------------
  
  for (j in seq_along(cohorts)) {
    
    cohort <- cohorts[j]
    
    dat <- cohort_data[[cohort]]
    
    # Randomly permute the gene labels
    shuffled_genes <- sample(
      dat$gene,
      size = nrow(dat),
      replace = FALSE
    )
    
    # Match the shuffled gene labels back to the fixed
    # SEA-AD discovery list.
    #
    # The logFC and P value remain paired with one another.
    match_index <- match(
      rep$gene,
      shuffled_genes
    )
    
    available <- !is.na(match_index)
    
    # Direction of the shuffled logFC
    shuffled_logfc <- dat$logFC
    
    # Determine whether the shuffled gene's direction
    # agrees with the fixed SEA-AD discovery direction.
    #
    # Positive logFC = UP
    # Negative logFC = DOWN
    perm_direction <- ifelse(
      shuffled_logfc > 0,
      "UP",
      "DOWN"
    )
    
    perm_same[available, j] <-
      perm_direction[match_index[available]] ==
      discovery_direction[available]
    
    perm_p[available, j] <-
      dat$P.Value[match_index[available]]
  }
  
  # --------------------------------------------------------
  # Apply the exact primary replication rules
  # --------------------------------------------------------
  
  n_same <- rowSums(
    perm_same == TRUE,
    na.rm = TRUE
  )
  
  n_p <- rowSums(
    perm_p < 0.05,
    na.rm = TRUE
  )
  
  high <- sum(
    n_same >= 3 &
      n_p >= 2,
    na.rm = TRUE
  )
  
  combined <- sum(
    n_same >= 3 |
      n_p >= 2,
    na.rm = TRUE
  )
  
  c(
    high = high,
    combined = combined
  )
}

# ------------------------------------------------------------
# 9. Run 10,000 permutations
# ------------------------------------------------------------

cat("\nRunning", n_perm, "permutations...\n")

perm_results <- matrix(
  NA_integer_,
  nrow = n_perm,
  ncol = 2
)

colnames(perm_results) <- c(
  "high",
  "combined"
)

pb <- txtProgressBar(
  min = 0,
  max = n_perm,
  style = 3
)

for (i in seq_len(n_perm)) {
  
  perm_results[i, ] <- run_one_permutation()
  
  if (i %% 100 == 0) {
    setTxtProgressBar(pb, i)
  }
}

close(pb)

perm_results <- as.data.frame(perm_results)

# ------------------------------------------------------------
# 10. Empirical P values
#
# +1 correction:
#
# (number of permutations >= observed + 1) /
# (number of permutations + 1)
# ------------------------------------------------------------

empirical_p_high <-
  (sum(perm_results$high >= observed_high) + 1) /
  (n_perm + 1)

empirical_p_combined <-
  (sum(
    perm_results$combined >= observed_combined
  ) + 1) /
  (n_perm + 1)

# ------------------------------------------------------------
# 11. Summary statistics
# ------------------------------------------------------------

high_mean <- mean(perm_results$high)
high_q95 <- quantile(
  perm_results$high,
  0.95
)

combined_mean <- mean(
  perm_results$combined
)

combined_q95 <- quantile(
  perm_results$combined,
  0.95
)

cat("\nPermutation results\n")
cat("----------------------------\n")

cat(
  "High-confidence observed:",
  observed_high,
  "\n"
)

cat(
  "High-confidence null mean:",
  round(high_mean, 2),
  "\n"
)

cat(
  "High-confidence null 95th percentile:",
  high_q95,
  "\n"
)

cat(
  "High-confidence empirical P:",
  empirical_p_high,
  "\n\n"
)

cat(
  "Combined observed:",
  observed_combined,
  "\n"
)

cat(
  "Combined null mean:",
  round(combined_mean, 2),
  "\n"
)

cat(
  "Combined null 95th percentile:",
  combined_q95,
  "\n"
)

cat(
  "Combined empirical P:",
  empirical_p_combined,
  "\n"
)

# ------------------------------------------------------------
# 12. Save permutation results
# ------------------------------------------------------------

saveRDS(
  perm_results,
  output_rds
)

# ------------------------------------------------------------
# 13. Create plotting data
# ------------------------------------------------------------

plot_high <- data.frame(
  count = perm_results$high
)

plot_combined <- data.frame(
  count = perm_results$combined
)

# ------------------------------------------------------------
# 14. Panel A: High-confidence genes
# ------------------------------------------------------------

pA <- ggplot(
  plot_high,
  aes(x = count)
) +
  geom_histogram(
    binwidth = 1,
    boundary = -0.5,
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = observed_high,
    linewidth = 1
  ) +
  annotate(
    "text",
    x = observed_high,
    y = Inf,
    label = "Observed = 36",
    vjust = 1.5,
    hjust = 1.05,
    size = 4
  ) +
  labs(
    title = "A",
    subtitle = "High-confidence replication",
    x = "Number of genes",
    y = "Permutation count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(
      size = 12
    )
  )

# ------------------------------------------------------------
# 15. Panel B: High + Moderate genes
# ------------------------------------------------------------

pB <- ggplot(
  plot_combined,
  aes(x = count)
) +
  geom_histogram(
    binwidth = 1,
    boundary = -0.5,
    linewidth = 0.2
  ) +
  geom_vline(
    xintercept = observed_combined,
    linewidth = 1
  ) +
  annotate(
    "text",
    x = observed_combined,
    y = Inf,
    label = "Observed = 125",
    vjust = 1.5,
    hjust = 1.05,
    size = 4
  ) +
  labs(
    title = "B",
    subtitle = "High + moderate replication",
    x = "Number of genes",
    y = "Permutation count"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(
      size = 12
    )
  )

# ------------------------------------------------------------
# 16. Combine panels
# ------------------------------------------------------------

figure7 <- pA + pB +
  plot_annotation(
    title = "Empirical permutation analysis of cross-cohort replication",
    subtitle = "10,000 independent gene-label permutations"
  )

# ------------------------------------------------------------
# 17. Save publication-quality JPG
# ------------------------------------------------------------

ggsave(
  filename = output_figure,
  plot = figure7,
  width = 10,
  height = 5.5,
  units = "in",
  dpi = 600,
  quality = 100
)

# ------------------------------------------------------------
# 18. Save text summary
# ------------------------------------------------------------

sink(output_summary)

cat("Figure 7 empirical permutation analysis\n")
cat("========================================\n\n")

cat("Input:\n")
cat(input_file, "\n\n")

cat("Number of permutations:\n")
cat(n_perm, "\n\n")

cat("Observed:\n")
cat("High-confidence =", observed_high, "\n")
cat("Moderate =", observed_moderate, "\n")
cat("Combined =", observed_combined, "\n\n")

cat("High-confidence null:\n")
cat("Mean =", high_mean, "\n")
cat("95th percentile =", high_q95, "\n")
cat("Empirical P =", empirical_p_high, "\n\n")

cat("Combined null:\n")
cat("Mean =", combined_mean, "\n")
cat("95th percentile =", combined_q95, "\n")
cat("Empirical P =", empirical_p_combined, "\n")

sink()

cat("\nDone.\n")
cat("Created:", output_figure, "\n")
cat("Created:", output_rds, "\n")
cat("Created:", output_summary, "\n")