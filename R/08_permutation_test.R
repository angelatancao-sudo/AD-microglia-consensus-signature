# ============================================================
# 08_permutation_test.R
#
# Reviewer-requested empirical null analysis
#
# Purpose:
#   Test whether the observed cross-cohort replication
#   concordance is greater than expected by chance.
#
# The analysis preserves:
#   1. The 504 SEA-AD discovery genes
#   2. SEA-AD discovery directions
#   3. The observed gene availability pattern in each
#      replication cohort
#   4. The observed logFC/P-value relationship within
#      each replication cohort
#
# Gene identities are permuted independently within
# each replication cohort.
#
# Original replication criteria:
#   High:
#       >=3 same-direction cohorts AND
#       >=2 cohorts with nominal P < 0.05
#
#   Moderate:
#       not High AND
#       >=3 same-direction OR
#       >=2 nominal P < 0.05
#
# ============================================================

# ------------------------------------------------------------
# Project folder. Set the AD_PROJECT_ROOT environment variable to the
# folder that contains Data/, RDS/, Output/ and Supplementary_Tables/,
# or start R with that folder as the working directory.
# ------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("AD_PROJECT_ROOT", unset = getwd())

set.seed(20260829)

# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(dplyr)

# ------------------------------------------------------------
# 2. Load original replication assessment
# ------------------------------------------------------------

RDS_DIR <- file.path(PROJECT_ROOT, "RDS")

rep <- readRDS(
  file.path(RDS_DIR, "replication_assessment_PRIMARY.rds")
)

# ------------------------------------------------------------
# 3. Confirm the expected structure
# ------------------------------------------------------------

expected_genes <- 504

if (nrow(rep) != expected_genes) {
  stop(
    paste0(
      "Expected 504 discovery genes, but found ",
      nrow(rep),
      "."
    )
  )
}

replication_datasets <- c(
  "GSE174367",
  "GSE157827",
  "GSE160936",
  "GSE188545"
)

# ------------------------------------------------------------
# 4. Extract discovery information
# ------------------------------------------------------------

discovery_genes <- rep$gene

discovery_direction <- rep$discovery_direction

names(discovery_direction) <- discovery_genes

# ------------------------------------------------------------
# 5. Build cohort-specific data frames
#
# Each cohort contains:
#   gene
#   logFC
#   P value
#
# Only genes actually available in that cohort are included.
# ------------------------------------------------------------

cohort_data <- list()

for (ds in replication_datasets) {
  
  logfc_col <- paste0(ds, "_logFC")
  p_col     <- paste0(ds, "_P.Value")
  
  tmp <- rep %>%
    select(
      gene,
      all_of(logfc_col),
      all_of(p_col)
    ) %>%
    rename(
      logFC = all_of(logfc_col),
      P.Value = all_of(p_col)
    ) %>%
    filter(
      !is.na(logFC),
      !is.na(P.Value)
    )
  
  cohort_data[[ds]] <- tmp
}

# ------------------------------------------------------------
# 6. Check availability
# ------------------------------------------------------------

cat("============================================\n")
cat("Replication dataset availability\n")
cat("============================================\n\n")

for (ds in replication_datasets) {
  
  cat(
    ds,
    ":",
    nrow(cohort_data[[ds]]),
    "of",
    length(discovery_genes),
    "genes available\n"
  )
}

# ------------------------------------------------------------
# 7. Function to classify one permutation
# ------------------------------------------------------------

calculate_permutation <- function() {
  
  # Initialize counts for each of the 504 discovery genes
  n_same <- setNames(
    integer(length(discovery_genes)),
    discovery_genes
  )
  
  n_nominal <- setNames(
    integer(length(discovery_genes)),
    discovery_genes
  )
  
  # ----------------------------------------------------------
  # Independently permute gene labels in each cohort
  # ----------------------------------------------------------
  
  for (ds in replication_datasets) {
    
    dat <- cohort_data[[ds]]
    
    # Randomly permute the gene labels while keeping
    # logFC and P value paired.
    shuffled_genes <- sample(dat$gene)
    
    # Create lookup table:
    # original discovery gene <- randomly assigned
    # replication result
    permuted <- data.frame(
      gene = shuffled_genes,
      logFC = dat$logFC,
      P.Value = dat$P.Value,
      stringsAsFactors = FALSE
    )
    
    # Keep only genes that belong to the 504 discovery set
    permuted <- permuted[
      permuted$gene %in% discovery_genes,
      ,
      drop = FALSE
    ]
    
    # Determine expected discovery direction
    expected_dir <- discovery_direction[permuted$gene]
    
    # Determine replication direction
    replication_dir <- ifelse(
      permuted$logFC > 0,
      "UP",
      ifelse(
        permuted$logFC < 0,
        "DOWN",
        NA_character_
      )
    )
    
    same <- !is.na(replication_dir) &
      replication_dir == expected_dir
    
    nominal <- !is.na(permuted$P.Value) &
      permuted$P.Value < 0.05
    
    # Add counts
    n_same[permuted$gene] <-
      n_same[permuted$gene] + as.integer(same)
    
    n_nominal[permuted$gene] <-
      n_nominal[permuted$gene] + as.integer(nominal)
  }
  
  # ----------------------------------------------------------
  # Apply the ORIGINAL classification rules
  # ----------------------------------------------------------
  
  high <- n_same >= 3 & n_nominal >= 2
  
  moderate <- !high &
    (n_same >= 3 | n_nominal >= 2)
  
  data.frame(
    n_high = sum(high),
    n_moderate = sum(moderate),
    n_high_moderate = sum(high | moderate),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# 8. Run permutations
# ------------------------------------------------------------

n_permutations <- 10000

cat("\n============================================\n")
cat("Running", n_permutations, "permutations\n")
cat("============================================\n\n")

permutation_results <- vector(
  "list",
  n_permutations
)

for (i in seq_len(n_permutations)) {
  
  permutation_results[[i]] <-
    calculate_permutation()
  
  if (i %% 1000 == 0) {
    cat(
      "Completed:",
      i,
      "/",
      n_permutations,
      "\n"
    )
  }
}

permutation_results <- bind_rows(
  permutation_results
)

# ------------------------------------------------------------
# 9. Observed values
# ------------------------------------------------------------

observed_high <- sum(
  rep$confidence == "High",
  na.rm = TRUE
)

observed_moderate <- sum(
  rep$confidence == "Moderate",
  na.rm = TRUE
)

observed_high_moderate <-
  observed_high + observed_moderate

# ------------------------------------------------------------
# 10. Empirical P values
#
# +1 correction avoids zero empirical P values.
# ------------------------------------------------------------

empirical_p_high <-
  (sum(
    permutation_results$n_high >= observed_high
  ) + 1) /
  (n_permutations + 1)

empirical_p_high_moderate <-
  (sum(
    permutation_results$n_high_moderate >=
      observed_high_moderate
  ) + 1) /
  (n_permutations + 1)

# ------------------------------------------------------------
# 11. Null distribution summaries
# ------------------------------------------------------------

summary_results <- data.frame(
  
  metric = c(
    "Observed_high_confidence",
    "Null_mean_high_confidence",
    "Null_SD_high_confidence",
    "Null_95th_percentile_high_confidence",
    "Empirical_P_high_confidence",
    
    "Observed_high_plus_moderate",
    "Null_mean_high_plus_moderate",
    "Null_SD_high_plus_moderate",
    "Null_95th_percentile_high_plus_moderate",
    "Empirical_P_high_plus_moderate"
  ),
  
  value = c(
    
    observed_high,
    
    mean(permutation_results$n_high),
    sd(permutation_results$n_high),
    quantile(
      permutation_results$n_high,
      0.95
    ),
    empirical_p_high,
    
    observed_high_moderate,
    
    mean(permutation_results$n_high_moderate),
    sd(permutation_results$n_high_moderate),
    quantile(
      permutation_results$n_high_moderate,
      0.95
    ),
    empirical_p_high_moderate
  )
)

# ------------------------------------------------------------
# 12. Print results
# ------------------------------------------------------------

cat("\n============================================\n")
cat("OBSERVED RESULTS\n")
cat("============================================\n\n")

cat(
  "High-confidence genes:",
  observed_high,
  "\n"
)

cat(
  "High + Moderate genes:",
  observed_high_moderate,
  "\n"
)

cat("\n============================================\n")
cat("EMPIRICAL NULL RESULTS\n")
cat("============================================\n\n")

cat(
  "Mean null High:",
  round(mean(permutation_results$n_high), 2),
  "\n"
)

cat(
  "95th percentile null High:",
  quantile(
    permutation_results$n_high,
    0.95
  ),
  "\n"
)

cat(
  "Empirical P for >= observed High:",
  empirical_p_high,
  "\n\n"
)

cat(
  "Mean null High + Moderate:",
  round(
    mean(permutation_results$n_high_moderate),
    2
  ),
  "\n"
)

cat(
  "95th percentile null High + Moderate:",
  quantile(
    permutation_results$n_high_moderate,
    0.95
  ),
  "\n"
)

cat(
  "Empirical P for >= observed High + Moderate:",
  empirical_p_high_moderate,
  "\n\n"
)

# ------------------------------------------------------------
# 13. Save outputs
# ------------------------------------------------------------

OUTPUT_DIR <- file.path(
  dirname(RDS_DIR),
  "Permutation"
)

dir.create(
  OUTPUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  permutation_results,
  file.path(
    OUTPUT_DIR,
    "primary_replication_permutation_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  summary_results,
  file.path(
    OUTPUT_DIR,
    "primary_replication_permutation_summary.csv"
  ),
  row.names = FALSE
)

saveRDS(
  permutation_results,
  file.path(
    OUTPUT_DIR,
    "primary_replication_permutation_results.rds"
  )
)

saveRDS(
  summary_results,
  file.path(
    OUTPUT_DIR,
    "primary_replication_permutation_summary.rds"
  )
)

cat("\n============================================\n")
cat("Permutation analysis complete.\n")
cat("Outputs saved to:\n")
cat(OUTPUT_DIR, "\n")
cat("============================================\n")