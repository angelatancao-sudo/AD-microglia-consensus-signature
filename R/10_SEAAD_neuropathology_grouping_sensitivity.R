# ============================================================
# 10_SEAAD_neuropathology_grouping_sensitivity.R
# ============================================================
# Purpose:
#   The primary SEA-AD analysis groups donors by cognitive status
#   (Dementia = "AD", No dementia = "Control"). This script checks
#   whether the 48-gene consensus signature holds when SEA-AD donors
#   are grouped by AD neuropathology instead.
#
#   Part A: cross-tabulate cognitive status against the SEA-AD
#           "Overall AD neuropathological Change" (ADNC) score.
#   Part B: rerun SEA-AD DE with a pathology-based grouping:
#             AD pathology = ADNC Intermediate or High
#             Control      = ADNC Not AD or Low
#           using the same pipeline as 04_DE_analysis
#           (filterByExpr, TMM, voom, lmFit, eBayes robust = TRUE),
#           then compare the 48 consensus genes with the primary model.
#
# Inputs (already created by earlier scripts, nothing new to download):
#   RDS/SEAAD_pseudobulk_counts.rds
#   RDS/SEAAD_pseudobulk_meta.rds
#   Supplementary_Table_2_DE_Consensus_Signatures.xlsx (sheet "Consensus48")
#
# Outputs (written to Output/Pathology_check/):
#   A_crosstab_cognitive_vs_ADNC.csv
#   B_pathology_DE_all_genes.csv
#   B_consensus48_primary_vs_pathology.csv
#   B_summary.txt
#
# HOW TO RUN: open in RStudio and click Source.
# Results are written to Output/Pathology_check/.
# ============================================================

# ------------------------------------------------------------
# Project folder. Set the AD_PROJECT_ROOT environment variable to the
# folder that contains Data/, RDS/, Output/ and Supplementary_Tables/,
# or start R with that folder as the working directory.
# ------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("AD_PROJECT_ROOT", unset = getwd())

ROOT_DIR <- PROJECT_ROOT

supp2_file <- file.path(PROJECT_ROOT, "Supplementary_Tables", "Supplementary_Table_2_DE_Consensus_Signatures.xlsx")

library(edgeR)
library(limma)
library(dplyr)
library(readxl)

rds_dir <- file.path(ROOT_DIR, "RDS")
out_dir <- file.path(ROOT_DIR, "Output", "Pathology_check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- readRDS(file.path(rds_dir, "SEAAD_pseudobulk_counts.rds"))
meta   <- readRDS(file.path(rds_dir, "SEAAD_pseudobulk_meta.rds"))
meta   <- as.data.frame(meta)

# Align metadata to the count matrix (same as SEAAD sensitivity.R)
meta <- meta[match(colnames(counts), meta[["Donor.ID"]]), ]
stopifnot(all(meta[["Donor.ID"]] == colnames(counts)))
cat("Donors:", ncol(counts), " Genes:", nrow(counts), "\n\n")

# ------------------------------------------------------------
# Find the ADNC column. Its exact name depends on how the H5AD and
# the Excel metadata were merged, so search by name.
# ------------------------------------------------------------
adnc_candidates <- grep("neuropath|ADNC", names(meta), ignore.case = TRUE, value = TRUE)
cat("Possible ADNC columns found:\n"); print(adnc_candidates)

if (length(adnc_candidates) == 0) {
  cat("\nAll metadata columns:\n"); print(names(meta))
  stop("No ADNC column found. Check the metadata column names printed above.")
}

# Prefer the column whose values look like Not AD / Low / Intermediate / High
adnc_col <- NA
for (cc in adnc_candidates) {
  vals <- unique(trimws(as.character(meta[[cc]])))
  if (any(vals %in% c("High", "Intermediate")) && any(vals %in% c("Not AD", "Low"))) {
    adnc_col <- cc; break
  }
}
if (is.na(adnc_col)) {
  for (cc in adnc_candidates) { cat("\n", cc, ":\n"); print(table(meta[[cc]], useNA = "ifany")) }
  stop("Could not recognise ADNC values. Check the value tables printed above.")
}
cat("\nUsing ADNC column:", adnc_col, "\n")

adnc <- factor(trimws(as.character(meta[[adnc_col]])),
               levels = c("Not AD", "Low", "Intermediate", "High"))

# ------------------------------------------------------------
# PART A: cognitive status x ADNC
# ------------------------------------------------------------
tab <- table(Cognitive = meta[["diagnosis_std"]], ADNC = adnc, useNA = "ifany")
cat("\n=== PART A: cognitive status (rows) x ADNC (columns) ===\n")
print(tab)
write.csv(as.data.frame.matrix(tab), file.path(out_dir, "A_crosstab_cognitive_vs_ADNC.csv"))

# ------------------------------------------------------------
# PART B: DE with pathology grouping
# ------------------------------------------------------------
path_group <- ifelse(adnc %in% c("Intermediate", "High"), "AD",
              ifelse(adnc %in% c("Not AD", "Low"), "Control", NA))
keep_donor <- !is.na(path_group)
cat("\n=== PART B: pathology grouping ===\n")
print(table(path_group, useNA = "ifany"))

run_de <- function(cnt, group) {
  dx <- factor(group, levels = c("Control", "AD"))
  design <- model.matrix(~ 0 + dx); colnames(design) <- levels(dx)
  dge <- DGEList(cnt)
  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  fit <- contrasts.fit(fit, makeContrasts(AD - Control, levels = design))
  fit <- eBayes(fit, robust = TRUE)
  res <- topTable(fit, number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  res
}

# Primary model, recomputed here so both use the same code
primary   <- run_de(counts, as.character(meta[["diagnosis_std"]]))
pathology <- run_de(counts[, keep_donor], path_group[keep_donor])
write.csv(pathology, file.path(out_dir, "B_pathology_DE_all_genes.csv"), row.names = FALSE)

# Sanity check: recomputed primary should match the published discovery logFC
genes48_tab <- read_excel(supp2_file, sheet = "Consensus48")
genes48 <- genes48_tab$gene
stopifnot(length(genes48) == 48)

cmp <- data.frame(gene = genes48) %>%
  left_join(genes48_tab %>% select(gene, published_logFC = discovery_logFC), by = "gene") %>%
  left_join(primary   %>% select(gene, primary_logFC = logFC, primary_FDR = adj.P.Val), by = "gene") %>%
  left_join(pathology %>% select(gene, pathology_logFC = logFC, pathology_P = P.Value,
                                 pathology_FDR = adj.P.Val), by = "gene") %>%
  mutate(same_direction = sign(primary_logFC) == sign(pathology_logFC),
         retained_fraction = pathology_logFC / primary_logFC)
write.csv(cmp, file.path(out_dir, "B_consensus48_primary_vs_pathology.csv"), row.names = FALSE)

max_diff <- max(abs(cmp$published_logFC - cmp$primary_logFC), na.rm = TRUE)
n_det    <- sum(!is.na(cmp$pathology_logFC))
n_same   <- sum(cmp$same_direction, na.rm = TRUE)
n_p05    <- sum(cmp$same_direction & cmp$pathology_P < 0.05, na.rm = TRUE)
n_half   <- sum(cmp$retained_fraction >= 0.5, na.rm = TRUE)
r_48     <- cor(cmp$primary_logFC, cmp$pathology_logFC, use = "complete.obs")
common   <- intersect(primary$gene, pathology$gene)
r_all    <- cor(primary[match(common, primary$gene), "logFC"],
                pathology[match(common, pathology$gene), "logFC"])

summary_lines <- c(
  paste("ADNC column used:", adnc_col),
  "Cognitive status x ADNC:", capture.output(print(tab)), "",
  paste("Pathology grouping: AD =", sum(path_group == "AD", na.rm = TRUE),
        " Control =", sum(path_group == "Control", na.rm = TRUE),
        " excluded (missing ADNC) =", sum(!keep_donor)),
  paste("Check: max |published - recomputed primary logFC| among 48 genes =", signif(max_diff, 3),
        "(should be ~0; a large value means the primary model was not reproduced)"),
  paste("Consensus genes detected in pathology model:", n_det, "/ 48"),
  paste("Same direction as primary:", n_same, "/", n_det),
  paste("Same direction and P < 0.05:", n_p05, "/", n_det),
  paste("Retained >= 50% of primary logFC:", n_half, "/", n_det),
  paste("Correlation of logFC, 48 genes:", round(r_48, 3)),
  paste("Correlation of logFC, all genes:", round(r_all, 3))
)
writeLines(summary_lines, file.path(out_dir, "B_summary.txt"))
cat("\n", paste(summary_lines, collapse = "\n"), "\n")
cat("\nDone. Results written to:", out_dir, "\n")
