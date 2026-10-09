# ============================================================
# 11_GSE188545_without_cluster18.R
# ============================================================
# Purpose:
#   GSE188545 cluster 18 contains 1,098 nuclei from 7 donors, 94% from
#   AD donors. This sensitivity analysis removes cluster 18, keeps only
#   cluster 7, and checks whether the 48-gene consensus signature still
#   holds.
#
# What it does (nothing in the primary results is overwritten):
#   Step 0  Reproduce the published replication with the ORIGINAL
#           GSE188545 results, as a check that this script matches
#           04_DE_analysis (should give 36 high-confidence genes,
#           38 RankProd genes, and the same 48 genes).
#   Step 1  Rebuild GSE188545 donor pseudobulk from cluster 7 only.
#   Step 2  Rerun GSE188545 DE with the same pipeline as 04_DE_analysis.
#   Step 3  Redo vote-counting and RankProd with the new GSE188545 results.
#   Step 4  Compare with the published 48 genes.
#
# Inputs (already created by 03E and 04):
#   RDS/GSE188545_microglia_PRIMARY_seurat.rds
#   RDS/all_DE_results_PRIMARY.rds
#   RDS/discovery_genes.rds
#   Final tables and figures/Supplementary_Table_2_DE_Consensus_Signatures.xlsx
#
# Outputs (written to Output/Cluster18_check/):
#   GSE188545_cluster7_only_DE.csv
#   GSE188545_cells_per_donor_before_after.csv
#   consensus48_with_and_without_cluster18.csv
#   new_consensus_without_cluster18.csv
#   summary.txt
#
# HOW TO RUN: open in RStudio and click Source. It may take 5-15 minutes
# (loading the Seurat object and RankProd are the slow steps).
# Results are written to Output/Cluster18_check/.
# ============================================================

# ------------------------------------------------------------
# Project folder. Set the AD_PROJECT_ROOT environment variable to the
# folder that contains Data/, RDS/, Output/ and Supplementary_Tables/,
# or start R with that folder as the working directory.
# ------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("AD_PROJECT_ROOT", unset = getwd())

ROOT_DIR <- PROJECT_ROOT
supp2_file <- file.path(PROJECT_ROOT, "Supplementary_Tables", "Supplementary_Table_2_DE_Consensus_Signatures.xlsx")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(readxl)
  library(RankProd)
})

rds_dir <- file.path(ROOT_DIR, "RDS")
out_dir <- file.path(ROOT_DIR, "Output", "Cluster18_check")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_lines <- c()
say <- function(...) { s <- paste0(...); cat(s, "\n"); log_lines <<- c(log_lines, s) }

# ------------------------------------------------------------
# Shared functions (same rules as 04_DE_analysis_AD61026.R)
# ------------------------------------------------------------
run_de <- function(counts, dx_values) {
  dx <- factor(as.character(dx_values), levels = c("Control", "AD"))
  stopifnot(!any(is.na(dx)))
  dge <- DGEList(counts = counts, group = dx)
  design <- model.matrix(~ 0 + dx); colnames(design) <- c("Control", "AD")
  keep <- filterByExpr(dge, design = design)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  fit <- contrasts.fit(fit, makeContrasts(AD_vs_Control = AD - Control, levels = design))
  fit <- eBayes(fit, robust = TRUE)
  res <- topTable(fit, coef = "AD_vs_Control", number = Inf, sort.by = "P", adjust.method = "BH")
  res$gene <- rownames(res)
  res
}

vote_count <- function(disc, rep_results) {
  bind_rows(lapply(seq_len(nrow(disc)), function(i) {
    g <- disc$gene[i]; up <- disc$logFC[i] > 0
    same <- c(); nom <- c()
    for (r in rep_results) {
      row <- r[r$gene == g, ]
      if (nrow(row) > 0) {
        same <- c(same, (row$logFC[1] > 0) == up)
        nom  <- c(nom, row$P.Value[1] < 0.05)
      }
    }
    n_same <- sum(same); n_nom <- sum(nom)
    high <- n_same >= 3 && n_nom >= 2
    mod  <- !high && (n_same >= 3 || n_nom >= 2)
    data.frame(gene = g, n_same_direction = n_same, n_nominal_p_lt_0.05 = n_nom,
               confidence = ifelse(high, "High", ifelse(mod, "Moderate", "Not_replicated")))
  }))
}

rankprod_same_dir <- function(disc, rep_results) {
  common <- sort(Reduce(intersect, lapply(rep_results, function(r) r$gene)))
  fc <- sapply(rep_results, function(r) r$logFC[match(common, r$gene)])
  rownames(fc) <- common
  set.seed(42)
  rp <- RankProducts(data = fc, cl = rep(1, ncol(fc)), logged = TRUE, na.rm = TRUE,
                     plot = FALSE, rand = 42, gene.names = common)
  down <- common[rp$pfp[, 1] < 0.05]
  up   <- common[rp$pfp[, 2] < 0.05]
  disc$gene[(disc$logFC > 0 & disc$gene %in% up) | (disc$logFC < 0 & disc$gene %in% down)]
}

consensus <- function(disc, rep_results) {
  vc <- vote_count(disc, rep_results)
  rp <- rankprod_same_dir(disc, rep_results)
  list(vc = vc, rp = rp,
       high = vc$gene[vc$confidence == "High"],
       moderate = vc$gene[vc$confidence == "Moderate"],
       genes = union(vc$gene[vc$confidence == "High"], rp))
}

# ------------------------------------------------------------
# Load published inputs
# ------------------------------------------------------------
all_de <- readRDS(file.path(rds_dir, "all_DE_results_PRIMARY.rds"))
say("DE result sets found: ", paste(names(all_de), collapse = ", "))
rep_ids <- c("GSE174367", "GSE157827", "GSE160936", "GSE188545")
find_set <- function(id) {
  hit <- grep(id, names(all_de), value = TRUE)
  if (length(hit) != 1) stop("Could not find DE results for ", id, " in all_DE_results_PRIMARY.rds")
  all_de[[hit]]
}
rep_orig <- setNames(lapply(rep_ids, find_set), rep_ids)

disc <- readRDS(file.path(rds_dir, "discovery_genes.rds"))
disc <- as.data.frame(disc)[, c("gene", "logFC")]
say("Discovery genes: ", nrow(disc), " (expected 504)")

published48 <- read_excel(supp2_file, sheet = "Consensus48")$gene
stopifnot(length(published48) == 48)

# ------------------------------------------------------------
# Step 0: reproduce the published result
# ------------------------------------------------------------
say("\nStep 0: reproducing published replication (takes a few minutes)...")
orig <- consensus(disc, rep_orig)
say("  High-confidence: ", length(orig$high), " (published 36)")
say("  Moderate: ", length(orig$moderate), " (published 89)")
say("  RankProd same-direction: ", length(orig$rp), " (published 38)")
say("  Consensus genes: ", length(orig$genes), " (published 48); identical to published list: ",
    setequal(orig$genes, published48))
if (!setequal(orig$genes, published48)) {
  say("  WARNING: could not reproduce the published 48. Results below may not be comparable.")
}

# ------------------------------------------------------------
# Step 1: GSE188545 pseudobulk from cluster 7 only
# ------------------------------------------------------------
say("\nStep 1: loading GSE188545 microglia (clusters 7 + 18)...")
seu <- readRDS(file.path(rds_dir, "GSE188545_microglia_PRIMARY_seurat.rds"))
seu <- tryCatch(JoinLayers(seu, assay = "RNA"), error = function(e) seu)
meta <- seu@meta.data
meta$cell_id <- rownames(meta)
meta$cluster_id <- as.character(meta$seurat_clusters)
say("  Nuclei by cluster: ", paste(names(table(meta$cluster_id)), table(meta$cluster_id), sep = "=", collapse = ", "))

counts_all <- GetAssayData(seu, assay = "RNA", layer = "counts")
keep_cells <- meta$cell_id[meta$cluster_id == "7"]

per_donor <- meta %>%
  group_by(donor_id, diagnosis_std) %>%
  summarise(nuclei_clusters_7_18 = n(), nuclei_cluster_7 = sum(cluster_id == "7"),
            nuclei_cluster_18 = sum(cluster_id == "18"), .groups = "drop")
write.csv(per_donor, file.path(out_dir, "GSE188545_cells_per_donor_before_after.csv"), row.names = FALSE)

donors <- per_donor$donor_id[per_donor$nuclei_cluster_7 >= 10]    # same >= 10 nuclei rule as the paper
dropped <- setdiff(per_donor$donor_id, donors)
say("  Donors kept (>= 10 cluster-7 nuclei): ", length(donors),
    if (length(dropped)) paste0("; dropped: ", paste(dropped, collapse = ", ")) else "")

pb <- sapply(donors, function(d) {
  cells <- intersect(meta$cell_id[meta$donor_id == d], keep_cells)
  Matrix::rowSums(counts_all[, cells, drop = FALSE])
})
dx <- per_donor$diagnosis_std[match(donors, per_donor$donor_id)]
say("  AD donors: ", sum(dx == "AD"), "  Control donors: ", sum(dx == "Control"))
rm(seu, counts_all); gc(verbose = FALSE)

# ------------------------------------------------------------
# Step 2: GSE188545 DE, cluster 7 only
# ------------------------------------------------------------
say("\nStep 2: GSE188545 DE (cluster 7 only)...")
new188545 <- run_de(pb, dx)
write.csv(new188545, file.path(out_dir, "GSE188545_cluster7_only_DE.csv"), row.names = FALSE)

old <- rep_orig$GSE188545
both <- intersect(old$gene, new188545$gene)
say("  Correlation of GSE188545 logFC, with vs without cluster 18 (all genes): ",
    round(cor(old$logFC[match(both, old$gene)], new188545$logFC[match(both, new188545$gene)]), 3))

# ------------------------------------------------------------
# Step 3: replication with the new GSE188545 results
# ------------------------------------------------------------
say("\nStep 3: redoing vote-counting and RankProd (takes a few minutes)...")
rep_new <- rep_orig; rep_new$GSE188545 <- new188545
newc <- consensus(disc, rep_new)
say("  High-confidence: ", length(newc$high), " (was ", length(orig$high), ")")
say("  Moderate: ", length(newc$moderate), " (was ", length(orig$moderate), ")")
say("  RankProd same-direction: ", length(newc$rp), " (was ", length(orig$rp), ")")
say("  Consensus genes: ", length(newc$genes), " (was ", length(orig$genes), ")")

# ------------------------------------------------------------
# Step 4: compare with the published 48
# ------------------------------------------------------------
cmp <- data.frame(gene = published48) %>%
  left_join(disc %>% rename(SEAAD_logFC = logFC), by = "gene") %>%
  left_join(old %>% select(gene, GSE188545_with18_logFC = logFC, GSE188545_with18_P = P.Value), by = "gene") %>%
  left_join(new188545 %>% select(gene, GSE188545_without18_logFC = logFC, GSE188545_without18_P = P.Value), by = "gene") %>%
  left_join(newc$vc %>% select(gene, new_confidence = confidence, new_n_same = n_same_direction,
                               new_n_nominal = n_nominal_p_lt_0.05), by = "gene") %>%
  mutate(
    GSE188545_same_direction_without18 = sign(GSE188545_without18_logFC) == sign(SEAAD_logFC),
    new_rankprod_same_direction = gene %in% newc$rp,
    still_in_consensus = gene %in% newc$genes
  )
write.csv(cmp, file.path(out_dir, "consensus48_with_and_without_cluster18.csv"), row.names = FALSE)

newlist <- data.frame(gene = newc$genes) %>%
  left_join(disc, by = "gene") %>%
  mutate(in_published48 = gene %in% published48,
         high_confidence = gene %in% newc$high,
         rankprod_same_direction = gene %in% newc$rp)
write.csv(newlist, file.path(out_dir, "new_consensus_without_cluster18.csv"), row.names = FALSE)

say("\nStep 4: the published 48 genes without cluster 18")
say("  GSE188545 same direction as SEA-AD: ",
    sum(cmp$GSE188545_same_direction_without18, na.rm = TRUE), " / ",
    sum(!is.na(cmp$GSE188545_without18_logFC)), " detected (with cluster 18: ",
    sum(sign(cmp$GSE188545_with18_logFC) == sign(cmp$SEAAD_logFC), na.rm = TRUE), ")")
say("  GSE188545 same direction and P < 0.05: ",
    sum(cmp$GSE188545_same_direction_without18 & cmp$GSE188545_without18_P < 0.05, na.rm = TRUE),
    " (with cluster 18: ",
    sum(sign(cmp$GSE188545_with18_logFC) == sign(cmp$SEAAD_logFC) & cmp$GSE188545_with18_P < 0.05, na.rm = TRUE), ")")
say("  Still meet consensus criteria: ", sum(cmp$still_in_consensus), " / 48")
say("  Still same direction in >= 3 of 4 cohorts: ", sum(cmp$new_n_same >= 3, na.rm = TRUE), " / 48")
say("  New genes entering the consensus: ", sum(!newlist$in_published48),
    if (any(!newlist$in_published48)) paste0(" (", paste(newlist$gene[!newlist$in_published48], collapse = ", "), ")") else "")

writeLines(log_lines, file.path(out_dir, "summary.txt"))
cat("\nDone. Results written to:", out_dir, "\n")
