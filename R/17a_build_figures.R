# ============================================================
# 17a_build_figures.R
# ============================================================
# Project:
#   AD61026 Alzheimer disease microglia multi-dataset pseudobulk study
#
# Purpose:
#   Build final manuscript figures as publication-ready JPG files.
#
# Final figure terminology:
#   Use "consensus signature" for the final 48-gene signature.
#   Use "replication-supported" for broader high + moderate evidence.
#   Do not call drug-pair supporting analyses validation.
#
# Important updates:
#   1. Forces the correct Output folder.
#   2. Removes titles/subtitles from all figure images.
#   3. Rebuilds Figure 1 as a cleaner workflow diagram with fewer words.
#   4. Uses Figure 2 as the SEA-AD volcano plot.
#   5. Fixes Figure 3 dataset labels, including SEA-AD and GSE labels.
#   6. Fixes Figure 5 pathway direction colors.
#   7. Builds Figure 6 using flexible column detection.
#   8. Shortens Figure 7 PK/DDI labels for readability.
#
# Output folder:
#   <project folder>/Output/Manuscript_Final_Figures_JPG_AD61026
#
# Final JPG figures:
#   Figure_1__study_workflow.jpg
#   Figure_2__SEAAD_discovery_volcano.jpg
#   Figure_3__strict_consensus_signature_heatmap.jpg
#   Figure_4__strict_48_consensus_signature_direction.jpg
#   Figure_5__pathway_enrichment_summary.jpg
#   Figure_6__single_drug_evidence_matrix.jpg
#   Figure_7__drug_pair_support_matrix.jpg
#
# Recommended run:
#   source("R/11_build_manuscript_figures_JPG_AD61026.R")
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

# Sections 4, 6, 7, 8 and 10 build earlier versions of the workflow, heatmap,
# signature-direction and pathway figures and a drug-pair figure. None of these
# is in the published paper (see 17b-17d for the published Figs 3-5; Fig 1 is a
# diagram). They are kept for transparency and skipped by default.
RUN_SUPERSEDED_FIGURES <- FALSE

BASE_DIR <- PROJECT_ROOT
OUT_DIR  <- file.path(BASE_DIR, "Output")

FINAL_FIG_DIR <- file.path(OUT_DIR, "Manuscript_Final_Figures_JPG_AD61026")
FINAL_FIG_JPG_DIR <- file.path(FINAL_FIG_DIR, "JPG")
FINAL_FIG_PDF_DIR <- file.path(FINAL_FIG_DIR, "PDF_internal_check")
FINAL_FIG_QC_DIR <- file.path(FINAL_FIG_DIR, "QC")
FINAL_FIG_SOURCE_DIR <- file.path(FINAL_FIG_DIR, "Figure_Source_Data_Used")

dir.create(FINAL_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FINAL_FIG_JPG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FINAL_FIG_PDF_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FINAL_FIG_QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FINAL_FIG_SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

required_pkgs <- c(
  "dplyr", "tidyr", "stringr", "tibble", "readr",
  "ggplot2", "forcats", "scales", "purrr"
)

missing_pkgs <- required_pkgs[
  !sapply(required_pkgs, requireNamespace, quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them before running this script, for example:\n",
    "install.packages(c(",
    paste0("'", missing_pkgs, "'", collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(scales)
  library(purrr)
})

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

cat("============================================================\n")
cat("  11: Build manuscript figures JPG\n")
cat("============================================================\n\n")

cat("Project folder:\n")
cat("  ", BASE_DIR, "\n\n")

cat("Figure output folder:\n")
cat("  ", FINAL_FIG_DIR, "\n\n")


# ============================================================
# SECTION 2: Global figure options
# ============================================================

FIG_DPI <- 600

BASE_FONT_SIZE <- 12
AXIS_TEXT_SIZE <- 10
AXIS_TITLE_SIZE <- 11
LEGEND_TEXT_SIZE <- 9

DISCOVERY_FDR_CUTOFF <- 0.05
DISCOVERY_ABS_LOGFC_CUTOFF <- 0.25


# ============================================================
# SECTION 3: Helper functions
# ============================================================

theme_ad <- function(base_size = BASE_FONT_SIZE) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = AXIS_TITLE_SIZE),
      axis.text = element_text(size = AXIS_TEXT_SIZE),
      strip.text = element_text(size = AXIS_TEXT_SIZE, face = "bold"),
      legend.title = element_text(size = LEGEND_TEXT_SIZE, face = "bold"),
      legend.text = element_text(size = LEGEND_TEXT_SIZE),
      panel.grid.minor = element_blank(),
      plot.margin = margin(8, 10, 8, 10)
    )
}

clean_label <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "_", " ")
  x <- str_to_title(x)
  
  # Dataset labels
  x <- str_replace_all(x, "\\bSeaad\\b", "SEA-AD")
  x <- str_replace_all(x, "\\bSea Ad\\b", "SEA-AD")
  x <- str_replace_all(x, "\\bGse", "GSE")
  x <- str_replace_all(x, "\\bMean Replication\\b", "Mean replication")
  
  # Common biology labels
  x <- str_replace_all(x, "\\bUv\\b", "UV")
  x <- str_replace_all(x, "\\bMyc\\b", "MYC")
  x <- str_replace_all(x, "\\bMtorc1\\b", "mTORC1")
  x <- str_replace_all(x, "\\bMapk\\b", "MAPK")
  x <- str_replace_all(x, "\\bOxphos\\b", "OXPHOS")
  x <- str_replace_all(x, "\\bDgidb\\b", "DGIdb")
  x <- str_replace_all(x, "\\bCmap\\b", "CMap")
  x <- str_replace_all(x, "\\bIlincs\\b", "iLINCS")
  x <- str_replace_all(x, "\\bDreimt\\b", "DREIMT")
  x <- str_replace_all(x, "\\bPk Ddi\\b", "PK/DDI")
  x <- str_replace_all(x, "\\bParkinsons\\b", "Parkinson's")
  x <- str_replace_all(x, "\\bAlzheimers\\b", "Alzheimer's")
  x <- str_replace_all(x, "\\bHuntingtons\\b", "Huntington's")
  
  x
}

drug_display <- function(x) {
  x <- as.character(x)
  case_when(
    x == "alda-1" ~ "ALDA-1",
    x == "vx-745" ~ "VX-745",
    x == "valproic acid" ~ "Valproic acid",
    x == "atorvastatin" ~ "Atorvastatin",
    x == "cediranib" ~ "Cediranib",
    x == "cilostazol" ~ "Cilostazol",
    x == "memantine" ~ "Memantine",
    x == "metformin" ~ "Metformin",
    x == "prednisone" ~ "Prednisone",
    x == "tipifarnib" ~ "Tipifarnib",
    TRUE ~ str_to_title(x)
  )
}

flatten_candidates <- function(file_candidates) {
  out <- unlist(file_candidates, recursive = TRUE, use.names = FALSE)
  out <- as.character(out)
  out <- out[!is.na(out) & out != ""]
  unique(out)
}

get_all_search_files <- function() {
  search_roots <- c(
    file.path(OUT_DIR, "Manuscript_Display_Index_AD61026"),
    file.path(OUT_DIR, "Manuscript_Display_Index_AD61026", "Figures"),
    file.path(OUT_DIR, "Manuscript_Source_Data"),
    file.path(OUT_DIR, "Manuscript_Table_Figure_Source_Record_AD61026"),
    file.path(OUT_DIR, "Manuscript_Final_Tables_AD61026"),
    file.path(OUT_DIR, "DE"),
    file.path(OUT_DIR, "Pathway"),
    file.path(OUT_DIR, "Drug_Prioritization"),
    file.path(OUT_DIR, "Drug_Prioritization", "processed_drug_results"),
    file.path(OUT_DIR, "QC"),
    OUT_DIR,
    BASE_DIR
  )
  
  search_roots <- unique(search_roots[file.exists(search_roots)])
  
  all_files <- unlist(lapply(search_roots, function(root) {
    list.files(root, recursive = TRUE, full.names = TRUE)
  }))
  
  all_files <- unique(all_files)
  all_files <- all_files[file.exists(all_files)]
  all_files <- all_files[file.info(all_files)$isdir == FALSE]
  
  all_files
}

ALL_SEARCH_FILES <- get_all_search_files()
ALL_SEARCH_BASENAMES <- basename(ALL_SEARCH_FILES)

cat("Searchable files indexed:", length(ALL_SEARCH_FILES), "\n\n")

find_source_file <- function(file_candidates) {
  candidates <- flatten_candidates(file_candidates)
  
  if (length(candidates) == 0 || length(ALL_SEARCH_FILES) == 0) {
    return(NA_character_)
  }
  
  for (candidate in candidates) {
    
    exact_hits <- ALL_SEARCH_FILES[ALL_SEARCH_BASENAMES == candidate]
    if (length(exact_hits) > 0) {
      return(exact_hits[order(file.info(exact_hits)$mtime, decreasing = TRUE)][1])
    }
    
    suffix <- paste0("__", candidate)
    suffix_hits <- ALL_SEARCH_FILES[endsWith(ALL_SEARCH_BASENAMES, suffix)]
    if (length(suffix_hits) > 0) {
      return(suffix_hits[order(file.info(suffix_hits)$mtime, decreasing = TRUE)][1])
    }
    
    loose_hits <- ALL_SEARCH_FILES[grepl(candidate, ALL_SEARCH_BASENAMES, fixed = TRUE)]
    if (length(loose_hits) > 0) {
      return(loose_hits[order(file.info(loose_hits)$mtime, decreasing = TRUE)][1])
    }
  }
  
  NA_character_
}

read_source <- function(file_candidates, required = FALSE) {
  path <- find_source_file(file_candidates)
  
  if (is.na(path) || !file.exists(path)) {
    if (required) {
      warning("Required source not found: ", paste(flatten_candidates(file_candidates), collapse = "; "))
    }
    return(list(data = tibble(), path = NA_character_))
  }
  
  ext <- tolower(tools::file_ext(path))
  
  df <- tryCatch({
    if (ext == "csv") {
      readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
    } else if (ext %in% c("tsv", "txt")) {
      readr::read_tsv(path, show_col_types = FALSE, guess_max = 100000)
    } else {
      tibble(unsupported_file_type = basename(path))
    }
  }, error = function(e) {
    tibble(read_error = as.character(e$message), source_file = basename(path))
  })
  
  list(data = as_tibble(df), path = path)
}

copy_source_used <- function(path, figure_id) {
  if (!is.na(path) && file.exists(path)) {
    out <- file.path(FINAL_FIG_SOURCE_DIR, paste0(figure_id, "__", basename(path)))
    file.copy(path, out, overwrite = TRUE)
  }
}

save_figure <- function(plot, figure_id, width, height) {
  # Figure titles and subtitles are removed from image files.
  # Use manuscript captions for titles.
  
  plot <- plot +
    labs(title = NULL, subtitle = NULL) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )
  
  jpg_path <- file.path(FINAL_FIG_JPG_DIR, paste0(figure_id, ".jpg"))
  pdf_path <- file.path(FINAL_FIG_PDF_DIR, paste0(figure_id, ".pdf"))
  
  ggplot2::ggsave(
    filename = jpg_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = FIG_DPI,
    device = "jpeg",
    bg = "white"
  )
  
  ggplot2::ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = "pdf",
    bg = "white"
  )
  
  tibble(
    figure_id = figure_id,
    jpg_path = normalizePath(jpg_path, winslash = "/", mustWork = FALSE),
    pdf_path = normalizePath(pdf_path, winslash = "/", mustWork = FALSE),
    width_in = width,
    height_in = height,
    dpi = FIG_DPI,
    status = "CREATED"
  )
}

pick_col <- function(df, candidates) {
  candidates <- flatten_candidates(candidates)
  
  exact <- candidates[candidates %in% colnames(df)]
  if (length(exact) > 0) return(exact[1])
  
  lower_cols <- tolower(colnames(df))
  lower_cand <- tolower(candidates)
  idx <- match(lower_cand, lower_cols)
  idx <- idx[!is.na(idx)]
  
  if (length(idx) > 0) return(colnames(df)[idx[1]])
  
  NA_character_
}

p_to_neglog10 <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  if (all(is.na(p))) return(rep(NA_real_, length(p)))
  if (any(p <= 0, na.rm = TRUE)) {
    min_positive <- min(p[p > 0], na.rm = TRUE)
    p[p <= 0] <- min_positive / 10
  }
  -log10(p)
}

figure_log <- list()

add_log <- function(figure_id, status, note, source_path = NA_character_) {
  figure_log[[length(figure_log) + 1]] <<- tibble(
    figure_id = figure_id,
    status = status,
    note = note,
    source_path = source_path
  )
}


# ============================================================
# SECTION 4: Figure 1 — Study workflow
# ============================================================
if (RUN_SUPERSEDED_FIGURES) {

cat("Building Figure 1...\n")

# ------------------------------------------------------------
# Clean workflow figure.
# Uses large font and fewer words.
# Current final numbers:
#   504 SEA-AD discovery genes
#   48 strict consensus genes
#   117 expanded non-MT drug-query genes
# ------------------------------------------------------------

workflow_boxes <- tibble(
  id = c(
    "Discovery",
    "Replication",
    "Step1",
    "DiscoveryGenes",
    "Step2",
    "StrictSignature",
    "Step3",
    "Programs",
    "Step4",
    "Pairs"
  ),
  label = c(
    "Discovery\nSEA-AD\n84 donors\n42 AD / 42 control",
    "Replication cohorts\n4 primary datasets\n62 donors\n+ 1 supplementary dataset",
    "Step 1: Pseudobulk differential expression\nedgeR / limma-voom  |  donor-level analysis",
    "504 discovery genes\nFDR < 0.05\n|logFC| > 0.25",
    "Step 2: Consensus signature construction\nVote-counting replication  |  RankProd meta-analysis",
    "Strict consensus\n48 genes\n24 up / 24 down\nExpanded 117 non-MT query",
    "Step 3: Pathway enrichment\nHallmark  |  KEGG  |  GO BP",
    "5 biological programs\nimmune  |  mitochondrial  |  stress/proteostasis\ncytoskeletal remodeling  |  lipid/receptor regulation",
    "Step 4: Drug prioritization\nCLUE/CMap  |  iLINCS  |  DREIMT  |  DGIdb",
    "Drug-pair hypotheses\nstrict rule-based prioritization\npair-level supporting analyses"
  ),
  xmin = c(0.5, 3.0, 0.5, 0.5, 3.0, 6.7, 0.5, 0.8, 5.1, 5.4),
  xmax = c(2.4, 8.9, 8.9, 2.4, 8.9, 8.9, 4.4, 4.1, 8.9, 8.6),
  ymin = c(7.4, 7.4, 5.9, 4.6, 4.2, 2.7, 1.4, 0.3, 1.4, 0.3),
  ymax = c(8.7, 8.7, 6.7, 5.4, 5.2, 3.7, 2.3, 1.0, 2.3, 1.0),
  fill = c(
    "#DDEBF7",
    "#E2F0D9",
    "#F2F2F2",
    "#DDEBF7",
    "#E2F0D9",
    "#E2F0D9",
    "#FCE4D6",
    "#FCE4D6",
    "#F4CCCC",
    "#F4CCCC"
  ),
  border = c(
    "#2F75B5",
    "#70AD47",
    "#7F7F7F",
    "#2F75B5",
    "#70AD47",
    "#70AD47",
    "#ED7D31",
    "#ED7D31",
    "#C00000",
    "#C00000"
  ),
  text_size = c(4.3, 4.0, 3.8, 3.3, 3.8, 3.3, 3.6, 2.9, 3.6, 2.9),
  fontface = c(
    "bold", "bold", "bold", "bold", "bold",
    "bold", "bold", "plain", "bold", "plain"
  )
)

replication_boxes <- tibble(
  cohort = c("GSE174367", "GSE157827", "GSE160936", "GSE188545", "GSE243292"),
  detail = c("18 donors", "20 donors", "12 donors", "12 donors", "10 donors\nsupplementary"),
  xmin = c(3.25, 4.35, 5.45, 6.55, 7.65),
  xmax = c(4.15, 5.25, 6.35, 7.45, 8.55),
  ymin = 7.62,
  ymax = 8.25
)

workflow_arrows <- tibble(
  x = c(
    1.45, 5.95, 1.45, 2.40, 5.95,
    7.80, 4.95, 2.45, 7.00
  ),
  xend = c(
    1.45, 5.95, 1.45, 3.00, 5.95,
    7.80, 2.50, 2.45, 7.00
  ),
  y = c(
    7.40, 7.40, 5.90, 5.00, 4.20,
    2.70, 2.70, 1.40, 1.40
  ),
  yend = c(
    6.70, 6.70, 5.40, 5.00, 3.70,
    2.30, 2.30, 1.00, 1.00
  )
)

p1 <- ggplot() +
  geom_rect(
    data = workflow_boxes,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = fill,
      color = border
    ),
    linewidth = 0.85
  ) +
  geom_text(
    data = workflow_boxes,
    aes(
      x = (xmin + xmax) / 2,
      y = (ymin + ymax) / 2,
      label = label,
      size = text_size,
      fontface = fontface
    ),
    lineheight = 0.88
  ) +
  geom_rect(
    data = replication_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white",
    color = "#70AD47",
    linewidth = 0.55
  ) +
  geom_text(
    data = replication_boxes,
    aes(
      x = (xmin + xmax) / 2,
      y = (ymin + ymax) / 2 + 0.08,
      label = cohort
    ),
    size = 2.55,
    fontface = "bold",
    lineheight = 0.82
  ) +
  geom_text(
    data = replication_boxes,
    aes(
      x = (xmin + xmax) / 2,
      y = (ymin + ymax) / 2 - 0.13,
      label = detail
    ),
    size = 2.30,
    lineheight = 0.80
  ) +
  geom_segment(
    data = workflow_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.13, "in")),
    linewidth = 0.55,
    color = "grey35"
  ) +
  coord_cartesian(xlim = c(0.25, 9.15), ylim = c(0.05, 8.95), clip = "off") +
  scale_fill_identity() +
  scale_color_identity() +
  scale_size_identity() +
  theme_void(base_size = BASE_FONT_SIZE) +
  theme(
    legend.position = "none",
    plot.margin = margin(8, 8, 8, 8)
  )

figure_log[[length(figure_log) + 1]] <- save_figure(
  p1,
  "Figure_1__study_workflow",
  10.5,
  7.2
)

add_log(
  "Figure_1",
  "CREATED",
  "Workflow diagram redesigned with larger font, fewer words, and consensus signature terminology.",
  NA_character_
)


}

# ============================================================
# SECTION 5: Figure 2 — SEA-AD discovery volcano
# ============================================================

cat("Building Figure 2...\n")

vol_src <- read_source(
  c(
    "Figure_2__SEAAD_discovery_volcano_source.csv",
    "Figure3_SEAAD_volcano_source.csv",
    "SEAAD_DE_results_mito_adjusted_AD61026_v2.csv",
    "SuppTable3_SEAAD_DE_full_source.csv"
  ),
  required = TRUE
)

vol_df <- vol_src$data
copy_source_used(vol_src$path, "Figure_2")

if (nrow(vol_df) == 0 || "read_error" %in% colnames(vol_df)) {
  
  add_log("Figure_2", "SKIPPED", "Volcano source could not be read.", vol_src$path)
  
} else {
  
  gene_col <- pick_col(vol_df, c("gene", "Gene", "symbol", "SYMBOL"))
  logfc_col <- pick_col(vol_df, c("logFC", "discovery_logFC", "avg_log2FC", "log2FoldChange"))
  fdr_col <- pick_col(vol_df, c("adj.P.Val", "adj_p_val", "FDR", "padj", "qvalue", "adj.P.Val_original"))
  p_col <- pick_col(vol_df, c("P.Value", "pvalue", "p_val", "P.Value_original"))
  
  if (is.na(gene_col) || is.na(logfc_col) || is.na(fdr_col)) {
    
    add_log(
      "Figure_2",
      "SKIPPED",
      paste(
        "Could not detect required volcano columns. Columns were:",
        paste(colnames(vol_df), collapse = ", ")
      ),
      vol_src$path
    )
    
  } else {
    
    plot_vol <- vol_df %>%
      mutate(
        gene_plot = as.character(.data[[gene_col]]),
        logFC_plot = suppressWarnings(as.numeric(.data[[logfc_col]])),
        FDR_plot = suppressWarnings(as.numeric(.data[[fdr_col]])),
        P_plot = if (!is.na(p_col)) suppressWarnings(as.numeric(.data[[p_col]])) else FDR_plot,
        neglog10_plot = p_to_neglog10(FDR_plot),
        discovery_class = case_when(
          FDR_plot < DISCOVERY_FDR_CUTOFF & logFC_plot >= DISCOVERY_ABS_LOGFC_CUTOFF ~ "AD-upregulated",
          FDR_plot < DISCOVERY_FDR_CUTOFF & logFC_plot <= -DISCOVERY_ABS_LOGFC_CUTOFF ~ "AD-downregulated",
          TRUE ~ "Not discovery DEG"
        )
      ) %>%
      filter(!is.na(logFC_plot), !is.na(FDR_plot))
    
    label_df <- bind_rows(
      plot_vol %>%
        filter(discovery_class == "AD-upregulated") %>%
        arrange(FDR_plot) %>%
        slice_head(n = 6),
      plot_vol %>%
        filter(discovery_class == "AD-downregulated") %>%
        arrange(FDR_plot) %>%
        slice_head(n = 6)
    )
    
    p2 <- ggplot(plot_vol, aes(x = logFC_plot, y = neglog10_plot)) +
      geom_point(aes(color = discovery_class), alpha = 0.70, size = 1.1) +
      geom_vline(
        xintercept = c(-DISCOVERY_ABS_LOGFC_CUTOFF, DISCOVERY_ABS_LOGFC_CUTOFF),
        linetype = "dashed",
        linewidth = 0.35
      ) +
      geom_hline(
        yintercept = -log10(DISCOVERY_FDR_CUTOFF),
        linetype = "dashed",
        linewidth = 0.35
      ) +
      scale_color_manual(
        values = c(
          "AD-upregulated" = "#B2182B",
          "AD-downregulated" = "#2166AC",
          "Not discovery DEG" = "grey75"
        )
      ) +
      labs(
        x = "logFC, AD versus control",
        y = expression(-log[10]("FDR")),
        color = NULL
      ) +
      theme_ad()
    
    if (has_ggrepel && nrow(label_df) > 0) {
      p2 <- p2 +
        ggrepel::geom_text_repel(
          data = label_df,
          aes(label = gene_plot),
          size = 2.8,
          max.overlaps = 20,
          box.padding = 0.25,
          min.segment.length = 0.05,
          show.legend = FALSE
        )
    } else if (nrow(label_df) > 0) {
      p2 <- p2 +
        geom_text(
          data = label_df,
          aes(label = gene_plot),
          size = 2.5,
          vjust = -0.5,
          show.legend = FALSE
        )
    }
    
    figure_log[[length(figure_log) + 1]] <- save_figure(
      p2,
      "Figure_2__SEAAD_discovery_volcano",
      7.2,
      5.8
    )
    
    add_log(
      "Figure_2",
      "CREATED",
      paste0("Volcano created using columns: gene=", gene_col, ", logFC=", logfc_col, ", FDR=", fdr_col),
      vol_src$path
    )
  }
}


# ============================================================
# SECTION 6: Figure 3 — Consensus signature heatmap
# ============================================================
if (RUN_SUPERSEDED_FIGURES) {

cat("Building Figure 3...\n")

heat_src <- read_source(
  c(
    "Figure_3__strict_replication_heatmap_source.csv",
    "Figure4_replication_heatmap_strict_source.csv",
    "replication_assessment_PRIMARY.csv"
  ),
  required = TRUE
)

heat_df <- heat_src$data
copy_source_used(heat_src$path, "Figure_3")

if (nrow(heat_df) == 0 || "read_error" %in% colnames(heat_df)) {
  
  add_log("Figure_3", "SKIPPED", "Heatmap source could not be read.", heat_src$path)
  
} else {
  
  gene_col <- pick_col(heat_df, c("gene", "Gene", "symbol", "SYMBOL"))
  dataset_col <- pick_col(heat_df, c("dataset", "Dataset", "cohort", "Cohort", "study"))
  logfc_col <- pick_col(heat_df, c("logFC", "replication_logFC", "mean_replication_logFC", "value"))
  
  if (!is.na(gene_col) && !is.na(dataset_col) && !is.na(logfc_col)) {
    
    heat_long <- heat_df %>%
      transmute(
        gene = as.character(.data[[gene_col]]),
        dataset = as.character(.data[[dataset_col]]),
        logFC = suppressWarnings(as.numeric(.data[[logfc_col]]))
      )
    
  } else if (!is.na(gene_col)) {
    
    summary_cols <- c(
      gene_col,
      "discovery_logFC",
      "discovery_adj.P.Val",
      "n_same_direction",
      "n_nominal_p_lt_0.05",
      "mean_replication_logFC",
      "RankProd_pfp",
      "RankProd_pval",
      "RankProd_AveFC"
    )
    
    numeric_cols <- colnames(heat_df)[sapply(heat_df, is.numeric)]
    numeric_cols <- setdiff(numeric_cols, summary_cols)
    
    heat_long <- heat_df %>%
      select(all_of(gene_col), all_of(numeric_cols)) %>%
      rename(gene = all_of(gene_col)) %>%
      pivot_longer(
        cols = -gene,
        names_to = "dataset",
        values_to = "logFC"
      ) %>%
      mutate(logFC = suppressWarnings(as.numeric(logFC)))
    
  } else {
    
    heat_long <- tibble()
  }
  
  strict_src <- read_source(
    c(
      "final_strict_consensus_high_or_rankprod_AD61026.csv",
      "Main_Table_2__strict_48_consensus_signature.csv"
    ),
    required = FALSE
  )
  
  if (nrow(strict_src$data) > 0 && "gene" %in% colnames(strict_src$data)) {
    strict_genes <- strict_src$data$gene
    heat_long <- heat_long %>% filter(gene %in% strict_genes)
  }
  
  # Remove duplicate SEA-AD column if both Discovery and SEAAD are present.
  heat_long <- heat_long %>%
    mutate(dataset_lower = tolower(dataset)) %>%
    filter(!(dataset_lower %in% c("seead", "seaad", "sea-ad", "sea_ad") & "Discovery" %in% dataset)) %>%
    select(-dataset_lower)
  
  if (nrow(heat_long) == 0) {
    
    add_log(
      "Figure_3",
      "SKIPPED",
      "Could not convert heatmap source to gene/dataset/logFC long format.",
      heat_src$path
    )
    
  } else {
    
    gene_order <- heat_long %>%
      group_by(gene) %>%
      summarize(mean_logFC = mean(logFC, na.rm = TRUE), .groups = "drop") %>%
      arrange(mean_logFC) %>%
      pull(gene)
    
    heat_long <- heat_long %>%
      mutate(
        gene = factor(gene, levels = gene_order),
        dataset = clean_label(dataset)
      )
    
    p3 <- ggplot(heat_long, aes(x = dataset, y = gene, fill = logFC)) +
      geom_tile(color = "white", linewidth = 0.20) +
      scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        name = "logFC"
      ) +
      labs(
        x = NULL,
        y = NULL
      ) +
      theme_ad() +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(size = 6.5),
        panel.grid = element_blank()
      )
    
    figure_log[[length(figure_log) + 1]] <- save_figure(
      p3,
      "Figure_3__strict_consensus_signature_heatmap",
      7.8,
      9.2
    )
    
    add_log(
      "Figure_3",
      "CREATED",
      paste0(
        "Consensus signature heatmap created with ",
        n_distinct(heat_long$gene),
        " genes and ",
        n_distinct(heat_long$dataset),
        " datasets."
      ),
      heat_src$path
    )
  }
}


}

# ============================================================
# SECTION 7: Figure 4 — Strict 48-gene consensus signature
# ============================================================
if (RUN_SUPERSEDED_FIGURES) {

cat("Building Figure 4...\n")

strict_src <- read_source(
  c(
    "Figure_4__strict_signature_direction_source.csv",
    "final_strict_consensus_high_or_rankprod_AD61026.csv",
    "MainTable2_strict_consensus_signature_source.csv",
    "Main_Table_2__strict_48_consensus_signature.csv"
  ),
  required = TRUE
)

strict_df <- strict_src$data
copy_source_used(strict_src$path, "Figure_4")

if (nrow(strict_df) == 0 || "read_error" %in% colnames(strict_df)) {
  
  add_log("Figure_4", "SKIPPED", "Strict consensus signature source could not be read.", strict_src$path)
  
} else {
  
  gene_col <- pick_col(strict_df, c("gene", "Gene"))
  logfc_col <- pick_col(strict_df, c("discovery_logFC", "SEA-AD logFC", "logFC"))
  conf_col <- pick_col(strict_df, c("confidence", "Confidence"))
  
  if (is.na(gene_col) || is.na(logfc_col)) {
    
    add_log("Figure_4", "SKIPPED", "Could not detect gene/logFC columns for strict signature plot.", strict_src$path)
    
  } else {
    
    plot_strict <- strict_df %>%
      mutate(
        gene = as.character(.data[[gene_col]]),
        logFC = suppressWarnings(as.numeric(.data[[logfc_col]])),
        confidence = if (!is.na(conf_col)) as.character(.data[[conf_col]]) else "Strict consensus",
        Direction = ifelse(logFC >= 0, "AD-upregulated", "AD-downregulated")
      ) %>%
      filter(!is.na(logFC)) %>%
      arrange(logFC) %>%
      mutate(gene = factor(gene, levels = gene))
    
    p4 <- ggplot(plot_strict, aes(x = logFC, y = gene)) +
      geom_vline(xintercept = 0, linewidth = 0.35, color = "grey40") +
      geom_segment(aes(x = 0, xend = logFC, yend = gene, color = Direction), linewidth = 0.55) +
      geom_point(aes(color = Direction, shape = confidence), size = 2.1) +
      scale_color_manual(
        values = c(
          "AD-upregulated" = "#B2182B",
          "AD-downregulated" = "#2166AC"
        )
      ) +
      labs(
        x = "SEA-AD discovery logFC, AD versus control",
        y = NULL,
        color = NULL,
        shape = "Evidence class"
      ) +
      theme_ad() +
      theme(
        axis.text.y = element_text(size = 6.8),
        panel.grid.major.y = element_blank()
      )
    
    figure_log[[length(figure_log) + 1]] <- save_figure(
      p4,
      "Figure_4__strict_48_consensus_signature_direction",
      7.4,
      9.2
    )
    
    add_log(
      "Figure_4",
      "CREATED",
      paste0("Strict 48-gene consensus direction plot created with ", nrow(plot_strict), " genes."),
      strict_src$path
    )
  }
}


}

# ============================================================
# SECTION 8: Figure 5 — Pathway enrichment summary
# ============================================================
if (RUN_SUPERSEDED_FIGURES) {

cat("Building Figure 5...\n")

hallmark <- read_source(
  c("GSEA_SEAAD_Hallmark.csv", "Figure_5__GSEA_dotplot_source.csv", "Figure6_GSEA_dotplot_source.csv"),
  required = TRUE
)

kegg <- read_source(c("GSEA_SEAAD_KEGG.csv"), required = FALSE)
gobp <- read_source(c("GSEA_SEAAD_GO_BP.csv"), required = FALSE)

copy_source_used(hallmark$path, "Figure_5_Hallmark")
copy_source_used(kegg$path, "Figure_5_KEGG")
copy_source_used(gobp$path, "Figure_5_GOBP")

make_gsea_subset <- function(df, collection_label, n_keep) {
  
  if (nrow(df) == 0 || "read_error" %in% colnames(df)) {
    return(tibble())
  }
  
  pathway_col <- pick_col(df, c("pathway_label", "Pathway", "pathway"))
  nes_col <- pick_col(df, c("NES", "nes"))
  fdr_col <- pick_col(df, c("padj", "FDR", "qvalue"))
  
  if (is.na(pathway_col) || is.na(nes_col) || is.na(fdr_col)) {
    return(tibble())
  }
  
  df %>%
    mutate(
      Collection = collection_label,
      Pathway = clean_label(.data[[pathway_col]]),
      NES = suppressWarnings(as.numeric(.data[[nes_col]])),
      FDR = suppressWarnings(as.numeric(.data[[fdr_col]])),
      Direction = ifelse(NES >= 0, "Enriched in AD", "Enriched in control")
    ) %>%
    filter(!is.na(NES), !is.na(FDR), FDR < 0.05) %>%
    arrange(FDR, desc(abs(NES))) %>%
    slice_head(n = n_keep)
}

path_df <- bind_rows(
  make_gsea_subset(hallmark$data, "Hallmark", 10),
  make_gsea_subset(kegg$data, "KEGG", 6),
  make_gsea_subset(gobp$data, "GO BP", 8)
) %>%
  mutate(
    neglog10FDR = p_to_neglog10(FDR),
    Pathway_wrapped = str_wrap(Pathway, width = 42)
  )

if (nrow(path_df) == 0) {
  
  add_log("Figure_5", "SKIPPED", "No pathway source data available or no significant pathways detected.", NA_character_)
  
} else {
  
  pathway_order <- path_df %>%
    arrange(Collection, NES) %>%
    pull(Pathway_wrapped)
  
  path_df <- path_df %>%
    mutate(
      Pathway_wrapped = factor(Pathway_wrapped, levels = unique(pathway_order)),
      Collection = factor(Collection, levels = c("Hallmark", "KEGG", "GO BP"))
    )
  
  p5 <- ggplot(path_df, aes(x = NES, y = fct_reorder(Pathway_wrapped, NES))) +
    geom_vline(xintercept = 0, linewidth = 0.30, color = "grey45") +
    geom_point(aes(size = neglog10FDR, color = Direction), alpha = 0.88) +
    facet_grid(Collection ~ ., scales = "free_y", space = "free_y") +
    scale_color_manual(
      values = c(
        "Enriched in AD" = "#B2182B",
        "Enriched in control" = "#2166AC"
      )
    ) +
    scale_size_continuous(
      name = expression(-log[10]("FDR")),
      range = c(2.0, 6.0)
    ) +
    labs(
      x = "Normalized enrichment score",
      y = NULL,
      color = NULL
    ) +
    theme_ad() +
    theme(
      axis.text.y = element_text(size = 7.5),
      strip.text.y = element_text(angle = 0),
      panel.grid.major.y = element_blank()
    )
  
  figure_log[[length(figure_log) + 1]] <- save_figure(
    p5,
    "Figure_5__pathway_enrichment_summary",
    8.2,
    9.0
  )
  
  add_log(
    "Figure_5",
    "CREATED",
    paste0("Pathway dotplot created with ", nrow(path_df), " selected significant pathways."),
    NA_character_
  )
}


}

# ============================================================
# SECTION 9: Figure 6 — Single-drug evidence matrix
# ============================================================

cat("Building Figure 6...\n")

drug_src <- read_source(
  c(
    "final_single_drug_STRICT_PAIR_INPUT_SUMMARY_AD61026.csv",
    "final_single_drug_FINAL_CURATION_SUMMARY_AD61026.csv",
    "Figure_6__single_drug_evidence_matrix_source.csv"
  ),
  required = TRUE
)

drug_df <- drug_src$data
copy_source_used(drug_src$path, "Figure_6")

if (nrow(drug_df) == 0 || "read_error" %in% colnames(drug_df)) {
  
  add_log("Figure_6", "SKIPPED", "Single-drug evidence source could not be read.", drug_src$path)
  
} else {
  
  drug_col <- pick_col(drug_df, c("drug_clean", "drug", "Drug", "candidate_drug"))
  pool_col <- pick_col(drug_df, c("manual_pool", "pair_pool", "candidate_pool", "Candidate-use category"))
  clue_col <- pick_col(drug_df, c("clue_strong_support", "CLUE_support", "clue_support"))
  ilincs_col <- pick_col(drug_df, c("ilincs_strong_support", "iLINCS_support", "ilincs_support"))
  dreimt_col <- pick_col(drug_df, c("dreimt_quality_label", "DREIMT", "dreimt_context", "DREIMT context"))
  dgidb_col <- pick_col(drug_df, c("dgidb_present", "DGIdb_present", "dgidb_support", "DGIdb target-level annotation"))
  tool_col <- pick_col(drug_df, c("transcriptomic_tool_count", "tool_count", "Transcriptomic tool count"))
  
  if (is.na(drug_col)) {
    
    add_log(
      "Figure_6",
      "SKIPPED",
      paste("Could not detect drug column. Columns were:", paste(colnames(drug_df), collapse = ", ")),
      drug_src$path
    )
    
  } else {
    
    plot_drug <- drug_df %>%
      mutate(
        drug_clean_plot = as.character(.data[[drug_col]]),
        tool_count_plot = if (!is.na(tool_col)) suppressWarnings(as.numeric(.data[[tool_col]])) else 0,
        pool_plot = if (!is.na(pool_col)) as.character(.data[[pool_col]]) else "candidate",
        clue_plot = if (!is.na(clue_col)) as.character(.data[[clue_col]]) else "unknown",
        ilincs_plot = if (!is.na(ilincs_col)) as.character(.data[[ilincs_col]]) else "unknown",
        dreimt_plot = if (!is.na(dreimt_col)) as.character(.data[[dreimt_col]]) else "unknown",
        dgidb_plot = if (!is.na(dgidb_col)) as.character(.data[[dgidb_col]]) else "unknown"
      ) %>%
      mutate(
        Drug = drug_display(drug_clean_plot),
        Candidate_class = case_when(
          str_detect(pool_plot, regex("primary", ignore_case = TRUE)) ~ "Primary",
          str_detect(pool_plot, regex("context", ignore_case = TRUE)) ~ "Context",
          str_detect(pool_plot, regex("explor", ignore_case = TRUE)) ~ "Explor.",
          TRUE ~ "Candidate"
        ),
        CLUE = case_when(
          clue_plot %in% c("TRUE", "True", "true", "1", "Supported") ~ "Yes",
          clue_plot %in% c("FALSE", "False", "false", "0", "Not supported") ~ "No",
          TRUE ~ "NA"
        ),
        iLINCS = case_when(
          ilincs_plot %in% c("TRUE", "True", "true", "1", "Supported") ~ "Yes",
          ilincs_plot %in% c("FALSE", "False", "false", "0", "Not supported") ~ "No",
          TRUE ~ "NA"
        ),
        DREIMT = case_when(
          str_detect(dreimt_plot, regex("FDR", ignore_case = TRUE)) ~ "FDR",
          str_detect(dreimt_plot, regex("tau", ignore_case = TRUE)) ~ "Tau",
          str_detect(dreimt_plot, regex("mixed", ignore_case = TRUE)) ~ "Mixed",
          str_detect(dreimt_plot, regex("clean", ignore_case = TRUE)) ~ "Clean",
          TRUE ~ "No"
        ),
        DGIdb = case_when(
          dgidb_plot %in% c("TRUE", "True", "true", "1") ~ "Target",
          str_detect(dgidb_plot, regex("target|DGIdb|present|yes", ignore_case = TRUE)) ~ "Target",
          dgidb_plot %in% c("FALSE", "False", "false", "0") ~ "No",
          TRUE ~ "No"
        )
      )
    
    drug_order <- plot_drug %>%
      mutate(
        pool_rank = case_when(
          Candidate_class == "Primary" ~ 1,
          Candidate_class == "Context" ~ 2,
          Candidate_class == "Explor." ~ 3,
          TRUE ~ 4
        )
      ) %>%
      arrange(pool_rank, desc(tool_count_plot), Drug) %>%
      pull(Drug)
    
    drug_long <- plot_drug %>%
      mutate(Drug = factor(Drug, levels = rev(unique(drug_order)))) %>%
      select(Drug, Candidate_class, CLUE, iLINCS, DREIMT, DGIdb) %>%
      pivot_longer(
        cols = c("Candidate_class", "CLUE", "iLINCS", "DREIMT", "DGIdb"),
        names_to = "Evidence_layer",
        values_to = "Label"
      ) %>%
      mutate(
        Evidence_layer = factor(
          Evidence_layer,
          levels = c("Candidate_class", "CLUE", "iLINCS", "DREIMT", "DGIdb"),
          labels = c("Candidate class", "CLUE/CMap", "iLINCS", "DREIMT", "DGIdb")
        ),
        Evidence_score = case_when(
          Label %in% c("Primary", "Yes", "FDR", "Target", "Clean") ~ 3,
          Label %in% c("Context", "Tau") ~ 2,
          Label %in% c("Explor.", "Mixed") ~ 1,
          TRUE ~ 0
        )
      )
    
    p6 <- ggplot(drug_long, aes(x = Evidence_layer, y = Drug, fill = Evidence_score)) +
      geom_tile(color = "white", linewidth = 0.35) +
      geom_text(aes(label = Label), size = 3.0) +
      scale_fill_gradient(
        low = "grey92",
        high = "#2B8CBE",
        limits = c(0, 3),
        breaks = c(0, 1, 2, 3),
        name = "Support\nlevel"
      ) +
      labs(
        x = NULL,
        y = NULL
      ) +
      theme_ad() +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1),
        panel.grid = element_blank()
      )
    
    figure_log[[length(figure_log) + 1]] <- save_figure(
      p6,
      "Figure_6__single_drug_evidence_matrix",
      8.2,
      6.0
    )
    
    add_log(
      "Figure_6",
      "CREATED",
      paste0("Single-drug matrix created with ", n_distinct(drug_long$Drug), " drugs."),
      drug_src$path
    )
  }
}


# ============================================================
# SECTION 10: Figure 7 — Drug-pair supporting evidence matrix
# ============================================================
if (RUN_SUPERSEDED_FIGURES) {
# ------------------------------------------------------------
# Force the correct project folder.
# This prevents OUT_DIR from being inherited from another script.
# ------------------------------------------------------------


cat("Building Figure 7...\n")

pair_src <- read_source(
  c(
    "Figure_7__drug_pair_program_complementarity_source.csv",
    "final_pair_supporting_evidence_table_PKDDI_REFINED_v2_AD61026.csv"
  ),
  required = TRUE
)

pair_df <- pair_src$data
copy_source_used(pair_src$path, "Figure_7")

if (nrow(pair_df) == 0 || "read_error" %in% colnames(pair_df)) {
  
  add_log("Figure_7", "SKIPPED", "Pair supporting evidence source could not be read.", pair_src$path)
  
} else {
  
  expected_pair_cols <- c(
    "Drug_1",
    "Drug_2",
    "Combined_transcriptomic_tool_count",
    "Program_union_n",
    "Pair_has_any_neuro_trial_context",
    "n_STRING_cross_drug_edges",
    "refined_PK_DDI_feasibility"
  )
  
  if (!all(expected_pair_cols %in% colnames(pair_df))) {
    
    add_log(
      "Figure_7",
      "SKIPPED",
      paste("Pair source missing expected columns. Columns were:", paste(colnames(pair_df), collapse = ", ")),
      pair_src$path
    )
    
  } else {
    
    plot_pair <- pair_df %>%
      mutate(
        Pair = paste0(drug_display(Drug_1), " + ", drug_display(Drug_2)),
        Pair = stringr::str_replace_all(Pair, "VX-745 / neflamapimod", "VX-745"),
        Pair = stringr::str_replace_all(Pair, "VX-745/neflamapimod", "VX-745"),
        Pair = stringr::str_replace_all(Pair, "VX-745 / Neflamapimod", "VX-745"),
        pair_rank = row_number(),
        Transcriptomic = case_when(
          Combined_transcriptomic_tool_count >= 4 ~ "Strong",
          Combined_transcriptomic_tool_count >= 3 ~ "Moderate",
          TRUE ~ "Limited"
        ),
        Program_coverage = case_when(
          Program_union_n >= 4 ~ "Broad",
          Program_union_n == 3 ~ "Moderate",
          TRUE ~ "Limited"
        ),
        Clinical_context = ifelse(Pair_has_any_neuro_trial_context, "Neuro context", "No neuro context"),
        STRING = ifelse(n_STRING_cross_drug_edges > 0, "Direct edge", "No direct edge"),
        PK_DDI = as.character(refined_PK_DDI_feasibility)
      ) %>%
      arrange(pair_rank)
    
    pair_order <- rev(plot_pair$Pair)
    
    pair_long <- plot_pair %>%
      select(Pair, Transcriptomic, Program_coverage, Clinical_context, STRING, PK_DDI) %>%
      pivot_longer(
        cols = -Pair,
        names_to = "Evidence_layer",
        values_to = "Evidence"
      ) %>%
      mutate(
        Pair = factor(Pair, levels = pair_order),
        Evidence_layer = factor(
          Evidence_layer,
          levels = c("Transcriptomic", "Program_coverage", "Clinical_context", "STRING", "PK_DDI"),
          labels = c(
            "Transcriptomic\nsupport",
            "Program\ncoverage",
            "Clinical\ncontext",
            "STRING\ntargets",
            "PK/DDI\nfeasibility"
          )
        ),
        Evidence_score = case_when(
          Evidence %in% c("Strong", "Broad", "Neuro context", "Direct edge", "Low", "Low-to-moderate") ~ 3,
          Evidence %in% c("Moderate") ~ 2,
          Evidence %in% c("Unknown / experimental", "Unknown/experimental", "Limited", "No direct edge") ~ 1,
          Evidence %in% c("No neuro context") ~ 0,
          TRUE ~ 1
        ),
        Label = case_when(
          Evidence == "Strong" ~ "Strong",
          Evidence == "Moderate" ~ "Mod.",
          Evidence == "Broad" ~ "Broad",
          Evidence == "Limited" ~ "Limited",
          Evidence == "Neuro context" ~ "Yes",
          Evidence == "No neuro context" ~ "No",
          Evidence == "Direct edge" ~ "Edge",
          Evidence == "No direct edge" ~ "No edge",
          Evidence %in% c("Unknown / experimental", "Unknown/experimental") ~ "Experimental",
          TRUE ~ Evidence
        )
      )
    
    p7 <- ggplot(pair_long, aes(x = Evidence_layer, y = Pair, fill = Evidence_score)) +
      geom_tile(color = "white", linewidth = 0.35) +
      geom_text(aes(label = Label), size = 2.9) +
      scale_fill_gradient(
        low = "grey92",
        high = "#238B45",
        limits = c(0, 3),
        breaks = c(0, 1, 2, 3),
        name = "Support\nlevel"
      ) +
      labs(
        x = NULL,
        y = NULL
      ) +
      theme_ad() +
      theme(
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid = element_blank()
      )
    
    figure_log[[length(figure_log) + 1]] <- save_figure(
      p7,
      "Figure_7__drug_pair_support_matrix",
      9.0,
      6.8
    )
    
    add_log(
      "Figure_7",
      "CREATED",
      paste0("Pair support matrix created with ", n_distinct(pair_long$Pair), " pairs."),
      pair_src$path
    )
  }
}


}

# ============================================================
# SECTION 11: Save figure build QC logs
# ============================================================

figure_log_df <- bind_rows(figure_log)

readr::write_csv(
  figure_log_df,
  file.path(FINAL_FIG_QC_DIR, "figure_build_QC_log_AD61026.csv")
)

jpg_manifest <- tibble(
  jpg_file = list.files(FINAL_FIG_JPG_DIR, pattern = "\\.jpg$", full.names = FALSE),
  jpg_path = list.files(FINAL_FIG_JPG_DIR, pattern = "\\.jpg$", full.names = TRUE)
) %>%
  mutate(
    jpg_path = normalizePath(jpg_path, winslash = "/", mustWork = FALSE),
    file_size_kb = round(file.info(jpg_path)$size / 1024, 1)
  )

readr::write_csv(
  jpg_manifest,
  file.path(FINAL_FIG_QC_DIR, "generated_JPG_manifest_AD61026.csv")
)

source_used <- figure_log_df %>%
  filter(!is.na(source_path), source_path != "") %>%
  distinct(figure_id, source_path)

readr::write_csv(
  source_used,
  file.path(FINAL_FIG_QC_DIR, "figure_source_files_used_AD61026.csv")
)


# ============================================================
# SECTION 12: Print final summary
# ============================================================

cat("\n============================================================\n")
cat("  Figure build complete\n")
cat("============================================================\n\n")

cat("JPG figure folder:\n")
cat("  ", FINAL_FIG_JPG_DIR, "\n\n")

cat("PDF internal-check folder:\n")
cat("  ", FINAL_FIG_PDF_DIR, "\n\n")

cat("QC folder:\n")
cat("  ", FINAL_FIG_QC_DIR, "\n\n")

cat("Generated JPG files:\n")
print(jpg_manifest)

cat("\nFigure build QC log:\n")
print(figure_log_df)

cat("\nImportant check:\n")
cat("  Open every JPG file and confirm axis labels are readable before manuscript submission.\n")
cat("  Figure titles are intentionally removed from the images.\n")
cat("  Use manuscript captions for figure titles and descriptions.\n")
cat("============================================================\n")
