# ============================================================
# Fig 4: Cross-cohort heatmap of the 48-gene consensus signature
# ============================================================
# Replaces SECTION 6 of 11_build_manuscript_figures_JPG_AD61026.R.
# Fixes: duplicate SEA-AD/Discovery column, asymmetric color scale,
# unexplained grey cells, gene order different from Table 2.
#
# Inputs (the supplementary workbooks):
#   S1 Table, sheet "Replication_Assessment" (per-cohort logFC and P values)
#   S2 Table, sheet "Consensus48"             (the 48 consensus genes)
# Output: Fig4.tif (7.5 in wide, 300 dpi)

# ------------------------------------------------------------
# Project folder. Set the AD_PROJECT_ROOT environment variable to the
# folder that contains Data/, RDS/, Output/ and Supplementary_Tables/,
# or start R with that folder as the working directory.
# ------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("AD_PROJECT_ROOT", unset = getwd())

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

dir.create(file.path(PROJECT_ROOT, "Output", "Figures"), recursive = TRUE, showWarnings = FALSE)

s1_file <- file.path(PROJECT_ROOT, "Supplementary_Tables", "Supplementary_Table_1_Cohort_QC_Replication.xlsx")
s2_file <- file.path(PROJECT_ROOT, "Supplementary_Tables", "Supplementary_Table_2_DE_Consensus_Signatures.xlsx")

cohorts <- c("GSE174367", "GSE157827", "GSE160936", "GSE188545")
col_labels <- c(
  "SEA-AD"    = "SEA-AD\nMTG+DLPFC\n84 donors",
  "GSE174367" = "GSE174367\nPFC\n18 donors",
  "GSE157827" = "GSE157827\nPFC\n20 donors",
  "GSE160936" = "GSE160936\nEC+SSC\n12 donors",
  "GSE188545" = "GSE188545\nMTG\n12 donors",
  "Mean"      = "Mean\nreplication"
)

strict <- read_excel(s2_file, sheet = "Consensus48")
rep    <- read_excel(s1_file, sheet = "Replication_Assessment") %>%
  filter(gene %in% strict$gene)

# Long format: one row per gene x column
long <- bind_rows(
  rep %>% transmute(gene, column = "SEA-AD", logFC = discovery_logFC, P = NA_real_),
  lapply(cohorts, function(cc) {
    rep %>% transmute(
      gene, column = cc,
      logFC = as.numeric(.data[[paste0(cc, "_logFC")]]),
      P     = as.numeric(.data[[paste0(cc, "_P.Value")]])
    )
  }) %>% bind_rows(),
  rep %>% transmute(gene, column = "Mean", logFC = mean_replication_logFC, P = NA_real_)
)

# Row order: up genes (strongest first), then down genes (strongest last),
# matching Table 2
gene_info <- rep %>%
  transmute(gene, sea = discovery_logFC,
            group = ifelse(sea > 0, "Up in AD (24)", "Down in AD (24)"))
# ggplot draws the first factor level at the bottom
gene_levels <- c(
  gene_info %>% filter(sea < 0) %>% arrange(sea) %>% pull(gene),
  gene_info %>% filter(sea > 0) %>% arrange(sea) %>% pull(gene)
)

long <- long %>%
  left_join(gene_info, by = "gene") %>%
  mutate(
    gene   = factor(gene, levels = gene_levels),
    column = factor(column, levels = names(col_labels)),
    group  = factor(group, levels = c("Up in AD (24)", "Down in AD (24)")),
    fill   = pmax(pmin(logFC, 2), -2),             # cap at +/-2
    star   = ifelse(!is.na(P) & P < 0.05, "*", ""),
    star_col = ifelse(abs(fill) > 1.4, "white", "black")
  )

p <- ggplot(long, aes(x = column, y = gene)) +
  geom_tile(aes(fill = fill), color = "white", linewidth = 0.3) +
  geom_text(aes(label = star, color = star_col), size = 4, vjust = 0.75) +
  geom_text(data = filter(long, is.na(logFC)), label = "n.d.",
            size = 2.6, color = "grey30") +
  scale_color_identity() +
  scale_fill_gradientn(
    colours = c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B"),
    values  = scales::rescale(c(-2, -1, 0, 1, 2)),
    limits  = c(-2, 2),
    breaks  = c(-2, -1, 0, 1, 2),
    labels  = c("≤ −2", "−1", "0", "1", "≥ 2"),
    na.value = "grey75",
    name = "logFC\n(AD vs control)"
  ) +
  scale_x_discrete(labels = col_labels, position = "top") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL,
       caption = "* nominal P < 0.05 in that cohort; n.d. = not detected") +
  theme_minimal(base_size = 9, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(face = "italic", size = 7.5),
    axis.text.x.top = element_text(size = 7.5, lineheight = 0.9),
    strip.placement = "outside",
    strip.text.y.left = element_text(face = "bold", angle = 90),
    panel.spacing.y = unit(6, "pt"),
    plot.caption = element_text(hjust = 0, size = 7.5)
  )

ggsave(file.path(PROJECT_ROOT, "Output", "Figures", "Fig4.tif"), p, width = 7.5, height = 8.3, dpi = 300,
       compression = "lzw", bg = "white")
