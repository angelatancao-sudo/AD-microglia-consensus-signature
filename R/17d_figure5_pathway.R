# ============================================================
# Fig 5: Pathway enrichment dot plot (AD- and control-enriched)
# ============================================================
# Replaces the pathway section of 11_build_manuscript_figures_JPG_AD61026.R.
# Shows the selected AD-enriched pathways (as before) plus the top
# control-enriched pathways by FDR, so both directions are visible.
# Input: S3 Table workbook. Output: Fig5.tif (7.5 in wide, 300 dpi)

# ------------------------------------------------------------
# Project folder. Set the AD_PROJECT_ROOT environment variable to the
# folder that contains Data/, RDS/, Output/ and Supplementary_Tables/,
# or start R with that folder as the working directory.
# ------------------------------------------------------------
PROJECT_ROOT <- Sys.getenv("AD_PROJECT_ROOT", unset = getwd())

library(readxl)
library(dplyr)
library(ggplot2)

dir.create(file.path(PROJECT_ROOT, "Output", "Figures"), recursive = TRUE, showWarnings = FALSE)

s3_file <- file.path(PROJECT_ROOT, "Supplementary_Tables", "Supplementary_Table_3_Pathway_Analysis.xlsx")
norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))

shown <- list(
  "Hallmark" = c("Hypoxia", "UV Response Up", "MYC Targets V1", "Androgen Response",
                 "Complement", "mTORC1 Signaling", "P53 Pathway", "Apoptosis",
                 "Epithelial Mesenchymal Transition", "Myogenesis"),
  "KEGG"     = c("Oxidative Phosphorylation", "Ribosome", "Parkinson's Disease",
                 "Regulation Of Actin Cytoskeleton", "Focal Adhesion",
                 "MAPK Signaling Pathway"),
  "GO BP"    = c("Proton Transmembrane Transport", "Alpha Beta T Cell Activation",
                 "Regulation Of Apoptotic Signaling Pathway",
                 "Negative Regulation Of Cell Adhesion", "Viral Process",
                 "Negative Regulation Of Immune System Process",
                 "Positive Regulation Of Programmed Cell Death",
                 "Lymphocyte Differentiation")
)
sheets <- c("Hallmark" = "GSEA_Hallmark_all", "KEGG" = "GSEA_KEGG_all",
            "GO BP" = "GSEA_GO_BP_significant")
n_neg <- 5

plot_df <- bind_rows(lapply(names(sheets), function(k) {
  d <- read_excel(s3_file, sheet = sheets[[k]])
  pos <- d %>% filter(norm(pathway_label) %in% norm(shown[[k]])) %>%
    mutate(label = shown[[k]][match(norm(pathway_label), norm(shown[[k]]))])
  neg <- d %>% filter(NES < 0, padj < 0.05) %>% arrange(padj) %>%
    slice_head(n = n_neg) %>% mutate(label = pathway_label)
  bind_rows(pos, neg) %>% mutate(collection = k)
})) %>%
  mutate(
    collection = factor(collection, levels = names(sheets)),
    direction  = ifelse(NES > 0, "AD", "Control"),
    neglog     = pmin(-log10(padj), 10)
  ) %>%
  group_by(collection) %>%
  arrange(NES, .by_group = TRUE) %>%
  mutate(label = factor(paste(collection, label, sep = "__"),
                        levels = unique(paste(collection, label, sep = "__")))) %>%
  ungroup()

p <- ggplot(plot_df, aes(x = NES, y = label)) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
  geom_point(aes(colour = direction, size = neglog), alpha = 0.9) +
  scale_y_discrete(labels = function(x) sub(".*__", "", x)) +
  scale_colour_manual(values = c(AD = "#B2182B", Control = "#2166AC"),
                      name = "Enriched in") +
  scale_size_continuous(name = "−log10(FDR)", range = c(2, 7),
                        breaks = c(2, 4, 6, 8)) +
  scale_x_continuous(limits = c(-2.6, 2.75), breaks = -2:2) +
  facet_grid(collection ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Normalized enrichment score (NES)", y = NULL) +
  theme_bw(base_size = 9, base_family = "Arial") +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        strip.text.y = element_text(face = "bold", angle = 0))

ggsave(file.path(PROJECT_ROOT, "Output", "Figures", "Fig5.tif"), p, width = 7.5, height = 6.5, dpi = 300,
       compression = "lzw", bg = "white")
