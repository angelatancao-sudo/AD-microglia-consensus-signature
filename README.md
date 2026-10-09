# Cross-cohort transcriptomic profiling of Alzheimer's disease microglia

Analysis code for the manuscript *"Cross-cohort transcriptomic profiling identifies a
reproducible Alzheimer's disease microglial consensus signature."*

The code identifies AD-associated microglial genes in the SEA-AD discovery cohort, tests
their replication in four independent single-nucleus RNA-seq cohorts with donor-level
pseudobulk analysis, and runs the sensitivity, orthogonal-support, pathway, and exploratory
drug-prioritization analyses reported in the paper.

## Requirements

- R 4.6.0 (scripts were also syntax-checked with R 4.3).
- CRAN packages: `dplyr`, `tidyr`, `tibble`, `readr`, `readxl`, `stringr`, `purrr`, `forcats`,
  `ggplot2`, `ggrepel`, `patchwork`, `pheatmap`, `RColorBrewer`, `scales`, `Matrix`, `Seurat`,
  `hdf5r`, `R.utils`, `writexl`, `openxlsx`, `janitor`, `httr`, `jsonlite`, `msigdbr`.
- Bioconductor packages: `edgeR`, `limma`, `fgsea`, `RankProd`.

## Data

All data are public. Place the files in a `Data/` folder inside the project folder, using the
file names listed in `R/01_setup.R`.

| Dataset | Role | Source |
|---|---|---|
| SEA-AD (MTG + DLPFC) | Discovery | SEA-AD Microglia-and-Immune multi-regional h5ad file and donor metadata from the Allen Institute SEA-AD release |
| GSE174367 | Replication | GEO |
| GSE157827 | Replication | GEO (raw 10x files) |
| GSE160936 | Replication | GEO (raw 10x files and series metadata) |
| GSE188545 | Replication | GEO (raw 10x files) |
| GSE216999 | Orthogonal support (xenograft) | GEO; downloaded automatically by `13_GSE216999_xenograft_pseudobulk.R` |
| GSE243292 | Supplementary dataset processed by scripts 03e and 04; not used in the reported results | GEO |

Several later scripts read the published supplementary tables (S1–S3 Tables). Place them in a
`Supplementary_Tables/` folder inside the project folder.

## Setup

Set the project folder once, either with an environment variable:

```r
Sys.setenv(AD_PROJECT_ROOT = "path/to/project")
```

or by starting R with the project folder as the working directory. Every script reads this
location; no other paths need to be edited. `01_setup.R` creates the `RDS/` and `Output/`
folders and saves a project settings file used by the later scripts.

## Run order

| Script | Purpose | Paper |
|---|---|---|
| `01_setup.R` | Folders, dataset registry, analysis parameters | Methods |
| `02_load_SEAAD.R` | SEA-AD microglia, donor-level pseudobulk | Methods; Table 1 |
| `03a_parse_GSE160936_metadata.R` | Parse GSE160936 sample metadata | Methods |
| `03b_load_GEO_and_cluster.R` | Load raw GEO datasets and cluster all nuclei | Methods |
| `03c_validate_microglia_clusters.R` | Marker-based review of microglial clusters | Methods; S1 Table |
| `03d_extract_microglia_pseudobulk.R` | Extract microglial clusters, donor-level pseudobulk | Methods; S1 Table |
| `03e_build_preannotated_GEO_pseudobulk.R` | Pseudobulk for author-annotated datasets | Methods |
| `03f_check_pseudobulk_inputs.R` | Checks before differential expression | S1 Table |
| `04_differential_expression_and_replication.R` | edgeR/limma-voom DE, vote-counting, RankProd | Results; Fig 2; S1–S2 Tables |
| `05_pathway_analysis.R` | GSEA (Hallmark, KEGG, GO BP) and ORA | Results; Table 3; S3 Table |
| `06_consensus_signature_and_source_data.R` | 48-gene consensus signature and source tables | Results; Table 2 |
| `07_SEAAD_mitochondrial_sensitivity.R` | Mitochondrial-fraction sensitivity analysis | Results; S2 Table |
| `08_permutation_test.R` | Gene-label permutation test of replication | Results; Fig 3 |
| `09_SEAAD_covariate_sensitivity.R` | Age, sex, and PMI sensitivity analysis | Results; S2 Table |
| `10_SEAAD_neuropathology_grouping_sensitivity.R` | SEA-AD regrouped by AD neuropathologic change | Results; S2 Table |
| `11_GSE188545_without_cluster18.R` | Replication repeated without GSE188545 cluster 18 | Results; S1 Table |
| `12_consensus_genes_by_cell_type.R` | Expression of consensus genes across brain cell types | Results; S2 Table |
| `13_GSE216999_xenograft_pseudobulk.R` | Mouse-level pseudobulk of xenografted human microglia | Results; S4 Table |
| `14_GSE216999_supplementary_table.R` | Assembles the xenograft supplementary table | S4 Table |
| `15a_drug_query_inputs.R` | Gene lists for perturbational drug queries | Methods; S5 Table |
| `15b_drug_prioritization_interim.R` | Combine CLUE/CMap, iLINCS, and DGIdb results | S5 Table |
| `15c_drug_prioritization_final.R` | Add DREIMT results and check signature direction | S5 Table |
| `15d_drug_name_QC.R` | Harmonize and collapse drug names | S5 Table |
| `15e_drug_curation_DREIMT_context.R` | DREIMT context classification and final curation | S5 Table |
| `15f_single_drug_summary.R` | Final single-drug candidate table | Table 4; Fig 6 |
| `16_build_main_tables.R` | Main tables | Tables 1–4 |
| `17a_build_figures.R` | Volcano plot and drug evidence matrix | Figs 2 and 6 |
| `17b_figure3_permutation.R` | Permutation figure | Fig 3 |
| `17c_figure4_heatmap.R` | Cross-cohort heatmap | Fig 4 |
| `17d_figure5_pathway.R` | Pathway enrichment dot plot | Fig 5 |

Fig 1 (study workflow) is a diagram and is not generated by code.

## Notes

- The drug-prioritization step (15a–15e) uses web tools. The query gene lists written by
  `15a` were submitted to CLUE/CMap, iLINCS, and DREIMT through their websites, and the
  exported results were saved to `Output/Drug_Prioritization/` before running `15b`–`15e`.
  DGIdb is queried through its public API.
- DREIMT threshold: `15c` and `15d` flag DREIMT reversal at tau < -75 as an interim screen.
  The final classification in `15e`, which is the one reported in the paper and S5 Table, uses
  tau < -80. All final candidates with a DREIMT reversal call have tau below -92, so the two
  cutoffs give the same result for them.
- Drug-level summaries pool all CLUE/CMap rows for a drug (all compound IDs and name variants):
  the lowest median score and the highest number of negative signatures are kept. The
  `CLUE_CMap_filtered` sheet of S5 Table shows one row per drug.
- An earlier version of this project included a drug-pair analysis that is not in the published
  paper. Its scripts are not included. Code that remains from it is switched off by default:
  `RUN_SUPERSEDED_FIGURES <- FALSE` in `17a` (earlier figure versions and a drug-pair figure) and
  `RUN_LEGACY_SUPPLEMENT <- FALSE` in `16` (an earlier ten-table supplement). Files written by
  `15e` and `15f` whose names contain "pair" hold the single-drug candidate categories used for
  Table 4; no drug pairs are generated.
- With the default settings, `17a` builds only Figs 2 and 6; the published Figs 3-5 are built
  by `17b`-`17d`.
- Some intermediate file and object names contain the internal project code `AD61026` or the
  label `Strict48`, an earlier name for the 48-gene consensus signature (sheet `Consensus48`
  in S2 Table).
