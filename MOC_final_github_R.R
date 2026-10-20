# Proteomic Profiling of Mucinous Ovarian Cancer

This repository contains the code used to reproduce the analyses and figures presented in:
  
  ## Comprehensive Proteomic Profiling Distinguishes Mucinous Ovarian Cancer from Extra-Ovarian Tumors
  
  All figures are generated using the original analysis and plotting code used in the manuscript.

All analyses were performed in **R version 4.5.1**.

---
  
  ## Analysis Workflow
  
  Run scripts in the following order:
  
1. `00_prepare_data.R`
2. `01_umap_compute_plots.R`
3. `02_consensus_clustering.R`
4. `03_DEP_analysis.R`
5. `04_GSEA_analysis.R`
6. `05_celltype_enrichment.R`
7. `06_Fig1.R` (including Fig 1e and Fig 1f scripts)
8. `07_Fig3.R`
9. `08_Fig4.R`
10. `09_Fig5.R`
11. `10_Fig6.R`
12. `11_survival.R`
13. `12_lasso.R`
14. `RNA dataset processing - see separate section`

---
  
  ## Input Data
  
  The primary input dataset for all proteomic analyses is:
  
  | File | Purpose |
  |--------|--------|
  | `combat_Omicodf5626.csv` | Combined batch-corrected proteomic cohort used for clustering, differential expression, pathway enrichment, ESTIMATE, therapeutic target and survival analyses. |
  
  Metadata and clinical annotations are available within the PRIDE repository and are merged during analysis using consistent sample identifiers.
  The integrated dataset (`combat_Omicodf5626.csv`) is generated from Cohorts 1-2 (`E0200_P15_working_matrix_UniID_2_.tsv`) and Cohort 3 (`cohort3_proteinmatrix.csv`), both of which are available through PRIDE.
  
  ## External Analyses
  
  Some analyses were performed using external tools:
  
  - **STRING-db**: protein-protein interaction network analyses.
- **WebGestalt**: pathway and cell-type enrichment analyses.
- **Werner et al. dataset**: independent validation of selected protein biomarkers.

Relevant input and output files are provided in the `results/` directory where applicable.

---
  
  ## Analytical Methods
  
  ### Differential Protein Abundance
  
- limma
- empirical Bayes moderation
- duplicateCorrelation modelling for repeated sampling

### Pathway Analysis

- GSVA / ssGSEA
- fgsea
- WebGestalt

### Clustering

- ConsensusClusterPlus
- PAM clustering
- Spearman distance metric

### Machine Learning

- Support Vector Machine (SVM)
- LASSO feature selection
- SHAP feature importance analysis

---
  
  ## Differential Abundance Comparisons
  
  Multiple differential abundance analyses were performed depending on the biological comparison of interest.

Examples include:
  
- Figure 1: MOC versus extra-ovarian tumour types (GEJ, PDAC, colorectal cancer, etc.)
- Figure 3: infiltrative versus expansile MOC
- Figure 5: ERBB2-amplified versus non-amplified MOC
- Figure 6: immune-high versus immune-low MOC
- Supplementary analyses: MOC versus high-grade serous ovarian carcinoma

These comparisons were analysed independently and subsequently used for downstream pathway analyses, cell-type enrichment, machine learning, and figure generation.

### Enrichment Analyses

Different enrichment methods were used depending on the biological question:
  
- **GSVA (ssGSEA)**: pathway scoring across samples and consensus clusters.
- **fgsea**: ranked enrichment visualisation.
- **WebGestalt**: selected pathway and cell-type enrichment analyses.

---
  
  ## Reproducibility
  
  All analyses were performed using fixed random seeds where applicable.

---
  
  ## Data Availability
  
  Proteomic datasets will be made available through PRIDE upon publication.

## Code Availability

All code required to reproduce the analyses presented in this study is provided across the repositories described below.

### Main R Repository

This repository contains all proteomic data processing, differential abundance analyses, clustering, pathway enrichment analyses, figure generation, survival analyses, and LASSO feature selection workflows.

### Machine Learning Repository

The machine-learning workflow is implemented separately in the Python script `MOC_classifier.py`, which reproduces classifier training, validation, SHAP interpretation, and clinical application analyses.

The final classifier uses a predefined 37-protein diagnostic panel reported in Supplementary Table 4 and implemented directly within the Python classifier script.

# Workflow Overview
|
  +--> Differential abundance analyses (03_DEP_analysis.R)
|
  +--> Consensus clustering (02_consensus_clustering.R)
|
  +--> Pathway enrichment (04_GSEA_analysis.R)
|
  +--> Figure generation (06-10)
|
  +--> LASSO feature selection (12_lasso.R)
|
  v
Supplementary Table 4
|
  v
MOC_classifier.py
|
  +--> Internal validation
+--> Independent validation
+--> External validation
+--> SHAP interpretation

## Figure-to-Script Mapping

| Figure | Script |
  |----------|----------|
  | Figure 1 | 06_Fig1.R |
  | Figure 2 | MOC_classifier.py |
  | Figure 3 | 07_Fig3.R |
  | Figure 4 | 08_Fig4.R |
  | Figure 5 | 09_Fig5.R |
  | Figure 6 | 10_Fig6.R |
  | Supplementary Figure 4G | 11_survival.R |
  | Supplementary Figure 3H | 07_Fig3.R |
  | LASSO feature selection | 12_lasso.R |

# ============================================================
# 00_prepare_data.R
# ============================================================
#
# Purpose:
# Prepare the analysis-ready proteomic matrix for all downstream
# analyses and figure generation.
#
# Inputs:
#   ../data/combat_Omicodf5626.csv
#       - batch-corrected and imputed protein abundance matrix
#
#   ../data/normal_liver_uniprot.csv
#       - Human Protein Atlas-derived list of liver-enriched proteins
#
# Processing steps:
#   1. Filter samples passing tumour content quality control
#   2. Remove proteins enriched in normal liver tissue
#
# Output:
#   combat_df


## Preprocessing Workflow

The integrated proteomic dataset (`combat_Omicodf5626.csv`) was generated as follows:
  
1. Discovery cohorts (E0060, E0147, E0153) were processed separately from the OMICO validation cohort.
2. Missing protein abundance values were imputed within each dataset using the MinProb method (imputeLCMD package).
3. The imputed datasets were merged using shared protein identifiers.
4. Batch effects between acquisition cohorts were corrected using ComBat (sva package).
5. The resulting batch-corrected matrix (`combat_Omicodf5626.csv`) was used for all downstream analyses.
6. Additional tumour-content and liver-protein filtering generated as described below.

# ============================================================
# 00_preprocessing.R
# ============================================================
#
# Purpose:
# Generate the integrated proteomic dataset used for
# downstream analyses.
#
# Workflow:
# 1. Process Cohorts 1-2
# 2. Process OMICO (Cohort 3)
# 3. Merge datasets
# 4. Batch correction using ComBat
#
# Output:
# combat_Omicodf5626.csv
#
# ============================================================

library(tidyverse)
library(MSnbase)
library(imputeLCMD)
library(sva)

# ============================================================
# PART 1: PROCESS COHORTS 1-2
# ============================================================

cohort12 <- read.csv("../data/E0200_P15_working_matrix_UniID_2_.tsv")
cohort12 <- E0200_P15_working_matrix_UniID_2_

# Convert to data frame if tibble
cohort12 <- as.data.frame(cohort12)

# Use Protein column as row names
rownames(cohort12) <- cohort12$Protein
cohort12$Protein <- NULL

# Transpose
cohort12_t <- as.data.frame(t(cohort12))

# Add Sample_ID column
cohort12_t$Sample_ID <- rownames(cohort12_t)

# Move Sample_ID to first column
cohort12_t <- cohort12_t[, c("Sample_ID", setdiff(names(cohort12_t), "Sample_ID"))]

rownames(cohort12_t) <- NULL

# Add metadata
names(Sample_metadata_clean)[names(Sample_metadata_clean) == "Sample ID"] <- "Sample_ID"

metadata12 <- Sample_metadata_clean %>%
  select(
    Sample_ID,
    Label,
    sample_cohort,
    Include_in_study,
    subject_collaborator_patient_id
  ) 

cohort12_t <- left_join(
  cohort12_t,
  metadata12,
  by = "Sample_ID"
)

# Keep only tumour types included in the study
cohort12_t <- cohort12_t %>%
  filter(
    Label %in% c(
      "Colorectal_adenocarcinoma",
      "Pancreas_adenocarcinoma",
      "High_grade_serous_ovary",
      "Mucinous_ovarian",
      "UpperGI_adenocarcinoma",
      "Small_bowel_adenocarcinoma"
    )
  )

cohort12_t[, 2:10634] <- lapply(
  cohort12_t[, 2:10634],
  as.numeric
)

cohort12_t <- cohort12_t %>%
  group_by(subject_collaborator_patient_id) %>%
  summarise(
    Sample_ID = Sample_ID[1],
    Label = Label[1],
    sample_cohort = sample_cohort[1],
    Include_in_study = Include_in_study[1],
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# Expression matrix
expr12 <- as.matrix(cohort12_t[, 6:10638])
rownames(expr12) <- cohort12_t$Sample_ID
expr12 <- t(expr12)

# MinProb imputation
msn12 <- MSnSet(exprs = expr12)
msn12_imp <- impute(msn12, method = "MinProb")

cohort12_imp <- as.data.frame(t(exprs(msn12_imp)))
cohort12_imp$Sample_ID <- rownames(cohort12_imp)

metadata12 <- Sample_metadata_clean %>%
  select(Sample_ID, Label, sample_cohort, Include_in_study)

cohort12_final <- left_join(
  cohort12_imp,
  metadata12,
  by = "Sample_ID"
)

colnames(cohort12_final) <- sub("#.*", "", colnames(cohort12_final))

# ============================================================
# PART 2: PROCESS OMICO (COHORT 3)
# ============================================================

cohort3_proteinmatrix <- read.csv("../data/cohort3_proteinmatrix.csv")

# Cohort 3 protein matrix
omico <- cohort3_proteinmatrix
names(omico)[names(omico) == "LIMS.ID1"] <- "Sample_ID"

# Metadata
metadata3 <- Sample_metadata_clean %>%
  select(
    Sample_ID,
    Label,
    sample_cohort,
    Include_in_study
  )

# Add metadata and remove excluded samples BEFORE imputation
omico <- left_join(
  omico,
  metadata3,
  by = "Sample_ID"
) 

# Identify protein columns
protein_cols <- setdiff(
  colnames(omico),
  c("Sample_ID", "Label", "sample_cohort", "Include_in_study")
)

# Convert protein columns to numeric
omico[protein_cols] <- lapply(
  omico[protein_cols],
  as.numeric
)

# Create expression matrix
expr3 <- as.matrix(omico[, protein_cols])
rownames(expr3) <- omico$Sample_ID
expr3 <- t(expr3)

# MinProb imputation
msn3 <- MSnSet(exprs = expr3)
msn3_imp <- impute(msn3, method = "MinProb")

# Back to sample x protein format
omico_imp <- as.data.frame(t(exprs(msn3_imp)))
omico_imp$Sample_ID <- rownames(omico_imp)

omico_final <- left_join(
  omico_imp,
  metadata3,
  by = "Sample_ID"
)

# ============================================================
# PART 3: MERGE DATASETS
# ============================================================

common_cols <- intersect(
  colnames(cohort12_final),
  colnames(omico_final)
)

cohort12_final <- cohort12_final[, common_cols]
omico_final <- omico_final[, common_cols]

merged_df <- bind_rows(
  cohort12_final,
  omico_final
)

merged_df <- merged_df %>%
  left_join(
    Sample_metadata_clean %>%
      select(Sample_ID, Include_in_study),
    by = "Sample_ID"
  )

# ============================================================
# PART 4: COMBAT BATCH CORRECTION
# ============================================================

protein_cols <- setdiff(
  colnames(merged_df),
  c("Sample_ID", "Label", "sample_cohort", "Include_in_study")
)

expr <- as.matrix(merged_df[, protein_cols])
expr <- t(expr)

mod <- model.matrix(~ Label, data = merged_df)

combat_matrix <- ComBat(
  dat = expr,
  batch = merged_df$sample_cohort,
  mod = mod,
  par.prior = TRUE,
  prior.plots = FALSE
)

colnames(combat_matrix) <- merged_df$Sample_ID

combat_df <- as.data.frame(t(combat_matrix))
combat_df$Sample_ID <- rownames(combat_df)

combat_df <- left_join(
  combat_df,
  merged_df %>%
    select(Sample_ID, Label, sample_cohort, Include_in_study),
  by = "Sample_ID"
)

# ============================================================
# SAVE OUTPUT
# ============================================================

library(openxlsx)

write.xlsx(
  combat_df,
  "../data/combat_Omicodf5626.xlsx",
  rowNames = FALSE
)

write.csv(
  combat_df,
  "../data/combat_Omicodf5626.csv",
  row.names = FALSE
)

#
# ============================================================

library(tidyverse)

# ============================================================
# Load analysis matrix
# ============================================================

combat_df <- read.csv("../data/combat_Omicodf5626.csv")

# ============================================================
# Filter samples
# Retain samples meeting tumour content criteria
# ============================================================

data_filtered <- combat_df %>%
  filter(Include_in_study == "Y")

## Note: Sample 03U0T was originally labelled as mucinous ovarian carcinoma in the metadata file but was subsequently confirmed to be a small bowel adenocarcinoma. 
## This sample was excluded from all MOC-specific analyses.

# Correct sample annotation based on clinical review
combat_df <- combat_df %>%
  mutate(
    Label = case_when(
      Sample_ID == "03U0T" ~ "Small_bowel_adenocarcinoma",
      TRUE ~ Label
    )
  )

# ============================================================
# Remove normal liver-associated proteins
#
# Liver-enriched proteins were obtained from the Human Protein
# Atlas and removed to minimise clustering driven by
# background liver tissue in metastatic biopsy samples.
# ============================================================
# This corresponds to the tissue category rna liver Tissue file downloaded from the Human Protein Atlas 
liver_genes <- read.csv("../data/normal_liver_uniprot.csv")$Uniprot
liver_genes <- unique(na.omit(liver_genes))

data_filtered <- data_filtered %>%
  dplyr::select(-any_of(liver_genes))

# ============================================================
# Save analysis-ready dataset
# ============================================================

write.csv(
  data_filtered,
  "../data/analysis_ready.csv",
  row.names = FALSE
)

# ============================================================
# 01_umap_compute_plots.R
# ============================================================
#
# Purpose:
# Generate UMAP visualisations of the proteomic cohort and
# associated supplementary figures.
#
# Inputs:
#   ../data/analysis_ready.csv
#
# Outputs:
#   ../data/UMAP_coordinates.csv
#
# Figures:
#   Figure 1A
#   Supplementary Figure 1B-E
#
# Notes:
#   UMAP was performed using Euclidean distance on the full
#   protein abundance matrix with fixed parameters:
#     n_neighbors = 15
#     min_dist = 0.1
#
# ============================================================

library(tidyverse)
library(umap)
library(ggplot2)

set.seed(1234)

# Load clean data
read.csv("../data/combat_Omicodf5626.csv")

# Correct sample annotation based on clinical review
Sample_metadata_clean <- Sample_metadata_clean %>%
  mutate(
    Label = case_when(
      Sample_ID == "03U0T" ~ "Small_bowel_adenocarcinoma",
      TRUE ~ Label
    )
  )

data <- data_filtered %>%
  left_join(
    Sample_metadata_clean %>%
      dplyr::select(Sample_ID, sample_site),
    by = "Sample_ID"
  )

# UMAP parameters (used consistently across all panels)
# Parameters selected through grid search and fixed for all UMAP visualizations.
config <- umap.defaults
config$n_neighbors <- 15
config$min_dist <- 0.1

# Run UMAP
umap_fit <- data %>%
  dplyr::select(where(is.numeric)) %>% # Select protein abundance matrix only
  umap(config = config)

# Store coordinates
umap_df <- as.data.frame(umap_fit$layout) %>%
  setNames(c("UMAP1", "UMAP2")) %>%
  bind_cols(
    data %>%
      dplyr::select(Sample_ID, Label, sample_cohort, sample_site)
  )

# Save for reuse
write.csv(umap_df, "../data/UMAP_coordinates.csv", row.names = FALSE)

# Load precomputed UMAP
umap_df <- read.csv("../data/UMAP_coordinates.csv")

# =========================
# Common colors used across all UMAP panels
# =========================

diag_colors <- c(
  "Appendix_adenocarcinoma" = "pink",
  "Colorectal_adenocarcinoma" = "red",
  "High_grade_serous_ovary" = "lightblue",
  "Mucinous_ovarian" = "black",
  "Pancreas_adenocarcinoma" = "#66a61e",
  "Small_bowel_adenocarcinoma" = "purple",
  "UpperGI_adenocarcinoma" = "orange"
)

clean_theme <- theme_classic(base_size = 16)

# =========================
# FIG 1A: Diagnosis UMAP
# =========================

p_A <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Label)) +
  geom_point(size = 3) +
  scale_color_manual(values = diag_colors) +
  labs(color = "Diagnosis") +
  clean_theme

ggsave("figures/Fig1A_UMAP_diagnosis.pdf", p_A, width = 6, height = 5)

# =========================
# SUPP FIG 1B: Cohort
# =========================

# Collapse Omico cohort identifiers into a single label
df_cohort <- umap_df %>%
  mutate(cohort_group = ifelse(grepl("^EC", sample_cohort),
                               "Omico", sample_cohort))

p_B <- ggplot(df_cohort,
              aes(UMAP1, UMAP2,
                  color = Label,
                  shape = cohort_group)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = diag_colors) +
  labs(color = "Diagnosis", shape = "Cohort") +
  clean_theme

ggsave("figures/SuppFig1B_cohort_umap.pdf", p_B, width = 6, height = 5)

# =========================
# SUPP FIG 1C: Sample site
# =========================

p_C <- ggplot(umap_df,
              aes(UMAP1, UMAP2,
                  color = sample_site,
                  shape = Label)) +
  geom_point(size = 2.5) +
  labs(color = "Sample site", shape = "Diagnosis") +
  clean_theme

ggsave("figures/SuppFig1C_sample_site_umap.pdf", p_C, width = 6, height = 5)

# =========================
# SUPP FIG 1D: Ovary only
# =========================

df_ovary <- umap_df %>%
  filter(sample_site == "Ovary")

p_D <- ggplot(df_ovary,
              aes(UMAP1, UMAP2, color = Label)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = diag_colors) +
  clean_theme

ggsave("figures/SuppFig1D_ovary_umap.pdf", p_D, width = 6, height = 5)

# =========================
# SUPP FIG 1E: Mucinous subset
# =========================

df_mucinous <- umap_df %>%
  filter(Label %in% c("Appendix_adenocarcinoma", "Mucinous_ovarian"))

p_E <- ggplot(df_mucinous,
              aes(UMAP1, UMAP2, color = Label)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = diag_colors) +
  clean_theme

ggsave("figures/SuppFig1E_mucinous_umap.pdf", p_E, width = 6, height = 5)


# ============================================================
# 02_consensus_clustering.R
# ============================================================

# ConsensusClusterPlus automatically generates:
# - consensus matrices
# - CDF plots
# - delta area plots
# These outputs are provided as Supplementary Figure 4.

# Input:
#   ../data/analysis_ready.csv
#
# Output:
#   ../results/MOC_clusters_k5.csv
#
# Downstream usage:
#   Figure 4A
#   Figure 4B
#   Figure 5
#   Figure 6
#   Survival analyses
#
# Clustering method:
#   ConsensusClusterPlus
#   PAM clustering
#   Spearman distance
#   1000 resampling iterations
#
# The optimal number of clusters (k = 5) was determined using:
#   - cumulative distribution function (CDF) plots
#   - delta area plots
#   - consensus matrices
#   - biological interpretability
#
# ============================================================

library(ConsensusClusterPlus)
library(tidyverse)

set.seed(1234)

# ============================================================
# Load MOC cohort
# ============================================================

combat_df <- read.csv("../data/combat_Omicodf5626.csv")

data_MOC <- combat_df %>%
  filter(Label == "Mucinous_ovarian")

# ============================================================
# Generate protein expression matrix
# ============================================================

# Select highly variable proteins for consensus clustering
# Protein variability was quantified using the median absolute deviation (MAD)
# across samples. The top 25% most variable proteins were retained for
# consensus clustering

# Select numeric protein expression columns only
protein_data <- data_MOC %>%
  dplyr::select(where(is.numeric))

# Convert to protein x sample matrix
data <- t(as.matrix(protein_data))

# Calculate MAD for each protein
gene_mads <- apply(data, 1, mad)

# Retain top 25% most variable proteins
top_25_percent <- ceiling(0.25 * length(gene_mads))
selected_genes <- names(sort(gene_mads, decreasing = TRUE)[1:top_25_percent])
filtered_data <- data[selected_genes, ]

# Create clustering input matrix
expr <- as.matrix(filtered_data)
colnames(expr) <- data_MOC$Sample_ID

# ============================================================
# Consensus clustering
# ============================================================

results <- ConsensusClusterPlus(
  expr,
  maxK = 8,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "pam",
  distance = "spearman",
  seed = 1234,
  plot = "pdf"
)

# ============================================================
# Selected solution: k = 5
# ============================================================

clusters <- results[[5]]$consensusClass

cluster_df <- data.frame(
  Sample_ID = names(clusters),
  clusterspearman = clusters
)

# ============================================================
# Save cluster assignments
# ============================================================

write.csv(
  cluster_df,
  "../results/MOC_clusters_k5.csv",
  row.names = FALSE
)

# ============================================================
# 03_DEP_analysis.R
# ============================================================
#
# Purpose:
# Perform differential protein abundance analysis using limma.
#
# Outputs from this script are used throughout the repository
# for downstream pathway enrichment, cell-type enrichment,
# machine learning feature selection and figure generation.
#
# Figures using these results:
#   Figure 1
#   Figure 3
#   Figure 5
#   Figure 6
#   Supplementary Figures
#
# Input:
#   ../data/analysis_ready.csv
#
# Output:
#   Multiple differential abundance result tables (*.csv)
#
# ============================================================

# ============================================================
# run_limma_dep()
#
# Generic wrapper for differential protein abundance analysis.
#
# Arguments:
#   data       : input dataframe
#   group_col  : grouping variable
#   group1     : comparison group
#   group2     : reference group
#   output_file: output csv filename
#
# Returns:
#   limma differential abundance results table
#
# ============================================================

library(tidyverse)
library(limma)
library(org.Hs.eg.db)
library(AnnotationDbi)

run_limma_dep <- function(data, group_col, group1, group2, output_file) {
  
  # Subset data
  data_subset <- data %>%
    filter(.data[[group_col]] %in% c(group1, group2)) %>%
    mutate(Group = ifelse(.data[[group_col]] == group1, 1, 0))
  
  # Expression matrix
  expr <- data_subset %>%
    select(where(is.numeric)) %>% # Expression matrix consists exclusively of protein abundance measurements after preprocessing.
    as.matrix() %>%
    t()
  
  # Design matrix
  design <- model.matrix(~ Group, data = data_subset)
  
  # Fit model
  fit <- lmFit(expr, design)
  fit <- eBayes(fit)
  
  deg <- topTable(fit, coef = 2, number = Inf)
  
  # Map UniProt identifiers to gene symbols for interpretation and downstream figure labeling.
  deg$symbol <- mapIds(
    org.Hs.eg.db,
    keys = rownames(deg),
    keytype = "UNIPROT",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  # Save results
  write.csv(deg, output_file)
  
  return(deg)
}

# 8. Volcano plots
# This code was used for the following figures 
# Figure 1 b-d
# Figure 3b
# Figure 5d
# Figure 6d
# Supplementary Figure 1f
# Supplementary Figure 2d-g
# ============================================================

# Significance thresholds used throughout the manuscript:
# log2FC >= 1
# adjusted P-value <= 0.05

plot_volcano<- function(df,
                        logFC_cutoff = 1,
                        pval_cutoff = 0.05,
                        label_n_up = 6,
                        label_n_down = 5,
                        title = "",
                        legend_title = "Expression") {
  
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  
  # classify
  df <- df %>%
    mutate(Expression = case_when(
      logFC >= logFC_cutoff & adj.P.Val <= pval_cutoff ~ "Up-regulated",
      logFC <= -logFC_cutoff & adj.P.Val <= pval_cutoff ~ "Down-regulated",
      TRUE ~ "Unchanged"
    ))
  
  # Select top positively and negatively associated proteins for figure annotation.
  top_up <- df %>%
    filter(Expression == "Up-regulated") %>%
    top_n(label_n_up, logFC)
  
  top_down <- df %>%
    filter(Expression == "Down-regulated") %>%
    top_n(label_n_down, -logFC)
  
  top_genes <- bind_rows(top_up, top_down)
  
  # plot
  p <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = Expression)) +
    
    geom_point() +
    
    scale_color_manual(
      values = c(
        "Up-regulated" = "red",
        "Down-regulated" = "blue",
        "Unchanged" = "grey"
      ),
      name = legend_title
    ) +
    
    geom_hline(
      yintercept = -log10(pval_cutoff),
      linetype = "dashed",
      color = "black"
    ) +
    
    geom_vline(
      xintercept = c(-logFC_cutoff, logFC_cutoff),
      linetype = "dashed",
      color = "black"
    ) +
    
    coord_cartesian(xlim = c(-5, 5)) +
    
    geom_text_repel(
      data = top_genes,
      aes(label = symbol),
      size = 5,
      box.padding = 0.35,
      point.padding = 0.3,
      segment.color = "grey50",
      show.legend = FALSE
    ) +
    
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10 Adjusted P-Value"
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18),
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 18)
    )
  
  return(p)
}

# ============================================================
# 04_GSEA_analysis.R
# ============================================================

library(tidyverse)
library(GSVA)
library(msigdbr)
library(org.Hs.eg.db)
library(AnnotationDbi)

set.seed(1234)

# ============================================================
# NOTE:
#
# This script performs pathway enrichment using GSVA/ssGSEA,
# which is the primary method used to derive sample-level and
# cluster-level pathway scores in this study (e.g. Fig 4, Fig 5).
#
# Other enrichment approaches are implemented separately:
# - fgsea: used for ranked enrichment plots (Fig 1e, Fig 6E)
# - WebGestalt: used for selected analyses (Fig 3, Supp Fig 1g)
# - ORA (clusterProfiler): used for lipid enrichment (Fig 5)
#
# ============================================================

# ============================================================
# 1. Load expression data
# ============================================================

data <- read.csv("../data/combat_df.csv")
data_MOC <- combat_df %>% filter(Label == "Mucinous_ovarian")

# Exclude re-labelled MOC sample 
data_MOC <- data_MOC %>% filter(Sample_ID !="03U0T")

# ============================================================
# 2. Build expression matrix
# ============================================================

expr <- data %>%
  select(where(is.numeric)) %>%
  as.matrix()

rownames(expr) <- data$Sample_ID
expr <- t(expr)

# ============================================================
# 3. Map UniProt → Entrez required by MSigDB Hallmark gene sets
# ============================================================

entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = rownames(expr),
  keytype = "UNIPROT",
  column = "ENTREZID",
  multiVals = "first"
)

expr <- expr[!is.na(entrez_ids), ]
rownames(expr) <- entrez_ids[!is.na(entrez_ids)]
expr <- expr[!duplicated(rownames(expr)), ]

# ============================================================
# 4. Load Hallmark gene sets for pathway-level analyses
# ============================================================

hallmark <- msigdbr("Homo sapiens", category = "H")
gene_sets <- split(hallmark$entrez_gene, hallmark$gs_name)

# ============================================================
# 5. Run ssGSEA (GSVA)
# ============================================================

params <- ssgseaParam(
  exprData = expr,
  geneSets = gene_sets,
  normalize = FALSE,
  alpha = 0.25
)

ssgsea_res <- gsva(params)

# ============================================================
# 6. Save output for downstream cluster analyses
# ============================================================

write.csv(
  ssgsea_res,
  "../results/GSVA_ssGSEA_scores.csv"
)

# ============================================================
# 05_celltype_enrichment.R
# Perform cell type enrichment analysis using Human Cell Landscape reference and fgsea
## Figure 1f
# ============================================================
library(GSEABase)
library(AnnotationDbi)
library(tidyverse)
library(fgsea)

dep <- read.csv("../results/DEP_MOC_vs_HGSOC.csv")

dep$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys = dep$RowNames, #RowNames = Uniprot identifiers
  keytype = "UNIPROT",
  column = "SYMBOL",
  multiVals = "first"
)

dep <- dep %>%
  filter(!is.na(SYMBOL))

ranked <- dep %>%
  group_by(SYMBOL) %>%
  summarise(logFC = logFC[which.max(abs(logFC))]) %>%
  deframe() %>%
  sort(decreasing = TRUE)

# Human Cell Landscape reference signatures
gmt <- gmtPathways("../data/HCL_celltype_reference.gmt")

res <- fgsea(gmt, ranked, nperm = 1000)

write.csv(res, "../results/celltype_fgsea.csv")


# ============================================================
# 06_Fig1.R
# ============================================================
#
# Purpose:
# Generate Figure 1 and associated supplementary figures.
#
# Figures:
#   Figure 1B-D: Differential abundance analyses
#   Figure 1E: Hallmark pathway enrichment
#   Figure 1F: Cell-type enrichment
#   

library(tidyverse)
library(limma)
library(Biobase)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(fgsea)
library(msigdbr)
library(data.table)
library(pheatmap)
library(ggrepel)


# ============================================================
# Filter samples
# Retain samples meeting tumour content criteria
# ============================================================
combat_df <- read.csv("../data/combat_Omicodf5625.csv")
data_filtered <- combat_Omicodf5625 %>%
  filter(Include_tumour_content == "Y")
data_filtered <- data_filtered %>%
  dplyr::select(-any_of(liver_genes))


# Differential expression analyses for Figure 1
# ============================================================

# Create output directories
dir.create("../results", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

run_DEP <- function(data_filtered, group_vector, output_prefix, title_text) {
  
  Biopsy_dummies <- data.frame(data_filtered, Group = group_vector)
  
  # Expression matrix
  X_matrix <- Biopsy_dummies %>% dplyr::select(where(is.numeric), -Group)
  X_matrix <- t(as.matrix(X_matrix))
  
  # Metadata
  p_risk <- Biopsy_dummies %>% dplyr::select(-where(is.numeric))
  pdf_risk <- as.data.frame(p_risk)
  
  colnames(X_matrix) <- rownames(pdf_risk)
  
  eset <- ExpressionSet(
    assayData = X_matrix,
    phenoData = AnnotatedDataFrame(pdf_risk)
  )
  
  # Design
  design <- model.matrix(~ Group, data = pdf_risk)
  
  # Standard limma model
  fit <- lmFit(eset, design)
  fit <- eBayes(fit)
  
  # Results
  res <- topTable(fit, coef = 2, number = Inf)
  
  res$RowNames <- rownames(res)
  
  # Mapping
  res$symbol <- mapIds(
    org.Hs.eg.db,
    keys = res$RowNames,
    keytype = "UNIPROT",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  # Classification
  res <- res %>%
    mutate(Expression = case_when(
      logFC >= 1 & adj.P.Val <= 0.05 ~ "Up-regulated",
      logFC <= -1 & adj.P.Val <= 0.05 ~ "Down-regulated",
      TRUE ~ "Unchanged"
    ))
  
  # Save table
  write.csv(res, paste0("../results/", output_prefix, ".csv"), row.names = FALSE)
  
  # Volcano (your style)
  top_up <- res %>% filter(Expression == "Up-regulated") %>% top_n(6, logFC)
  top_down <- res %>% filter(Expression == "Down-regulated") %>% top_n(5, -logFC)
  top_genes <- bind_rows(top_up, top_down)
  
  p <- ggplot(res, aes(x = logFC, y = -log10(adj.P.Val), color = Expression)) +
    geom_point() +
    scale_color_manual(
      values = c(
        "Up-regulated" = "red",
        "Down-regulated" = "blue",
        "Unchanged" = "grey"
      )
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    coord_cartesian(xlim = c(-5,5)) +
    geom_text_repel(
      data = top_genes,
      aes(label = symbol),
      size = 5
    ) +
    theme_minimal()
  
  ggsave(paste0("figures/", output_prefix, ".pdf"), p, width = 10, height = 8)
  
  return(res)
}

# Figure 1b: MOC versus all other samples
df <- data_filtered %>%
  filter(Label %in% c(
    "Mucinous_ovarian",
    "Appendix_adenocarcinoma",
    "Colorectal_adenocarcinoma",
    "Pancreas_adenocarcinoma",
    "UpperGI_adenocarcinoma",
    "Small_bowel_adenocarcinoma"
  ))

Group <- ifelse(df$Label == "Mucinous_ovarian", 1, 0)

res_Fig1B <- run_DEP(
  df,
  Group,
  "Fig1B_MOC_vs_all",
  "MOC vs all extra-ovarian tumours"
)

# Figure 1C: MOC vs. GEJ
df <- data_filtered %>%
  filter(Label %in% c("Mucinous_ovarian", "UpperGI_adenocarcinoma"))

Group <- ifelse(df$Label == "Mucinous_ovarian", 1, 0)

res_Fig1C <- run_DEP(
  df,
  Group,
  "Fig1C_MOC_vs_GEJ",
  "MOC vs gastro-esophageal adenocarcinoma"
)

# Figure 1d: MOC vs. PDAC
df <- data_filtered %>%
  filter(Label %in% c("Mucinous_ovarian", "Pancreas_adenocarcinoma"))

Group <- ifelse(df$Label == "Mucinous_ovarian", 1, 0)

res_Fig1D <- run_DEP(
  df,
  Group,
  "Fig1D_MOC_vs_PDAC",
  "MOC vs pancreatic adenocarcinoma"
)

# ============================================================
# Fig 1E: GSVA pathway enrichment across cancer diagnoses
# ============================================================

set.seed(1234)

data_filtered <- read.csv("../data/combat_Omicodf5625.csv")
data_filtered <- combat_Omicodf5625 %>% filter(Include_tumour_content=="Y")

# ============================================================
# 1. Prepare expression matrix
# ============================================================

labels <- data_filtered$Label

expression_data <- data_filtered %>%
  dplyr::select(where(is.numeric)) %>%
  as.matrix()

expression_data <- t(expression_data)
colnames(expression_data) <- data_filtered$Sample_ID

# ============================================================
# 2. Compute mean expression per diagnosis
# ============================================================

label_levels <- unique(labels)

label_means <- sapply(label_levels, function(label) {
  idx <- which(labels == label)
  if (length(idx) > 1) {
    rowMeans(expression_data[, idx, drop = FALSE])
  } else {
    expression_data[, idx]
  }
})

colnames(label_means) <- label_levels

# ============================================================
# 3. Map UniProt → Entrez
# ============================================================

entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = rownames(label_means),
  keytype = "UNIPROT",
  column = "ENTREZID",
  multiVals = "first"
)

label_means <- label_means[!is.na(entrez_ids), ]
rownames(label_means) <- entrez_ids[!is.na(entrez_ids)]
label_means <- label_means[!duplicated(rownames(label_means)), ]

# ============================================================
# 4. Load Hallmark gene sets
# ============================================================

hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
gene_sets <- split(hallmark_sets$entrez_gene, hallmark_sets$gs_name)

# ============================================================
# 5. Run GSVA (ssGSEA)
# ============================================================

params <- ssgseaParam(
  exprData = label_means,
  geneSets = gene_sets,
  normalize = FALSE,
  alpha = 0.25
)

ssgsea_results <- gsva(params)

# ============================================================
# 6. Heatmap visualisation
# ============================================================

ssgsea_results_sub <- ssgsea_results[, c(
  "Mucinous_ovarian",
  "Pancreas_adenocarcinoma",
  "UpperGI_adenocarcinoma"
)]

annotation_col <- data.frame(
  Diagnosis = colnames(ssgsea_results_sub)
)
rownames(annotation_col) <- colnames(ssgsea_results_sub)

pdf("figures/Fig1E_GSVA_heatmap.pdf", width = 8, height = 8)

pheatmap(
  ssgsea_results_sub,
  scale = "row",
  clustering_distance_rows = "correlation",
  cluster_cols = FALSE,
  clustering_method = "complete",
  annotation_col = annotation_col,
  main = "GSVA enrichment across cancer diagnoses",
  fontsize_row = 6,
  fontsize_col = 8
)

dev.off()

# =========================
# Figure 1f_celltype_enrichment.R
# =========================

# =========================
# 1. Load DEP results
# =========================

dep_list <- list(
  MOC_vs_all = read.csv("../results/Fig1B_MOC_vs_all.csv"),
  MOC_vs_PDAC = read.csv("../results/Fig1D_MOC_vs_PDAC.csv"),
  MOC_vs_GEJ = read.csv("../results/Fig1C_MOC_vs_GEJ.csv")
)

# =========================
# 2. Map UniProt → gene SYMBOL
# =========================

map_uniprot_to_symbol <- function(df) {
  
  mapping <- mapIds(
    org.Hs.eg.db,
    keys = df$RowNames, #RowNames = Uniprot identifier
    keytype = "UNIPROT",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  df$SYMBOL <- mapping
  
  df %>%
    filter(!is.na(SYMBOL))
}

dep_list <- lapply(dep_list, map_uniprot_to_symbol)

# =========================
# 3. Create ranked lists
# =========================

ranked_lists <- lapply(dep_list, function(df) {
  
  ranked <- df %>%
    group_by(SYMBOL) %>%
    summarise(logFC = logFC[which.max(abs(logFC))]) %>%
    deframe()
  
  sort(ranked, decreasing = TRUE)
})

# =========================
# 4. Load HCL gene sets (GMT)
# =========================

gmt <- gmtPathways("../data/HCL_celltype_reference.gmt")

# =========================
# 5. Run fgsea
# =========================

gsea_results <- lapply(names(ranked_lists), function(name) {
  
  res <- fgsea(
    pathways = gmt,
    stats = ranked_lists[[name]],
    minSize = 10,
    maxSize = 500,
    nperm = 1000
  )
  
  res$CancerType <- name
  res
})

gsea_df <- bind_rows(gsea_results)

# =========================
# 6. Filter + format
# =========================

plot_df <- gsea_df %>%
  filter(padj < 0.05) %>%
  group_by(CancerType) %>%
  slice_max(NES, n = 10) %>%
  ungroup()

# clean labels
plot_df$pathway <- gsub("_", " ", plot_df$pathway)

# =========================
# 7. Plot (Fig 1f)
# =========================

p <- ggplot(plot_df, aes(x = CancerType, y = pathway)) +
  
  geom_point(aes(size = -log10(padj), color = NES)) +
  
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  
  scale_size_continuous(range = c(3, 8)) +
  
  theme_classic(base_size = 16) +
  
  labs(
    x = "Cancer type",
    y = "Cell type signature",
    size = "-log10(adj p-value)",
    color = "NES"
  ) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("figures/Fig1f_celltype_enrichment.pdf", p, width = 8, height = 6)

# ============================================================
# Supplementary Fig 1g
# ============================================================
# NOTE:
# GSEA results were generated using WebGestalt (online tool).
# Pathways shown here represent manually curated significant
# results exported from WebGestalt output and plotted in R
# for visualization.
# ============================================================

library(tidyverse)
library(ggplot2)
library(stringr)

# ============================================================
# 1. Load WebGestalt results
# ============================================================

reactome_gsea <- read.csv("../results/webgestalt_reactome.csv")
hallmark_gsea <- read.csv("../results/webgestalt_hallmark.csv")

# ============================================================
# 2. Select pathways of interest
# ============================================================

hallmark_keep <- c(
  "HALLMARK_COMPLEMENT",
  "HALLMARK_COAGULATION",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_PROTEIN_SECRETION"
)

reactome_keep <- c(
  "Antimicrobial peptides",
  "O-linked glycosylation of mucins",
  "O-linked glycosylation",
  "Neutrophil degranulation",
  "Activation of Matrix Metalloproteinases",
  "Degradation of the extracellular matrix",
  "Interferon alpha/beta signaling",
  "Nitric oxide stimulates guanylate cyclase"
)

# ============================================================
# 3. Prepare plotting data
# ============================================================

hallmark_sel <- hallmark_gsea %>%
  filter(GeneSet %in% hallmark_keep) %>%
  transmute(
    pathway = Description,
    NES = NES,
    FDR = FDR,
    source = "Hallmark"
  )

reactome_sel <- reactome_gsea %>%
  filter(Description %in% reactome_keep) %>%
  transmute(
    pathway = Description,
    NES = NES,
    FDR = FDR,
    source = "Reactome"
  )

gsea_plot_df <- bind_rows(hallmark_sel, reactome_sel)

# ============================================================
# 4. Clean pathway labels
# ============================================================

gsea_plot_df <- gsea_plot_df %>%
  mutate(
    pathway_label = case_when(
      source == "Hallmark" ~ str_to_sentence(pathway),
      source == "Reactome" ~ str_wrap(pathway, width = 30)
    )
  )

# ============================================================
# 5. Plot
# ============================================================

p <- ggplot(
  gsea_plot_df,
  aes(
    x = NES,
    y = reorder(pathway_label, NES),
    size = -log10(FDR),
    colour = NES > 0
  )
) +
  
  geom_point() +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  
  facet_wrap(~ source) +
  
  scale_colour_manual(
    values = c(
      "TRUE"  = "#D55E00",  # enriched
      "FALSE" = "#0072B2"   # depleted
    ),
    guide = "none"
  ) +
  
  scale_size_continuous(
    name  = "-log10(adj p-value)",
    range = c(5, 12)
  ) +
  
  scale_x_continuous(
    expand = expansion(mult = c(0.15, 0.15))
  ) +
  
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  
  theme_classic(base_size = 16) +
  theme(
    text = element_text(color = "black"),
    legend.position = "bottom",
    strip.text = element_text(size = 14, face = "bold"),
    plot.margin = margin(15, 25, 15, 25)
  )

# ============================================================
# 6. Save
# ============================================================

ggsave(
  "figures/SuppFig1G_GSEA_dotplot.pdf",
  plot = p,
  width = 12,
  height = 10
)

# ============================================================
# Supplementary Figure 2a–c: STRING network analysis
# ------------------------------------------------------------

# Purpose:
# Generate Supplementary Figure 2 and export protein lists
# used for STRING network analysis.
#
# Panels:
#   Supplementary Figure 2A-C
#       STRING network inputs
#
#   Supplementary Figure 2D-F
#       Volcano plots generated using duplicateCorrelation
#
#   Supplementary Figure 2G
#       MOC versus HGSOC differential abundance analysis
#
# Inputs:
#   ../data/analysis_ready.csv
#   Figure 1 differential abundance results
#
# Outputs:
#   STRING input files
#   Volcano plots
#   HGSOC comparison results
#
# ============================================================

# NOTE:
# STRING networks were generated using STRING-db (online tool)
# based on top 150 differentially expressed proteins derived
# from the following comparisons:
#
# a) MOC vs PDAC          → tumour-specific pancreatic biology
# b) MOC vs GEJ           → tumour-specific gastro-esophageal biology
# c) MOC vs UGI           → combined GEJ + PDAC (shared upper GI biology)
#
# The UGI comparison was performed specifically to identify
# shared proteomic features across upper gastrointestinal tumors.
# ============================================================

# NOTE:
# run_DEP() is defined in 06_Fig1.R and should be sourced prior to execution of this script.

data <- read.csv("../data/analysis_ready.csv")

library(tidyverse)

# ============================================================
# 1. Load DEP results
# ============================================================

deg_pdac <- read.csv("../results/Fig1D_MOC_vs_PDAC.csv")   # a)
deg_gej  <- read.csv("../results/Fig1C_MOC_vs_GEJ.csv")    # b)

# ============================================================
# 2. Generate UGI comparison (c)
# ============================================================

# Combine GEJ + PDAC as upper GI group
df_ugi <- data %>%
  filter(Label %in% c(
    "Mucinous_ovarian",
    "UpperGI_adenocarcinoma",
    "Pancreas_adenocarcinoma"
  ))

group_ugi <- ifelse(df_ugi$Label == "Mucinous_ovarian", 1, 0)

deg_ugi <- run_DEP(
  df_ugi,
  group_ugi,
  "STRING_MOC_vs_UGI",
  "MOC vs upper gastrointestinal tumours (combined GEJ + PDAC)"
)

# ============================================================
# 3. Function: extract top 150 proteins
# ============================================================

get_top150 <- function(df, direction = "down") {
  
  df %>%
    filter(adj.P.Val <= 0.05) %>%
    arrange(if (direction == "down") logFC else desc(logFC)) %>%
    slice_head(n = 150)
}

# ============================================================
# 4. Extract protein sets for STRING
# ============================================================

# a) PDAC-driven biology (proteins enriched in PDAC vs MOC)
top_pdac <- get_top150(deg_pdac, "down")

# b) GEJ-driven biology (proteins enriched in GEJ vs MOC)
top_gej  <- get_top150(deg_gej, "down")

# c) Shared upper GI biology (proteins enriched in UGI vs MOC)
top_ugi  <- get_top150(deg_ugi, "down")

# ============================================================
# 5. Export for STRING input
# ============================================================

write.csv(top_pdac$symbol,
          "../results/SuppFig2A_PDAC_top150.csv",
          row.names = FALSE)

write.csv(top_gej$symbol,
          "../results/SuppFig2B_GEJ_top150.csv",
          row.names = FALSE)

write.csv(top_ugi$symbol,
          "../results/SuppFig2C_UGI_top150.csv",
          row.names = FALSE)


# ============================================================
# duplicateCorrelation sensitivity analyses (Supp Fig 2d-f)
# ============================================================

combat_df <- read.csv("../data/combat_Omicodf5625.csv")
data_filtered <- combat_Omicodf5625 %>%
  filter(Include_tumour_content == "Y")
data_filtered <- data_filtered %>%
  dplyr::select(-any_of(liver_genes))

run_DEP_dupCor <- function(data, group_vector, output_prefix, title_text) {
  
  Biopsy_dummies <- data.frame(data, Group = group_vector)
  
  # Expression matrix
  X_matrix <- Biopsy_dummies %>%
    dplyr::select(where(is.numeric), -Group)
  
  X_matrix <- data.matrix(X_matrix)
  X_matrix <- t(X_matrix)
  
  # Metadata
  p_risk <- Biopsy_dummies %>%
    dplyr::select(-where(is.numeric))
  
  pdf_risk <- as.data.frame(p_risk)
  
  colnames(X_matrix) <- rownames(pdf_risk)
  
  eset <- ExpressionSet(
    assayData = X_matrix,
    phenoData = AnnotatedDataFrame(pdf_risk)
  )
  
  # Design matrix
  design <- model.matrix(~ Group, data = pdf_risk)
  
  # duplicateCorrelation
  block <- pdf_risk$sample_site
  
  corfit <- duplicateCorrelation(
    eset,
    design = design,
    block = block
  )
  
  fit <- lmFit(
    eset,
    design = design,
    block = block,
    correlation = corfit$consensus.correlation
  )
  
  fit <- eBayes(fit)
  
  # Results
  res <- topTable(fit, coef = 2, number = Inf)
  
  res$RowNames <- rownames(res)
  
  # Annotation
  res$symbol <- mapIds(
    org.Hs.eg.db,
    keys = res$RowNames,
    keytype = "UNIPROT",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  # Classification
  res <- res %>%
    mutate(
      Expression = case_when(
        logFC >= 1 & adj.P.Val <= 0.05 ~ "Up-regulated",
        logFC <= -1 & adj.P.Val <= 0.05 ~ "Down-regulated",
        TRUE ~ "Unchanged"
      )
    )
  
  # Save results
  write.csv(
    res,
    paste0("../results/", output_prefix, ".csv"),
    row.names = FALSE
  )
  
  return(res)
}
    
# Supp Fig 2d - MOC vs all

df <- data_filtered %>%
  filter(Label %in% c(
    "Mucinous_ovarian",
    "Appendix_adenocarcinoma",
    "Colorectal_adenocarcinoma",
    "Pancreas_adenocarcinoma",
    "UpperGI_adenocarcinoma",
    "Small_bowel_adenocarcinoma"
  ))

Group <- ifelse(df$Label == "Mucinous_ovarian", 1, 0)

deg_all_dupCor <- run_DEP_dupCor(
  df,
  Group,
  "DEP_dupCor_Fig1B_MOC_vs_all",
  "MOC vs all extra-ovarian tumours"
)

# Supp Fig 2e - MOC vs PDAC

df <- data_filtered %>%
  filter(Label %in% c(
    "Mucinous_ovarian",
    "Pancreas_adenocarcinoma"
  ))

Group <- ifelse(df$Label == "Mucinous_ovarian", 1, 0)

deg_pdac_dupCor <- run_DEP_dupCor(
  df,
  Group,
  "DEP_dupCor_Fig1D_MOC_vs_PDAC",
  "MOC vs PDAC"
)

# Supp Fig 2f - MOC vs GEJ

df <- data_filtered %>%
  filter(Label %in% c(
    "Mucinous_ovarian",
    "UpperGI_adenocarcinoma"
  ))

Group <- ifelse(df$Label == "Mucinous_ovarian", 1, 0)

deg_gej_dupCor <- run_DEP_dupCor(
  df,
  Group,
  "DEP_dupCor_Fig1C_MOC_vs_GEJ",
  "MOC vs GEJ"
)

# Supplementary Figures 2d-f
p_all <- plot_volcano(deg_all_dupCor,
                      title = "MOC vs all extra-ovarian tumours (duplicateCorrelation)")

p_pdac <- plot_volcano(deg_pdac_dupCor,
                       title = "MOC vs PDAC (duplicateCorrelation)")

p_gej <- plot_volcano(deg_gej_dupCor,
                      title = "MOC vs GEJ (duplicateCorrelation)")

ggsave("figures/SuppFig2d_MOC_vs_all.pdf", p_all, width=10, height=8)
ggsave("figures/SuppFig2e_MOC_vs_PDAC.pdf", p_pdac, width=10, height=8)
ggsave("figures/SuppFig2f_MOC_vs_GEJ.pdf", p_gej, width=10, height=8)


# Supplementary Figure 2g: MOC vs. HGSOC
df <- data %>%
  filter(Label %in% c("Mucinous_ovarian", "High_grade_serous_ovary"))

Group <- ifelse(df$Label == "Mucinous_ovarian", 1, 0)

res_Supp2g <- run_DEP(
  df,
  Group,
  "DEP_MOC_vs_HGSOC",
  "MOC vs high-grade serous ovarian cancer"
)

# ============================================================
# 07_Fig3: Clinical and Histological Features of MOC
# ============================================================
#
# Figures:
#   Figure 3A-E : UMAP visualisation of clinical variables
#   Figure 3B   : Infiltrative versus expansile differential abundance
#   Figure 3C   : Cell-type enrichment (WebGestalt-derived)
#   Supplementary Figure 3H : Werner dataset validation
#
# Inputs:
#   ../data/MOC_cluster_input.csv
#   ../results/DEP_MOC_pattern.csv
#   ../data/Werner_MC_vs_HGSC_logFC.csv
#
# ============================================================

library(tidyverse)
library(umap)
library(patchwork)
library(ggrepel)
library(fgsea)
library(GSEABase)

# ============================================================
# 1. Load data
# ============================================================

# NOTE:
# This dataset corresponds to the same MOC cohort used for clustering (Fig 4)
# uses only MOC samples with > 20% cancer content

# Generate Figure 4 clustering input
# Derived from merged analysis matrix by:
# - retaining MOC samples only
# - retaining samples with >20% tumour content

data <- combat_df %>%
  filter(
    Label == "Mucinous_ovarian",
    Include_in_study == "Y"
  )

MOC_clinicaldata <- MOC_clinicaldata %>% filter(Tumour_percent != "0")

# ============================================================
# 2. Run UMAP (fixed parameters)
# ============================================================
set.seed(1234)
config <- umap.defaults
config$n_neighbors <- 5
config$min_dist <- 0.001

umap_fit <- data %>%
  dplyr::select(where(is.numeric)) %>%
  umap(config = config)

umap_df <- as.data.frame(umap_fit$layout)

umap_df <- as.data.frame(umap_fit$layout) %>%
  bind_cols(
    MOC_clinicaldata %>%
      dplyr::select(
        Sample_ID,
        ER,
        KRAS,
        ERBB2.any,
        TP53.any,
        Grade,
        Stage,
        Histological_pattern
      )
  )

colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")

# ============================================================
# 3. Clean variables
# ============================================================

umap_df <- umap_df %>%
  mutate(
    Stage_display = ifelse(Stage %in% c("1", "2", "3"),
                           Stage, "Unknown"),
    Grade_display = ifelse(Grade %in% c("1", "2", "3"),
                           Grade, "Unknown"),
    Pattern_display = ifelse(
      Histological_pattern %in% c("Infiltrative", "Expansile", "Mixed"),
      Histological_pattern, "Unknown"
    )
  )

umap_df <- umap_df %>%
  mutate(
    ER_Status = case_when(
      ER %in% c("Positive", "Negative") ~ ER,
      TRUE ~ "Unknown"
    ),
    ERBB2_Status = case_when(
      ERBB2.any %in% c("Amplified", "HER2+", "AMP") ~ "Amplified",
      ERBB2.any %in% c("Not amplified", "HER2-", "NOT AMP") ~ "Not amplified",
      TRUE ~ "Unknown"
    )
  )

# ============================================================
# 4. Plot helper
# ============================================================

plot_umap <- function(df, var) {
  
  ggplot(
    df %>% filter(.data[[var]] != "Unknown"),
    aes(UMAP1, UMAP2, color = .data[[var]])
  ) +
    geom_point(size = 2, alpha = 0.8) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white"),
      axis.line = element_line(color = "grey")
    )
}

# ============================================================
# 5. Generate panels (Fig 3A–E)
# ============================================================

p1 <- plot_umap(umap_df, "ERBB2_Status")
p2 <- plot_umap(umap_df, "Stage_display")
p3 <- plot_umap(umap_df, "Grade_display")
p4 <- plot_umap(umap_df, "Pattern_display")
p5 <- plot_umap(umap_df, "ER_Status")

library(patchwork)

final_umap <- (p1 | p2) / (p3 | p4) / (p5 | plot_spacer())

ggsave(
  "figures/Fig3A-UMAP_clinical.pdf",
  final_umap,
  width = 10,
  height = 12
)

# ============================================================
# -----------------------
# Fig 3B: Volcano
# Expansile vs Infiltrative
# -----------------------
# ============================================================

library(dplyr)
library(limma)
library(Biobase)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ============================================================
# Expansile vs Infiltrative MOC DEP analysis
# ============================================================

Infiltrative <- combat_Omicodf5625 %>%
  filter(
    Label == "Mucinous_ovarian",
    Include_tumour_content == "Y"
  ) %>%
  left_join(
    MOC_clinicaldata %>%
      dplyr::select(Sample_ID, Histological_pattern),
    by = "Sample_ID"
  ) %>%
  filter(
    Histological_pattern %in% c("Expansile", "Infiltrative")
  )

Infiltrative_column <- ifelse(
  Infiltrative$Histological_pattern == "Infiltrative",
  1,
  0
)

Biopsy_dummies <- data.frame(
  Infiltrative,
  Infiltrative_column
)

# Protein matrix
X_matrix <- Biopsy_dummies[, 1:10416]
X_matrix <- data.matrix(X_matrix)
X_matrix <- t(X_matrix)

# Metadata
p_risk <- Biopsy_dummies[, 10417:ncol(Biopsy_dummies)]
pdf_risk <- as.data.frame(p_risk)

colnames(X_matrix) <- rownames(pdf_risk)

eset_risk <- ExpressionSet(
  assayData = X_matrix,
  phenoData = AnnotatedDataFrame(pdf_risk)
)

design_lymph <- model.matrix(
  ~ Infiltrative_column,
  data = pdf_risk
)

fit_risk <- lmFit(eset_risk, design_lymph)
fit_risk <- eBayes(fit_risk)

deg_pattern <- topTable(
  fit_risk,
  coef = 2,
  number = Inf
)

write.csv(
  deg_pattern,
  "../results/DEP_MOC_pattern.csv",
  row.names = TRUE
)

deg_pattern$UNIPROT <- rownames(deg_pattern)

deg_pattern$symbol <- mapIds(
  org.Hs.eg.db,
  keys = deg_pattern$UNIPROT,
  keytype = "UNIPROT",
  column = "SYMBOL",
  multiVals = "first"
)


deg_pattern <- deg_pattern %>%
  mutate(
    status = case_when(
      logFC >= 1 & adj.P.Val <= 0.05 ~ "Up-regulated",
      logFC <= -1 & adj.P.Val <= 0.05 ~ "Down-regulated",
      TRUE ~ "Unchanged"
    )
  ) %>%
  filter(!is.na(logFC), !is.na(adj.P.Val))

top_genes <- deg_pattern %>%
  filter(status != "Unchanged") %>%
  slice_max(abs(logFC), n = 10)

p_volcano <- ggplot(
  deg_pattern,
  aes(logFC, -log10(adj.P.Val), color = status)
) +
  geom_point() +
  scale_color_manual(
    values = c(
      "Up-regulated" = "red",
      "Down-regulated" = "blue",
      "Unchanged" = "grey"
    )
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  geom_text_repel(
    data = top_genes,
    aes(label = symbol),
    size = 4
  ) +
  theme_classic(base_size = 16)

ggsave(
  "figures/Fig3B_volcano_pattern.pdf",
  p_volcano,
  width = 6,
  height = 5
)

# ============================================================
# Fig 3C: Cell-type enrichment (WebGestalt-derived)
# ============================================================

library(tidyverse)
library(ggplot2)

# NOTE:
# Cell-type enrichment was performed using WebGestalt.
#
# Significant pathways and enrichment statistics were
# exported from WebGestalt and manually imported here
# for figure generation.
#
# This analysis cannot be reproduced entirely within R
# because enrichment was performed using the WebGestalt
# web interface.

celltype_df <- data.frame(
  GeneSet = c("Fetal stromal cell","Fetal Neuron","Astrocyte","Epithelial cell (intermediated)",
              "Epithelial cell","AT2 cell","Intercalated cell","Hepatocyte/Endodermal cell",
              "Gastric chief cell","Fetal enterocyte","Enterocyte progenitor","Goblet cell",
              "Enterocyte","Fetal mesenchymal progenitor","Fibroblast","Macrophage","Stromal cell",
              "Smooth muscle cell","Immature sertoli cell","Dendritic cell"),
  NES = c(2.0779,1.9337,1.9303,-1.9526,-2.0492,-2.1410,-2.1572,-2.3587,-2.5641,-2.5921,
          -2.7367,-3.0816,-3.4387,1.8593,1.8747,1.8205,1.5202,1.4632,1.4108,1.4292),
  FDR = c(rep(2.2e-16,13),0.00011399,0.00014248,0.00018998,0.0096075,0.019948,0.026929,0.027800)
)

celltype_df <- celltype_df %>%
  mutate(
    Direction = ifelse(NES > 0, "Infiltrative", "Expansile"),
    logFDR = -log10(FDR)
  ) %>%
  arrange(NES)

p_ct <- ggplot(
  celltype_df,
  aes(NES, reorder(GeneSet, NES), color = Direction, size = logFDR)
) +
  geom_point() +
  scale_color_manual(values = c("Infiltrative" = "red", "Expansile" = "blue")) +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    color = "Histological growth pattern",
    size = "-log10(FDR)"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )

ggsave("figures/Fig3C_celltype_enrichment.pdf", p_ct, width = 7, height = 5)


# ============================================================
# Supplementary Fig 3h: Werner dataset validation (MC vs HGSC)
# ============================================================

library(tidyverse)
library(ggplot2)

# NOTE:
# Differential abundance values are derived from the published
# Werner dataset (MC vs HGSC comparison).
# logFC values were extracted from supplementary material of this published dataset.

# ============================================================
# 1. Load Werner dataset results
# ============================================================

werner <- read.csv("../data/Werner_MC_vs_HGSC_logFC.csv")

# ============================================================
# 2. Define ML-selected proteins
# ============================================================

ml_proteins <- c("CEACAM5", "TSPAN8", "S100P", "CDH17", "GOLM1", "B4GALT4", "CA9", 
                 "CLCA1", "WT1", "BCAT1", "SPON1")  # Machine-learning selected proteins evaluated in the independent Werner cohort

plot_df <- werner %>%
  filter(Protein %in% ml_proteins)

# ============================================================
# 3. Plot
# ============================================================

p <- ggplot(plot_df,
            aes(x = logFC,
                y = reorder(Protein, logFC),
                color = logFC > 0)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  labs(
    x = "log2 Fold Change (MC vs HGSC)",
    y = "Protein",
    color = "Direction"
  ) +
  theme_minimal(base_size = 16)

ggsave(
  "figures/SuppFig3H_Werner_validation.pdf",
  p,
  width = 6,
  height = 5
)

# ============================================================
# 08_Fig4.R

# ============================================================
# Supplementary Fig 4: Consensus clustering diagnostics
# ============================================================

# NOTE:
# These plots are automatically generated by ConsensusClusterPlus
# and saved as PDF during clustering (02_consensus_clustering.R)
# They include:
# - CDF plots
# - Delta area plots
# - Consensus matrices (k = 2–8)

# ============================================================
# Consensus clustering of MOC samples
# ------------------------------------------------------------
# Clustering performed using ConsensusClusterPlus
# with PAM clustering and Spearman distance.
#
# Results used for:
# - Fig 4A (cluster heatmap)
# - Fig 4B (pathway enrichment)
# - Fig 5D–E (cluster analyses)
# - Fig 6A (ADC targets)
# ============================================================

# =========================
## Fig 4A Heatmap
# =========================

library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(stringr)
library(grid)
library(tidyverse)
library(GSVA)
library(msigdbr)
library(org.Hs.eg.db)
library(AnnotationDbi)

# ============================================================
# 1. Load data
# ============================================================
# data filtered to include only mucinous ovary

data_heatmap <- combat_Omicodf5625 %>%
  filter(
    Label == "Mucinous_ovarian",
    Include_tumour_content == "Y"
  ) %>%
  left_join(
    MOC_clinicaldata %>%
      dplyr::select(
        Sample_ID,
        Unsupervised_cluster,
        Grade,
        Stage.y,
        Histological_pattern,
        ERBB2.any,
        TP53.any,
        CDKN2A.any,
        KRAS,
        ER,
        Cancer_percent
      ),
    by = "Sample_ID"
  )

table(data_heatmap$Unsupervised_cluster)

expr <- data_heatmap %>%
  dplyr::select(where(is.numeric)) %>%
  as.matrix()

rownames(expr) <- data_heatmap$Sample_ID
expr <- t(expr)

# =========================
# 2. Build annotation table
# =========================
# Clinical, molecular and proteomic annotations displayed above the consensus clustering heatmap.
annotation_col <- data.frame(
  Cluster = data_heatmap$Unsupervised_cluster,
  Grade = data_heatmap$Grade,
  Stage = data_heatmap$Stage.y,
  `Histological Pattern` = data_heatmap$Histological_pattern,
  ERBB2 = scale(data_heatmap$P04626),
  TP53 = scale(data_heatmap$Q12888),
  SMAD4 = scale(data_heatmap$Q13485),
  PRGR = scale(data_heatmap$P06401),
  TGFBI = scale(data_heatmap$Q15582),
  `ERBB2 gene` = data_heatmap$ERBB2.any,
  `TP53 gene` = data_heatmap$TP53.any,
  `CDKN2A gene` = data_heatmap$CDKN2A.any,
  `KRAS gene` = data_heatmap$KRAS,
  `ER IHC` = data_heatmap$ER,
  `Cancer percent` = data_heatmap$Cancer_percent,
  check.names = FALSE
)

rownames(annotation_col) <- data_heatmap$Sample_ID

# =========================
# 3. Clean categorical variables
# =========================

clean_factor <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "NA"] <- "Unknown"
  factor(x)
}

annotation_col <- annotation_col %>%
  mutate(
    across(c(
      Cluster, Grade, Stage,
      `Histological Pattern`,
      `TP53 gene`, `CDKN2A gene`,
      `ERBB2 gene`, `KRAS gene`,
      `ER IHC`, `Cancer percent`
    ), clean_factor)
  )

# =========================
# 4. Stage simplification (I/II/III/IV)
# =========================

annotation_col$Stage <- case_when(
  str_detect(annotation_col$Stage, "^I[^I]") ~ "1",
  str_detect(annotation_col$Stage, "^II[^I]") ~ "2",
  str_detect(annotation_col$Stage, "^III") ~ "3",
  str_detect(annotation_col$Stage, "^IV") ~ "4",
  TRUE ~ "Unknown"
)

annotation_col$Stage <- factor(annotation_col$Stage)

# =========================
# 5. Order samples by cluster
# =========================

annotation_col$Cluster <- factor(annotation_col$Cluster)
ordered_samples <- rownames(annotation_col)[order(annotation_col$Cluster)]

expr_matrix <- expr[, ordered_samples]
annotation_col <- annotation_col[ordered_samples, ]

# =========================
# 6. Select high-variance proteins
# =========================
# Restrict heatmap to the top 5% most variable proteins to improve visualization of cluster-specific structure.

gene_mads <- apply(expr_matrix, 1, mad, na.rm = TRUE)
cutoff <- quantile(gene_mads, 0.95)

expr_matrix <- expr_matrix[gene_mads > cutoff, ]

# =========================
# 7. Scale expression
# =========================

expr_matrix_scaled <- t(scale(t(expr_matrix)))

# =========================
# 8. Annotation colours
# =========================

ann_colors <- list(
  
  Cluster = c(
    "1" = "#66C2A5",
    "2" = "#FC8D62",
    "3" = "#8DA0CB",
    "4" = "#E78AC3",
    "5" = "#A6D854"
  ),
  
  Grade = c("1"="#66c2a5","2"="#fc8d62","3"="#8da0cb","Unknown"="#CCCCCC"),
  
  Stage = c("1"="#e78ac3","2"="#a6d854","3"="#ffd92f","4"="#e5c494","Unknown"="#CCCCCC"),
  
  `Histological Pattern` = c(
    "Infiltrative"="#b3b",
    "Expansile"="#beb",
    "Unknown"="#CCCCCC"
  ),
  
  `Cancer percent` = c(
    "<20"="white","20-50"="yellow","50-80"="orange",">80"="red","Unknown"="#CCCCCC"
  ),
  
  `ERBB2 gene` = c(
    "Amplified"="navy","Not amplified"="pink","Unknown"="#CCCCCC"
  ),
  
  `TP53 gene` = c(
    "No"="blue","Yes"="red","Unknown"="#CCCCCC"
  ),
  
  `CDKN2A gene` = c(
    "No"="orange","Yes"="lightgreen","Unknown"="#CCCCCC"
  ),
  
  `KRAS gene` = c(
    "MUT"="lightblue","WT"="black","Unknown"="#CCCCCC"
  ),
  
  `ER IHC` = c(
    "Negative"="yellow","Positive"="purple","Unknown"="#CCCCCC"
  )
)

# continuous variables
ann_colors$ERBB2 <- colorRamp2(c(-2, 0, 2), c("blue","white","red"))
ann_colors$TP53 <- colorRamp2(c(-2, 0, 2), c("blue","white","red"))
ann_colors$SMAD4 <- colorRamp2(c(-2, 0, 2), c("blue","white","red"))
ann_colors$PRGR <- colorRamp2(c(-2, 0, 2), c("blue","white","red"))
ann_colors$TGFBI <- colorRamp2(c(-2, 0, 2), c("blue","white","red"))

# =========================
# 9. Heatmap annotation
# =========================

ha <- HeatmapAnnotation(
  df = annotation_col,
  col = ann_colors,
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 12, fontface = "bold")
)

# =========================
# 10. Heatmap
# =========================

ht <- Heatmap(
  expr_matrix_scaled,
  name = "Expression",
  col = colorRamp2(c(-4, 0, 4), c("blue","white","red")),
  top_annotation = ha,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_split = annotation_col$Cluster,
  heatmap_legend_param = list(
    title = "Relative Protein Expression"
  )
)

# =========================
# 11. Save
# =========================

pdf("figures/Fig4A_heatmap.pdf", width = 15, height = 15)
draw(ht)
dev.off()

# =========================
# Fig 4B: Hallmark pathways (ssGSEA per cluster)
# =========================

# =========================
# 1. Load data
# =========================

data_hallmark_heatmap <- combat_Omicodf5625 %>%
  filter(
    Label == "Mucinous_ovarian",
    Include_tumour_content == "Y"
  ) %>%
  left_join(
    MOC_clinicaldata %>%
      dplyr::select(
        Sample_ID,
        Unsupervised_cluster
      ),
    by = "Sample_ID")

data_hallmark_heatmap$Unsupervised_cluster <- factor(
  data_hallmark_heatmap$Unsupervised_cluster,
  levels = c(1,2,3,4,5)
)

# =========================
# 2. Expression matrix
# =========================

expr <- data_hallmark_heatmap %>%
  dplyr::select(where(is.numeric)) %>%
  as.matrix()

rownames(expr) <- data_hallmark_heatmap$Sample_ID
expr <- t(expr)

# =========================
# 3. Map UniProt → Entrez
# =========================

entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = rownames(expr),
  keytype = "UNIPROT",
  column = "ENTREZID",
  multiVals = "first"
)

expr <- expr[!is.na(entrez_ids), ]
rownames(expr) <- entrez_ids[!is.na(entrez_ids)]
expr <- expr[!duplicated(rownames(expr)), ]

# =========================
# 4. Average per cluster
# =========================
# Average protein abundance across samples within each consensus cluster prior to pathway enrichment analysis.

cluster_vec <- data_hallmark_heatmap$Unsupervised_cluster
# clusterspearman = consensus cluster assignment (k = 5)

valid_clusters <- unique(cluster_vec[!is.na(cluster_vec)])

cluster_means <- sapply(valid_clusters, function(cl) {
  idx <- which(cluster_vec == cl)
  if (length(idx) > 1) {
    rowMeans(expr[, idx, drop = FALSE])
  } else {
    expr[, idx]
  }
})
colnames(cluster_means) <- valid_clusters

# =========================
# 5. Hallmark gene sets
# =========================

hallmark <- msigdbr(species = "Homo sapiens", category = "H")
gene_sets <- split(hallmark$entrez_gene, hallmark$gs_name)

# =========================
# 6. ssGSEA
# =========================

params <- ssgseaParam(
  exprData = cluster_means,
  geneSets = gene_sets,
  normalize = FALSE,
  alpha = 0.25
)

ssgsea <- gsva(params)

# =========================
# 7. Scale (row-wise)
# =========================

ssgsea_scaled <- t(scale(t(ssgsea)))

# =========================
# 8. Plot heatmap
# =========================

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

ht <- Heatmap(
  ssgsea_scaled,
  name = "Enrichment",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE
)

pdf("figures/Fig4B_ssGSEA_hallmark_clusters.pdf", width = 8, height = 10)
draw(ht)
dev.off()

# ============================================================
# 09_Fig5.R
# ============================================================
#
# Purpose:
# Generate Figure 5 and associated supplementary analyses.
#
# Panels:
#   Figure 5A : ERBB2 protein expression by amplification status
#   Figure 5B : Differential abundance (ERBB2 amplified vs non-amplified)
#   Figure 5C : Lipid-related GO enrichment
#   Figure 5D : Proliferation markers across clusters
#   Figure 5E : KRAS pathway enrichment across clusters
#
# Supplementary:
#   TCGA-BRCA validation
#   TCGA-STAD validation
#
# Inputs:
#   ../data/MOC_cluster_input.csv
#   ../results/MOC_clusters_k5.csv
#   ../results/DEP_Fig5_ERBB2.csv
#
# ============================================================

library(TCGAbiolinks)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(openxlsx)
library(ggpubr)


# ============================================================
# Load MOC cohort
# ============================================================

data <- combat_Omicodf5625 %>%
  filter(
    Label == "Mucinous_ovarian",
    Include_tumour_content == "Y"
  ) %>%
  left_join(
    MOC_clinicaldata %>%
      dplyr::select(
        Sample_ID,
        ERBB2.any
      ),
    by = "Sample_ID"
  )

# ============================================================
# -----------------------
# FIG 5A: ERBB2 protein
# -----------------------
# ============================================================

data <- data %>%
  mutate(
    ERBB2_status = case_when(
      ERBB2.any == "Amplified" ~ "Amplified",
      ERBB2.any == "Not amplified" ~ "Not amplified",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(ERBB2_status))

data$ERBB2_status <- factor(
  data$ERBB2_status,
  levels = c("Not amplified", "Amplified")
)

p_5A <- ggplot(data, aes(ERBB2_status, P04626)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(method = "t.test", size = 5) +
  labs(
    x = "ERBB2 gene status",
    y = "ERBB2 protein expression"
  ) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

ggsave("figures/Fig5A_ERBB2_protein.pdf", p_5A, width = 5, height = 5)

# ============================================================
# -----------------------
# FIG 5B: Volcano (DEP)
# -----------------------
# ============================================================

# NOTE: MOC_data is DEP result ERBB2 amp vs. non-amp
# DEP_Fig5_ERBB2.csv contains limma differential expression results
# comparing ERBB2-amplified versus ERBB2-non-amplified MOC samples.
# Results generated in the differential expression analysis script
MOC_data <- read.csv("../results/DEP_Fig5_ERBB2.csv")

df <- MOC_data %>%
  mutate(
    status = case_when(
      logFC >= 1 & adj.P.Val <= 0.05 ~ "Up",
      logFC <= -1 & adj.P.Val <= 0.05 ~ "Down",
      TRUE ~ "NS"
    )
  )

top_genes <- df %>%
  filter(status != "NS") %>%
  slice_max(abs(logFC), n = 10)

p_5B <- ggplot(df, aes(logFC, -log10(adj.P.Val), color = status)) +
  geom_point() +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = top_genes, aes(label = symbol), size = 4) +
  labs(
    x = "log2 Fold Change",
    y = "-log10 adj P-value"
  ) +
  theme_classic(base_size = 16)

ggsave("figures/Fig5B_volcano_ERBB2.pdf", p_5B, width = 6, height = 5)

# ============================================================
# -----------------------
# FIG 5C: Lipid GO enrichment
# -----------------------
# ============================================================

dep_genes <- MOC_data %>%
  filter(adj.P.Val <= 0.05 & logFC > 0) %>%
  pull(symbol) %>%
  unique()

go_res <- enrichGO(
  gene = dep_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pvalueCutoff = 0.05
)

go_lipid <- as.data.frame(go_res) %>%
  filter(grepl("lipid|fatty|cholesterol|phospholipid",
               Description,
               ignore.case = TRUE)) %>%
  mutate(logP = -log10(p.adjust)) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10)

p_5C <- ggplot(go_lipid,
               aes(Count,
                   reorder(Description, Count),
                   color = logP,
                   size = Count)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    x = "Gene count",
    y = "GO term"
  ) +
  theme_minimal(base_size = 14)

ggsave("figures/Fig5C_GO_lipid.pdf", p_5C, width = 7, height = 5)

# ============================================================
# -----------------------
# SUPPLEMENTARY FIGURE 5
# TCGA validation
# -----------------------
# ============================================================

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)
tcga <- GDCprepare(query)

expr <- assay(tcga)

erbb2 <- expr["ENSG00000141736.14", ]

# Define ERBB2-high tumors as the upper quartile of ERBB2 transcript abundance.
group <- ifelse(
  erbb2 > quantile(erbb2, 0.75),
  "Amplified",
  "Not amplified"
)

group <- factor(group)

design <- model.matrix(~ group)
fit <- eBayes(lmFit(expr, design))

deg <- topTable(fit, coef = 2, number = Inf)

genes <- rownames(deg)[deg$logFC > 0 & deg$adj.P.Val < 0.05]
genes <- sub("\\..*", "", genes)

go_tcga <- enrichGO(
  gene = genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP"
)

go_tcga_lipid <- as.data.frame(go_tcga) %>%
  filter(grepl("lipid|fatty|cholesterol", Description, ignore.case = TRUE)) %>%
  mutate(logP = -log10(p.adjust)) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10)

p_sup <- ggplot(go_tcga_lipid,
                aes(Count,
                    reorder(Description, Count),
                    color = logP,
                    size = Count)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal(base_size = 14)

ggsave("figures/SuppFig5_TCGA_BRCA_lipid.pdf", p_sup, width = 8, height = 5)


# ============================================================
# -----------------------
# SUPPLEMENTARY FIGURE 5
# TCGA gastric validation
# -----------------------
# ============================================================

query <- GDCquery(
  project = "TCGA-STAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)
tcga <- GDCprepare(query)

expr <- assay(tcga)

# ERBB2 Ensembl ID consistent with your BRCA script
erbb2 <- expr["ENSG00000141736.14", ]

group <- ifelse(
  erbb2 > quantile(erbb2, 0.75),
  "Amplified",
  "Not amplified"
)

group <- factor(group)

design <- model.matrix(~ group)

fit <- eBayes(lmFit(expr, design))
deg <- topTable(fit, coef = 2, number = Inf)

genes <- rownames(deg)[deg$logFC > 0 & deg$adj.P.Val < 0.05]
genes <- sub("\\..*", "", genes)

go_tcga <- enrichGO(
  gene = genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP"
)

go_tcga_lipid <- as.data.frame(go_tcga) %>%
  filter(grepl("lipid|fatty|cholesterol", Description, ignore.case = TRUE)) %>%
  mutate(logP = -log10(p.adjust)) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10)

p_sup <- ggplot(
  go_tcga_lipid,
  aes(
    Count,
    reorder(Description, Count),
    color = logP,
    size = Count
  )
) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal(base_size = 14)

ggsave("figures/SuppFig5_TCGA_STAD_lipid.pdf", p_sup, width = 8, height = 5)

# =========================
# Fig 5D
# Proliferation markers (Cluster 2 vs others)
# =========================

data <- combat_Omicodf5625 %>%
  filter(
    Label == "Mucinous_ovarian",
    Include_tumour_content == "Y"
  ) %>%
  left_join(
    MOC_clinicaldata %>%
      dplyr::select(
        Sample_ID,
        Unsupervised_cluster
      ),
    by = "Sample_ID"
  )


data <- data %>%
  mutate(
    clusterspearman = Unsupervised_cluster
  )

data$clusterspearman <- factor(
  data$clusterspearman,
  levels = c(1,2,3,4,5)
)

# =========================
# 2. Define proteins
# =========================

protein_map <- c(
  "P11388" = "TOP2A",
  "P46013" = "MKI67",
  "P24821" = "PCNA"
)

# =========================
# 3. Reshape
# =========================

df_long <- data %>%
  dplyr::select(clusterspearman, all_of(names(protein_map))) %>%
  pivot_longer(
    cols = all_of(names(protein_map)),
    names_to = "UniProt",
    values_to = "intensity"
  ) %>%
  mutate(
    Gene = recode(UniProt, !!!protein_map),
    scaled_intensity = scale(intensity)[,1]
  )

# =========================
# 4. Define comparison group
# =========================

data <- data %>%
  mutate(group = ifelse(clusterspearman == 2, "Cluster 2", "Other"))

# =========================
# 5. Statistical testing
# =========================

pvals <- sapply(names(protein_map), function(pid) {
  wilcox.test(data[[pid]] ~ data$group)$p.value
})

pvals_adj <- p.adjust(pvals, method = "BH")

pval_df <- data.frame(
  UniProt = names(protein_map),
  p_adj = pvals_adj
)

pval_labels <- data.frame(
  Gene = unname(protein_map),
  p_adj = pvals_adj
)

pval_labels <- pval_labels %>%
  left_join(
    df_long %>%
      group_by(Gene) %>%
      summarise(y = max(scaled_intensity) * 1.15),
    by = "Gene"
  ) %>%
  mutate(
    label = paste0("BH p = ", signif(p_adj, 3)),
    x = "3"
  )  

ggplot(df_long, aes(x = as.factor(clusterspearman), y = scaled_intensity, fill = as.factor(clusterspearman))) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ Gene, scales = "free_y") +
  geom_text(data = pval_labels, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "",
       x = "Cluster",
       y = "Relative Protein Expression",
       fill = "Cluster") +
  theme_minimal(base_size = 18) +
  theme(
    strip.text = element_text(face = "bold", size = 18),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/Fig5D_proliferation_markers.pdf",
  width = 8,
  height = 6
)

# =========================
# Fig 5E KRAS related pathways

# ============================================================
# Fig 5E: KRAS pathway enrichment (cluster-level heatmap)
# ============================================================

set.seed(1234)

# ============================================================
# 1. Load data + clusters
# ============================================================

data <- combat_Omicodf5625 %>%
  filter(
    Label == "Mucinous_ovarian",
    Include_tumour_content == "Y"
  ) %>%
  left_join(
    MOC_clinicaldata %>%
      dplyr::select(
        Sample_ID,
        Unsupervised_cluster
      ),
    by = "Sample_ID"
  )

data$clusterspearman <- factor(
  data$Unsupervised_cluster,
  levels = c(1,2,3,4,5)
)

# ============================================================
# 2. Build expression matrix
# ============================================================

expr <- data %>%
  select(where(is.numeric)) %>%
  as.matrix()

rownames(expr) <- data$Sample_ID
expr <- t(expr)

# ============================================================
# 3. Map UniProt → Entrez
# ============================================================

entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = rownames(expr),
  keytype = "UNIPROT",
  column = "ENTREZID",
  multiVals = "first"
)

expr <- expr[!is.na(entrez_ids), ]
rownames(expr) <- entrez_ids[!is.na(entrez_ids)]
expr <- expr[!duplicated(rownames(expr)), ]

# ============================================================
# 4. MSigDB C6 (oncogenic signatures)
# ============================================================

c6_sets <- msigdbr(species = "Homo sapiens", category = "C6")
gene_sets <- split(c6_sets$entrez_gene, c6_sets$gs_name)

# ============================================================
# 5. Run ssGSEA
# ============================================================

params <- ssgseaParam(
  exprData = expr,
  geneSets = gene_sets,
  normalize = FALSE,
  alpha = 0.25
)

ssgsea_res <- gsva(params)

# ============================================================
# 6. Subset KRAS pathways
# ============================================================

kras_paths <- rownames(ssgsea_res)[grepl("KRAS", rownames(ssgsea_res))]
ssgsea_kras <- ssgsea_res[kras_paths, , drop = FALSE]

# ensure numeric
ssgsea_kras <- as.matrix(ssgsea_kras)
storage.mode(ssgsea_kras) <- "numeric"

# ============================================================
# 7. Match annotation
# ============================================================

annotation <- data %>%
  select(Sample_ID, clusterspearman) %>%
  distinct()

common_ids <- intersect(colnames(ssgsea_kras), annotation$Sample_ID)

ssgsea_kras <- ssgsea_kras[, common_ids, drop = FALSE]

annotation <- annotation[
  match(common_ids, annotation$Sample_ID),
]

rownames(annotation) <- annotation$Sample_ID

# ============================================================
# 8. AVERAGE BY CLUSTER (final biological version)
# ============================================================

clusters_vec <- annotation$clusterspearman
cluster_levels <- levels(clusters_vec)

cluster_means <- sapply(cluster_levels, function(cl) {
  idx <- which(clusters_vec == cl)
  
  if (length(idx) > 1) {
    rowMeans(ssgsea_kras[, idx, drop = FALSE])
  } else {
    ssgsea_kras[, idx]
  }
})

colnames(cluster_means) <- cluster_levels

# ============================================================
# 9. Scale to [-1, 1]
# ============================================================

cluster_scaled <- t(scale(t(cluster_means)))
cluster_scaled[!is.finite(cluster_scaled)] <- 0

cluster_scaled[cluster_scaled > 1] <- 1
cluster_scaled[cluster_scaled < -1] <- -1

# clean labels
rownames(cluster_scaled) <- gsub("_", " ", rownames(cluster_scaled))
rownames(cluster_scaled) <- gsub("\\.", " ", rownames(cluster_scaled))

# ============================================================
# 10. Plot (ComplexHeatmap)
# ============================================================

cluster_colors <- c(
  "1" = "#66C2A5",
  "2" = "#FC8D62",
  "3" = "#8DA0CB",
  "4" = "#E78AC3",
  "5" = "#A6D854"
)

ha <- HeatmapAnnotation(
  Cluster = factor(colnames(cluster_scaled), levels = cluster_levels),
  col = list(Cluster = cluster_colors)
)

ht <- Heatmap(
  cluster_scaled,
  name = "GSEA\nscore",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  top_annotation = ha,
  width = unit(8, "cm"), # adjust down
  show_row_names = TRUE,
  row_names_max_width = unit(6, "cm"),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = FALSE,
)

pdf("figures/Fig5E_KRAS_clusters.pdf", width = 16, height = 10)

draw(
  ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legends = TRUE,
  padding = unit(c(5, 5, 5, 30), "mm")
)
dev.off()

# ============================================================
# 10_Fig6.R
# ============================================================

library(tidyverse)
library(pheatmap)
library(tibble)

# ============================================================
# 1. Load MOC data and cluster assignments
# ============================================================

data <- combat_Omicodf5625 %>%
  filter(
    Label == "Mucinous_ovarian",
    Include_tumour_content == "Y"
  ) %>%
  left_join(
    MOC_clinicaldata %>%
      select(
        Sample_ID,
        Unsupervised_cluster
      ),
    by = "Sample_ID"
  )

data$Cluster <- factor(
  data$Unsupervised_cluster,
  levels = c(1, 2, 3, 4, 5)
)

# ============================================================
# 2. Define markers (UniProt → Gene)
# ============================================================

uniprot_to_gene <- c(
  "P09758" = "TROP2",
  "Q96NY8" = "NECTIN4",
  "P04626" = "ERBB2",
  "Q5ZPR3" = "CD276",
  "P56856" = "CLDN18",
  "Q7Z7D3" = "B7H4",
  "P06731" = "CEACAM5"
)

colnames(data)[colnames(data) %in% names(uniprot_to_gene)] <-
  uniprot_to_gene[
    colnames(data)[colnames(data) %in% names(uniprot_to_gene)]
  ]

markers <- unname(uniprot_to_gene)

# ============================================================
# 3. Build matrix
# ============================================================

heatmap_mat <- data %>%
  select(Sample_ID, all_of(markers)) %>%
  column_to_rownames("Sample_ID") %>%
  as.matrix()

# Z-score per marker
heatmap_scaled <- scale(
  heatmap_mat,
  center = TRUE,
  scale = TRUE
)

heatmap_scaled <- t(heatmap_scaled)

# ============================================================
# 4. Colours
# ============================================================

cluster_colors <- c(
  "1" = "#66C2A5",
  "2" = "#FC8D62",
  "3" = "#8DA0CB",
  "4" = "#E78AC3",
  "5" = "#A6D854"
)

# ============================================================
# 5. Cluster annotation
# ============================================================

annotation_cluster <- data %>%
  select(
    Sample_ID,
    Cluster
  ) %>%
  column_to_rownames("Sample_ID")

ord_cluster <- order(annotation_cluster$Cluster)

mat_cluster <- heatmap_scaled[, ord_cluster]

ann_cluster <- annotation_cluster[
  ord_cluster,
  ,
  drop = FALSE
]

ann_colors_cluster <- list(
  Cluster = cluster_colors
)

# ============================================================
# 6. Plot heatmap
# ============================================================

pdf(
  "figures/Fig6A_drug_targets_clusters.pdf",
  width = 10,
  height = 8
)

pheatmap(
  mat_cluster,
  annotation_col = ann_cluster,
  annotation_colors = ann_colors_cluster,
  color = colorRampPalette(
    c("blue", "white", "red")
  )(100),
  breaks = seq(-3, 3, length.out = 101),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize = 16
)

dev.off()

# ============================================================
# ============================================================
# FIG 6B: Drug targets by diagnosis
# ============================================================

data <- combat_Omicodf5625 %>%
  filter(
    Include_tumour_content == "Y"
  )


# UniProt -> gene names
uniprot_to_gene <- c(
  "P09758" = "TROP2",
  "Q96NY8" = "NECTIN4",
  "P04626" = "ERBB2",
  "Q5ZPR3" = "CD276",
  "P56856" = "CLDN18",
  "Q7Z7D3" = "B7H4",
  "P06731" = "CEACAM5"
)

colnames(data)[colnames(data) %in% names(uniprot_to_gene)] <-
  uniprot_to_gene[
    colnames(data)[colnames(data) %in% names(uniprot_to_gene)]
  ]

markers <- unname(uniprot_to_gene)

# Build matrix
heatmap_mat <- data %>%
  select(Sample_ID, all_of(markers)) %>%
  column_to_rownames("Sample_ID") %>%
  as.matrix()

# Z-score per marker
heatmap_scaled <- scale(
  heatmap_mat,
  center = TRUE,
  scale = TRUE
)

heatmap_scaled <- t(heatmap_scaled)

# Diagnosis annotation
annotation_diag <- data %>%
  mutate(
    Diagnosis = gsub("_", " ", Label)
  ) %>%
  select(
    Sample_ID,
    Diagnosis
  ) %>%
  column_to_rownames("Sample_ID")

annotation_diag$Diagnosis <- factor(
  annotation_diag$Diagnosis
)
annotation_diag <- annotation_diag %>%
  mutate(
    Diagnosis = recode(
      Diagnosis,
      "UpperGI adenocarcinoma" = "Gastro-esophageal adenocarcinoma"
    )
  )

# Order by diagnosis
ord_diag <- order(annotation_diag$Diagnosis)

mat_diag <- heatmap_scaled[, ord_diag]

ann_diag <- annotation_diag[
  ord_diag,
  ,
  drop = FALSE
]

ann_colors_diag <- list(
  Diagnosis = diagnosis_colors
)

pdf(
  "figures/Fig6B_drug_targets_diagnosis.pdf",
  width = 10,
  height = 8
)

pheatmap(
  mat_diag,
  annotation_col = ann_diag,
  annotation_colors = ann_colors_diag,
  color = colorRampPalette(
    c("blue", "white", "red")
  )(100),
  breaks = seq(-3, 3, length.out = 101),
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize = 16
)

dev.off()

# ============================================================
# Fig 6C
# ESTIMATE immune score across cancer types
# ============================================================

library(tidyverse)
library(ggplot2)


# Generate scores 

# ============================================================
# ESTIMATE immune and stromal scores
# Calculated from the ComBat-corrected proteomics matrix
# ============================================================

library(tidyverse)
library(tidyestimate)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Create protein expression matrix
expr_matrix <- combat_Omicodf5625[, 1:10416]

expr_matrix <- t(as.matrix(expr_matrix))
colnames(expr_matrix) <- combat_Omicodf5625$Sample_ID

# Convert UniProt IDs to HGNC symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(expr_matrix),
  keytype = "UNIPROT",
  column = "SYMBOL",
  multiVals = "first"
)

rownames(expr_matrix) <- gene_symbols

# Remove unmapped and duplicated genes
expr_matrix <- expr_matrix[!is.na(rownames(expr_matrix)), ]
expr_matrix <- expr_matrix[!duplicated(rownames(expr_matrix)), ]

# Filter to genes recognised by ESTIMATE
filtered <- filter_common_genes(
  expr_matrix,
  id = "hgnc_symbol",
  find_alias = TRUE
)

# Calculate ESTIMATE scores
estimate_results <- estimate_score(
  filtered,
  is_affymetrix = FALSE
)

estimate_results$Sample_ID <- estimate_results$sample

# Add metadata
estimate_results <- merge(
  estimate_results,
  combat_Omicodf5625 %>%
    select(
      Sample_ID,
      Label,
      sample_cohort
    ),
  by = "Sample_ID",
  all.x = TRUE
)


estimate_results <- estimate_results %>%
  left_join(
    combat_Omicodf5625 %>%
      select(
        Sample_ID,
        Label,
        sample_cohort
      ),
    by = "Sample_ID"
  ) %>%
  left_join(
    Sample_metadata_clean %>%
      select(
        Sample_ID,
        Include_in_study
      ),
    by = "Sample_ID"
  )

# Save results
write.csv(
  estimate_results,
  "../results/ESTIMATE_scores.csv",
  row.names = FALSE
)

# ============================================================
# 1. Load data
# ============================================================

data <- read.csv("../data/ESTIMATE_scores.csv")
data<- estimate_results

# ============================================================
# 2. Clean labels
# ============================================================

data <- data %>% filter(Include_in_study == "Y")

data <- data %>%
  mutate(
    Diagnosis = gsub("_", " ", Label.x)
  )

# ============================================================
# 3. Define comparison (MOC vs others)
# ============================================================

data <- data %>%
  mutate(
    Group = ifelse(Diagnosis == "Mucinous ovarian",
                   "MOC",
                   "Other")
  )

# ============================================================
# 4. Statistical test
# ============================================================

pval <- wilcox.test(immune ~ Group, data = data)$p.value

p_label <- paste0("Wilcoxon p = ", signif(pval, 3))

# ============================================================
# 5. Plot
# ============================================================

p <- ggplot(data, aes(x = Diagnosis, y = immune, fill = Diagnosis)) +
  
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  
  geom_jitter(
    width = 0.2,
    size = 2,
    alpha = 0.6,
    color = "black"
  ) +
  
  annotate(
    "text",
    x = 2,
    y = max(data$immune, na.rm = TRUE) * 1.05,
    label = p_label,
    size = 5
  ) +
  
  scale_x_discrete(labels = c(
    "Appendix adenocarcinoma" = "AA",
    "Colorectal adenocarcinoma" = "CRC",
    "High grade serous ovary" = "HGSOC",
    "Mucinous ovarian" = "MOC",
    "Pancreas adenocarcinoma" = "PDAC",
    "Small bowel adenocarcinoma" = "SB",
    "UpperGI adenocarcinoma" = "GEJ"
  )) +
  
  scale_y_continuous(
    limits = c(-1000, 1900),
    breaks = seq(-1000, 1900, by = 500)
  ) +
  
  labs(
    x = "",
    y = "ESTIMATE immune score"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )

# ============================================================
# 6. Save
# ============================================================

ggsave(
  "figures/Fig6C_immune_scores.pdf",
  p,
  width = 7,
  height = 5
)

# ============================================================
# Fig 6E: GSEA enrichment plots
# Interferon + Complement
# ============================================================

library(tidyverse)
library(fgsea)
library(msigdbr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(patchwork)

# ============================================================
# 1. Load DEP results
# ============================================================

# ============================================================
# Immune-high vs immune-low DEP
# Top 25% vs bottom 75% ESTIMATE immune score
# ============================================================

# ============================================================
# Immune-high vs immune-low DEP
# Top 25% vs bottom 75% ESTIMATE immune score
# ============================================================

df <- combat_Omicodf5625 %>%
  filter(
    Label == "Mucinous_ovarian",
    Include_tumour_content == "Y"
  ) %>%
  left_join(
    estimate_results %>%
      select(
        Sample_ID,
        immune
      ),
    by = "Sample_ID"
  )

df <- df %>%
  mutate(
    Immune_group = ifelse(
      rank(-immune) <= 0.25 * n(),
      "High",
      "Low"
    )
  )

table(df$Immune_group)


summary(df$immune)

df %>%
  arrange(desc(immune)) %>%
  select(Sample_ID, immune, Immune_group) %>%
  head(20)

Group <- ifelse(
  df$Immune_group == "High",
  1,
  0
)

Biopsy_dummies <- data.frame(df, Group = Group)

X_matrix <- Biopsy_dummies[, 1:10416]
X_matrix <- data.matrix(X_matrix)
X_matrix <- t(X_matrix)

p_risk <- Biopsy_dummies[, 10417:ncol(Biopsy_dummies)]
pdf_risk <- as.data.frame(p_risk)

colnames(X_matrix) <- rownames(pdf_risk)

eset_risk <- ExpressionSet(
  assayData = X_matrix,
  phenoData = AnnotatedDataFrame(pdf_risk)
)

# Add grouping variable to metadata
pdf_risk$Group <- Group

# Design matrix
design <- model.matrix(
  ~ Group,
  data = pdf_risk
)

# Differential expression
fit_risk <- lmFit(
  eset_risk,
  design
)

fit_risk <- eBayes(fit_risk)

# Extract results
MOC_DEP_results <- topTable(
  fit_risk,
  coef = 2,
  number = Inf
)

# Add UniProt IDs
MOC_DEP_results$RowNames <- rownames(MOC_DEP_results)

# Map UniProt IDs to gene symbols
MOC_DEP_results$symbol <- mapIds(
  org.Hs.eg.db,
  keys = MOC_DEP_results$RowNames,
  keytype = "UNIPROT",
  column = "SYMBOL",
  multiVals = "first"
)

# Classify proteins
MOC_DEP_results <- MOC_DEP_results %>%
  mutate(
    Expression = case_when(
      logFC >= 1 & adj.P.Val <= 0.05 ~ "Up-regulated",
      logFC <= -1 & adj.P.Val <= 0.05 ~ "Down-regulated",
      TRUE ~ "Unchanged"
    )
  )

# Save results for FGSEA
write.csv(
  MOC_DEP_results,
  "../results/MOC_DEP_results.csv",
  row.names = FALSE
)

# QC checks
table(MOC_DEP_results$Expression)

MOC_DEP_results %>%
  arrange(adj.P.Val) %>%
  select(symbol, logFC, adj.P.Val) %>%
  head(20)

# ============================================================
# Fig 6D: Immune-high vs immune-low volcano plot
# ============================================================

library(tidyverse)
library(ggrepel)

MOC_data <- read.csv("../results/MOC_DEP_results.csv")

MOC_data <- MOC_data %>%
  mutate(
    status = case_when(
      logFC >= 1 & adj.P.Val <= 0.05 ~ "Up-regulated",
      logFC <= -1 & adj.P.Val <= 0.05 ~ "Down-regulated",
      TRUE ~ "Unchanged"
    )
  ) %>%
  filter(
    !is.na(logFC),
    !is.na(adj.P.Val)
  )

# DEG counts
table(MOC_data$status)

top_up <- MOC_data %>%
  filter(status == "Up-regulated") %>%
  slice_max(logFC, n = 5)

top_down <- MOC_data %>%
  filter(status == "Down-regulated") %>%
  slice_min(logFC, n = 5)

top_genes <- bind_rows(top_up, top_down)

p_volcano <- ggplot(
  MOC_data,
  aes(
    x = logFC,
    y = -log10(adj.P.Val),
    color = status
  )
) +
  geom_point() +
  scale_color_manual(
    values = c(
      "Up-regulated" = "red",
      "Down-regulated" = "blue",
      "Unchanged" = "grey"
    )
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  geom_text_repel(
    data = top_genes,
    aes(label = symbol),
    size = 4
  ) +
  labs(
    x = "log2 Fold Change",
    y = "-log10 adjusted P-value"
  ) +
  theme_classic(base_size = 16)

ggsave(
  "figures/Fig6D_Immune_high_vs_low_volcano.pdf",
  p_volcano,
  width = 10,
  height = 8
)    

# ============================================================
# 2. Map UniProt → Entrez
# ============================================================

entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = MOC_data$RowNames,
  keytype = "UNIPROT",
  column = "ENTREZID",
  multiVals = "first"
)

MOC_data$entrez <- entrez_ids

# ============================================================
# 3. Ranked list
# ============================================================

ranked_genes <- setNames(MOC_data$logFC, MOC_data$entrez)
ranked_genes <- ranked_genes[!is.na(names(ranked_genes))]
ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# ============================================================
# 4. Gene sets
# ============================================================

hallmark <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

reactome <- msigdbr(
  species = "Homo sapiens",
  category = "C2"
) %>%
  filter(gs_subcollection == "CP:REACTOME")

kegg <- msigdbr(
  species = "Homo sapiens",
  category = "C2"
) %>%
  filter(gs_subcollection == "CP:KEGG_LEGACY")

pathways <- c(
  split(hallmark$entrez_gene, hallmark$gs_name),
  split(reactome$entrez_gene, reactome$gs_name),
  split(kegg$entrez_gene, kegg$gs_name)
)

# ============================================================
# 5. Run FGSEA
# ============================================================

fgsea_res <- fgsea(
  pathways = pathways,
  stats = ranked_genes,
  nperm = 1000
) %>%
  as.data.frame()

# ============================================================
# 6. Pathways to plot
# ============================================================

ifn_paths <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "REACTOME_INTERFERON_GAMMA_SIGNALING"
)

comp_paths <- c(
  "HALLMARK_COMPLEMENT",
  "REACTOME_COMPLEMENT_CASCADE",
  "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES"
)

# ============================================================
# 7. Build enrichment curves
# ============================================================

ifn_df <- bind_rows(
  lapply(ifn_paths, function(path) {
    d <- plotEnrichment(
      pathways[[path]],
      ranked_genes
    )$data
    d$Pathway <- path
    d
  })
)

comp_df <- bind_rows(
  lapply(comp_paths, function(path) {
    d <- plotEnrichment(
      pathways[[path]],
      ranked_genes
    )$data
    d$Pathway <- path
    d
  })
)

# ============================================================
# 8. Plot
# ============================================================

p_ifn <- ggplot(
  ifn_df,
  aes(x = rank, y = ES, colour = Pathway)
) +
  geom_line(linewidth = 1.2) +
  theme_classic(base_size = 16) +
  labs(
    title = "Interferon Signalling",
    x = "Ranked Genes",
    y = "Enrichment Score"
  )

p_comp <- ggplot(
  comp_df,
  aes(x = rank, y = ES, colour = Pathway)
) +
  geom_line(linewidth = 1.2) +
  theme_classic(base_size = 16) +
  labs(
    title = "Complement Pathways",
    x = "Ranked Genes",
    y = "Enrichment Score"
  )

# ============================================================
# 9. Combine
# ============================================================

final_plot <- p_ifn + p_comp

# ============================================================
# 10. Save
# ============================================================

ggsave(
  "figures/Fig6E_GSEA_IFN_Complement.pdf",
  final_plot,
  width = 24,
  height = 5
)


# ============================================================
# Fig 6F: Antigen-processing pathways in immune-high MOC
# ============================================================

# Perform ssGSEA for MOC samples 

library(tidyverse)
library(GSVA)
library(msigdbr)

# ============================================================
# MOC expression matrix
# ============================================================

moc_data <- combat_Omicodf5625 %>%
  filter(
    Label == "Mucinous_ovarian",
    Include_tumour_content == "Y"
  )

expr <- moc_data[, 1:10416]

expr <- t(as.matrix(expr))
colnames(expr) <- moc_data$Sample_ID

# ============================================================
# Pathway databases
# ============================================================

hallmark <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

reactome <- msigdbr(
  species = "Homo sapiens",
  category = "C2"
) %>%
  filter(gs_subcollection == "CP:REACTOME")

kegg <- msigdbr(
  species = "Homo sapiens",
  category = "C2"
) %>%
  filter(gs_subcollection == "CP:KEGG_LEGACY")

Hallmark_list <- split(
  hallmark$gene_symbol,
  hallmark$gs_name
)

Reactome_list <- split(
  reactome$gene_symbol,
  reactome$gs_name
)

KEGG_list <- split(
  kegg$gene_symbol,
  kegg$gs_name
)

gene_sets <- c(
  Hallmark_list,
  Reactome_list,
  KEGG_list
)

# ============================================================
# UniProt -> Gene Symbol
# ============================================================

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(expr),
  keytype = "UNIPROT",
  column = "SYMBOL",
  multiVals = "first"
)

rownames(expr) <- gene_symbols

expr <- expr[!is.na(rownames(expr)), ]
expr <- expr[!duplicated(rownames(expr)), ]

# ============================================================
# ssGSEA
# ============================================================

params <- ssgseaParam(
  exprData = expr,
  geneSets = gene_sets,
  normalize = FALSE,
  alpha = 0.25
)

ssgsea_results <- gsva(params)


library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Define Top 3 immune samples
df <- df %>%
  mutate(
    immune_group = ifelse(
      rank(-immune) <= 3,
      "High",
      "Low"
    )
  )

df <- df %>%
  left_join(
    MOC_clinicaldata %>%
      select(
        Sample_ID,
        Unsupervised_cluster
      ),
    by = "Sample_ID"
  )

# Extract antigen-processing pathways
reactome_antigen_paths <- grep(
  "ANTIGEN_PROCESSING",
  unique(reactome$gs_name),
  value = TRUE,
  ignore.case = TRUE
)

kegg_antigen_paths <- grep(
  "ANTIGEN",
  unique(kegg$gs_name),
  value = TRUE,
  ignore.case = TRUE
)
antigen_processing_paths <- c(
  reactome_antigen_paths,
  kegg_antigen_paths
)

# ssGSEA scores
antigen_scores <- ssgsea_results[
  antigen_processing_paths,
  ,
  drop = FALSE
]

# Sample annotations
annotation_col <- data.frame(
  `Immune Group` = df$immune_group,
  Cluster = factor(
    df$Unsupervised_cluster.x,
    levels = c(1, 2, 3, 4, 5)
  )
)

rownames(annotation_col) <- colnames(antigen_scores)

# Order by immune group
ord <- order(annotation_col$`Immune.Group`)

antigen_scores <- antigen_scores[, ord]
annotation_col <- annotation_col[ord, ]

# Z-score pathways
scaled_matrix <- t(
  scale(
    t(antigen_scores)
  )
)

# Clean pathway names
rownames(scaled_matrix) <- gsub(
  "_",
  " ",
  rownames(scaled_matrix)
)

# Colours
cluster_colors <- c(
  "1" = "#1b9e77",
  "2" = "#d95f02",
  "3" = "#7570b3",
  "4" = "#e7298a",
  "5" = "#66a61e"
)

ann_colors <- list(
  Cluster = cluster_colors,
  `Immune Group` = c(
    "High" = "#e41a1c",
    "Low" = "#377eb8"
  )
)

top_annot <- HeatmapAnnotation(
  df = annotation_col,
  col = ann_colors
)

ht <- Heatmap(
  scaled_matrix,
  name = "ssGSEA",
  col = colorRamp2(
    c(-2, 0, 2),
    c("navy", "white", "firebrick")
  ),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_names = FALSE,
  top_annotation = top_annot
)

pdf(
  "figures/Fig6F_Antigen_processing_heatmap.pdf",
  width = 12,
  height = 8
)

draw(
  ht,
  heatmap_legend_side = "top",
  annotation_legend_side = "top"
)

dev.off()

# ============================================================
# Fig 6G: Antigen-presentation proteins
# Top 3 immune versus remaining MOCs
# ============================================================

library(tidyverse)
library(ggpubr)

proteins_of_interest <- c(
  "P01732",  # CD8A
  "P01909",  # HLA-DQA1
  "P04233",  # CD74
  "Q03518"   # TAP1
)

protein_labels <- c(
  "P01732" = "CD8A",
  "P01909" = "HLA-DQA1",
  "P04233" = "CD74",
  "Q03518" = "TAP1"
)

plot_data <- df %>%
  mutate(
    immune_group = ifelse(
      rank(-immune) <= 3,
      "Top 3",
      "Low"
    )
  ) %>%
  select(
    Sample_ID,
    immune_group,
    all_of(proteins_of_interest)
  ) %>%
  pivot_longer(
    cols = all_of(proteins_of_interest),
    names_to = "Protein",
    values_to = "Expression"
  )

plot_data$Protein <- factor(
  plot_data$Protein,
  levels = names(protein_labels),
  labels = protein_labels
)

plot_data$immune_group <- factor(
  plot_data$immune_group,
  levels = c("Low", "Top 3")
)

p <- ggplot(
  plot_data,
  aes(
    x = Protein,
    y = Expression,
    fill = immune_group
  )
) +
  geom_boxplot(
    outlier.colour = "black"
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.6,
    size = 0.8
  ) +
  stat_compare_means(
    method = "kruskal.test",
    label = "p.format",
    label.y = 20
  ) +
  labs(
    x = "",
    y = "Relative Protein Expression",
    fill = "Immune Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  "figures/Fig6G_Immune_markers.pdf",
  p,
  width = 8,
  height = 5
)

# ============================================================
# 11_survival.R
# ============================================================
# Supplementary Fig 4g and Supplementary Fig 6a

library(tidyverse)
library(survival)
library(survminer)

data <- MOC_clinicaldata %>%
  select(
    Sample_ID,
    Unsupervised_cluster,
    OS.days,
    OS.Status
  )

data$clusterspearman <- factor(
  data$Unsupervised_cluster,
  levels = c(1, 2, 3, 4, 5)
)

data <- data %>%
  filter(
    !is.na(OS.days),
    !is.na(OS.Status)
  )

data <- data %>%
  mutate(
    OS.days = as.numeric(OS.days),
    OS.Status = as.numeric(OS.Status)
  ) %>%
  filter(
    !is.na(OS.days),
    !is.na(OS.Status)
  )

fit <- survfit(
  Surv(OS.days, OS.Status) ~ clusterspearman,
  data = data
)

surv_plot <- ggsurvplot(
  fit,
  data = data,
  pval = TRUE,
  risk.table = TRUE
)

pdf(
  "../figures/Fig6_survival_clusters.pdf",
  width = 8,
  height = 6
)

print(surv_plot)

dev.off()


library(tidyverse)
library(survival)
library(survminer)

data <- MOC_clinicaldata %>%
  left_join(
    estimate_results %>%
      select(
        Sample_ID,
        immune
      ),
    by = "Sample_ID"
  ) %>%
  mutate(
    immune_group = ifelse(
      rank(-immune) <= 0.25 * n(),
      "High",
      "Low"
    ),
    OS.days = as.numeric(OS.days),
    OS.Status = as.numeric(OS.Status)
  ) %>%
  filter(
    !is.na(OS.days),
    !is.na(OS.Status)
  )

fit <- survfit(
  Surv(OS.days, OS.Status) ~ immune_group,
  data = data
)

surv_plot <- ggsurvplot(
  fit,
  data = data,
  pval = TRUE,
  risk.table = TRUE,
  legend.title = "Immune Group",
  legend.labs = c("Low", "High")
)

pdf(
  "figures/Fig6_OS_immune_group.pdf",
  width = 8,
  height = 6
)

print(surv_plot)

dev.off()

# ================================
# Fig 6f_ADC SRM Cluster Plot (MOC only)
# ================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)


# ================================
# Fig 6f_ADC SRM Cluster Plot (MOC only)
# ================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# ============================================================
# SRM validation
#
# This section:
# 1. Calculates SRM log2 fold-changes between MOC and CRC.
# 2. Calculates corresponding DIA log2 fold-changes.
# 3. Assesses cross-platform concordance.
# 4. Evaluates expression of candidate therapeutic targets.
# 5. Compares target expression across MOC subtypes.
# 6. Assesses association between target expression and tumour content.
#
# Inputs:
# - SRM_MOC.csv
# - combat_Omicodf5626.csv
#
# Outputs:
# - SRM_DIA_concordance.csv
# - Cross-platform concordance plots
# - Target protein expression plots
# - Cluster-specific expression plots
# ============================================================

# ------------------------------------------------------------
# Load SRM validation dataset (Supplementary Table 5)
# ------------------------------------------------------------

SRM_MOC <- read.csv("../data/SRM_MOC.csv")

# ------------------------------------------------------------
# Target protein expression across MOC clusters
# ------------------------------------------------------------

MOC <- SRM_MOC %>%
  filter(
    Diagnosis == "Mucinous_ovarian",
    !is.na(Cluster),
    Cluster != "NA"
  )

plot_long <- MOC %>%
  pivot_longer(
    cols = c(
      `sp|P06731|CEAM5_HUMAN`,
      `sp|P56856|CLD18_HUMAN`,
      `sp|P09758|TACD2_HUMAN`,
      `sp|P04626|ERBB2_HUMAN`,
      `sp|Q5ZPR3|CD276_HUMAN`,
      `sp|Q96NY8|NECT4_HUMAN`,
      `sp|Q7Z7D3|VTCN1_HUMAN`
    ),
    names_to = "Protein",
    values_to = "Expression"
  ) %>%
  mutate(
    Protein = recode(
      Protein,
      "sp|P06731|CEAM5_HUMAN" = "CEACAM5",
      "sp|P56856|CLD18_HUMAN" = "CLDN18",
      "sp|P09758|TACD2_HUMAN" = "TROP2",
      "sp|P04626|ERBB2_HUMAN" = "ERBB2",
      "sp|Q5ZPR3|CD276_HUMAN" = "CD276",
      "sp|Q96NY8|NECT4_HUMAN" = "NECTIN4",
      "sp|Q7Z7D3|VTCN1_HUMAN" = "B7H4"
    ),
    log2_expr = log2(Expression),
    Cluster = factor(as.character(Cluster),
                     levels = c("1","2","3","4","5"))
  )

cluster_cols <- c(
  "1" = "#66C2A5",
  "2" = "#FC8D62",
  "3" = "#8DA0CB",
  "4" = "#E78AC3",
  "5" = "#A6D854"
)

ggplot(plot_long,
       aes(x = Cluster,
           y = log2_expr,
           fill = Cluster)) +
  geom_violin(
    trim = FALSE,
    width = 0.6,
    scale = "width",
    alpha = 0.5,
    colour = "black"
  ) +
  geom_jitter(
    aes(colour = Cluster),
    width = 0.1,
    size = 1.3,
    alpha = 0.7
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.35,
    linewidth = 0.4,
    colour = "black"
  ) +
  facet_wrap(~Protein,
             scales = "free_y",
             ncol = 3) +
  scale_fill_manual(values = cluster_cols) +
  scale_color_manual(values = cluster_cols) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 13,
                              face = "bold")
  ) +
  labs(
    x = "Cluster",
    y = "Protein Expression"
  )

# ------------------------------------------------------------
# Calculate SRM log2FC (MOC vs CRC)
# ------------------------------------------------------------

SRM_fc <- SRM_MOC %>%
  filter(Diagnosis %in% c("Mucinous_ovarian", "Colorectal")) %>%
  group_by(Diagnosis) %>%
  summarise(across(where(is.numeric), median, na.rm = TRUE)) %>%
  pivot_longer(
    cols = -Diagnosis,
    names_to = "Protein",
    values_to = "Value"
  ) %>%
  pivot_wider(names_from = Diagnosis, values_from = Value) %>%
  mutate(
    SRM_log2FC = log2(Mucinous_ovarian / Colorectal),
    Protein = sub("sp\\|([^|]+)\\|.*", "\\1", Protein)
  ) %>%
  select(Protein, SRM_log2FC)

# ------------------------------------------------------------
# Calculate DIA log2FC (MOC vs CRC)
# ------------------------------------------------------------

DIA_df <- combat_Omicodf5626
colnames(DIA_df) <- sub("sp\\|([^|]+)\\|.*", "\\1", colnames(DIA_df))

DIA_fc <- DIA_df %>%
  filter(Label %in% c("Mucinous_ovarian", "Colorectal_adenocarcinoma")) %>%
  group_by(Label) %>%
  summarise(across(where(is.numeric), median, na.rm = TRUE)) %>%
  pivot_longer(
    cols = -Label,
    names_to = "Protein",
    values_to = "Value"
  ) %>%
  pivot_wider(names_from = Label, values_from = Value) %>%
  mutate(
    DIA_log2FC = log2(Mucinous_ovarian /
                        Colorectal_adenocarcinoma)
  ) %>%
  select(Protein, DIA_log2FC)

# ------------------------------------------------------------
# Merge SRM and DIA results
# ------------------------------------------------------------

combined_fc <- inner_join(SRM_fc, DIA_fc, by = "Protein")

protein_labels <- c(
  "P00533" = "EGFR",
  "P04626" = "ERBB2",
  "P06731" = "CEACAM5",
  "P08581" = "MET",
  "P09758" = "TROP2",
  "P12004" = "PCNA",
  "P13726" = "F3",
  "P15328" = "FOLR1",
  "P15391" = "CD19",
  "P16070" = "CD44",
  "P20138" = "CD33",
  "P20273" = "CD22",
  "P21860" = "ERBB3",
  "P30530" = "AXL",
  "P40259" = "CD79B",
  "P43121" = "MUC18",
  "P46013" = "MKI67",
  "P51654" = "GPC3",
  "P56747" = "CLDN6",
  "P56856" = "CLDN18",
  "Q5ZPR3" = "CD276",
  "Q7Z7D3" = "VTCN1",
  "Q96NY8" = "NECTIN4"
)

combined_fc <- combined_fc %>%
  mutate(
    Label = protein_labels[Protein]
  )

write.csv(
  combined_fc,
  "../results/SRM_DIA_concordance.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# Concordance plot
# ------------------------------------------------------------

key_proteins <- c(
  "ERBB2",
  "CLDN18",
  "TROP2",
  "NECTIN4",
  "CEACAM5"
)

combined_fc <- combined_fc %>%
  mutate(
    Highlight = ifelse(Label %in% key_proteins,
                       "Key", "Other")
  )

ggplot(combined_fc,
       aes(x = DIA_log2FC,
           y = SRM_log2FC)) +
  geom_point(aes(colour = Highlight),
             size = 3) +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = "dashed") +
  geom_vline(xintercept = 0,
             linetype = "dotted",
             colour = "grey70") +
  geom_hline(yintercept = 0,
             linetype = "dotted",
             colour = "grey70") +
  geom_text_repel(
    data = subset(combined_fc,
                  Highlight == "Key"),
    aes(label = Label),
    colour = "#d73027",
    fontface = "bold"
  ) +
  scale_colour_manual(
    values = c("Key" = "#d73027",
               "Other" = "black")
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = "DIA log2FC (MOC vs CRC)",
    y = "SRM log2FC (MOC vs CRC)"
  )

# ============================================================
# 12_lasso.R LASSO feature selection for proteomic classifier
# ============================================================

library(tidyverse)
library(glmnet)

set.seed(123)

## Machine learning feature selection

#Feature selection was performed using LASSO regression (glmnet, alpha = 1) applied to differentially expressed proteins across tumour types.

#- Differentially expressed proteins from multiple comparisons (MOC vs HGSOC, CRC, PDAC, UGI) were used as input features.
#- Multinomial logistic regression with L1 regularization was performed using 10-fold cross-validation.
#- The optimal model was selected using lambda.min.

#Outputs:
#- `/results/LASSO_selected_features.csv` (global feature set)
#- `/results/LASSO_selected_features_per_class.csv` (class-specific features)

#Selected proteins were used for downstream validation (e.g. Werner dataset).

## 12_lasso.R performs feature selection analyses used to identify candidate diagnostic proteins.

## The final 37-protein classifier panel is provided in Supplementary Table 4 and is implemented directly within the Python machine-learning workflow.

# ============================================================
# 1. Load differential expression results
# ============================================================

# Differential abundance results generated using run_DEP().
#
# Some files originate from figure-specific analyses
# (e.g. Fig1B_MOC_vs_all.csv, Fig1D_MOC_vs_PDAC.csv),
# whereas others were generated as standalone biological
# comparisons (e.g. DEP_MOC_vs_CRC.csv, DEP_MOC_vs_HGSOC.csv).
#
# All files contain differential abundance results used for
# feature pre-filtering prior to LASSO feature selection.

HGSOC_data <- read.csv("../results/DEP_MOC_vs_HGSOC.csv")
MOC_data   <- read.csv("../results/Fig1B_MOC_vs_all.csv")
CRC_data   <- read.csv("../results/DEP_MOC_vs_CRC.csv")
UGI_data   <- read.csv("../results/STRING_MOC_vs_UGI.csv")
PDAC_data  <- read.csv("../results/Fig1D_MOC_vs_PDAC.csv")

# ============================================================
# 2. Select DE proteins (feature pre-filtering)
# ============================================================

get_de_features <- function(df) {
  df$RowNames[df$Expression %in% c("Up-regulated", "Down-regulated")]
}

selected_features <- unique(c(
  get_de_features(HGSOC_data),
  get_de_features(MOC_data),
  get_de_features(CRC_data),
  get_de_features(UGI_data),
  get_de_features(PDAC_data)
))

# ============================================================
# 3. Load full dataset
# ============================================================

df <- read.csv("../data/analysis_ready.csv")

# Keep only selected features present in dataset
features_in_df <- intersect(selected_features, colnames(df))

df_subset <- df[, c(features_in_df, "Label", "Sample_ID")]
rownames(df_subset) <- df_subset$Sample_ID

# Small bowel adenocarcinoma removed due to limited sample size.
# Remove rare class if needed
df_subset <- df_subset %>%
  filter(Label != "Small_bowel_adenocarcinoma")

# ============================================================
# 4. Prepare matrices
# ============================================================

df_subset$Label <- as.factor(df_subset$Label)

X <- df_subset %>%
  select(-Label, -Sample_ID) %>%
  as.matrix()

# glmnet requires numeric labels
y <- as.numeric(df_subset$Label) - 1  # 0-based

# ============================================================
# 5. Train/test split
# ============================================================

train_idx <- sample(seq_len(nrow(df_subset)), size = 0.8 * nrow(df_subset))

X_train <- X[train_idx, ]
y_train <- y[train_idx]

X_test  <- X[-train_idx, ]
y_test  <- y[-train_idx]

# ============================================================
# 6. LASSO (multinomial)
# ============================================================

cvfit <- cv.glmnet(
  X_train, y_train,
  family = "multinomial",
  alpha = 1,
  nfolds = 10
)

best_lambda <- cvfit$lambda.min

# ============================================================
# 7. Extract selected features
# ============================================================

coef_list <- coef(cvfit, s = best_lambda)

selected_features_per_class <- lapply(coef_list, function(mat) {
  feats <- rownames(mat)[mat[,1] != 0]
  setdiff(feats, "(Intercept)")
})

selected_features_all <- unique(unlist(selected_features_per_class))

# ============================================================
# 8. Evaluate model
# ============================================================

pred_test <- predict(cvfit, newx = X_test, s = best_lambda, type = "class")
accuracy <- mean(as.numeric(pred_test) == y_test)

print(paste("Test set accuracy:", round(accuracy, 3)))

# ============================================================
# 9. Save outputs (CRITICAL FOR REPRODUCIBILITY)
# ============================================================

write.csv(
  selected_features_all,
  "../results/LASSO_selected_features.csv",
  row.names = FALSE
)

write.csv(
  stack(selected_features_per_class),
  "../results/LASSO_selected_features_per_class.csv",
  row.names = FALSE
)

# ============================================================
# 10. Plot CV curve
# ============================================================

pdf("figures/LASSO_CV_plot.pdf", width = 6, height = 5)
plot(cvfit)
dev.off()



# Transcriptomic and Proteomic Profiling of Mucinous Ovarian Cancer

---
  ## Analysis Workflow
  
  Run scripts in the following order:
  
  - `14_prepare_MOC_proteomics_data.R`
- `15_prepare_MOC_transcriptomic_data.R`
- `16_Supplementary_Figure.7.R`

---
  
  ## Input Data
  
  The primary input dataset for all transcriptomic analyses is:
  
  | File | Purpose |
  |--------|--------|
  | `batch_corrected_matrices.rds` | Batch-corrected transcriptomic data from the GAMuT MOC cohort, used for gene expression visualization and mRNA-protein correlation analysis. |
  
  The primary input dataset for all proteomic analyses is:
  
  | File | Purpose |
  |--------|--------|
  | `expr_PC_Imp.csv` | Imputed proteomic data from the ProCan MOC cohort, used for mRNA-protein correlation analysis. |
  
  UniProt accessions were mapped to gene symbols using the UniProt REST API, and gene symbols to Ensembl gene IDs using org.Hs.eg.db and the Ensembl REST API. The mapping is provided as:
  
  | File | Purpose |
  |--------|--------|
  | `uniprot_gene_map.csv` | UniProt accession to gene symbol and Ensembl gene ID mapping, used in `16_Supplementary_Figure.7.R`. |
  | `uniprot_gene_map_info.txt` | Query date and database versions (UniProt, Ensembl, org.Hs.eg.db) used for the mapping. |
  
  Metadata are available within the PRIDE repository.

---
  
  ## Figure-to-Script Mapping
  
  | Figure | Script |
  |----------|----------|
  | Supplementary Figure 7 | `16_Supplementary_Figure.7.R` |
  
  # ============================================================
# 14_prepare_MOC_proteomics_data.R
# ============================================================
#
# Purpose:
#   Generate the imputed MOC proteomic dataset used for
#   downstream analyses.
#
# Input:
#   input/MOC_Matrix_forshare.xlsx
#
# Outputs:
#   data/expr_PC_Imp.csv
#   data/meta_ProCan.csv
#

suppressMessages({
  library(tidyverse)
  library(MSnbase)
  library(readxl)
})

# Corresponds to the combat_df dataset filtered for the MOC cases
raw <- as.data.frame(read_excel("input/MOC_Matrix_forshare.xlsx"))

# Protein abundance matrix (proteins x cases); replicates of the same case are averaged
expr_C <- raw %>%
  select(-Label, -Sample_ID) %>%
  group_by(Collaborator_ID) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  column_to_rownames("Collaborator_ID") %>%
  t()

# Proteins missing in all replicates of a case are averaged to NaN; convert them to NA
expr_C[is.nan(expr_C)] <- NA

# Sample metadata for the proteomic data
meta_C <- raw %>%
  select(Sample_ID = Collaborator_ID, Label) %>%
  distinct()

# Quantile normalisation, then MinProb imputation of missing values (if not already done as part of original workflow)
expr_PC_Imp <- MSnbase::MSnSet(exprs = expr_C) %>%
  MSnbase::normalise("quantiles") %>%
  MSnbase::impute(method = "MinProb", q = 0.01) %>%
  Biobase::exprs()

# ============================================================
# Save analysis-ready dataset
# ============================================================

write.csv(expr_PC_Imp, "data/expr_PC_Imp.csv")
write.csv(meta_C, "data/meta_ProCan.csv")

# ============================================================
# 15_prepare_MOC_transcriptomic_data.R
# ============================================================
#
# Purpose:
#   Generate the batch-corrected MOC transcriptomic dataset used
#   for downstream analyses.
#
# Inputs:
#   input/salmon.merged.gene_counts_length_scaled_PE.tsv
#   input/salmon.merged.gene_counts_length_scaled_SE.tsv
#   input/RNA_metadata.csv
#
# Output:
#   data/batch_corrected_matrices.rds
#

suppressMessages({
  library(tidyverse)
  library(limma)
  library(edgeR)
  library(DESeq2)
})

meta <- read_csv("input/RNA_metadata.csv", col_types = cols())

# Gene-level counts (Salmon, length-scaled TPM) from paired-end and single-end sequencing
counts_pe <- read_tsv("input/salmon.merged.gene_counts_length_scaled_PE.tsv", col_types = cols())
counts_se <- read_tsv("input/salmon.merged.gene_counts_length_scaled_SE.tsv", col_types = cols())

counts_merged <- dplyr::full_join(counts_pe, counts_se, by = c("gene_id", "gene_name"))

# Integer count matrix (genes x samples), with columns in metadata order
count_mat <- counts_merged %>%
  dplyr::select(-gene_name) %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix() %>%
  round()

storage.mode(count_mat) <- "integer"
count_mat <- count_mat[, meta$sample]

# Filter out lowly expressed genes, but keep the classifier genes KLHL40 and LDHAL6B
keep <- filterByExpr(count_mat, group = meta$class_pca)
keep[rownames(count_mat) %in% c("ENSG00000157119", "ENSG00000171989")] <- TRUE
count_mat <- count_mat[keep, ]

# Variance-stabilising transformation (VST)
vst_mat <- DESeq2::vst(count_mat, blind = TRUE)

# Remove batch effects while preserving differences between classes
# (corrected values are used for visualization and correlation only)
design <- model.matrix(~ class_pca, data = meta)
vst_removebatch_composite <- removeBatchEffect(vst_mat, batch = meta$batch, design = design)

# Gene ID to gene name map
gene_id_map <- counts_merged %>% dplyr::select(gene_id, gene_name)

# Sample annotations for heatmap visualization
sample_annotations <- meta[, c("sample", "original_id", "class_pca", "Classification", "GROUP")]

# Average replicate samples of the same case for mRNA-protein correlation analysis
# (matrix columns follow metadata order, so meta$`Sample ID` gives the case of each column)
vst_removebatch_composite_by_id <- avearrays(vst_removebatch_composite,
                                             ID = paste0("s_", meta$`Sample ID`))

sample_annotations_by_id <- meta %>%
  dplyr::distinct(`Sample ID`, class_pca, Classification, GROUP) %>%
  dplyr::mutate(sample = paste0("s_", `Sample ID`)) %>%
  dplyr::select(sample, `Sample ID`, class_pca, Classification, GROUP)

# ============================================================
# Save analysis-ready dataset
# ============================================================

saveRDS(
  list(
    vst_removebatch_composite = vst_removebatch_composite,
    gene_id_map = gene_id_map,
    sample_annotations = sample_annotations,
    vst_removebatch_composite_by_id = vst_removebatch_composite_by_id,
    sample_annotations_by_id = sample_annotations_by_id
  ),
  "data/batch_corrected_matrices.rds"
)

# ============================================================
# 16_Supplementary_Figure.7.R
# ============================================================
#
# Purpose:
#   Generate Supplementary Figure 7:
#   - Heatmaps of RNA expression for selected genes
#   - Scatter plots of mRNA-protein correlation for selected genes
#   - Distribution of global mRNA-protein correlation
#
# Inputs:
#   data/batch_corrected_matrices.rds
#   data/expr_PC_Imp.csv
#   input/id_mapping_protein_rna.csv
#   data/uniprot_gene_map.csv
#
# Outputs:
#   output/Heatmap_RNA_long_gene_list_vst_removebatch.pdf
#   output/Heatmap_RNA_short_gene_list_vst_removebatch.pdf
#   output/RNA_protein_correlation_long_gene_list_vst.csv
#   output/RNA_protein_correlation_short_gene_list_vst.csv
#   output/RNA_protein_scatter_long_gene_list_vst.pdf
#   output/RNA_protein_scatter_short_gene_list_vst.pdf
#   output/RNA_protein_correlation_all_expressed_genes_vst.csv
#   output/RNA_protein_correlation_distribution_vst.pdf
#

suppressMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
})

# ============================================================
# Load data
# ============================================================

corrected <- readRDS("data/batch_corrected_matrices.rds")
gene_id_map <- corrected$gene_id_map

# Remove sample ID prefixes ("sample_", "s_") and replicate suffixes ("-a", "-b", "-c")
clean_sample_id <- function(x) sub("-[abc]$", "", sub("^(sample_|s_)", "", x))

# Per-sample matrix (replicates kept), used for heatmaps
vst_mat <- corrected$vst_removebatch_composite
colnames(vst_mat) <- clean_sample_id(colnames(vst_mat))

sample_annotations <- corrected$sample_annotations
sample_annotations$sample <- clean_sample_id(sample_annotations$sample)

# Per-case matrix (replicates averaged), used for mRNA-protein correlation
vst_by_id <- corrected$vst_removebatch_composite_by_id
colnames(vst_by_id) <- clean_sample_id(colnames(vst_by_id))

# Machine learning classifier candidates
long_gene_list <- c("CEACAM20", "KLHL40", "ANXA10", "GCG", "CA9", "MACROD2", "LDHAL6B",
                    "ERBB2", "VSIG1", "CLDN18", "WNK3", "GOLM1", "AGTPBP1", "CEACAM5",
                    "PGC", "TMEM175", "ERC1", "CHI3L1", "KRT20", "CLCA1", "GALNT1", "SATB2",
                    "OCLN", "GPA33", "S100P", "B4GALT4", "BCAT1", "SPRR3", "CDH17", "KRT7",
                    "LIPF", "TSPAN8", "ACADL", "RPN2", "WT1", "SPON1", "KBTBD11")

# Antibody-drug conjugate targets
short_gene_list <- c("CD276", "CEACAM5", "VTCN1", "CLDN18", "ERBB2", "NECTIN4", "TACSTD2")

# Subset a matrix (rows = gene IDs) to the given gene names;
# returns a data frame with gene names as row names and a gene_id column
subset_by_gene_name <- function(mat, genes) {
  gene_map <- gene_id_map %>% filter(gene_name %in% genes, gene_id %in% rownames(mat))
  out <- mat[gene_map$gene_id, , drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene_id")
  rownames(out) <- gene_map$gene_name
  out
}

# ============================================================
# Heatmaps of RNA expression
# ============================================================

rna_mat_long_vst <- subset_by_gene_name(vst_mat, long_gene_list)
rna_mat_short_vst <- subset_by_gene_name(vst_mat, short_gene_list)

# Column annotation: classification of each sample
annotation_col_class <- data.frame(Classification = sample_annotations$class_pca,
                                   row.names = sample_annotations$sample)

class_levels <- sort(unique(annotation_col_class$Classification))
class_colors <- brewer.pal(max(3, length(class_levels)), "Set2")[seq_along(class_levels)]
names(class_colors) <- class_levels
annotation_colors_class <- list(Classification = class_colors)

# Row-scaled heatmap (z-scores capped at -3 to 3); rows and columns clustered by correlation distance
draw_heatmap <- function(mat, file, title) {
  pheatmap(
    mat,
    scale = "row",
    breaks = seq(-3, 3, length.out = 101),
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize_col = 5,
    cellwidth = 6,
    cellheight = 12,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "complete",
    border_color = NA,
    annotation_col = annotation_col_class,
    annotation_colors = annotation_colors_class,
    main = title,
    filename = file
  )
}

# Generate matrix for pheatmap()
long_heatmap_dt <- rna_mat_long_vst %>%
  select(-gene_id) %>%
  as.matrix()

short_heatmap_dt <- rna_mat_short_vst %>%
  select(-gene_id) %>%
  as.matrix()

draw_heatmap(long_heatmap_dt,
             "output/Heatmap_RNA_long_gene_list_vst_removebatch.pdf",
             "RNA expression (vst, batch-corrected) - long gene list")

draw_heatmap(short_heatmap_dt,
             "output/Heatmap_RNA_short_gene_list_vst_removebatch.pdf",
             "RNA expression (vst, batch-corrected) - short gene list")

# ============================================================
# Match RNA and protein data
# ============================================================

# UniProt accession to gene symbol and Ensembl gene ID mapping
uniprot_map <- read.csv("data/uniprot_gene_map.csv", colClasses = "character")

# Proteomic data annotated with gene names
# (proteins mapped to several Ensembl IDs are kept as several rows)
protein_data <- read.csv("data/expr_PC_Imp.csv", check.names = FALSE) %>%
  dplyr::rename(uniprot = 1) %>%
  dplyr::left_join(uniprot_map, by = "uniprot", relationship = "many-to-many") %>%
  dplyr::rename(gene_name = gene)

# Cases with matched transcriptomic and proteomic data
id_map <- read.csv("input/id_mapping_protein_rna.csv", colClasses = "character")
matched_map <- id_map %>% filter(Sample_ID_RNA != "")

rna_cols <- matched_map$Sample_ID_RNA
protein_cols <- matched_map$Sample_ID_protein

# Spearman correlation; returns NA if the test cannot be computed
spearman <- function(x, y) {
  test <- tryCatch(suppressWarnings(cor.test(x, y, method = "spearman")),
                   error = function(e) NULL)
  if (is.null(test)) return(tibble(spearman_r = NA_real_, p_value = NA_real_))
  tibble(spearman_r = unname(test$estimate), p_value = test$p.value)
}

# All genes with both transcriptomic and proteomic data, matched by Ensembl gene ID
rna_protein_pairs <- gene_id_map %>%
  filter(gene_id %in% rownames(vst_by_id)) %>%
  inner_join(protein_data, by = "gene_id", suffix = c("_rna", "_protein"))

# Ensembl gene IDs matched to multiple proteins: keep one protein per gene, prioritising
#   1) proteins whose gene symbol matches the RNA gene name
#   2) the protein with the highest mean abundance
rna_protein_pairs <- rna_protein_pairs %>%
  mutate(name_match = coalesce(gene_name_rna == gene_name_protein, FALSE),
         mean_abundance = rowMeans(across(where(is.numeric)))) %>%
  group_by(gene_id) %>%
  arrange(desc(name_match), desc(mean_abundance), .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  ungroup()

# ============================================================
# mRNA-protein correlation of selected genes
# ============================================================

rna_vst_long <- subset_by_gene_name(vst_by_id, long_gene_list)
rna_vst_short <- subset_by_gene_name(vst_by_id, short_gene_list)

# Keep genes with both transcriptomic and proteomic data
long_genes_both <- rownames(rna_vst_long)[rna_vst_long$gene_id %in% rna_protein_pairs$gene_id]
short_genes_both <- rownames(rna_vst_short)[rna_vst_short$gene_id %in% rna_protein_pairs$gene_id]

# RNA and protein values of matched cases, one row per gene per case
# (protein values are matched by the Ensembl gene ID of the RNA row)
get_pairs <- function(rna_mat, genes) {
  map_dfr(genes, ~ tibble(
    gene = .x,
    rna = unlist(rna_mat[.x, rna_cols]),
    protein = unlist(rna_protein_pairs[match(rna_mat[.x, "gene_id"], rna_protein_pairs$gene_id), protein_cols])
  ))
}

pairs_long <- get_pairs(rna_vst_long, long_genes_both)
pairs_short <- get_pairs(rna_vst_short, short_genes_both)

cor_long_vst <- pairs_long %>%
  group_by(gene) %>%
  summarise(n = n(), spearman(rna, protein)) %>%
  arrange(desc(spearman_r))

cor_short_vst <- pairs_short %>%
  group_by(gene) %>%
  summarise(n = n(), spearman(rna, protein)) %>%
  arrange(desc(spearman_r))

write.csv(cor_long_vst, "output/RNA_protein_correlation_long_gene_list_vst.csv", row.names = FALSE)
write.csv(cor_short_vst, "output/RNA_protein_correlation_short_gene_list_vst.csv", row.names = FALSE)

# Scatter plots faceted by gene, ordered by correlation coefficient
plot_marker_scatter <- function(pairs, cor_df, title) {
  labels <- cor_df %>%
    mutate(label = sprintf("%s (r=%.2f, p=%.3g)", gene, spearman_r, p_value))
  
  pairs %>%
    left_join(select(labels, gene, label), by = "gene") %>%
    mutate(label = factor(label, levels = labels$label)) %>%
    ggplot(aes(x = rna, y = protein)) +
    geom_point(size = 1.6, alpha = 0.8, color = "#2A3F8F") +
    geom_smooth(method = "lm", se = FALSE, color = "firebrick", linewidth = 0.5) +
    facet_wrap(~ label, scales = "free") +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(size = 7)) +
    labs(title = title, x = "RNA expression (vst)", y = "Protein abundance (log2, imputed)")
}

pdf("output/RNA_protein_scatter_long_gene_list_vst.pdf", width = 16, height = 14)
print(plot_marker_scatter(pairs_long, cor_long_vst,
                          "Long gene list: RNA (vst, batch-corrected) vs protein"))
dev.off()

pdf("output/RNA_protein_scatter_short_gene_list_vst.pdf", width = 10, height = 9)
print(plot_marker_scatter(pairs_short, cor_short_vst,
                          "Short gene list: RNA (vst, batch-corrected) vs protein"))
dev.off()

# ============================================================
# Global mRNA-protein correlation
# ============================================================

rna_all <- vst_by_id[rna_protein_pairs$gene_id, rna_cols, drop = FALSE]
protein_all <- as.matrix(rna_protein_pairs[, protein_cols])

all_cor_vst <- rna_protein_pairs %>%
  transmute(gene_id, gene_name = gene_name_rna, uniprot, n = length(rna_cols)) %>%
  bind_cols(map_dfr(seq_len(nrow(rna_all)), ~ spearman(rna_all[.x, ], protein_all[.x, ])))

write.csv(all_cor_vst, "output/RNA_protein_correlation_all_expressed_genes_vst.csv", row.names = FALSE)

# Summary statistics shown on the histogram
binwidth <- 0.01
n_genes <- nrow(all_cor_vst)
sig_pos_pct <- mean(all_cor_vst$p_value < 0.05 & all_cor_vst$spearman_r > 0, na.rm = TRUE) * 100
median_r <- median(all_cor_vst$spearman_r, na.rm = TRUE)
x_range <- range(all_cor_vst$spearman_r, na.rm = TRUE)

# Histogram of correlation coefficients with a density curve scaled to gene counts
p_all_cor_vst <- ggplot(all_cor_vst, aes(x = spearman_r)) +
  geom_histogram(binwidth = binwidth, fill = "#7B93DB", color = "white", linewidth = 0.05) +
  geom_density(aes(y = after_stat(density) * n_genes * binwidth), color = "#2A3F8F", linewidth = 1) +
  annotate("text", x = x_range[1], y = Inf, vjust = 1.5, hjust = 0, size = 4,
           label = sprintf("%.1f%% significant positive\ncorrelation (p<0.05)", sig_pos_pct)) +
  annotate("text", x = x_range[2], y = Inf, vjust = 1.5, hjust = 1, size = 4,
           label = sprintf("median = %.2f", median_r)) +
  theme_bw() +
  labs(
    title = "Tumor gene-wise mRNA-protein correlation (vst)",
    subtitle = sprintf("%d genes with matched RNA-protein data, n=%d cases", n_genes, length(rna_cols)),
    x = "Spearman's correlation (r)",
    y = "Number of genes"
  )

pdf("output/RNA_protein_correlation_distribution_vst.pdf", width = 6, height = 5)
print(p_all_cor_vst)
dev.off()
