################################################
# Title: Tofacitinib for Immune Skin Conditions in Down Syndrome: Analysis of whole blood Transcriptome 
# Author(s):
#   - Neetha Paul Eduthan
#   - Matthew Galbraith
# affiliation(s):
#   - Linda Crnic Institute for Down syndrome
#   - University of Colorado Anschutz
################################################

### Summary:  
# DESeq2 analysis of differential gene expression in whole blood transcriptome data. Paired-end 
# strand-specific globin-depleted polyA+ libraries generated and sequenced by Novogene. Assess differential 
# gene expression at baseline and 16 weeks of Tofacitinib (TOFA) treatment, and between individuals 
# with (T21) and without (D21) Down syndrome in Human Trisome Project (HTP) data set, followed by 
# comparison of treatment and T21 effects.
# See README.md for more details
# 

### Data type(s):
# Clinical trial (TOFA) datasets:
#    A. Participant-level metadata
#    B. Visit/Event-level metadata
#    C. PAXgene whole blood RNAseq data (Counts); DOI: 
#    D. PAXgene whole blood RNAseq data (RPKMs); DOI: 10.5281/zenodo.19954573
#      
# Human Trisome Project (HTP) datasets:
#    D. Participant-level metadata
#    E. Visit/Event-level metadata
#    F. PAXgene whole blood RNAseq data (Counts); DOI: 
#    G. PAXgene whole blood RNAseq data (RPKMs); DOI: 10.5281/zenodo.20044079

# 0 General Setup -----
# RUN THIS FIRST TIME - Initialize and install packages with renv:
# renv::init(bioconductor = TRUE)
#
# To install the exact versions of all R packages base on renv.lock file (requires matching R version):
# renv::restore()

## 0.1 Load required libraries ----
library("DESeq2") # differential expression analysis
library("edgeR") # for cpm() function (can also be used for differential expression analysis)
library("limma") # for removeBatchEffect() function (can also be used for differential expression analysis)
library("BiocParallel") # enables mutli-cpu for some of DEseq2 functions
library("apeglm") # used with DESeq2 to 'moderate' fold-changes
library("readxl") # reading Excel files
library("openxlsx") # for exporting results as Excel workbooks
# library("genefilter") # used for function rowVars?
library("tidyverse") # required for ggplot2, dplyr etc
library("ggforce") # used for sina plots
library("ggrastr") # required for rasterizing some layers of plots
library("ggrepel") # required for using geom_text and geom_text_repel() to make sample labels for PCA plot
library("ggsignif") # required for adding pvalue annotations
library("RColorBrewer") # color palettes
library("circlize") # color scale generation
library("tidyHeatmap") # tidy heatmaps
# library("factoextra") # extraction and visualization for PCA
library("conflicted") # force all conflicts to become errors
conflicts_prefer( # declare preferences in cases of conflict
  dplyr::filter,
  dplyr::select,
  dplyr::count,
  dplyr::rename,
  base::paste,
  matrixStats::rowVars,
  dplyr::bind_rows
)
library("here")
#

## 0.2 Set file name parameters ----
# TOFA files:
tofa_participant_meta_data_file <- here("data/TOFA_Participant_metadata_zenodo_v1.txt") # Available on request
tofa_visit_meta_data_file <- here("data/TOFA_Visit_metadata_zenodo_v1.txt") # Available on request
tofa_counts_file <- here("data/TOFA_PAXgene_RNAseq_Counts_zenodo_v1.txt.gz") # Source: 10.5281/zenodo.19954573
tofa_rpkms_file <- here("data/TOFA_PAXgene_RNAseq_RPKMs_zenodo_v1.txt.gz") # Source: 10.5281/zenodo.19954573
# HTP data files
htp_participant_metadata_file <- here("data/HTP_Participant_metadata_zenodo_v1.txt")
htp_visit_metadata_file <- here("data/HTP_Visit_metadata_zenodo_v1.txt")
htp_counts_file <- here("data/HTP_PAXgene_RNAseq_Counts_zenodo_v1.txt.gz") # Source: 10.5281/zenodo.19954573
htp_rpkms_file <- here("data/HTP_PAXgene_RNAseq_RPKMs_long_zenodo_v1.txt.gz") # Source: 10.5281/zenodo.19954573
#
min_cpm <- 0.5 # used for low count filtering; default is 0.5
min_samples <- "auto" # used for low count filtering; use a number, "all", or "auto" (sets to half number of samples)
#
standard_colors <- c("Control" = "grey30", "T21" = "#009b4e", "Baseline" = "#999999", "16 week" = "#6baed6")
out_file_prefix <- "TOFA_DESeq2.R_" # should match this script title
# End required parameters ###
source(here("helper_functions_DESeq.R")) # load helper functions
#

# 1 TOFA 16 weeks vs. baseline analysis ------
## 1.1 Read in TOFA files ------
### 1.1.1 Read in TOFA RNAseq data ----

# Read in TOFA Counts data 
tofa_counts_long <- tofa_counts_file |> 
  read_tsv() |> 
  rename(raw_count = Value)
#
tofa_counts_long # 11,707,766 rows
tofa_counts_long |> distinct(ParticipantID) # 42 Participants
tofa_counts_long |> distinct(VisitID) # 193 VisitIDs = PAXgene Whole Blood RNA samples
tofa_counts_long |> distinct(EnsemblID) # 60,662 EnsemblIDs
tofa_counts_long |> distinct(Sequencing_batch) # 3 sequencing batches
#

# Read in TOFA RPKMs data 
tofa_rpkms_long <- tofa_rpkms_file |> 
  read_tsv() |> 
  rename(RPKM = Value)
#
tofa_rpkms_long # 11,707,766 rows
tofa_rpkms_long |> distinct(ParticipantID) # 42 Participants
tofa_rpkms_long |> distinct(VisitID) # 193 VisitIDs = PAXgene Whole Blood RNA samples
tofa_rpkms_long |> distinct(EnsemblID) # 60,662 EnsemblIDs
#

### 1.1.2. Read in TOFA meta data ----
# Read in Participant level metadata
tofa_participant_meta_data <- tofa_participant_meta_data_file |> 
  read_tsv() |> 
  mutate(
    Sex = fct_relevel(Sex, c("Female", "Male")) # set factor levels
  )
#
tofa_participant_meta_data # 47 rows
tofa_participant_meta_data |> distinct(ParticipantID) # 47 participants in this table
tofa_participant_meta_data |> count(Endpoint_eligible, Participant_notes) # 5 participants not eligible for Endpoint analyses
#
# Read in Event/Visit level meta data 
tofa_visit_meta_data <- tofa_visit_meta_data_file |> 
  read_tsv() |> 
  # keep only the timepoints included in this analysis
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  mutate(
    Event_Name = fct_relevel(Event_Name, c("Baseline", "16 week"))
  ) 
#
tofa_visit_meta_data # 85 rows
tofa_visit_meta_data |> distinct(ParticipantID) # 43 Participants with samples
tofa_visit_meta_data |> distinct(VisitID) # 85 Visits/Samples
#

## 1.2 Join RNAseq data with metadata ----
# Join RNAseq samples with Participant and Visit metadata 
tofa_meta_data <- tofa_counts_long |> 
  distinct(ParticipantID, VisitID, Sequencing_batch) |>
  mutate(Sequencing_batch = fct_relevel(as.character(Sequencing_batch), c("1", "2", "3"))) |>
  inner_join(tofa_participant_meta_data) |> 
  filter(Endpoint_eligible == TRUE) |> # ensure only endpoint eligible participants
  inner_join(tofa_visit_meta_data) |> # keep only baseline and 16 week samples
  mutate(Event_Name = fct_relevel(Event_Name, c("Baseline", "16 week")),
         Sequencing_batch = fct_drop(Sequencing_batch)) |> 
  arrange(Sequencing_batch, ParticipantID, Event_Name) |> # set order for plots
  mutate(VisitID = fct_inorder(VisitID))
#
tofa_meta_data |> count(Event_Name, Sequencing_batch) 
# Event_Name Sequencing_batch       n
# 1 Baseline   1                   25
# 2 Baseline   2                   17
# 3 16 week    1                   18
# 4 16 week    2                   21

# generate a groups df to use with deseq containing select co-variates from metadata
tofa_groups <- tofa_meta_data |>
  select(ParticipantID, VisitID, Sequencing_batch, Event_Name, Age_years_at_visit, Sex)

# Join Counts data with metadata
tofa_counts_meta <- tofa_counts_long |> 
  mutate(Sequencing_batch = fct_relevel(as.character(Sequencing_batch), c("1", "2", "3"))) |>
  inner_join(tofa_meta_data) |> # returns 4,913,622 of 11,707,766 rows
  mutate(Sequencing_batch = fct_drop(Sequencing_batch)) 
tofa_counts_meta |> distinct(VisitID) # 81 samples
tofa_counts_meta |> distinct(VisitID, Event_Name) |> count(Event_Name) 
# Event_Name       n
# 1 Baseline      42
# 2 16 week       39

# Check raw read count distributions across samples (not normalized) 
tofa_counts_meta |>
  mutate(VisitID = fct_relevel(VisitID, tofa_groups |> arrange(Event_Name) |> pull(VisitID) %>%  as.character)) |> 
  ggplot(aes(VisitID, log2(raw_count + 0.1), color = Event_Name)) + # CUSTOMIZE AS NEEDED
  geom_sina(size = 0.01) +
  geom_boxplot(notch = FALSE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.2, color = "black", fill = "transparent", size = 0.75) +
  scale_color_manual(values = standard_colors) +
  facet_grid(~ Sequencing_batch, scale = "free_x", space = "free") +
  labs(title = "TOFA: Raw read count distributions across samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_TOFA_counts_unfiltered.png")), width = 30, height = 5, units = "in")

# Join RPKMs data with metadata
tofa_rpkms_meta <- tofa_rpkms_long |> 
  mutate(Sequencing_batch = fct_relevel(as.character(Sequencing_batch), c("1", "2", "3"))) |>
  inner_join(tofa_meta_data) |> # returns 4,913,622 of 11,707,766 rows
  mutate(Sequencing_batch = fct_drop(Sequencing_batch)) 
tofa_rpkms_meta |> distinct(VisitID) # 81 samples
tofa_rpkms_meta |> distinct(VisitID, Event_Name) |> count(Event_Name) 
# Event_Name       n
# 1 Baseline      42
# 2 16 week       39

## Convert Counts and RPKMs to wide format 
# Subselect counts data and rpkms to match meta_data
tofa_counts <- tofa_counts_long |>
  select(EnsemblID, VisitID, raw_count) |>
  pivot_wider(names_from = VisitID, values_from = raw_count) |>
  select(EnsemblID, tofa_meta_data |> pull(VisitID))
#

tofa_rpkms <- tofa_rpkms_long |>
  select(EnsemblID, VisitID, RPKM) |>
  pivot_wider(names_from = VisitID, values_from = RPKM) |>
  select(EnsemblID, tofa_meta_data |> pull(VisitID))
#

## 1.3 Filter/remove genes with low expression ----
#
## Summary of total reads per sample 
tofa_counts |> select(-EnsemblID) |> colSums() |> summary()
tofa_counts |> 
  select(-EnsemblID) |> 
  colSums() |>
  enframe(name = "VisitID", value = "Total_reads") |> 
  arrange(Total_reads)
#
# Keep only rows (transcripts / genes) with greater than `min_cpm` cpm in `min_samples`
# NOTE: 10 counts=0.5 cpm for 20 million reads,  15 counts=0.5 cpm for 30 million reads...
# from Michael Love: https://support.bioconductor.org/p/95840/
# The independent filtering is designed only to filter out low count genes to
# the extent that they are not enriched with small p-values. Here the problem is
# not independent filtering, but that these two genes get a small p-value rather
# than being filtered or having an insignificant p-value. Datasets can be
# different in many ways, and for whatever reason, these two genes survive the
# filtering and get a counterintuitive small p-value. I'd recommend you just use
# a more strict filter in the very beginning, e.g. at least three samples with
# counts greater than 10
## Filter by minimum counts per million 
# Check min_samples and calculate if needed (ie if number is not supplied)
if (min_samples == "all") {
  min_samples=ncol(tofa_counts) - 1
} else if (min_samples == "auto") {
  min_samples=(ncol(tofa_counts) - 1) / 2
}
tofa_before <- tofa_counts %>% 
  transmute(
    EnsemblID = EnsemblID, 
    row_sum = rowSums(select(., -EnsemblID))
  ) |> 
  filter(row_sum > 0) |> 
  nrow()
tofa_cpm_data <- tofa_counts |> 
  column_to_rownames("EnsemblID") |> 
  cpm()
tofa_keep <- tofa_cpm_data |> 
  as_tibble(rownames = "EnsemblID") |> 
  pivot_longer(-EnsemblID, names_to = "VisitID", values_to = "cpm") |> 
  mutate(cpm > min_cpm) |> # check against min_cpm
  filter(`cpm > min_cpm` == TRUE) |> # and filter
  dplyr::count(EnsemblID) |> # count samples remaining per EnsemblID
  filter(n >= min_samples) # filter against min_samples
tofa_counts_filtered <- tofa_counts |> 
  filter(EnsemblID %in% tofa_keep$EnsemblID)
## Summarize rows before and after filtering 
cat("Total number of rows: ", tofa_counts |> nrow())
cat("Number of rows with non-zero total read counts before filtering: ", tofa_before)
cat("Number of rows after filtering: ", tofa_counts_filtered |> nrow()) # 16720
#

# Gather to long format
tofa_counts_filtered_long <- tofa_counts_filtered |>
  pivot_longer(-EnsemblID, names_to = "VisitID", values_to = "raw_count") |>
  inner_join(tofa_meta_data)

# Check filtered read count distributions for each sample (not normalized) 
tofa_counts_filtered_long |>
  mutate(VisitID = fct_relevel(VisitID, tofa_groups |> arrange(Event_Name) |> pull(VisitID) %>% as.character)) |>
  ggplot(aes(VisitID, log2(raw_count + 0.1), color = Event_Name)) + # CUSTOMIZE AS NEEDED
  geom_sina(size = 0.01) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.2, color = "black", fill = "transparent", size = 0.75) +
  scale_color_manual(values = standard_colors) +
  facet_grid(~ Sequencing_batch, scale = "free_x", space = "free") +
  labs(title = "TOFA: CPM-filtered read count distributions across samples", x = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_TOFA_counts_filtered.png")), width = 30, height = 5, units = "in")
#

## 1.4. Generate DESeqDataSet object(s) and run DESeq2 analysis ----
tofa_formula <- as.formula(paste0("~", "Sequencing_batch +", "ParticipantID +", "Event_Name"))
#
dds_tofa <- DESeqDataSetFromMatrix(
  countData = tofa_counts_filtered |>
    select(EnsemblID, tofa_groups |> pull(VisitID)) |>  # ensures correct order of columns
    column_to_rownames("EnsemblID"), # must be converted to data frame from tibble
  colData = tofa_groups,
  design = tofa_formula
)
# Check meta data read in to DESeqDataSet (multivariable version):
colData(dds_tofa)
#

# Run DESeq2 analysis
dds_tofa <- DESeq(dds_tofa)
#

## 1.5. Get DESeq2 results ---- 
### 1.5.1 Define comparisons of interest ----
tofa_comparisons <- list( 
  c("Event_Name", "16.week", "Baseline")
)
#
### 1.5.2 Results summaries ----
dds_tofa |> get_results_sum(tofa_comparisons[[1]], show_ind_filt_off=FALSE)
#

### 1.5.3 Assemble DESeq2 results table(s) ----
# generate results table with adjusted log2FCs using lfcShrink()
res_multi_BatchStudyID_16.week_vs_Baseline <- dds_tofa |> # SET to correct dds object
  get_results_tbl(
    contrast = tofa_comparisons[[1]],
    shrink_type = "apeglm"
  )

#
res_multi_BatchStudyID_16.week_vs_Baseline

# volcano plot of results
res_multi_BatchStudyID_16.week_vs_Baseline %>%
  volcano_plot_lab(
    title = "TOFA: 16 week vs. Baseline",
    subtitle = paste0(
      paste("Model:", paste(tofa_formula, collapse = ""), "\n"), # CUSTOMIZE formula
      "[Down: ", (.) %>% filter(padj < 0.1 & FoldChange_adj <1) |> nrow(), "; Up: ", (.) |> filter(padj < 0.1 & FoldChange_adj >1) |> nrow(), "]"
    ),
    labels = TRUE,
    n_labels = 3,
    raster = TRUE
  ) +
  expand_limits(x = c(-1, 1) * 2)
ggsave(filename = here("plots", paste0(out_file_prefix, "TOFA_16week_vs_Baseline", "_Volcano.png")), width = 5, height = 5, units = "in")
#

## 1.6 GSEA Hallmarks analysis ----
# Download Human GSEA hallmarks from MSigDB: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# Read in Human GSEA hallmarks
hallmarks <- here("data/GSEA/human/h.all.v7.4.symbols.gmt") %>%
  fgsea::gmtPathways(gmt.file = .)

### 1.6.1 Generate ranks ----
ranks_tofa <- res_multi_BatchStudyID_16.week_vs_Baseline |>
  filter(!is.na(log2FoldChange_adj)) |> # need to remove NA rows that will break enplot
  select(ID = Gene_name, t = log2FoldChange_adj) |>
  arrange(-abs(t)) |>
  distinct(ID, .keep_all = TRUE) |> # to avoid duplicates
  tibble::deframe() # convert to named num vector
#
### 1.6.2 Run fgsea ----
set.seed(123) # set this for reproducibility
hallmarks_tofa <- run_fgsea2(geneset = hallmarks, ranks = ranks_tofa, weighted = FALSE)
hallmarks_tofa |> filter(padj < 0.1) # 36 sig. hallmarks
#

# 2 HTP T21 vs D21 analysis -----
## 2.1 Read in HTP files ----
### 2.1.1 Read in HTP RNAseq files -----
htp_counts_long <- htp_counts_file |> 
  read_tsv() |> 
  rename(raw_count = Value)
#
htp_rpkms_long <- htp_rpkms_file |> 
  read_tsv()

### 2.1.2 Read in HTP metadata ------
htp_metadata <- htp_participant_metadata_file |> 
  read_tsv() |>
  inner_join(
    htp_visit_metadata_file |> 
      read_tsv()
  )
#

## 2.2 Join RNAseq data with metadata  ----
htp_counts_meta <- htp_counts_long |> 
  inner_join(htp_metadata)

htp_rpkms_meta <- htp_rpkms_long |> 
  inner_join(htp_metadata)

# Generate a groups df with relevant metadata columns
# Subset to ages between 12-40 to match with age of TOFA cohort
htp_groups <- htp_counts_meta |> 
  distinct(ParticipantID, VisitID, Karyotype, Sex, Age_years_at_visit, Sample_source_code) |> # 400 rows
  filter(dplyr::between(Age_years_at_visit, 12, 40)) |> # 282 samples
  mutate(
    Karyotype = fct_relevel(Karyotype, "Control"),
    Sex = fct_relevel(Sex, "Female")
  )
htp_groups |> count(Karyotype) # 52 D21 + 230 T21
htp_groups |> count(Karyotype, Sex)
htp_groups |> count(Karyotype, Sample_source_code)

# Subselect counts data and rpkms to match meta_data
# (ensures same samples and order as in meta_data in wide format)
htp_counts_sub <- htp_counts_long |> 
  pivot_wider(id_cols = EnsemblID, names_from = VisitID, values_from = raw_count) |>
  select(EnsemblID, htp_groups |> pull(VisitID))
#

# Check raw read count distributions across samples (not normalized) 
htp_counts_meta |>
  filter(VisitID %in% htp_groups$VisitID) |>
  mutate(VisitID = fct_relevel(VisitID, htp_groups |> arrange(Karyotype) |> pull(VisitID) %>%  as.character)) |> 
  ggplot(aes(VisitID, log2(raw_count + 0.1), color = Karyotype)) + # CUSTOMIZE AS NEEDED
  geom_sina(size = 0.01) +
  geom_boxplot(notch = FALSE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.2, color = "black", fill = "transparent", size = 0.75) +
  scale_color_manual(values = standard_colors) +
  labs(title = "HTP: Raw read count distributions across samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_HTP_counts_unfiltered.png")), width = 30, height = 5, units = "in")


## 2.3 Filter by minimum counts per million ----
htp_before <- htp_counts_sub %>% 
  transmute(
    EnsemblID, 
    row_sum = rowSums(select(., -EnsemblID))
  ) |> 
  filter(row_sum > 0) |> 
  nrow()
htp_cpm_data <- htp_counts_sub |> 
  column_to_rownames("EnsemblID") |> 
  cpm()
htp_keep <- htp_cpm_data |> 
  as_tibble(rownames = "EnsemblID") |> 
  pivot_longer(-EnsemblID, names_to = "VisitID", values_to = "cpm") |> 
  mutate(cpm > min_cpm) |> # check against min_cpm
  filter(`cpm > min_cpm` == TRUE) |> # and filter
  dplyr::count(EnsemblID) |> # count samples remaining per EnsemblID
  filter(n >= min_samples) # filter against min_samples
htp_counts_filtered <- htp_counts_sub |> 
  filter(EnsemblID %in% htp_keep$EnsemblID)

# Summarize rows before and after filtering 
cat("Total number of rows: ", htp_counts_sub |> nrow())
cat("Number of rows with non-zero total read counts before filtering: ", htp_before)
cat("Number of rows after filtering: ", htp_counts_filtered |> nrow())
#
# Check filtered read count distributions for each sample (not normalized) 
htp_counts_filtered |>
  pivot_longer(-EnsemblID, names_to = "VisitID", values_to = "raw_count") |>
  inner_join(htp_groups) |>
  mutate(VisitID = fct_relevel(VisitID, htp_groups |> arrange(Karyotype) |> pull(VisitID) %>% as.character)) |>
  ggplot(aes(VisitID, log2(raw_count + 0.1), color = Karyotype)) + # CUSTOMIZE AS NEEDED
  geom_sina(size = 0.01) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.2, color = "black", fill = "transparent", size = 0.75) +
  scale_color_manual(values = standard_colors) +
  labs(title = "HTP: CPM-filtered read count distributions across samples", x = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_HTP_counts_filtered.png")), width = 30, height = 5, units = "in")
#

## 2.4. Generate DESeqDataSet object(s) and run DESeq2 analysis ----
htp_formula <- as.formula(paste0("~ ", "Sex + Age_years_at_visit + Sample_source_code +", "Karyotype"))
#
dds_htp <- DESeqDataSetFromMatrix(
  countData = htp_counts_filtered |> 
    select(EnsemblID, htp_groups |> pull(VisitID)) |>  # ensures correct order of columns
    column_to_rownames("EnsemblID"), # must be converted to data frame from tibble
  colData = htp_groups,
  design = htp_formula
)
# Check meta data read in to DESeqDataSet (multivariable version):
colData(dds_htp)
#

# Run DESeq2 analysis 
dds_htp <- DESeq(dds_htp)
#

## 2.5. Get DESeq2 results ---- 
### 2.5.1 Define comparisons of interest ----
htp_comparisons <- list( 
  c("Karyotype", "T21", "Control") 
)
#
### 2.5.2 Results summaries ----
dds_htp |> get_results_sum(htp_comparisons[[1]], show_ind_filt_off=FALSE)
#
### 2.5.3 Assemble DESeq2 results table(s) ----
# Initialize empty vector to store names of results objects for later reference
res_SexAgeSource_T21_vs_Control <- dds_htp |> 
  get_results_tbl(
    contrast = htp_comparisons[[1]],
    shrink_type = "apeglm"
  )
#

res_SexAgeSource_T21_vs_Control

# volcano plot of results:
res_SexAgeSource_T21_vs_Control %>%
  volcano_plot_lab(
    title = "HTP ages 12-40: T21 vs. Control",
    subtitle = paste0(
      paste("Model:", paste(htp_formula, collapse = ""), "\n"), # CUSTOMIZE formula
      "[Down: ", (.) %>% filter(padj < 0.1 & FoldChange_adj <1) |> nrow(), "; Up: ", (.) |> filter(padj < 0.1 & FoldChange_adj >1) |> nrow(), "]"
    ),
    labels = TRUE,
    n_labels = 3,
    raster = TRUE
  ) +
  expand_limits(x = c(-1, 1) * 4)
ggsave(filename = here("plots", paste0(out_file_prefix, "P4C_12_40_T21_vs_Control", "_Volcano.pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#

## 2.6. GSEA Hallmarks analysis ----
#
### 2.6.1 Generate ranks ----
ranks_htp <- res_SexAgeSource_T21_vs_Control |> 
  filter(!is.na(log2FoldChange_adj)) |>  # need to remove NA rows that will break enplot
  select(ID = Gene_name, t = log2FoldChange_adj) |> 
  arrange(-abs(t)) |> 
  distinct(ID, .keep_all = TRUE) |>  # to avoid duplicates
  tibble::deframe() # convert to named num vector
#
### 2.6.2 Run fgsea ----
set.seed(123) # set this for reproducibility of results
hallmarks_htp <- run_fgsea2(geneset = hallmarks, ranks = ranks_htp, weighted = FALSE)
hallmarks_htp |> filter(padj < 0.1) 
#

# 3. Compare TOFA and HTP results ------

## 3.1.  FC comparison --------

# Create combined results table for FC
log2fc_cutoff <- 0 
fdr_cutoff <- 0.1 

htp_tofa_comb_FC <- res_SexAgeSource_T21_vs_Control |> 
  select(-log2FoldChange) |> # to be able to rename log2FoldChange_adj as log2FoldChange
  select(Geneid, Gene_name, chr, log2FoldChange = log2FoldChange_adj, pvalue, qvalue = padj) |> 
  inner_join(
    res_multi_BatchStudyID_16.week_vs_Baseline |> 
      select(-log2FoldChange) |> # to be able to rename log2FoldChange_adj as log2FoldChange
      select(Geneid, Gene_name, chr, log2FoldChange = log2FoldChange_adj, pvalue, qvalue = padj),
    by = c("Geneid", "Gene_name", "chr"),
    suffix = c("_T21", "_TOFA")
  ) |> 
  mutate(
    effect_relationship = case_when(
      is.na(log2FoldChange_T21) | is.na(log2FoldChange_TOFA) ~ "Incomplete",
      log2FoldChange_T21 > log2fc_cutoff & log2FoldChange_TOFA < log2fc_cutoff ~ "Opposite: T21 up, TOFA down",
      log2FoldChange_T21 < log2fc_cutoff & log2FoldChange_TOFA > log2fc_cutoff ~ "Opposite: T21 down, TOFA up",
      log2FoldChange_T21 > log2fc_cutoff & log2FoldChange_TOFA > log2fc_cutoff ~ "Same: both up",
      log2FoldChange_T21 < log2fc_cutoff & log2FoldChange_TOFA < log2fc_cutoff ~ "Same: both down",
      TRUE ~ "Near zero"
    ),
    both_significant = !is.na(qvalue_T21) & qvalue_T21 < fdr_cutoff & !is.na(qvalue_TOFA) & qvalue_TOFA < fdr_cutoff,
    T21_significant_only = !is.na(qvalue_T21) & qvalue_T21 < fdr_cutoff & (is.na(qvalue_TOFA) | qvalue_TOFA >= fdr_cutoff),
    TOFA_significant_only = !is.na(qvalue_TOFA) & qvalue_TOFA < fdr_cutoff & (is.na(qvalue_T21) | qvalue_T21 >= fdr_cutoff),
  ) |> 
  mutate(
    sig_TOFA_effect_on_T21  = case_when(
      both_significant == TRUE & str_detect(effect_relationship, "Opposite") ~ "reversed",
      both_significant == TRUE & str_detect(effect_relationship, "Same") ~ "exacerbated",
      both_significant == FALSE & qvalue_T21 < fdr_cutoff ~ "unchanged",
      TRUE ~ "unaffected by T21" # ignoring effects other than T21
    ),
    # reversal score captures strength of exacerbation/reversal: sign is inverse of product 
    # and magnitude is minimum abs(NES) between T21 and TOFA such that
    # high positive values => strong reversal; high neg => strong exacerbation
    reversal_score = -sign(log2FoldChange_T21 * log2FoldChange_TOFA) * 
      pmin(abs(log2FoldChange_T21), abs(log2FoldChange_TOFA))
  )
#

htp_tofa_comb_FC |> count(sig_TOFA_effect_on_T21)
# 1 exacerbated              291
# 2 reversed                 486
# 3 unaffected by T21       8022
# 4 unchanged               7903

### ranked plot ---------
htp_tofa_comb_FC |> 
  filter(sig_TOFA_effect_on_T21 != "unaffected by T21") |>  
  arrange(-log2FoldChange_T21) |> 
  mutate(rank = 1:length(Geneid)) |> 
  ggplot(aes(rank, log2FoldChange_T21)) + 
  geom_point_rast(data = . %>% filter(sig_TOFA_effect_on_T21 == "unchanged"),
                  aes(y = log2FoldChange_TOFA, color = sig_TOFA_effect_on_T21), alpha = 0.5, raster.dpi = 600) +
  geom_hline(yintercept = 0, color = "grey") + 
  geom_point_rast(data = . %>% filter(sig_TOFA_effect_on_T21 != "unchanged"),
                  aes(y = log2FoldChange_TOFA, color = sig_TOFA_effect_on_T21), alpha = 0.75, raster.dpi = 600) +
  geom_point_rast(aes(y = log2FoldChange_T21, color = "T21 DEGs"), alpha = 0.5, raster.dpi = 600) + 
  scale_x_continuous(breaks = c(1, 4340, 8680), expand = 0.010, labels = c("0%", "50%", "100%")) +
  scale_color_manual(values = c("reversed" = "#d73027", "exacerbated"= "#4575b4", "unchanged" = "grey90", "T21 DEGs" = "grey30")) + 
  labs(
    title = "Transcriptome: Effect of Tofa on T21 DEGs",
    subtitle = paste0("threshold: FDR ", fdr_cutoff, "; log2FC ", log2fc_cutoff),
    y = "log2FC"
  ) + 
  annotate(
    "text",
    x = Inf, y = Inf,
    label = htp_tofa_comb_FC |>  
      filter(sig_TOFA_effect_on_T21 != "unaffected by T21") |> 
      count(sig_TOFA_effect_on_T21) |> 
      pivot_wider(names_from = sig_TOFA_effect_on_T21, values_from = n, values_fill = 0) %>% 
      { paste0("Reversed: ", .$reversed, "\nExacerbated: ", .$exacerbated, "\nUnchanged: ", .$unchanged) },
    hjust = 1, vjust = 1,
    size = 4,
    color = "black"
  ) + 
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        aspect.ratio = 0.55)
ggsave(filename = here("plots", paste0(out_file_prefix, "T21_ranked_vs_tofa16wk_effects", ".pdf")), device = cairo_pdf, width = 8, height = 4, units = "in")
#

### barplot of all up/down combinations ----
htp_tofa_comb_FC |> 
  filter(both_significant == TRUE) |> 
  count(effect_relationship, sig_TOFA_effect_on_T21) |> 
  mutate(effect_relationship = str_remove(effect_relationship, "^[^:]*: ")) |> 
  mutate(effect_relationship = fct_relevel(
    effect_relationship, c("T21 up, TOFA down", "T21 down, TOFA up", "both up", "both down") |> rev())) |> 
  ggplot(aes(n, effect_relationship, fill = sig_TOFA_effect_on_T21)) + 
  geom_col(width = 0.5) +
  geom_text(aes(label = n), hjust = -0.15, size = 3) +
  labs(
    title = "Gene direction",
    subtitle = paste0("threshold: FDR ", fdr_cutoff, "; log2FC ", log2fc_cutoff),
    x = "Genes",
    y = NULL
  ) +
  scale_fill_manual(values = c("reversed" = "#d73027", "exacerbated"= "#4575b4", "unchanged" = "grey90", "T21 DEGs" = "grey30")) + 
  coord_cartesian(xlim = c(0, 320)) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        aspect.ratio = 1.5)
ggsave(filename = here("plots", paste0(out_file_prefix, "T21_vs_Tofa_gene_direction_counts", ".pdf")), device = cairo_pdf, width = 6, height = 6, units = "in")
#

## 3.2 GSEA comparison ----------
htp_tofa_comb_gsea <- hallmarks_htp |> 
  select(pathway, pvalue = pval, qvalue = padj, ES, NES, size, leadingEdge) |> 
  inner_join(
    hallmarks_tofa |> 
      select(pathway, pvalue = pval, qvalue = padj, ES, NES, size, leadingEdge),
    by = join_by(pathway),
    suffix = c("_T21", "_TOFA")
  ) |> 
  mutate(
    reversal_score = -sign(NES_T21 * NES_TOFA) * pmin(abs(NES_T21), abs(NES_TOFA)), # same logic as FC above
    combined_significance = case_when(
      qvalue_T21 < 0.1 & qvalue_TOFA < 0.1 ~ "both",
      qvalue_T21 >= 0.1 & qvalue_TOFA < 0.1 ~ "T21 only",
      qvalue_T21 < 0.1 & qvalue_TOFA >= 0.1 ~ "Tofa only",
      TRUE ~ "neither"
    ),
    reversal_class = case_when(
      is.na(NES_T21) | is.na(NES_TOFA) ~ "Incomplete",
      NES_T21 > 0 & NES_TOFA < 0 ~ "Reversal",
      NES_T21 < 0 & NES_TOFA > 0 ~ "Reversal",
      NES_T21 > 0 & NES_TOFA > 0 ~ "Same direction",
      NES_T21 < 0 & NES_TOFA < 0 ~ "Same direction",
      TRUE ~ "Near zero"
    )
  ) |> 
  mutate(pathway = str_remove(pathway, "HALLMARK_") |> str_replace_all("_", " ") |> str_to_title()) 

### pathway reversal bar plot ------
htp_tofa_comb_gsea |> 
  filter(combined_significance == "both") |> # important to do in addition to reversal score
  slice_max(abs(reversal_score), n =10) |> 
  arrange(reversal_score) |> 
  mutate(pathway = fct_inorder(pathway)) |> 
  ggplot(aes(x = reversal_score, y = pathway, fill = reversal_class)) +
  geom_vline(xintercept = 0, color = "grey", linewidth = 0.5) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c("Reversal" = "#d73027", "Same direction" = "#4575b4")) +
  labs(
    title = "Transcriptome: Pathway reversal rank",
    subtitle = "sig in both; top 10 pathways by reversal score",
    x = "Reversal score",
    y = NULL
  ) + 
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        aspect.ratio = 1.7)
ggsave(filename = here("plots", paste0(out_file_prefix, "T21_vs_Tofa_GSEA_reversal_bar", ".pdf")), device = cairo_pdf, width = 8, height = 4.8, units = "in")
#

### pathway scatter ---------
htp_tofa_comb_gsea |>
  ggplot(aes(x = NES_T21, y = NES_TOFA, color = reversal_class)) +
  geom_hline(yintercept = 0, color = "grey", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "grey", linewidth = 0.5) +
  geom_point(alpha = 0.8) +
  geom_text_repel(
    data = . %>% filter(combined_significance == "both") %>% 
      slice_max(abs(reversal_score), n =6),
    aes(label = pathway),
    size = 3,
    min.segment.length = 0,
    max.overlaps = 50,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("Reversal" = "#d73027", "Same direction" = "#4575b4")) +
  labs(
    title = "Transcritome: Hallmarks",
    x = "T21 NES",
    y = "TOFA NES",
    color = NULL
  ) + 
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        aspect.ratio = 1)
ggsave(filename = here("plots", paste0(out_file_prefix, "T21_vs_Tofa_GSEA_scatter", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#

### log2FC scatter for top pathways ---------

# get the list of top reversed pathways
htp_tofa_comb_gsea |> 
  filter(combined_significance == "both") |>  # important to do in addition to reversal score
  slice_max(abs(reversal_score), n =10) |> 
  pull(pathway)

# create hallmark genes df
top_reversal_hallmark_list <- list(
  HALLMARK_HEME_METABOLISM = hallmarks$HALLMARK_HEME_METABOLISM,
  HALLMARK_INTERFERON_GAMMA_RESPONSE = hallmarks$HALLMARK_INTERFERON_GAMMA_RESPONSE,
  HALLMARK_INFLAMMATORY_RESPONSE = hallmarks$HALLMARK_INFLAMMATORY_RESPONSE,
  HALLMARK_OXIDATIVE_PHOSPHORYLATION = hallmarks$HALLMARK_OXIDATIVE_PHOSPHORYLATION,
  HALLMARK_INTERFERON_ALPHA_RESPONSE = hallmarks$HALLMARK_INTERFERON_ALPHA_RESPONSE,
  HALLMARK_IL2_STAT5_SIGNALING = hallmarks$HALLMARK_IL2_STAT5_SIGNALING,
  HALLMARK_ALLOGRAFT_REJECTION = hallmarks$HALLMARK_ALLOGRAFT_REJECTION,
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = hallmarks$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,
  HALLMARK_MYC_TARGETS_V1 = hallmarks$HALLMARK_MYC_TARGETS_V1,
  HALLMARK_IL6_JAK_STAT3_SIGNALING = hallmarks$HALLMARK_IL6_JAK_STAT3_SIGNALING
)

top_reversal_hallmark_df <-enframe(top_reversal_hallmark_list, 
                                   name = "pathway", 
                                   value = "Gene_name") |>
  unnest(Gene_name)

top_reversal_hallmark_FC_df <- htp_tofa_comb_FC |> 
  inner_join(top_reversal_hallmark_df, by = "Gene_name") |> 
  select(Gene_name, pathway, contains("log2FoldChange"), reversal_score) |> 
  pivot_longer(cols = contains("log2FoldChange"), names_to = "comparison", values_to = "log2FC") |> 
  mutate(
    comparison = str_remove(comparison, "log2FoldChange_") |> str_to_upper(),
    pathway = str_remove(pathway, "HALLMARK_") |> str_replace_all("_", " ") |>  str_to_title()
  ) |> 
  mutate(pathway = fct_relevel(pathway, (htp_tofa_comb_gsea |> 
                                           filter(combined_significance == "both") |> # important to do in addition to reversal score
                                           slice_max(abs(reversal_score), n =10) |> 
                                           mutate(pathway = str_remove(pathway, "HALLMARK_") |> str_replace_all("_", " ") |> str_to_title()) |> 
                                           pull(pathway))))
#

htp_tofa_comb_FC |> 
  inner_join(top_reversal_hallmark_df, by = "Gene_name") |> 
  mutate(
    pathway = str_remove(pathway, "HALLMARK_") |> str_replace_all("_", " ") |> str_to_title(),
    pathway = fct_relevel(pathway, (htp_tofa_comb_gsea |> 
                                      filter(combined_significance == "both") |> # important to do in addition to reversal score
                                      slice_max(abs(reversal_score), n =10) |> 
                                      mutate(pathway = str_remove(pathway, "HALLMARK_") |> str_replace_all("_", " ") |> str_to_title()) |> 
                                      pull(pathway))),
    sig_TOFA_effect_on_T21 = fct_relevel(sig_TOFA_effect_on_T21, c("reversed", "exacerbated", "unchanged", "unaffected by T21")) 
  ) |> 
  ggplot(aes(log2FoldChange_T21, log2FoldChange_TOFA)) + 
  geom_point(aes(color = sig_TOFA_effect_on_T21), alpha = 0.75) + 
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) + 
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) + 
  theme(aspect.ratio = 1,
        legend.position = "bottom") +
  scale_color_manual(values = c("reversed" = "#d73027", "exacerbated"= "#4575b4", "unchanged" = "grey90", "unaffected by T21" = "grey90")) + 
  facet_wrap(~pathway, nrow = 1, scales = "free") +
  labs(title = "Transcriptome FCs: Top hallmarks")
ggsave(filename = here("plots", paste0(out_file_prefix, "T21_vs_Tofa_GSEA_FC_scatter", ".pdf")), device = cairo_pdf, width = 50, height = 3.2, units = "in", limitsize = F)
#

## 3.3 HTP + TOFA sina plot of example genes -----

### Combined batch effect removal ----
rpkms_combined <- bind_rows(
  tofa_rpkms_meta |> 
    mutate(Karyotype = "T21", batch = "TOFA", Sample_source_code = "TOFA") |> 
    select(ParticipantID, VisitID, Karyotype, Sex, Age = Age_years_at_visit, Event_Name, batch, Sample_source_code, EnsemblID, Gene_name, RPKM),
  htp_rpkms_meta |> 
    mutate(batch = "HTP", Event_Name = "Baseline") |> 
    select(ParticipantID, VisitID, Karyotype, Sex, Age = Age_years_at_visit, Event_Name, batch, Sample_source_code, EnsemblID, Gene_name, RPKM = Value)
) |> 
  # filter to genes passing cpm filter:
  filter(EnsemblID %in% tofa_keep$EnsemblID) |> 
  # combined group and batch:
  mutate(
    batch_combined = paste0(Sample_source_code, "_", batch),
    group = paste0(Karyotype, "_", Event_Name)
  )
rpkms_combined |> distinct(batch_combined)
rpkms_combined |> distinct(group)
rpkms_combined |> distinct(group, batch_combined)
#
sample_data_combined <- rpkms_combined |> 
  select(-c(RPKM))
unadj_data_combined <- rpkms_combined |> 
  select(VisitID, EnsemblID, RPKM) |> 
  mutate(RPKM = log2(RPKM)) |> # need to log2 transform for batch correction
  pivot_wider(names_from = VisitID, values_from = RPKM) |> 
  column_to_rownames(var = "EnsemblID")
#
rpkm_combined_BatchSource_adj <- unadj_data_combined |> 
  limma::removeBatchEffect(
    # adjust for batch_combined:
    batch = unadj_data_combined |> colnames() |> enframe(name = NULL, value = "VisitID") |> 
      inner_join(sample_data_combined |> select(VisitID, batch_combined) |> distinct()) |> pull(batch_combined),
    design = unadj_data_combined |> colnames() |> enframe(name = NULL, value = "VisitID") |> 
      inner_join(sample_data_combined |> select(VisitID, group) |> distinct()) %>% 
      model.matrix(~ group, data = .)
  ) |> 
  as_tibble(rownames = "EnsemblID") |> # convert back to tibble
  pivot_longer(-EnsemblID, names_to = "VisitID", values_to = "RPKM_adj") |> 
  mutate(RPKM_adj = 2^RPKM_adj) |> # remove log2 transformation
  inner_join(sample_data_combined)
#

### most reversed genes -----
reversed_examples <- htp_tofa_comb_FC |> 
  filter(both_significant == TRUE) |> 
  slice_max(reversal_score, n = 20) |> 
  arrange(-log2FoldChange_T21) |> 
  pull(Gene_name)
#
# create annotation df for ggsignif:
reversed_examples_signif_df <- rpkm_combined_BatchSource_adj |> 
  filter(Gene_name %in% reversed_examples) |> 
  # calculate y position per Gene_name and group
  mutate(extreme = rstatix::is_extreme((RPKM_adj)), .by = c(Gene_name, batch, group)) |> 
  filter(extreme != TRUE) |> 
  summarise(y_position = max(RPKM_adj, na.rm = TRUE) * 1.01, 
            .by = c(Gene_name, batch)) |> 
  mutate(
    xmin = if_else(batch == "HTP", "HTP:Control_Baseline", "TOFA:T21_Baseline"),
    xmax = if_else(batch == "HTP", "HTP:T21_Baseline", "TOFA:T21_16 week")
  ) |> 
  # add qvalue
  inner_join(htp_tofa_comb_FC |> 
               filter(Gene_name %in% reversed_examples) |> 
               select(Gene_name, contains("qvalue")) |> 
               pivot_longer(cols = -Gene_name, names_to = "batch", values_to = "qvalue") |> 
               mutate(batch = if_else(str_detect(batch, "T21"), "HTP", "TOFA"))) |> 
  mutate(annotations = if_else(
    qvalue < 0.01,
    scales::label_scientific(digits = 2)(qvalue),
    scales::number(qvalue, accuracy = 0.01)
  )) |> 
  mutate(Gene_name = fct_relevel(Gene_name, reversed_examples)) # control plotting order
#
reversed_examples_sina <- rpkm_combined_BatchSource_adj |> 
  filter(Gene_name %in% reversed_examples) |> 
  mutate(Gene_name = fct_relevel(Gene_name, reversed_examples)) |> # control plotting order
  mutate(group = paste0(batch, ":", group) %>%  
           fct_relevel(., c("HTP:Control_Baseline", "HTP:T21_Baseline", 
                            "TOFA:T21_Baseline", "TOFA:T21_16 week"))) |> 
  group_by(Gene_name, group) |> 
  mutate(extreme = rstatix::is_extreme((RPKM_adj))) |> 
  filter(extreme != TRUE) |> 
  ungroup() |> 
  ggplot(aes(group, RPKM_adj, color = group)) + # CUSTOMIZE
  ggrastr::rasterise(geom_sina(maxwidth = 0.5), dpi= 600) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.5, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~ Gene_name, scales = "free", nrow = 1) +
  scale_color_manual(values = c("TOFA:T21_Baseline" = "#999999", "TOFA:T21_16 week" = "#6baed6", "HTP:Control_Baseline" = "grey30", "HTP:T21_Baseline" = "#009b4e")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    aspect.ratio = 1.2,
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5)
  ) +
  labs(
    title = "Transcriptome: Top reversed by TOFA",
    subtitle = "batch+source adj.; Extreme outliers removed",
    x = NULL
  )  + 
  ggsignif::geom_signif(
    data = reversed_examples_signif_df,
    aes(xmin = xmin, xmax = xmax,
        annotations = annotations,
        y_position = y_position),
    inherit.aes = FALSE, manual = TRUE, tip_length = 0, color = "black", textsize = 3) 
ggsave(reversed_examples_sina, filename = here("plots", paste0(out_file_prefix, "tofa_reversed_examples_sina", ".pdf")), device = cairo_pdf, width = 80, height = 4, units = "in", limitsize = F)
#

### most exacerbated genes -----
exacerbated_examples <- htp_tofa_comb_FC |> 
  filter(both_significant == TRUE) |> 
  slice_min(reversal_score, n = 20) |> 
  arrange(-log2FoldChange_T21) |> 
  pull(Gene_name)

# create annotation df for ggsignif:
exacerbated_examples_signif_df <- rpkm_combined_BatchSource_adj |> 
  filter(Gene_name %in% exacerbated_examples) |> 
  # calculate y position per Gene_name and group
  mutate(extreme = rstatix::is_extreme((RPKM_adj)), .by = c(Gene_name, batch, group)) |> 
  filter(extreme != TRUE) |> 
  summarise(y_position = max(RPKM_adj, na.rm = TRUE) * 1.01, 
            .by = c(Gene_name, batch)) |> 
  mutate(
    xmin = if_else(batch == "HTP", "HTP:Control_Baseline", "TOFA:T21_Baseline"),
    xmax = if_else(batch == "HTP", "HTP:T21_Baseline", "TOFA:T21_16 week")
  ) |> 
  # add qvalue
  inner_join(htp_tofa_comb_FC |> 
               filter(Gene_name %in% exacerbated_examples) |> 
               select(Gene_name, contains("qvalue")) |> 
               pivot_longer(cols = -Gene_name, names_to = "batch", values_to = "qvalue") |> 
               mutate(batch = if_else(str_detect(batch, "t21"), "HTP", "TOFA"))) |> 
  mutate(annotations = if_else(
    qvalue < 0.01,
    scales::label_scientific(digits = 2)(qvalue),
    scales::number(qvalue, accuracy = 0.01)
  )) |> 
  mutate(Gene_name = fct_relevel(Gene_name, exacerbated_examples)) # control plotting order
#
exacerbated_examples_sina <- rpkm_combined_BatchSource_adj |> 
  filter(Gene_name %in% exacerbated_examples) |> 
  mutate(Gene_name = fct_relevel(Gene_name, exacerbated_examples)) |> # control plotting order
  mutate(group = paste0(batch, ":", group) %>%  
           fct_relevel(., c("HTP:Control_Baseline", "HTP:T21_Baseline", 
                            "TOFA:T21_Baseline", "TOFA:T21_16 week"))) |> 
  group_by(Gene_name, group) |> 
  mutate(extreme = rstatix::is_extreme((RPKM_adj))) |> 
  filter(extreme != TRUE) |> 
  ungroup() |> 
  ggplot(aes(group, RPKM_adj, color = group)) + # CUSTOMIZE
  ggrastr::rasterise(geom_sina(maxwidth = 0.5), dpi= 600) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.5, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~ Gene_name, scales = "free", nrow = 1) +
  scale_color_manual(values = c("TOFA:T21_Baseline" = "#999999", "TOFA:T21_16 week" = "#6baed6", "HTP:Control_Baseline" = "grey30", "HTP:T21_Baseline" = "#009b4e")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    aspect.ratio = 1.2,
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5)
  ) +
  labs(
    title = "Transcriptome: Top exacerbated by TOFA",
    subtitle = "batch+source adj.; Extreme outliers removed",
    x = NULL
  )  + 
  ggsignif::geom_signif(
    data = exacerbated_examples_signif_df,
    aes(xmin = xmin, xmax = xmax,
        annotations = annotations,
        y_position = y_position),
    inherit.aes = FALSE, manual = TRUE, tip_length = 0, color = "black", textsize = 3) 
ggsave(reversed_examples_sina, filename = here("plots", paste0(out_file_prefix, "tofa_exacerbated_examples_sina", ".pdf")), device = cairo_pdf, width = 80, height = 4, units = "in", limitsize = F)
#

### most neutral genes --------
neutral_examples <- htp_tofa_comb_FC |> 
  filter(both_significant == TRUE | T21_significant_only == TRUE) |> 
  # select some chr21 and non-chr21 examples
  mutate(chr21 = if_else(chr == "chr21", "chr21", "non-chr21")) |> 
  slice_min(abs(reversal_score), n = 2, by = chr21) |> 
  arrange(-log2FoldChange_T21) |> 
  pull(Gene_name)
# add interferon receptors:
neutral_examples <- c(neutral_examples, "IFNGR2", "IL10RB")
#

# create annotation df for ggsignif:
neutral_examples_signif_df <- rpkm_combined_BatchSource_adj |> 
  filter(Gene_name %in% neutral_examples) |> 
  # calculate y position per Gene_name and group
  mutate(extreme = rstatix::is_extreme((RPKM_adj)), .by = c(Gene_name, batch, group)) |> 
  filter(extreme != TRUE) |> 
  summarise(y_position = max(RPKM_adj, na.rm = TRUE) * 1.01, 
            .by = c(Gene_name, batch)) |> 
  mutate(
    xmin = if_else(batch == "HTP", "HTP:Control_Baseline", "TOFA:T21_Baseline"),
    xmax = if_else(batch == "HTP", "HTP:T21_Baseline", "TOFA:T21_16 week")
  ) |> 
  # add qvalue
  inner_join(htp_tofa_comb_FC |> 
               filter(Gene_name %in% neutral_examples) |> 
               select(Gene_name, contains("qvalue")) |> 
               pivot_longer(cols = -Gene_name, names_to = "batch", values_to = "qvalue") |> 
               mutate(batch = if_else(str_detect(batch, "t21"), "HTP", "TOFA"))) |> 
  mutate(annotations = if_else(
    qvalue < 0.01,
    scales::label_scientific(digits = 2)(qvalue),
    scales::number(qvalue, accuracy = 0.01)
  )) |> 
  mutate(Gene_name = fct_relevel(Gene_name, neutral_examples)) # control plotting order
#
neutral_examples_sina <- rpkm_combined_BatchSource_adj |> 
  filter(Gene_name %in% neutral_examples) |> 
  mutate(Gene_name = fct_relevel(Gene_name, neutral_examples)) |> # control plotting order
  mutate(group = paste0(batch, ":", group) %>%  
           fct_relevel(., c("HTP:Control_Baseline", "HTP:T21_Baseline", 
                            "TOFA:T21_Baseline", "TOFA:T21_16 week"))) |> 
  group_by(Gene_name, group) |> 
  mutate(extreme = rstatix::is_extreme((RPKM_adj))) |> 
  filter(extreme != TRUE) |> 
  ungroup() |> 
  ggplot(aes(group, RPKM_adj, color = group)) + # CUSTOMIZE
  ggrastr::rasterise(geom_sina(maxwidth = 0.5), dpi= 600) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.5, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~ Gene_name, scales = "free", nrow = 1) +
  scale_color_manual(values = c("TOFA:T21_Baseline" = "#999999", "TOFA:T21_16 week" = "#6baed6", "HTP:Control_Baseline" = "grey30", "HTP:T21_Baseline" = "#009b4e")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    aspect.ratio = 1.2,
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5)
  ) +
  labs(
    title = "Transcriptome: Most neutral by TOFA",
    subtitle = "batch+source adj.; Extreme outliers removed",
    x = NULL
  )  + 
  ggsignif::geom_signif(
    data = neutral_examples_signif_df,
    aes(xmin = xmin, xmax = xmax,
        annotations = annotations,
        y_position = y_position),
    inherit.aes = FALSE, manual = TRUE, tip_length = 0, color = "black", textsize = 3) 
ggsave(reversed_examples_sina, filename = here("plots", paste0(out_file_prefix, "tofa_neutral_examples_sina", ".pdf")), device = cairo_pdf, width = 80, height = 4, units = "in", limitsize = F)
#

## 3.4 Export results -----
# combine all results into a single excel file
list(
  "Transcriptome TOFA 16wk" = res_multi_BatchStudyID_16.week_vs_Baseline |>
    select(-c(FoldChange, log2FoldChange)) |> # remove non-lfcShrink() FCs
    transmute(Gene_name, chr, Geneid, comparison = "16 week vs Baseline", FoldChange = FoldChange_adj, 
              log2FoldChange = log2FoldChange_adj, pvalue, qvalue = padj, Model),
  "Transcriptome T21 vs D21" = res_SexAgeSource_T21_vs_Control |>
    select(-c(FoldChange, log2FoldChange)) |> # remove non-lfcShrink() FCs
    transmute(Gene_name, chr, Geneid, comparison = "T21 vs D21", FoldChange = FoldChange_adj, 
              log2FoldChange = log2FoldChange_adj, pvalue, qvalue = padj, Model),
  "Transcriptome T21 vs TOFA 16wk" = htp_tofa_comb_FC,
  "Transcriptome T21 vs TOFA GSEA" = htp_tofa_comb_gsea
) |> 
  export_excel(filename = "results")

########## End of Script #############

##### Save workspace    ----
################################################################################
save.image(file = here("rdata", paste0(out_file_prefix, ".RData")), compress = TRUE, safe = TRUE) # saves entire workspace (can be slow)

# Reload the workspace
# load(here("rdata", paste0(out_file_prefix, ".RData")))

# session_info ----
date()
sessionInfo()
################################################
