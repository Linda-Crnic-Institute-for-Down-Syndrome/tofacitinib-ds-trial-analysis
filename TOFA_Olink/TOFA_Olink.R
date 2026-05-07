################################################
# Title: Tofacitinib for Immune Skin Conditions in Down Syndrome: Analysis of Olink plasma proteins 
# Author(s):
#   - Neetha Paul Eduthan
#   - Matthew Galbraith
# affiliation(s):
#   - Linda Crnic Institute for Down syndrome
#   - University of Colorado Anschutz
################################################

### Summary:  
# Mixed-effects linear regression analysis to assess differential abundance of plasma proteins
# measured by the Olink assay at baseline and 16 weeks of Tofacitinib (TOFA) treatment,
# and linear regression analysis of differential protein abundance between individuals with (T21)
# and without (D21) Down syndrome in Human Trisome Project (HTP) data set, followed by comparison 
# of treatment and T21 effects.
# See README.md for more details
# 

### Data type(s):
# Clinical trial (TOFA) datasets:
#    A. Participant-level metadata; Available on request.
#    B. Visit/Event-level metadata; Available on request.
#    C. TOFA Olink data; DOI: 10.5281/zenodo.19962923
#      
# Human Trisome Project (HTP) datasets:
#    D. Participant-level metadata; DOI: 10.5281/zenodo.19962380
#    E. Visit/Event-level metadata; DOI: 10.5281/zenodo.19962380
#    F. HTP Olink data; DOI: 10.5281/zenodo.20046326

# 0 General Setup -----
# RUN THIS FIRST TIME - Initialize and install packages with renv:
# renv::init(bioconductor = TRUE)
#
# To install the exact versions of all R packages base on renv.lock file (requires matching R version):
# renv::restore()

## 0.1 Load required libraries ----
library("readxl") # used to read .xlsx files
library("openxlsx") # used for data export as Excel workbooks
library("tidyverse") # data manipulation, visualization, and tidy workflows
library("skimr") # data table summaries
library("rstatix") # for outlier detection 
library("fgsea") # to run GSEA 
library("limma")  # linear modeling for high-dimensional omics data
library("broom") # extract tidy model results
library("ggplot2") # to create plots
library("ggrepel") # required for labelling features
library("ggforce") # required for zooming and sina
library("ggrastr") # for rasterizing plots
library("conflicted") # force all conflicts to become errors
library("lme4") # for mixed models
library("lmerTest") # for mixed models
library("broom.mixed") # for mixed models
conflicts_prefer( # declare preferences in cases of conflict
  dplyr::filter,
  dplyr::count,
  dplyr::rename,
  dplyr::bind_rows
)
library("here") # generates path to current project directory
#

## 0.2 Set file name parameters ----
# TOFA files:
tofa_participant_meta_data_file <- here("data", "TOFA_Participant_metadata_zenodo_v1.txt") # Source: Available on request
tofa_visit_meta_data_file <- here("data", "TOFA_Visit_metadata_zenodo_v1.txt") # Source: Available on request
tofa_olink_data_file <- here("data", "TOFA_Olink_baseline_16wk_data_zenodo_v1.txt.gz") # Source: 10.5281/zenodo.19962923
# HTP data files
htp_participant_metadata_file <- here("data", "HTP_Participant_metadata_zenodo_v1.txt") # Source: 10.5281/zenodo.19962380
htp_visit_metadata_file <- here("data", "HTP_Visit_metadata_zenodo_v1.txt") # Source: 10.5281/zenodo.19962380
htp_olink_data_file <- here("data", "HTP_Olink_data_zenodo_v1.txt.gz") # Source: 10.5281/zenodo.20046326
#
standard_colors <- c("Control" = "grey30", "T21" = "#009b4e", "Baseline" = "#999999", "16 week" = "#6baed6")
out_file_prefix <- "TOFA_Olink.R_" # should match this script title
# End required parameters ###
source(here("helper_functions.R")) # load helper functions
#

# 1 TOFA 16 weeks vs. baseline analysis --------
## 1.1 Read in TOFA data ----
#### 1.1.1 Participant level meta data ----
tofa_participant_meta_data <- tofa_participant_meta_data_file |> 
  read_tsv() |> 
  mutate(Sex = fct_relevel(Sex, "Female"))
#
#### 1.1.2 Event/Visit level meta data ----
tofa_visit_meta_data <- tofa_visit_meta_data_file |> 
  read_tsv() |> 
  # keep only timepoints included in this analysis
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  mutate(
    Event_Name = fct_relevel(Event_Name, c("Baseline", "16 week"))
  ) 
#

#### 1.1.3 Read Olink proteomics data ----
tofa_olink_data <- tofa_olink_data_file |> 
  read_tsv() |>
  select(OlinkID = FeatureID, Assay = Feature_name, UniProt = Feature_UniProt_id, NPX = Value, everything())
#
tofa_olink_data |> distinct(ParticipantID) # 43
tofa_olink_data |> distinct(VisitID) # 84
tofa_olink_data |> distinct(Assay) # 5,401

#### 1.1.4. Join TOFA Olink data with metadata ------
# Excluding Endpoint eligible = FALSE samples
tofa_olink_data_meta <- tofa_olink_data |> 
  inner_join(tofa_participant_meta_data) |> 
  filter(Endpoint_eligible == TRUE) |> 
  inner_join(tofa_visit_meta_data) |> 
  arrange(ParticipantID, Event_Name) |> 
  select(ParticipantID, VisitID, Event_Name, OlinkID, UniProt, Assay, NPX, everything())
#

tofa_olink_data_meta |> distinct(VisitID, Event_Name) |> count(Event_Name)
# Event_Name     n
# Baseline      41
# 16 week       39

## 1.2 TOFA Mixed-effects linear regression -----
# ParticipantID as mixed effect and age at baseline and Sex as fixed effects

tofa_regressions_dat <- tofa_olink_data_meta |>
  # add a new column with age_at_baseline
  inner_join(
    tofa_olink_data_meta |> 
      filter(Event_Name == "Baseline") |> 
      distinct(ParticipantID, Age_years_at_visit) |> 
      select(ParticipantID, age_at_baseline = Age_years_at_visit)
  ) |> 
 nest(data = -c(OlinkID, UniProt, Assay))

## mixed effect model ~ Event_Name + Sex + age_at_baseline + (1|ParticipantID)  
tofa_regressions_multi_fixedSexAge_mixedStudyID <- tofa_regressions_dat |> 
  mutate(
    fit = map(data, ~ lmerTest::lmer(NPX ~ Event_Name + Sex + age_at_baseline + (1|ParticipantID), REML = FALSE, data = .x)),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?tidy.lm
    glanced = map(fit, broom.mixed::glance), # see ?glance.lm
    augmented = map(fit, broom.mixed::augment), # see ?augment.lm
    vifs = map(fit, ~car::vif(mod = .x) |> as_tibble(rownames = "term")) # check co-linearity of variables
  )
#

## 1.3 Extract results for TOFA 16weeks vs Baseline -----
tofa_lm_multi_fixedSexAge_mixedStudyID_res <- tofa_regressions_multi_fixedSexAge_mixedStudyID |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, "Event_Name")) |>
  select(OlinkID, Assay, UniProt, n_observations = nobs, log2FoldChange = estimate, conf.low, conf.high, statistic, pvalue = p.value) |> 
  mutate(qvalue = p.adjust(pvalue, method = "BH")) |> 
  mutate(
    comparison = "16 week vs Baseline",
    .after = UniProt) |> 
  mutate(model = "~Event_Name+(1|ParticipantID)+Age+Sex")

tofa_lm_multi_fixedSexAge_mixedStudyID_res |> filter(qvalue < 0.1) # 386 sig.

# volcano plot
tofa_lm_multi_fixedSexAge_mixedStudyID_res %>% 
  volcano_plot_lab_lm(
    title="Diff. abund. in 16 wks vs. baseline", 
    subtitle = paste0("~Event_Name+(1|StudyID)+Sex+Age","\n[Up: ", 
                      (.) |> filter(qvalue < 0.1 & log2FoldChange >0) |> nrow(), "; Down: ", 
                      (.) |> filter(qvalue < 0.1 & log2FoldChange <0) |> nrow(), "]", "\ndown                      up")
  )

## 1.4 GSEA -----
# Download Human GSEA hallmarks from MSigDB: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# Read in Human GSEA hallmarks
hallmarks <- here("data/GSEA/human/h.all.v7.4.symbols.gmt") %>% 
  fgsea::gmtPathways(gmt.file = .)
#

## Generate ranks 
ranks_tofa <- tofa_lm_multi_fixedSexAge_mixedStudyID_res |>
  mutate(t = -log10(pvalue) * log2FoldChange) |>
  select(ID = Assay, t) |>
  arrange(-t) |>
  tibble::deframe() # convert to named num vector
#
## Run fgsea 
set.seed(123) # for reproducibility
hallmarks_tofa <- run_fgsea2(geneset = hallmarks, ranks = ranks_tofa, weighted = FALSE)

hallmarks_tofa |> filter(padj < 0.1) # 8 sig. hallmarks

# 2 HTP T21 vs D21 analysis -----

## 2.1. Read in HTP data ------

#### 2.1.1 Read in HTP metadata ------
htp_metadata <- htp_participant_metadata_file |> 
  read_tsv() |> 
  inner_join(htp_visit_metadata_file |> 
               read_tsv()) 

#### 2.1.2 Read in HTP Olink data ------
htp_olink_data <- htp_olink_data_file |> 
  read_tsv() |>
  select(OlinkID = FeatureID, Assay = Feature_name, UniProt = Feature_UniProt_id, NPX = Value, everything())
#
htp_olink_data |> distinct(ParticipantID) # 935
htp_olink_data |> distinct(VisitID) # 935
htp_olink_data |> distinct(Assay) # 5,401

#### 2.1.3 Join HTP Olink data with metadata ------
# Subset to ages between 12-40 to match with age of TOFA cohort
htp_olink_data_meta <- htp_olink_data |> 
  inner_join(htp_metadata) |> 
  filter(dplyr::between(Age_years_at_visit, 12, 40)) 

htp_olink_data_meta |>  distinct(VisitID, Karyotype) |>  count(Karyotype)
# Karyotype     n
# Control     178
# T21         462

## 2.2 HTP linear regression -----
htp_regressions_dat <- htp_olink_data_meta |> 
  # remove extreme outliers
  mutate(extreme = rstatix::is_extreme(NPX), .by = c(Assay, Karyotype)) |> 
  filter(extreme == FALSE) |> 
  nest(data = -c(OlinkID, UniProt, Assay))

htp_regressions_multi_fixedSexAgeSource <- htp_regressions_dat |> 
  mutate(
    fit = map(data, ~ lm(NPX ~ Karyotype + Age_years_at_visit + Sex + Sample_source_code, data = .x)),
    tidied = map(fit, broom::tidy), # see ?tidy.lm
    glanced = map(fit, broom::glance), # see ?glance.lm
    augmented = map(fit, broom::augment), # see ?augment.lm
    vifs = map(fit, ~car::vif(mod = .x) |> as_tibble(rownames = "term")) # check co-linearity of variables
  )

## 2.3 Extract results for HTP T21 vs D21 -----
htp_multi_fixedSexAgeSource_res <- htp_regressions_multi_fixedSexAgeSource |> 
  unnest(tidied) |> 
  filter(str_detect(term, "Karyotype")) |>
  select(OlinkID, Assay, UniProt, log2FoldChange = estimate, pvalue = p.value) |> 
  mutate(qvalue = p.adjust(pvalue, method = "BH")) |> 
  mutate(
    comparison = "T21 vs D21",
    .after = UniProt) |> 
  mutate(model = "~Karyotype+Sex+Age+Sample_source")

htp_multi_fixedSexAgeSource_res |>  filter(qvalue < 0.1) # 2,655
#

# volcano plot
htp_multi_fixedSexAgeSource_res %>% 
  volcano_plot_lab_lm(
    title="Diff. abund. in T21 vs. D21", 
    subtitle = paste0("~Karyotype+Sex+Age+Sample_source","\n[Up: ", 
                      (.) |> filter(qvalue < 0.1 & log2FoldChange >0) |> nrow(), "; Down: ", 
                      (.) |> filter(qvalue < 0.1 & log2FoldChange <0) |> nrow(), "]", "\ndown                      up")
  )

## 2.4 GSEA -----

# Generate ranks 
ranks_htp <- htp_multi_fixedSexAgeSource_res |> 
  mutate(t = -log10(pvalue) * log2FoldChange) |> 
  select(ID = Assay, t) |> 
  arrange(-t) |> 
  tibble::deframe() # convert to named num vector
#
## Run fgsea 
set.seed(123) # for reproducibility
hallmarks_htp <- run_fgsea2(geneset = hallmarks, ranks = ranks_htp, weighted = FALSE)

hallmarks_htp |>  filter(padj < 0.1) # 27 sig. hallmarks

# 3. Compare TOFA and HTP results ------

## 3.1.  FC comparison --------

# Create combined results table for FC
log2fc_cutoff <- 0 
fdr_cutoff <- 0.1 

htp_tofa_comb_FC <- htp_multi_fixedSexAgeSource_res |> 
  select(Assay, log2FoldChange, pvalue, qvalue) |> 
  inner_join(
    tofa_lm_multi_fixedSexAge_mixedStudyID_res |> 
      select(Assay, log2FoldChange, pvalue, qvalue),
    by = c("Assay"),
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
# sig_TOFA_effect_on_T21     n
# exacerbated               88
# reversed                 208
# unaffected by T21       2746
# unchanged               2359


### ranked plot ---------
htp_tofa_comb_FC |> 
  filter(sig_TOFA_effect_on_T21 != "unaffected by T21") |>  
  arrange(-log2FoldChange_T21) |> 
  mutate(rank = 1:length(Assay)) |> 
  ggplot(aes(rank, log2FoldChange_T21)) + 
  geom_point_rast(data = . %>% filter(sig_TOFA_effect_on_T21 == "unchanged"),
                  aes(y = log2FoldChange_TOFA, color = sig_TOFA_effect_on_T21), alpha = 0.5, raster.dpi = 600) +
  geom_hline(yintercept = 0, color = "grey") + 
  geom_point_rast(data = . %>% filter(sig_TOFA_effect_on_T21 != "unchanged"),
                  aes(y = log2FoldChange_TOFA, color = sig_TOFA_effect_on_T21), alpha = 0.75, raster.dpi = 600) +
  geom_point_rast(aes(y = log2FoldChange_T21, color = "T21 DEGs"), alpha = 0.5, raster.dpi = 600) + 
  scale_x_continuous(breaks = c(1, 1327, 2655), expand = 0.010, labels = c("0%", "50%", "100%")) +
  scale_color_manual(values = c("reversed" = "#d73027", "exacerbated"= "#4575b4", "unchanged" = "grey90", "T21 DEGs" = "grey30")) + 
  labs(
    title = "Olink: Effect of Tofa on T21 DAPs",
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
    title = "Protein direction",
    subtitle = paste0("threshold: FDR ", fdr_cutoff, "; log2FC ", log2fc_cutoff),
    x = "Proteins",
    y = NULL
  ) +
  scale_fill_manual(values = c("reversed" = "#d73027", "exacerbated"= "#4575b4", "unchanged" = "grey90", "T21 DEGs" = "grey30")) + 
  coord_cartesian(xlim = c(0, 155)) +
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        aspect.ratio = 1.5)
ggsave(filename = here("plots", paste0(out_file_prefix, "T21_vs_Tofa_protein_direction_counts", ".pdf")), device = cairo_pdf, width = 6, height = 6, units = "in")
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
    title = "Olink: Pathway reversal rank",
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
    title = "Olink: Hallmarks",
    x = "T21 NES",
    y = "TOFA NES",
    color = NULL
  ) + 
  theme(panel.border = element_blank(),
        axis.line = element_line(),
        aspect.ratio = 1)
ggsave(filename = here("plots", paste0(out_file_prefix, "T21_vs_Tofa_GSEA_scatter", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#

### log2FC scatter for top hallmarks ---------

# get the list of top reversed pathways
htp_tofa_comb_gsea |> 
  filter(combined_significance == "both") |>  # important to do in addition to reversal score
  slice_max(abs(reversal_score), n =10) |> 
  pull(pathway)

# create hallmark genes df
top_reversal_hallmark_list <- list(
  HALLMARK_INFLAMMATORY_RESPONSE = hallmarks$HALLMARK_INFLAMMATORY_RESPONSE,
  HALLMARK_INTERFERON_GAMMA_RESPONSE = hallmarks$HALLMARK_INTERFERON_GAMMA_RESPONSE,
  HALLMARK_IL6_JAK_STAT3_SIGNALING = hallmarks$HALLMARK_IL6_JAK_STAT3_SIGNALING,
  HALLMARK_ALLOGRAFT_REJECTION = hallmarks$HALLMARK_ALLOGRAFT_REJECTION,
  HALLMARK_IL2_STAT5_SIGNALING = hallmarks$HALLMARK_IL2_STAT5_SIGNALING,
  HALLMARK_CHOLESTEROL_HOMEOSTASIS = hallmarks$HALLMARK_CHOLESTEROL_HOMEOSTASIS
)
top_reversal_hallmark_df <-enframe(top_reversal_hallmark_list, 
                                   name = "pathway", 
                                   value = "Assay") |>
  unnest(Assay)

top_reversal_hallmark_FC_df <- htp_tofa_comb_FC |> 
  inner_join(top_reversal_hallmark_df, by = "Assay") |> 
  select(Assay, pathway, contains("log2FoldChange"), reversal_score) |> 
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
  inner_join(top_reversal_hallmark_df, by = "Assay") |> 
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
  labs(title = "Olink FCs: Top hallmarks")
ggsave(filename = here("plots", paste0(out_file_prefix, "T21_vs_Tofa_GSEA_FC_scatter", ".pdf")), device = cairo_pdf, width = 50, height = 3.2, units = "in", limitsize = F)
#

## 3.3 HTP + TOFA sina plot of example proteins -----

### Combined batch effect removal ----
olink_combined <- bind_rows(
  tofa_olink_data_meta |> 
    mutate(Karyotype = "T21", batch = "TOFA", Sample_source_code = "TOFA") |> 
    select(ParticipantID, VisitID, Karyotype, Sex, Age = Age_years_at_visit, Event_Name, batch, Sample_source_code, Assay, NPX),
  htp_olink_data_meta |> 
    mutate(batch = "HTP", Event_Name = "Baseline") |> 
    select(ParticipantID, VisitID, Karyotype, Sex, Age = Age_years_at_visit, Event_Name, batch, Sample_source_code, Assay, NPX)
) |> 
  # combined group and batch:
  mutate(
    batch_combined = paste0(Sample_source_code, "_", batch),
    group = paste0(Karyotype, "_", Event_Name)
  )
olink_combined |> distinct(batch_combined)
olink_combined |> distinct(group)
olink_combined |> distinct(group, batch_combined)
#
sample_data_combined <- olink_combined |> 
  select(-c(NPX))
unadj_data_combined <- olink_combined |> 
  select(VisitID, Assay, NPX) |> 
  pivot_wider(names_from = VisitID, values_from = NPX) |> # NPX is already log-transformed
  column_to_rownames(var = "Assay")
#
olink_combined_BatchSource_adj <- unadj_data_combined |> 
  limma::removeBatchEffect(
    # adjust for batch_combined:
    batch = unadj_data_combined |> colnames() |> enframe(name = NULL, value = "VisitID") |> 
      inner_join(sample_data_combined |> select(VisitID, batch_combined) |> distinct()) |> pull(batch_combined),
    design = unadj_data_combined |> colnames() |> enframe(name = NULL, value = "VisitID") |> 
      inner_join(sample_data_combined |> select(VisitID, group) |> distinct()) %>% 
      model.matrix(~ group, data = .)
  ) |> 
  as_tibble(rownames = "Assay") |> # convert back to tibble
  pivot_longer(-Assay, names_to = "VisitID", values_to = "NPX_adj") |> 
  inner_join(sample_data_combined)
#

### most reversed proteins -----
reversed_examples <- htp_tofa_comb_FC |> 
  filter(both_significant == TRUE) |> 
  slice_max(reversal_score, n = 20) |> 
  arrange(-log2FoldChange_T21) |> 
  pull(Assay)
#
# create annotation df for ggsignif:
reversed_examples_signif_df <- olink_combined_BatchSource_adj |> 
  filter(Assay %in% reversed_examples) |> 
  # calculate y position per Assay and group
  mutate(extreme = rstatix::is_extreme((NPX_adj)), .by = c(Assay, batch, group)) |> 
  filter(extreme != TRUE) |> 
  summarise(y_position = max(NPX_adj, na.rm = TRUE) + 0.3, 
            .by = c(Assay, batch)) |> 
  mutate(
    xmin = if_else(batch == "HTP", "HTP:Control_Baseline", "TOFA:T21_Baseline"),
    xmax = if_else(batch == "HTP", "HTP:T21_Baseline", "TOFA:T21_16 week")
  ) |> 
  # add qvalue
  inner_join(htp_tofa_comb_FC |> 
               filter(Assay %in% reversed_examples) |> 
               select(Assay, contains("qvalue")) |> 
               pivot_longer(cols = -Assay, names_to = "batch", values_to = "qvalue") |> 
               mutate(batch = if_else(str_detect(batch, "T21"), "HTP", "TOFA"))) |> 
  mutate(annotations = if_else(
    qvalue < 0.01,
    scales::label_scientific(digits = 2)(qvalue),
    scales::number(qvalue, accuracy = 0.01)
  )) |> 
  mutate(Assay = fct_relevel(Assay, reversed_examples)) # control plotting order
#
reversed_examples_sina <- olink_combined_BatchSource_adj |> 
  filter(Assay %in% reversed_examples) |> 
  mutate(Assay = fct_relevel(Assay, reversed_examples)) |> # control plotting order
  mutate(group = paste0(batch, ":", group) %>%  
           fct_relevel(., c("HTP:Control_Baseline", "HTP:T21_Baseline", 
                            "TOFA:T21_Baseline", "TOFA:T21_16 week"))) |> 
  group_by(Assay, group) |> 
  mutate(extreme = rstatix::is_extreme((NPX_adj))) |> 
  filter(extreme != TRUE) |> 
  ungroup() |> 
  ggplot(aes(group, NPX_adj, color = group)) + # CUSTOMIZE
  ggrastr::rasterise(geom_sina(maxwidth = 0.5), dpi= 600) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.5, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~ Assay, scales = "free", nrow = 1) +
  scale_color_manual(values = c("TOFA:T21_Baseline" = "#999999", "TOFA:T21_16 week" = "#6baed6", "HTP:Control_Baseline" = "grey30", "HTP:T21_Baseline" = "#009b4e")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    aspect.ratio = 1.2,
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5)
  ) +
  labs(
    title = "Olink: Top reversed by TOFA",
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

### most exacerbated proteins -----
exacerbated_examples <- htp_tofa_comb_FC |> 
  filter(both_significant == TRUE) |> 
  slice_min(reversal_score, n = 20) |> 
  arrange(-log2FoldChange_T21) |> 
  pull(Assay)
#
# create annotation df for ggsignif:
exacerbated_examples_signif_df <- olink_combined_BatchSource_adj |> 
  filter(Assay %in% exacerbated_examples) |> 
  # calculate y position per Assay and group
  mutate(extreme = rstatix::is_extreme((NPX_adj)), .by = c(Assay, batch, group)) |> 
  filter(extreme != TRUE) |> 
  summarise(y_position = max(NPX_adj, na.rm = TRUE) + 0.3, 
            .by = c(Assay, batch)) |> 
  mutate(
    xmin = if_else(batch == "HTP", "HTP:Control_Baseline", "TOFA:T21_Baseline"),
    xmax = if_else(batch == "HTP", "HTP:T21_Baseline", "TOFA:T21_16 week")
  ) |> 
  # add qvalue
  inner_join(htp_tofa_comb_FC |> 
               filter(Assay %in% exacerbated_examples) |> 
               select(Assay, contains("qvalue")) |> 
               pivot_longer(cols = -Assay, names_to = "batch", values_to = "qvalue") |> 
               mutate(batch = if_else(str_detect(batch, "T21"), "HTP", "TOFA"))) |> 
  mutate(annotations = if_else(
    qvalue < 0.01,
    scales::label_scientific(digits = 2)(qvalue),
    scales::number(qvalue, accuracy = 0.01)
  )) |> 
  mutate(Assay = fct_relevel(Assay, exacerbated_examples)) # control plotting order
#
exacerbated_examples_sina <- olink_combined_BatchSource_adj |> 
  filter(Assay %in% exacerbated_examples) |> 
  mutate(Assay = fct_relevel(Assay, exacerbated_examples)) |> # control plotting order
  mutate(group = paste0(batch, ":", group) %>% 
           fct_relevel(., c("HTP:Control_Baseline", "HTP:T21_Baseline", 
                            "TOFA:T21_Baseline", "TOFA:T21_16 week"))) |> 
  group_by(Assay, group) |> 
  mutate(extreme = rstatix::is_extreme((NPX_adj))) |> 
  filter(extreme != TRUE) |> 
  ungroup() |> 
  ggplot(aes(group, NPX_adj, color = group)) + # CUSTOMIZE
  ggrastr::rasterise(geom_sina(maxwidth = 0.5), dpi= 600) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.5, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~ Assay, scales = "free", nrow = 1) +
  scale_color_manual(values = c("TOFA:T21_Baseline" = "#999999", "TOFA:T21_16 week" = "#6baed6", "HTP:Control_Baseline" = "grey30", "HTP:T21_Baseline" = "#009b4e")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    aspect.ratio = 1.2,
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5)
  ) +
  labs(
    title = "Olink: Top exacerbated by TOFA",
    subtitle = "batch+source adj.; Extreme outliers removed",
    x = NULL
  )  + 
  ggsignif::geom_signif(
    data = exacerbated_examples_signif_df,
    aes(xmin = xmin, xmax = xmax,
        annotations = annotations,
        y_position = y_position),
    inherit.aes = FALSE, manual = TRUE, tip_length = 0, color = "black", textsize = 3) 
ggsave(exacerbated_examples_sina, filename = here("plots", paste0(out_file_prefix, "tofa_exacerbated_examples_sina", ".pdf")), device = cairo_pdf, width = 80, height = 4, units = "in", limitsize = F)
#

### most neutral proteins ------
neutral_examples <- htp_tofa_comb_FC |> 
  filter(both_significant == TRUE | T21_significant_only == TRUE) |> 
  slice_min(abs(reversal_score), n = 20) |> 
  arrange(-log2FoldChange_T21) |> 
  pull(Assay)
#

# create annotation df for ggsignif:
neutral_examples_signif_df <- olink_combined_BatchSource_adj |> 
  filter(Assay %in% neutral_examples) |> 
  # calculate y position per Assay and group
  mutate(extreme = rstatix::is_extreme((NPX_adj)), .by = c(Assay, batch, group)) |> 
  filter(extreme != TRUE) |> 
  summarise(y_position = max(NPX_adj, na.rm = TRUE) + 0.3, 
            .by = c(Assay, batch)) |> 
  mutate(
    xmin = if_else(batch == "HTP", "HTP:Control_Baseline", "TOFA:T21_Baseline"),
    xmax = if_else(batch == "HTP", "HTP:T21_Baseline", "TOFA:T21_16 week")
  ) |> 
  # add qvalue
  inner_join(htp_tofa_comb_FC |> 
               filter(Assay %in% neutral_examples) |> 
               select(Assay, contains("qvalue")) |> 
               pivot_longer(cols = -Assay, names_to = "batch", values_to = "qvalue") |> 
               mutate(batch = if_else(str_detect(batch, "T21"), "HTP", "TOFA"))) |> 
  mutate(annotations = if_else(
    qvalue < 0.01,
    scales::label_scientific(digits = 2)(qvalue),
    scales::number(qvalue, accuracy = 0.01)
  )) |> 
  mutate(Assay = fct_relevel(Assay, neutral_examples)) # control plotting order
#
neutral_examples_sina <- olink_combined_BatchSource_adj |> 
  filter(Assay %in% neutral_examples) |> 
  mutate(Assay = fct_relevel(Assay, neutral_examples)) |> # control plotting order
  mutate(group = paste0(batch, ":", group) %>%  
           fct_relevel(., c("HTP:Control_Baseline", "HTP:T21_Baseline", 
                            "TOFA:T21_Baseline", "TOFA:T21_16 week"))) |> 
  group_by(Assay, group) |> 
  mutate(extreme = rstatix::is_extreme((NPX_adj))) |> 
  filter(extreme != TRUE) |> 
  ungroup() |> 
  ggplot(aes(group, NPX_adj, color = group)) + # CUSTOMIZE
  ggrastr::rasterise(geom_sina(maxwidth = 0.5), dpi= 600) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.5, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~ Assay, scales = "free", nrow = 1) +
  scale_color_manual(values = c("TOFA:T21_Baseline" = "#999999", "TOFA:T21_16 week" = "#6baed6", "HTP:Control_Baseline" = "grey30", "HTP:T21_Baseline" = "#009b4e")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    aspect.ratio = 1.2,
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5)
  ) +
  labs(
    title = "Olink: Most neutral by TOFA",
    subtitle = "batch+source adj.; Extreme outliers removed",
    x = NULL
  )  + 
  ggsignif::geom_signif(
    data = neutral_examples_signif_df,
    aes(xmin = xmin, xmax = xmax,
        annotations = annotations,
        y_position = y_position),
    inherit.aes = FALSE, manual = TRUE, tip_length = 0, color = "black", textsize = 3) 
ggsave(neutral_examples_sina, filename = here("plots", paste0(out_file_prefix, "tofa_neutral_examples_sina", ".pdf")), device = cairo_pdf, width = 80, height = 4, units = "in", limitsize = F)
#

## 3.4 Export results -----
# combine all results into a single excel file
list(
  "Olink TOFA 16wk" = tofa_lm_multi_fixedSexAge_mixedStudyID_res,
  "Olink T21 vs D21" = htp_multi_fixedSexAgeSource_res,
  "Olink T21 vs TOFA 16wk" = htp_tofa_comb_FC,
  "Olink T21 vs TOFA GSEA" = htp_tofa_comb_gsea
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
