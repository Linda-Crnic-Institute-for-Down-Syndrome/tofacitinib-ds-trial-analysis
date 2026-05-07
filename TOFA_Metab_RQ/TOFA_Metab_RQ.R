################################################
# Title: Tofacitinib for Immune Skin Conditions in Down Syndrome: Analysis of Relative Quantitation (RQ) plasma metabolite data
# Author(s):
#   - Neetha Paul Eduthan
#   - Matthew Galbraith
# affiliation(s):
#   - Linda Crnic Institute for Down syndrome
#   - University of Colorado Anschutz
################################################

### Summary:  
# Mixed-effects linear regression analysis of targeted plasma metabolite changes from baseline to 
# 16 weeks of TOFA treatment, and linear regression analysis of baseline T21 samples from the 
# trial and D21 samples from the HTP data set, followed by comparison of treatment and T21 effects.
# See README.md for more details
# 

### Data type(s):
# Clinical trial (TOFA) datasets:
#    A. Participant-level metadata; Available on request
#    B. Visit/Event-level metadata; Available on request
#    C. TOFA + HTP D21s metabolite data; DOI: 10.5281/zenodo.20043706
#      
# Human Trisome Project (HTP) datasets:
#    D. Participant-level metadata; DOI: 10.5281/zenodo.19962380
#    E. Visit/Event-level metadata; DOI: 10.5281/zenodo.19962380

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
# HTP data files
htp_participant_metadata_file <- here("data", "HTP_Participant_metadata_zenodo_v1.txt") # Source: 10.5281/zenodo.19962380
htp_visit_metadata_file <- here("data", "HTP_Visit_metadata_zenodo_v1.txt") # Source: 10.5281/zenodo.19962380
# HTP D21 + TOFA metabolite data file:
metab_RQ_data_file <- here("data", "TOFA_Metab_RQ_D21_baseline_16wk_data_zenodo_v1.txt.gz") # Source: 10.5281/zenodo.20043706
#
standard_colors <- c("Control" = "grey30", "T21" = "#009b4e", "Baseline" = "#999999", "16 week" = "#6baed6")
out_file_prefix <- "TOFA_Metab_RQ.R_" # should match this script title
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

#### 1.1.3 Read Metabolite data ----
metab_RQ_data <- metab_RQ_data_file |> 
  read_tsv() |>
  select(Analyte = FeatureID, CmpdID = Feature_cmpdID, Pathway = Feature_pathway, everything()) |> 
  mutate(log2_rel_abundance = log2(Value), .after = VisitID) # log-transform Value
#
metab_RQ_data |> distinct(ParticipantID) # 84
metab_RQ_data |> distinct(VisitID) # 125
metab_RQ_data |> distinct(Analyte) # 129

#### 1.1.4. Join Metabolite data with TOFA metadata ------
# Excluding Endpoint eligible = FALSE samples
tofa_metab_RQ_data_meta <- metab_RQ_data |> 
  inner_join(tofa_participant_meta_data) |> 
  filter(Endpoint_eligible == TRUE) |> 
  inner_join(tofa_visit_meta_data) |> 
  arrange(ParticipantID, Event_Name) |> 
  select(ParticipantID, VisitID, Event_Name, Analyte, CmpdID, Pathway, log2_rel_abundance, everything())
#

tofa_metab_RQ_data_meta |> distinct(VisitID, Event_Name) |> count(Event_Name)
# Event_Name     n
# Baseline      41
# 16 week       39

## 1.2 TOFA Mixed-effects linear regression -----
# ParticipantID as mixed effect and age at baseline and Sex as fixed effects

tofa_regressions_dat <- tofa_metab_RQ_data_meta |>
  select(ParticipantID, VisitID, Age_years_at_visit, Sex, Event_Name, Analyte, CmpdID, Pathway, log2_rel_abundance) |> 
  # add a new column with age_at_baseline
  inner_join(
    tofa_metab_RQ_data_meta |> 
      filter(Event_Name == "Baseline") |> 
      distinct(ParticipantID, Age_years_at_visit) |> 
      select(ParticipantID, age_at_baseline = Age_years_at_visit)
  ) |> 
 nest(data = -c(Analyte, CmpdID, Pathway))

## mixed effect model ~ Event_Name + Sex + age_at_baseline + (1|ParticipantID)  
tofa_regressions_multi_fixedSexAge_mixedStudyID <- tofa_regressions_dat |> 
  mutate(
    fit = map(data, ~ lmerTest::lmer(log2_rel_abundance ~ Event_Name + Sex + age_at_baseline + (1|ParticipantID), REML = FALSE, data = .x)),
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
  select(Analyte, CmpdID, Pathway, n_observations = nobs, log2FoldChange = estimate, conf.low, conf.high, statistic, pvalue = p.value) |> 
  mutate(qvalue = p.adjust(pvalue, method = "BH")) |> 
  mutate(
    comparison = "16 week vs Baseline",
    .after = Pathway) |> 
  mutate(model = "~Event_Name+(1|ParticipantID)+Age+Sex")

tofa_lm_multi_fixedSexAge_mixedStudyID_res |> filter(qvalue < 0.1) # 12 sig.

# volcano plot
tofa_lm_multi_fixedSexAge_mixedStudyID_res %>% 
  volcano_plot_lab_lm(
    title="Diff. abund. in 16 wks vs. baseline", 
    subtitle = paste0("~Event_Name+(1|StudyID)+Sex+Age","\n[Up: ", 
                      (.) |> filter(qvalue < 0.1 & log2FoldChange >0) |> nrow(), "; Down: ", 
                      (.) |> filter(qvalue < 0.1 & log2FoldChange <0) |> nrow(), "]", "\ndown                      up")
  )

# 2 HTP T21 vs D21 analysis -----

## 2.1. Read in HTP data ------

#### 2.1.1 Read in HTP metadata ------
htp_metadata <- htp_participant_metadata_file |> 
  read_tsv() |> 
  inner_join(htp_visit_metadata_file |> 
               read_tsv()) 

#### 2.1.3 Join Metabolite data with HTP metadata ------
# Subset to ages between 12-40 to match with age of TOFA cohort
htp_metab_RQ_data_meta <- metab_RQ_data |> 
  inner_join(htp_metadata) 

htp_metab_RQ_data_meta |>  distinct(VisitID, Karyotype) |>  count(Karyotype)
# Karyotype     n
# Control      40

## 2.3 HTP linear regression -----
T21vsD21_regressions_dat <- htp_metab_RQ_data_meta |> 
  select(ParticipantID, VisitID, Age_years_at_visit, Sex, Karyotype, Sample_source_code, Analyte, CmpdID, Pathway, log2_rel_abundance) |> 
  # add TOFA baseline samples
  bind_rows(
    tofa_metab_RQ_data_meta |> 
      filter(Event_Name == "Baseline") |> 
      mutate(Karyotype = "T21",
             Sample_source_code = "TOFA") |> 
      select(ParticipantID, VisitID, Age_years_at_visit, Sex, Karyotype, Sample_source_code, Analyte, CmpdID, Pathway, log2_rel_abundance)
  ) |> 
  # remove extreme outliers
  mutate(extreme = rstatix::is_extreme(log2_rel_abundance), .by = c(Analyte, Karyotype)) |> 
  filter(extreme == FALSE) |> 
  nest(data = -c(Analyte, CmpdID, Pathway))

# NOTE: Sample_source_code (TOFA vs HTP) is fully confounded with comparison groups
# (T21 vs D21), and therefore cannot be included as a separate covariate in the model.
T21vsD21_regressions_multi_fixedSexAge <- T21vsD21_regressions_dat |> 
  mutate(
    fit = map(data, ~ lm(log2_rel_abundance ~ Karyotype + Age_years_at_visit + Sex, data = .x)),
    tidied = map(fit, broom::tidy), # see ?tidy.lm
    glanced = map(fit, broom::glance), # see ?glance.lm
    augmented = map(fit, broom::augment), # see ?augment.lm
    vifs = map(fit, ~car::vif(mod = .x) |> as_tibble(rownames = "term")) # check co-linearity of variables
  )

## 2.4 Extract results for HTP T21 vs D21 -----
T21vsD21_multi_fixedSexAge_res <- T21vsD21_regressions_multi_fixedSexAge |> 
  unnest(tidied) |> 
  filter(str_detect(term, "Karyotype")) |>
  select(Analyte, CmpdID, Pathway, log2FoldChange = estimate, pvalue = p.value) |> 
  mutate(qvalue = p.adjust(pvalue, method = "BH")) |> 
  mutate(
    comparison = "T21 vs D21",
    .after = Pathway) |> 
  mutate(model = "~Karyotype+Sex+Age")

T21vsD21_multi_fixedSexAge_res |>  filter(qvalue < 0.1) # 55
#

# volcano plot
T21vsD21_multi_fixedSexAge_res %>% 
  volcano_plot_lab_lm(
    title="Diff. abund. in T21 vs. D21", 
    subtitle = paste0("~Karyotype+Sex+Age","\n[Up: ", 
                      (.) |> filter(qvalue < 0.1 & log2FoldChange >0) |> nrow(), "; Down: ", 
                      (.) |> filter(qvalue < 0.1 & log2FoldChange <0) |> nrow(), "]", "\ndown                      up")
  )


# 3. Compare TOFA and HTP results ------

## 3.1.  FC comparison --------

# Create combined results table for FC
log2fc_cutoff <- 0 
fdr_cutoff <- 0.1 

T21vsD21_tofa_comb_FC <- T21vsD21_multi_fixedSexAge_res |> 
  select(Analyte, log2FoldChange, pvalue, qvalue) |> 
  inner_join(
    tofa_lm_multi_fixedSexAge_mixedStudyID_res |> 
      select(Analyte, log2FoldChange, pvalue, qvalue),
    by = c("Analyte"),
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

T21vsD21_tofa_comb_FC |> count(sig_TOFA_effect_on_T21)
# sig_TOFA_effect_on_T21     n
# exacerbated                3
# reversed                   2
# unaffected by T21         74
# unchanged                 50

### ranked plot ---------
T21vsD21_tofa_comb_FC |> 
  filter(sig_TOFA_effect_on_T21 != "unaffected by T21") |>  
  arrange(-log2FoldChange_T21) |> 
  mutate(rank = 1:length(Analyte)) |> 
  ggplot(aes(rank, log2FoldChange_T21)) + 
  geom_point_rast(data = . %>% filter(sig_TOFA_effect_on_T21 == "unchanged"),
                  aes(y = log2FoldChange_TOFA, color = sig_TOFA_effect_on_T21), alpha = 0.5, raster.dpi = 600) +
  geom_hline(yintercept = 0, color = "grey") + 
  geom_point_rast(data = . %>% filter(sig_TOFA_effect_on_T21 != "unchanged"),
                  aes(y = log2FoldChange_TOFA, color = sig_TOFA_effect_on_T21), alpha = 0.75, raster.dpi = 600) +
  geom_point_rast(aes(y = log2FoldChange_T21, color = "T21 DEGs"), alpha = 0.5, raster.dpi = 600) + 
  scale_x_continuous(breaks = c(1, 28, 55), expand = 0.010, labels = c("0%", "50%", "100%")) +
  scale_color_manual(values = c("reversed" = "#d73027", "exacerbated"= "#4575b4", "unchanged" = "grey90", "T21 DEGs" = "grey30")) + 
  labs(
    title = "Metabolites: Effect of Tofa on T21 DAMs",
    subtitle = paste0("threshold: FDR ", fdr_cutoff, "; log2FC ", log2fc_cutoff),
    y = "log2FC"
  ) + 
  annotate(
    "text",
    x = Inf, y = Inf,
    label = T21vsD21_tofa_comb_FC |> 
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


## 3.2 HTP + TOFA sina plot of example metabolites -----

### Combined dataset ----
# No batch-effect correction (limma::removeBatchEffect) was applied prior to visualization,
# as all samples were generated in a single batch. Sample source is confounded with
# T21 vs D21 group status and cannot be adjusted without removing the biological signal.
metab_combined <- htp_metab_RQ_data_meta |>
  mutate(group = "D21") |>
  bind_rows(
    tofa_metab_RQ_data_meta |>
      mutate(group = Event_Name)
  ) |> 
  mutate(group = fct_relevel(group, c("D21", "Baseline", "16 week")))        

metab_combined |> distinct(group)
metab_combined |> distinct(VisitID, group) |> count(group)
#

### reversed metabolites -----
reversed_examples <- T21vsD21_tofa_comb_FC |> 
  filter(sig_TOFA_effect_on_T21 == "reversed") |> 
  arrange(-log2FoldChange_T21) |> 
  pull(Analyte)
#

# create annotation df for ggsignif:
reversed_examples_signif_df <- metab_combined |>
  filter(Analyte %in% reversed_examples) |> 
  mutate(extreme = rstatix::is_extreme(log2_rel_abundance), .by = c(Analyte, group)) |>
  filter(extreme != TRUE) |>
  summarise(y_position = max(log2_rel_abundance, na.rm = TRUE), .by = c(Analyte, group)) %>%
  {
    df <- .
    bind_rows(
      df |>
        filter(group %in% c("D21", "Baseline")) |>
        summarise(
          y_position = max(y_position) * 1.01,
          .by = Analyte
        ) |>
        mutate(xmin = "D21", xmax = "Baseline"),
      df |>
        filter(group %in% c("Baseline", "16 week")) |>
        summarise(
          y_position = max(y_position) * 1.01,
          .by = Analyte
        ) |>
        mutate(xmin = "Baseline", xmax = "16 week")
    )
  } |> 
  mutate(comparison = if_else(xmin == "D21", "T21vsD21", "TOFA")) |> 
  # add qvalue
  inner_join(T21vsD21_tofa_comb_FC |> 
               filter(Analyte %in% reversed_examples) |> 
               select(Analyte, contains("qvalue")) |> 
               pivot_longer(cols = -Analyte, names_to = "comparison", values_to = "qvalue") |> 
               mutate(comparison = if_else(str_detect(comparison, "T21"), "T21vsD21", "TOFA"))) |> 
  mutate(annotations = if_else(
    qvalue < 0.01,
    scales::label_scientific(digits = 2)(qvalue),
    scales::number(qvalue, accuracy = 0.01)
  )) |> 
  mutate(Analyte = fct_relevel(Analyte, reversed_examples)) # control plotting order
#
reversed_examples_sina <- metab_combined |>
  filter(Analyte %in% reversed_examples) |> 
  mutate(Analyte = fct_relevel(Analyte, reversed_examples)) |> # control plotting order
  group_by(Analyte, group) |>
  mutate(extreme = rstatix::is_extreme((log2_rel_abundance))) |>
  filter(extreme != TRUE) |>
  ungroup() |>
  ggplot(aes(group, log2_rel_abundance, color = group)) + # CUSTOMIZE
  ggrastr::rasterise(geom_sina(maxwidth = 0.5), dpi= 600) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.5, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~Analyte, scales = "free", nrow = 1) +
  scale_color_manual(values = c("Baseline" = "#999999", "16 week" = "#6baed6", "D21" = "grey30")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    aspect.ratio = 1.2,
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5)
  ) +
  labs(
    title = "Metabolites: Reversed by TOFA",
    subtitle = "unadjusted data; Extreme outliers removed",
    x = NULL
  )  + 
  ggsignif::geom_signif(
    data = reversed_examples_signif_df,
    aes(xmin = xmin, xmax = xmax,
        annotations = annotations,
        y_position = y_position),
    inherit.aes = FALSE, manual = TRUE, tip_length = 0, color = "black", textsize = 3) 
ggsave(reversed_examples_sina, filename = here("plots", paste0(out_file_prefix, "tofa_reversed_examples_sina", ".pdf")), device = cairo_pdf, width = 20, height = 4, units = "in", limitsize = F)
#

### exacerbated metabolites -----
exacerbated_examples <- T21vsD21_tofa_comb_FC |> 
  filter(sig_TOFA_effect_on_T21 == "exacerbated") |> 
  arrange(-log2FoldChange_T21) |> 
  pull(Analyte)
#
# create annotation df for ggsignif:
exacerbated_examples_signif_df <- metab_combined |>
  filter(Analyte %in% exacerbated_examples) |> 
  mutate(extreme = rstatix::is_extreme(log2_rel_abundance), .by = c(Analyte, group)) |>
  filter(extreme != TRUE) |>
  summarise(y_position = max(log2_rel_abundance, na.rm = TRUE), .by = c(Analyte, group)) %>%
  {
    df <- .
    bind_rows(
      df |>
        filter(group %in% c("D21", "Baseline")) |>
        summarise(
          y_position = max(y_position) * 1.01,
          .by = Analyte
        ) |>
        mutate(xmin = "D21", xmax = "Baseline"),
      df |>
        filter(group %in% c("Baseline", "16 week")) |>
        summarise(
          y_position = max(y_position) * 1.01,
          .by = Analyte
        ) |>
        mutate(xmin = "Baseline", xmax = "16 week")
    )
  } |> 
  mutate(comparison = if_else(xmin == "D21", "T21vsD21", "TOFA")) |> 
  # add qvalue
  inner_join(T21vsD21_tofa_comb_FC |> 
               filter(Analyte %in% exacerbated_examples) |> 
               select(Analyte, contains("qvalue")) |> 
               pivot_longer(cols = -Analyte, names_to = "comparison", values_to = "qvalue") |> 
               mutate(comparison = if_else(str_detect(comparison, "T21"), "T21vsD21", "TOFA"))) |> 
  mutate(annotations = if_else(
    qvalue < 0.01,
    scales::label_scientific(digits = 2)(qvalue),
    scales::number(qvalue, accuracy = 0.01)
  )) |> 
  mutate(Analyte = fct_relevel(Analyte, exacerbated_examples)) # control plotting order
#
exacerbated_examples_sina <- metab_combined |>
  filter(Analyte %in% exacerbated_examples) |> 
  mutate(Analyte = fct_relevel(Analyte, exacerbated_examples)) |> # control plotting order
  group_by(Analyte, group) |>
  mutate(extreme = rstatix::is_extreme((log2_rel_abundance))) |>
  filter(extreme != TRUE) |>
  ungroup() |>
  ggplot(aes(group, log2_rel_abundance, color = group)) + # CUSTOMIZE
  ggrastr::rasterise(geom_sina(maxwidth = 0.5), dpi= 600) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.5, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~Analyte, scales = "free", nrow = 1) +
  scale_color_manual(values = c("Baseline" = "#999999", "16 week" = "#6baed6", "D21" = "grey30")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    aspect.ratio = 1.2,
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5)
  ) +
  labs(
    title = "Metabolites: Exacerbated by TOFA",
    subtitle = "unadjusted data; Extreme outliers removed",
    x = NULL
  )  + 
  ggsignif::geom_signif(
    data = exacerbated_examples_signif_df,
    aes(xmin = xmin, xmax = xmax,
        annotations = annotations,
        y_position = y_position),
    inherit.aes = FALSE, manual = TRUE, tip_length = 0, color = "black", textsize = 3) 
ggsave(exacerbated_examples_sina, filename = here("plots", paste0(out_file_prefix, "tofa_exacerbated_examples_sina", ".pdf")), device = cairo_pdf, width = 20, height = 4, units = "in", limitsize = F)
#

## 3.3 Export results -----
# combine all results into a single excel file
list(
  "Metabolites TOFA 16wk" = tofa_lm_multi_fixedSexAge_mixedStudyID_res,
  "Metabolites T21 vs D21" = T21vsD21_multi_fixedSexAge_res,
  "Metabolites T21 vs TOFA 16wk" = T21vsD21_tofa_comb_FC
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
