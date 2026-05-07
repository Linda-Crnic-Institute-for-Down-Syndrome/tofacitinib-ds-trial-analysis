################################################
# Title: Tofacitinib for Immune Skin Conditions in Down Syndrome: Analysis of NULISA plasma proteins
# Author(s):
#   - Neetha Paul Eduthan
#   - Matthew Galbraith
# affiliation(s):
#   - Linda Crnic Institute for Down syndrome
#   - University of Colorado Anschutz
################################################

### Summary:  
# Mixed-effects linear regression analysis to assess differential abundance of plasma proteins
# measured by the NULISA assay at baseline and 16 weeks of Tofacitinib (TOFA) treatment,
# and linear regression analysis of differential protein abundance between individuals with (T21)
# and without (D21) Down syndrome in Human Trisome Project (HTP) data set, followed by comparison 
# of treatment and T21 effects.
# See README.md for more details
# 

### Data type(s):
# Clinical trial (TOFA) datasets:
#    A. Participant-level metadata; Available on request.
#    B. Visit/Event-level metadata; Available on request.
#    C. TOFA NULISA data; DOI: 10.5281/zenodo.20043773
#      
# Human Trisome Project (HTP) datasets:
#    D. Participant-level metadata; DOI: 10.5281/zenodo.19962380
#    E. Visit/Event-level metadata; DOI: 10.5281/zenodo.19962380
#    F. HTP NULISA data; DOI: 10.5281/zenodo.20043943

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
tofa_nulisa_data_file <- here("data", "TOFA_NULISA_baseline_16wk_data_zenodo_v1.txt.gz") # Source: 10.5281/zenodo.20043773
# HTP data files
htp_participant_metadata_file <- here("data", "HTP_Participant_metadata_zenodo_v1.txt") # Source: 10.5281/zenodo.19962380
htp_visit_metadata_file <- here("data", "HTP_Visit_metadata_zenodo_v1.txt") # Source: 10.5281/zenodo.19962380
htp_nulisa_data_file <- here("data", "HTP_NULISA_data_zenodo_v1.txt.gz") # Source: 10.5281/zenodo.20043943
#
standard_colors <- c("Control" = "grey30", "T21" = "#009b4e", "Baseline" = "#999999", "16 week" = "#6baed6")
out_file_prefix <- "TOFA_NULISA.R_" # should match this script title
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

#### 1.1.3 Read NULISA data ----
tofa_nulisa_data <- tofa_nulisa_data_file |> 
  read_tsv() |>
  select(Feature = Feature_name, AlamarTargetID = FeatureID, Category = Feature_category, Value, everything())
#
tofa_nulisa_data |> distinct(ParticipantID) # 42
tofa_nulisa_data |> distinct(VisitID) # 81
tofa_nulisa_data |> distinct(Feature) # 130

#### 1.1.4. Join TOFA NULISA data with metadata ------
# Excluding Endpoint eligible = FALSE samples
tofa_nulisa_data_meta <- tofa_nulisa_data |> 
  inner_join(tofa_participant_meta_data) |> 
  filter(Endpoint_eligible == TRUE) |> 
  inner_join(tofa_visit_meta_data) |> 
  arrange(ParticipantID, Event_Name) |> 
  select(ParticipantID, VisitID, Event_Name, Feature, AlamarTargetID, Category, Value, everything())
#

tofa_nulisa_data_meta |> distinct(VisitID, Event_Name) |> count(Event_Name)
# Event_Name     n
# Baseline      41
# 16 week       38

## 1.2 Calculate amyloid and Tau ratios ------

# prepare df with ratios
tofa_nulisa_ratios_data_meta <- tofa_nulisa_data_meta |> 
  filter(Feature %in% c("pTau-217", "Aβ42", "Aβ40")) |> 
  select(ParticipantID, VisitID, Age_years_at_visit, Sex, Event_Name, Feature, Value) |> 
  pivot_wider(names_from = Feature, values_from = Value) |> 
  # calculate log2(ratio) as difference since values are in log-scale
  mutate(ptau217_Aβ42_ratio = `pTau-217`-Aβ42,
         Aβ42_Aβ40_ratio = `Aβ42`-`Aβ40`) |> 
  pivot_longer(cols = -c(ParticipantID, VisitID, Age_years_at_visit, Sex, Event_Name), names_to = "Feature", values_to = "Value") |> 
  filter(Feature %in% c("ptau217_Aβ42_ratio", "Aβ42_Aβ40_ratio")) |> 
  mutate(Category = "ratio")

## 1.3 TOFA Mixed-effects linear regression -----
# ParticipantID as mixed effect and age at baseline and Sex as fixed effects

tofa_regressions_dat <- tofa_nulisa_data_meta |>
  select(ParticipantID, VisitID, Age_years_at_visit, Sex, Event_Name, Feature, AlamarTargetID, Category, Value) |> 
  # add ratios
  bind_rows(
    tofa_nulisa_ratios_data_meta
  ) |> 
  # add a new column with age_at_baseline
  inner_join(
    tofa_nulisa_data_meta |> 
      filter(Event_Name == "Baseline") |> 
      distinct(ParticipantID, Age_years_at_visit) |> 
      select(ParticipantID, age_at_baseline = Age_years_at_visit)
  ) |> 
 nest(data = -c(Feature, AlamarTargetID, Category))

## mixed effect model ~ Event_Name + Sex + age_at_baseline + (1|ParticipantID)  
tofa_regressions_multi_fixedSexAge_mixedStudyID <- tofa_regressions_dat |> 
  mutate(
    fit = map(data, ~ lmerTest::lmer(Value ~ Event_Name + Sex + age_at_baseline + (1|ParticipantID), REML = FALSE, data = .x)),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?tidy.lm
    glanced = map(fit, broom.mixed::glance), # see ?glance.lm
    augmented = map(fit, broom.mixed::augment), # see ?augment.lm
    vifs = map(fit, ~car::vif(mod = .x) |> as_tibble(rownames = "term")) # check co-linearity of variables
  )
#

## 1.4 Extract results for TOFA 16weeks vs Baseline -----
tofa_lm_multi_fixedSexAge_mixedStudyID_res <- tofa_regressions_multi_fixedSexAge_mixedStudyID |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, "Event_Name")) |>
  select(Feature, AlamarTargetID, Category, n_observations = nobs, log2FoldChange = estimate, conf.low, conf.high, statistic, pvalue = p.value) |> 
  mutate(qvalue = p.adjust(pvalue, method = "BH")) |> 
  mutate(
    comparison = "16 week vs Baseline",
    .after = Category) |> 
  mutate(model = "~Event_Name+(1|ParticipantID)+Age+Sex")

tofa_lm_multi_fixedSexAge_mixedStudyID_res |> filter(qvalue < 0.1) # 42 sig.

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

#### 2.1.2 Read in HTP NULISA data ------
htp_nulisa_data <- htp_nulisa_data_file |> 
  read_tsv() |>
  select(Feature = Feature_name, AlamarTargetID = FeatureID, Category = Feature_category, Value, everything())
#
htp_nulisa_data |> distinct(ParticipantID) # 224
htp_nulisa_data |> distinct(VisitID) # 224
htp_nulisa_data |> distinct(Feature) # 130

#### 2.1.3 Join HTP NULISA data with metadata ------
# Subset to ages between 12-40 to match with age of TOFA cohort
htp_nulisa_data_meta <- htp_nulisa_data |> 
  inner_join(htp_metadata) |> 
  filter(dplyr::between(Age_years_at_visit, 12, 40)) 

htp_nulisa_data_meta |>  distinct(VisitID, Karyotype) |>  count(Karyotype)
# Karyotype     n
# Control      78
# T21          64

## 2.2 Calculate amyloid and Tau ratios ------

# prepare df with ratios
htp_nulisa_ratios_data_meta <- htp_nulisa_data_meta |> 
  filter(Feature %in% c("pTau-217", "Aβ42", "Aβ40")) |> 
  select(ParticipantID, VisitID, Age_years_at_visit, Sex, Karyotype, Sample_batch, Sample_source_code, Feature, Value) |> 
  pivot_wider(names_from = Feature, values_from = Value) |> 
  # calculate log2(ratio) as difference since values are in log-scale
  mutate(ptau217_Aβ42_ratio = `pTau-217`-Aβ42,
         Aβ42_Aβ40_ratio = `Aβ42`-`Aβ40`) |> 
  pivot_longer(cols = -c(ParticipantID, VisitID, Age_years_at_visit, Sex, Karyotype, Sample_batch, Sample_source_code), names_to = "Feature", values_to = "Value") |> 
  filter(Feature %in% c("ptau217_Aβ42_ratio", "Aβ42_Aβ40_ratio")) |> 
  mutate(Category = "ratio")

## 2.3 HTP linear regression -----
htp_regressions_dat <- htp_nulisa_data_meta |> 
  select(ParticipantID, VisitID, Age_years_at_visit, Sex, Karyotype, Sample_batch, Sample_source_code, Feature, AlamarTargetID, Category, Value) |> 
  # add ratios
  bind_rows(
    htp_nulisa_ratios_data_meta
  ) |> 
  # remove extreme outliers
  mutate(extreme = rstatix::is_extreme(Value), .by = c(Feature, Karyotype)) |> 
  filter(extreme == FALSE) |> 
  nest(data = -c(Feature, AlamarTargetID, Category))

htp_regressions_multi_fixedSexAgeSourceBatch <- htp_regressions_dat |> 
  mutate(
    fit = map(data, ~ lm(Value ~ Karyotype + Age_years_at_visit + Sex + Sample_batch + Sample_source_code, data = .x)),
    tidied = map(fit, broom::tidy), # see ?tidy.lm
    glanced = map(fit, broom::glance), # see ?glance.lm
    augmented = map(fit, broom::augment), # see ?augment.lm
    vifs = map(fit, ~car::vif(mod = .x) |> as_tibble(rownames = "term")) # check co-linearity of variables
  )

## 2.4 Extract results for HTP T21 vs D21 -----
htp_multi_fixedSexAgeSourceBatch_res <- htp_regressions_multi_fixedSexAgeSourceBatch |> 
  unnest(tidied) |> 
  filter(str_detect(term, "Karyotype")) |>
  select(Feature, AlamarTargetID, Category, log2FoldChange = estimate, pvalue = p.value) |> 
  mutate(qvalue = p.adjust(pvalue, method = "BH")) |> 
  mutate(
    comparison = "T21 vs D21",
    .after = Category) |> 
  mutate(model = "~Karyotype+Sex+Age+Sample_source+Sample_batch")

htp_multi_fixedSexAgeSourceBatch_res |>  filter(qvalue < 0.1) # 69
#

# volcano plot
htp_multi_fixedSexAgeSourceBatch_res %>% 
  volcano_plot_lab_lm(
    title="Diff. abund. in T21 vs. D21", 
    subtitle = paste0("~Karyotype+Sex+Age+Sample_source+Sample_batch","\n[Up: ", 
                      (.) |> filter(qvalue < 0.1 & log2FoldChange >0) |> nrow(), "; Down: ", 
                      (.) |> filter(qvalue < 0.1 & log2FoldChange <0) |> nrow(), "]", "\ndown                      up")
  )


# 3. Compare TOFA and HTP results ------

## 3.1.  FC comparison --------

# Create combined results table for FC
log2fc_cutoff <- 0 
fdr_cutoff <- 0.1 

htp_tofa_comb_FC <- htp_multi_fixedSexAgeSourceBatch_res |> 
  select(Feature, log2FoldChange, pvalue, qvalue) |> 
  inner_join(
    tofa_lm_multi_fixedSexAge_mixedStudyID_res |> 
      select(Feature, log2FoldChange, pvalue, qvalue),
    by = c("Feature"),
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
# exacerbated               10
# reversed                  17
# unaffected by T21         63
# unchanged                 42

## 3.2 HTP + TOFA sina plot of example proteins -----

### Combined batch effect removal ----

nulisa_combined <- bind_rows(
  tofa_nulisa_data_meta |>
    mutate(Karyotype = "T21", Sample_source_code = "TOFA") |> 
    select(ParticipantID, VisitID, Karyotype, Sex, Age = Age_years_at_visit, Event_Name, Sample_batch, Sample_source_code, Feature, NPQ = Value),
  htp_nulisa_data_meta |> 
    mutate(Event_Name = "Baseline") |> 
    select(ParticipantID, VisitID, Karyotype, Sex, Age = Age_years_at_visit, Event_Name, Sample_batch, Sample_source_code, Feature, NPQ = Value)
) |> 
  # combined group and batch:
  mutate(
    batch_combined = paste0(Sample_source_code, "_", Sample_batch),
    group = paste0(Karyotype, "_", Event_Name)
  )
nulisa_combined |> distinct(batch_combined)
nulisa_combined |> distinct(group)
nulisa_combined |> distinct(group, batch_combined)
#
sample_data_combined <- nulisa_combined |> 
  select(-c(NPQ))
unadj_data_combined <- nulisa_combined |> 
  select(VisitID, Feature, NPQ) |> 
  pivot_wider(names_from = VisitID, values_from = NPQ) |> # NPX is already log-transformed
  column_to_rownames(var = "Feature")
#
nulisa_combined_BatchSource_adj <- unadj_data_combined |> 
  limma::removeBatchEffect(
    # adjust for batch_combined:
    batch = unadj_data_combined |> colnames() |> enframe(name = NULL, value = "VisitID") |> 
      inner_join(sample_data_combined |> select(VisitID, batch_combined) |> distinct()) |> pull(batch_combined),
    design = unadj_data_combined |> colnames() |> enframe(name = NULL, value = "VisitID") |> 
      inner_join(sample_data_combined |> select(VisitID, group) |> distinct()) %>% 
      model.matrix(~ group, data = .)
  ) |> 
  as_tibble(rownames = "Feature") |> # convert back to tibble
  pivot_longer(-Feature, names_to = "VisitID", values_to = "NPQ_adj") |> 
  inner_join(sample_data_combined)
#

### reversed proteins -----
reversed_examples <- htp_tofa_comb_FC |> 
  filter(sig_TOFA_effect_on_T21 == "reversed") |> 
  arrange(-log2FoldChange_T21) |> 
  pull(Feature)
#

# create annotation df for ggsignif:
reversed_examples_signif_df <- nulisa_combined_BatchSource_adj |> 
  filter(Feature %in% reversed_examples)|> 
  # re-code batch back to HTP/TOFA
  mutate(batch = if_else(Sample_source_code == "TOFA", "TOFA", "HTP")) |> 
  # calculate y position per gene and group
  mutate(extreme = rstatix::is_extreme((NPQ_adj)), .by = c(Feature, batch, group)) |> 
  filter(extreme != TRUE) |> 
  summarise(y_position = max(NPQ_adj, na.rm = TRUE) *1.01, 
            .by = c(Feature, batch)) |> 
  mutate(
    xmin = if_else(batch == "HTP", "HTP:Control_Baseline", "TOFA:T21_Baseline"),
    xmax = if_else(batch == "HTP", "HTP:T21_Baseline", "TOFA:T21_16 week")
  ) |> 
  # add qvalue
  inner_join(htp_tofa_comb_FC |> 
               filter(Feature %in% reversed_examples) |> 
               select(Feature, contains("qvalue")) |> 
               pivot_longer(cols = -Feature, names_to = "batch", values_to = "qvalue") |> 
               mutate(batch = if_else(str_detect(batch, "T21"), "HTP", "TOFA"))) |> 
  mutate(annotations = if_else(
    qvalue < 0.01,
    scales::label_scientific(digits = 2)(qvalue),
    scales::number(qvalue, accuracy = 0.01)
  )) |> 
  mutate(Feature = fct_relevel(Feature, reversed_examples)) # control plotting order
#
reversed_examples_sina <- nulisa_combined_BatchSource_adj |> 
  filter(Feature %in% reversed_examples) |> 
  mutate(Feature = fct_relevel(Feature, reversed_examples)) |>  # control plotting order
  # re-code batch back to HTP/TOFA
  mutate(batch = if_else(Sample_source_code == "TOFA", "TOFA", "HTP")) |> 
  mutate(group = paste0(batch, ":", group) %>% 
           fct_relevel(., c("HTP:Control_Baseline", "HTP:T21_Baseline", 
                            "TOFA:T21_Baseline", "TOFA:T21_16 week"))) |> 
  group_by(Feature, group) |> 
  mutate(extreme = rstatix::is_extreme((NPQ_adj))) |> 
  filter(extreme != TRUE) |> 
  ungroup() |> 
  ggplot(aes(group, NPQ_adj, color = group)) + # CUSTOMIZE
  ggrastr::rasterise(geom_sina(maxwidth = 0.5), dpi= 600) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.5, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~Feature, scales = "free", nrow = 1) +
  scale_color_manual(values = c("TOFA:T21_Baseline" = "#999999", "TOFA:T21_16 week" = "#6baed6", "HTP:Control_Baseline" = "grey30", "HTP:T21_Baseline" = "#009b4e")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    aspect.ratio = 1.2,
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5)
  ) +
  labs(
    title = "Nulisa: Reversed by TOFA",
    subtitle = "batch+source adj.; Extreme outliers removed",
    x = NULL
  )  + 
  ggsignif::geom_signif(
    data = reversed_examples_signif_df,
    aes(xmin = xmin, xmax = xmax,
        annotations = annotations,
        y_position = y_position),
    inherit.aes = FALSE, manual = TRUE, tip_length = 0, color = "black", textsize = 3) 
ggsave(reversed_examples_sina, filename = here("plots", paste0(out_file_prefix, "tofa_reversed_examples_sina", ".pdf")), device = cairo_pdf, width = 60, height = 4, units = "in", limitsize = F)
#

### exacerbated proteins -----
exacerbated_examples <- htp_tofa_comb_FC |> 
  filter(sig_TOFA_effect_on_T21 == "exacerbated") |> 
  arrange(-log2FoldChange_T21) |> 
  pull(Feature)
#
# create annotation df for ggsignif:
exacerbated_examples_signif_df <- nulisa_combined_BatchSource_adj |> 
  filter(Feature %in% exacerbated_examples) |> 
  # re-code batch back to HTP/TOFA
  mutate(batch = if_else(Sample_source_code == "TOFA", "TOFA", "HTP")) |> 
  # calculate y position per gene and group
  mutate(extreme = rstatix::is_extreme((NPQ_adj)), .by = c(Feature, batch, group)) |> 
  filter(extreme != TRUE) |> 
  summarise(y_position = max(NPQ_adj, na.rm = TRUE) *1.01, 
            .by = c(Feature, batch)) |> 
  mutate(
    xmin = if_else(batch == "HTP", "HTP:Control_Baseline", "TOFA:T21_Baseline"),
    xmax = if_else(batch == "HTP", "HTP:T21_Baseline", "TOFA:T21_16 week")
  ) |> 
  # add qvalue
  inner_join(htp_tofa_comb_FC |> 
               filter(Feature %in% exacerbated_examples) |> 
               select(Feature, contains("qvalue")) |> 
               pivot_longer(cols = -Feature, names_to = "batch", values_to = "qvalue") |> 
               mutate(batch = if_else(str_detect(batch, "T21"), "HTP", "TOFA"))) |> 
  mutate(annotations = if_else(
    qvalue < 0.01,
    scales::label_scientific(digits = 2)(qvalue),
    scales::number(qvalue, accuracy = 0.01)
  )) |> 
  mutate(Feature = fct_relevel(Feature, exacerbated_examples)) # control plotting order
#
exacerbated_examples_sina <- nulisa_combined_BatchSource_adj |> 
  filter(Feature %in% exacerbated_examples) |> 
  mutate(Feature = fct_relevel(Feature, exacerbated_examples)) |>  # control plotting order
  # re-code batch back to HTP/TOFA
  mutate(batch = if_else(Sample_source_code == "TOFA", "TOFA", "HTP")) |> 
  mutate(group = paste0(batch, ":", group) %>% 
           fct_relevel(., c("HTP:Control_Baseline", "HTP:T21_Baseline", 
                            "TOFA:T21_Baseline", "TOFA:T21_16 week"))) |> 
  group_by(Feature, group) |> 
  mutate(extreme = rstatix::is_extreme((NPQ_adj))) |> 
  filter(extreme != TRUE) |> 
  ungroup() |> 
  ggplot(aes(group, NPQ_adj, color = group)) + # CUSTOMIZE
  ggrastr::rasterise(geom_sina(maxwidth = 0.5), dpi= 600) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.5, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~Feature, scales = "free", nrow = 1) +
  scale_color_manual(values = c("TOFA:T21_Baseline" = "#999999", "TOFA:T21_16 week" = "#6baed6", "HTP:Control_Baseline" = "grey30", "HTP:T21_Baseline" = "#009b4e")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    aspect.ratio = 1.2,
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5)
  ) +
  labs(
    title = "Nulisa: Exacerbated by TOFA",
    subtitle = "batch+source adj.; Extreme outliers removed",
    x = NULL
  )  + 
  ggsignif::geom_signif(
    data = exacerbated_examples_signif_df,
    aes(xmin = xmin, xmax = xmax,
        annotations = annotations,
        y_position = y_position),
    inherit.aes = FALSE, manual = TRUE, tip_length = 0, color = "black", textsize = 3) 
ggsave(exacerbated_examples_sina, filename = here("plots", paste0(out_file_prefix, "tofa_exacerbated_examples_sina", ".pdf")), device = cairo_pdf, width = 60, height = 4, units = "in", limitsize = F)
#

### select Amyloid & Tau features ---------
select_features <- c("Aβ38", "Aβ40", "Aβ42", "NEFL", "GFAP", "pTau-181", "pTau-217", "pTau-231")
# create annotation df for ggsignif:
select_features_signif_df <- nulisa_combined_BatchSource_adj |> 
  filter(Feature %in% select_features) |> 
  # re-code batch back to HTP/TOFA
  mutate(batch = if_else(Sample_source_code == "TOFA", "TOFA", "HTP")) |> 
  # calculate y position per gene and group
  mutate(extreme = rstatix::is_extreme((NPQ_adj)), .by = c(Feature, batch, group)) |> 
  filter(extreme != TRUE) |> 
  summarise(y_position = max(NPQ_adj, na.rm = TRUE) *1.01, 
            .by = c(Feature, batch)) |> 
  mutate(
    xmin = if_else(batch == "HTP", "HTP:Control_Baseline", "TOFA:T21_Baseline"),
    xmax = if_else(batch == "HTP", "HTP:T21_Baseline", "TOFA:T21_16 week")
  ) |> 
  # add qvalue
  inner_join(htp_tofa_comb_FC |> 
               filter(Feature %in% select_features) |> 
               select(Feature, contains("qvalue")) |> 
               pivot_longer(cols = -Feature, names_to = "batch", values_to = "qvalue") |> 
               mutate(batch = if_else(str_detect(batch, "T21"), "HTP", "TOFA"))) |> 
  mutate(annotations = if_else(
    qvalue < 0.01,
    scales::label_scientific(digits = 2)(qvalue),
    scales::number(qvalue, accuracy = 0.01)
  )) |> 
  mutate(Feature = fct_relevel(Feature, select_features)) # control plotting order
#
select_features_sina <- nulisa_combined_BatchSource_adj |> 
  filter(Feature %in% select_features) |> 
  mutate(Feature = fct_relevel(Feature, select_features)) |>  # control plotting order
  # re-code batch back to HTP/TOFA
  mutate(batch = if_else(Sample_source_code == "TOFA", "TOFA", "HTP")) |> 
  mutate(group = paste0(batch, ":", group) %>% 
           fct_relevel(., c("HTP:Control_Baseline", "HTP:T21_Baseline", 
                            "TOFA:T21_Baseline", "TOFA:T21_16 week"))) |> 
  group_by(Feature, group) |> 
  mutate(extreme = rstatix::is_extreme((NPQ_adj))) |> 
  filter(extreme != TRUE) |> 
  ungroup() |> 
  ggplot(aes(group, NPQ_adj, color = group)) + # CUSTOMIZE
  ggrastr::rasterise(geom_sina(maxwidth = 0.5), dpi= 600) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.5, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~Feature, scales = "free", nrow = 1) +
  scale_color_manual(values = c("TOFA:T21_Baseline" = "#999999", "TOFA:T21_16 week" = "#6baed6", "HTP:Control_Baseline" = "grey30", "HTP:T21_Baseline" = "#009b4e")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    aspect.ratio = 1.2,
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5)
  ) +
  labs(
    title = "Nulisa: Amyloid & Tau features",
    subtitle = "batch+source adj.; Extreme outliers removed",
    x = NULL
  )  + 
  ggsignif::geom_signif(
    data = select_features_signif_df,
    aes(xmin = xmin, xmax = xmax,
        annotations = annotations,
        y_position = y_position),
    inherit.aes = FALSE, manual = TRUE, tip_length = 0, color = "black", textsize = 3) 
ggsave(select_features_sina, filename = here("plots", paste0(out_file_prefix, "tofa_amyloid_tau_sina", ".pdf")), device = cairo_pdf, width = 60, height = 4, units = "in", limitsize = F)
#

### Amyloid/Tau ratios ---------
# calculate ratios from combined adjusted data:
nulisa_combined_BatchSource_adj_ratios <- nulisa_combined_BatchSource_adj |> 
  filter(Feature %in% c("pTau-217", "Aβ42", "Aβ40")) |> 
  # re-code batch back to HTP/TOFA
  mutate(batch = if_else(Sample_source_code == "TOFA", "TOFA", "HTP")) |> 
  select(VisitID, group, batch, Feature, NPQ_adj)|> 
  pivot_wider(names_from = Feature, values_from = NPQ_adj) |> 
  # calculate log2(ratio) as difference since all values are in log-scale
  mutate(
    ptau217_Aβ42_ratio = `pTau-217`-Aβ42,
    Aβ42_Aβ40_ratio = `Aβ42`-`Aβ40`) |> 
  select(-c("pTau-217", "Aβ42", "Aβ40")) |> 
  pivot_longer(cols = -c(VisitID, group, batch), names_to = "Feature", values_to = "log2ratio") |> 
  mutate(ratio = 2^log2ratio) 

# create annotation df for ggsignif:
ratios_signif_df <- nulisa_combined_BatchSource_adj_ratios |> 
  # calculate y position per gene and group
  mutate(extreme = rstatix::is_extreme((ratio)), .by = c(Feature, batch, group)) |> 
  filter(extreme != TRUE) |> 
  summarise(y_position = max(ratio, na.rm = TRUE)*1.01, 
            .by = c(Feature, batch)) |> 
  mutate(
    xmin = if_else(batch == "HTP", "HTP:Control_Baseline", "TOFA:T21_Baseline"),
    xmax = if_else(batch == "HTP", "HTP:T21_Baseline", "TOFA:T21_16 week")
  ) |> 
  # add qvalue
  inner_join(htp_tofa_comb_FC |> 
               filter(str_detect(Feature, "ratio")) |> 
               select(Feature, contains("qvalue")) |> 
               pivot_longer(cols = -Feature, names_to = "batch", values_to = "qvalue") |> 
               mutate(batch = if_else(str_detect(batch, "T21"), "HTP", "TOFA"))) |> 
  mutate(annotations = if_else(
    qvalue < 0.01,
    scales::label_scientific(digits = 2)(qvalue),
    scales::number(qvalue, accuracy = 0.01)
  )) 
#

ratios_sina <- nulisa_combined_BatchSource_adj_ratios |>
  mutate(group = paste0(batch, ":", group) %>% 
           fct_relevel(., c("HTP:Control_Baseline", "HTP:T21_Baseline", 
                            "TOFA:T21_Baseline", "TOFA:T21_16 week"))) |> 
  group_by(Feature, group) |>
  mutate(extreme = rstatix::is_extreme((ratio))) |>
  filter(extreme != TRUE) |>
  ungroup() |>
  ggplot(aes(group, ratio, color = group)) + # CUSTOMIZE
  ggrastr::rasterise(geom_sina(maxwidth = 0.5), dpi= 600) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.5, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~Feature, scales = "free", nrow = 1) +
  scale_color_manual(values = c("TOFA:T21_Baseline" = "#999999", "TOFA:T21_16 week" = "#6baed6", "HTP:Control_Baseline" = "grey30", "HTP:T21_Baseline" = "#009b4e")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    aspect.ratio = 1.2,
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5)
  ) +
  labs(
    title = "Nulisa: Ratios",
    subtitle = "batch+source adj.; Extreme outliers removed",
    x = NULL
  )  + 
  ggsignif::geom_signif(
    data = ratios_signif_df,
    aes(xmin = xmin, xmax = xmax,
        annotations = annotations,
        y_position = y_position),
    inherit.aes = FALSE, manual = TRUE, tip_length = 0, color = "black", textsize = 3) 
ggsave(ratios_sina, filename = here("plots", paste0(out_file_prefix, "ratios_sina", ".pdf")), device = cairo_pdf, width = 20, height = 4, units = "in", limitsize = F)
#

## 3.3 Export results -----
# combine all results into a single excel file
list(
  "NULISA TOFA 16wk" = tofa_lm_multi_fixedSexAge_mixedStudyID_res,
  "NULISA T21 vs D21" = htp_multi_fixedSexAgeSourceBatch_res,
  "NULISA T21 vs TOFA 16wk" = htp_tofa_comb_FC
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
