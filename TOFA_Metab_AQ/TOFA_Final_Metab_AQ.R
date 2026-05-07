################################################
# Title: Tofacitinib for Immune Skin Conditions in Down Syndrome: Analysis of Plasma Endpoint Metabolites Absolute Quantitation (AQ)
# Author(s):
#   - Matthew Galbraith
# Affiliation(s):
#   - Linda Crnic Institute for Down Syndrome
#   - University of Colorado Anschutz
################################################

### Summary:  
# Analysis of endpoint metabolites absolute concentrations, measured by UHPLC-MS of plasma.
# See README.md for more details.
#  

### Data type(s):
# Clinical trial (TOFA) datasets:
#    * Participant-level metadata; Available on request.
#    * Visit/Event-level metadata; Available on request.
#    * Baseline obesity status; Available on request.
#    * COVID-19 history; Available on request.
#    * Plasma endpoint metabolites absolute quantitation; DOI: 10.5281/zenodo.20046361
#
# Human Trisome Project (HTP) datasets:
#    * Participant-level metadata; DOI: 10.5281/zenodo.19962380
#    * Visit/Event-level metadata; DOI: 10.5281/zenodo.19962380
#    * Plasma metabolites absolute quantitation; DOI: 10.5281/zenodo.20074289
# 


# 0 General Setup -----
# RUN THIS FIRST TIME - Initialize and install packages with renv:
# renv::init(bioconductor = TRUE)
#
# To install the exact versions of all R packages base on renv.lock file (requires matching R version):
# renv::restore()


## 0.1 Load required libraries ----
library("tidyverse") # required for ggplot2, dplyr etc
library("readxl") # reading Excel files
library("openxlsx") # for exporting results as Excel workbooks
library("ggforce") # used for sina plots
library("ggrastr") # required for rasterizing some layers of plots
library("ggrepel") # plot labels
library("patchwork") # assembling plots
library("rstatix") # stats functions including t tests
library("coin") # required for calculation of standardized effect sizes
library("lme4") # for mixed models
library("lmerTest") # for mixed models
library("broom.mixed") # for mixed models
library("conflicted") # force all conflicts to become errors
conflicts_prefer( # declare preferences in cases of conflict
  dplyr::filter,
  dplyr::select,
  dplyr::count
)
library("here") # generates path to current project directory
#

## 0.2 Set input files and other parameters ----
#
# Input data files
# Datasets used in this study can be obtained from the associated
# repositories (further details in README.md).
# Download each dataset to /data directory within this R project.
#
# Clinical trial datasets:
tofa_plasma_metab_AQ_data_file <- here("data/TOFA_Plasma_Metab_AQ_data_zenodo_v1.txt.gz") # Source: UPDATE
tofa_participant_meta_data_file <- here("data/TOFA_Participant_metadata_zenodo_v1.txt") # Source: Available on request
tofa_visit_meta_data_file <- here("data/TOFA_Visit_metadata_zenodo_v1.txt") # Source: Available on request
tofa_baseline_obesity_file <- here("data/TOFA_Baseline_Obesity_Status_zenodo_v1.txt") # Source: Available on request
tofa_covid_history_file <- here("data/TOFA_COVID_History_zenodo_v1.txt") # Source: Available on request
# Human Trisome Project datasets:
htp_participant_meta_data_file <- here("data/HTP_Participant_metadata_zenodo_v1.txt") # Source: 10.5281/zenodo.19962380
htp_visit_meta_data_file <- here("data/HTP_Visit_metadata_zenodo_v1.txt") # Source: 10.5281/zenodo.19962380
htp_plasma_metab_AQ_data_file <- here("data/HTP_Plasma_Metab_AQ_data_zenodo_v1.txt") # Source: UPDATE
#
# Other parameters:
standard_colors <- c("Baseline" = "#999999", "2 week" = "#c6dbef", "8 week" = "#9ecae1", "16 week" = "#6baed6", "40 week" = "#4292c6")
#
out_file_prefix <- "TOFA_Final_Metab_AQ_" # should match this script title
#
source(here("helper_functions.R")) # load helper functions
#


# 1 Read in data ----
## 1.1 Read in TOFA meta data ----
### 1.1.1 Participant level meta data ----
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
### 1.1.2 Event/Visit level meta data ----
tofa_visit_meta_data <- tofa_visit_meta_data_file |> 
  read_tsv() |> 
  mutate(
    Event_Name = fct_relevel(Event_Name, c("Baseline", "2 week", "8 week", "16 week","40 week")) # set factor levels
  )
#
tofa_visit_meta_data # 202 rows
tofa_visit_meta_data |> distinct(ParticipantID) # 43 Participants with samples
tofa_visit_meta_data |> distinct(VisitID) # 306 Visits/Samples
#

## 1.2 Read in TOFA Metab AQ data ----
tofa_plasma_metab_AQ_data <- tofa_plasma_metab_AQ_data_file |>
  read_tsv()
#
tofa_plasma_metab_AQ_data # 570 rows
tofa_plasma_metab_AQ_data |> distinct(ParticipantID) # 41 Participants
tofa_plasma_metab_AQ_data |> distinct(VisitID) # 190 VisitIDs = Plasma samples
tofa_plasma_metab_AQ_data |> distinct(Compound_Name) # 3 metabolites
#

## 1.3 Join with meta data, filter eligible etc, filter to endpoint cytokines ----
tofa_plasma_metab_AQ_data <- tofa_plasma_metab_AQ_data |> 
  inner_join(tofa_visit_meta_data) |> # returns 570 of 570 rows
  inner_join(tofa_participant_meta_data) |> # returns 570 of 570 rows
  filter(Endpoint_eligible == TRUE) # returns 570 of 570 rows
#
### 1.3.1 ParticipantID vs Event_Name summary ---- 
tofa_plasma_metab_AQ_data |> 
  distinct(ParticipantID, VisitID, Event_Name) |> 
  arrange(Event_Name) |> 
  pivot_wider(names_from = Event_Name, values_from = VisitID) |> 
  print(n = Inf)
# Note: Some Participants have missing samples and not all completed extension to 40 weeks.

## 1.4 Read in HTP data ----
### 1.4.1 Read in HTP Participant level meta data ----
htp_participant_meta_data <- htp_participant_meta_data_file |>
  read_tsv() |> 
  mutate(
    Karyotype = fct_relevel(Karyotype, c("Control", "T21")), # convert to factor and set order
    Sex = fct_relevel(Sex, "Female"), # convert to factor and set order
  )
# inspect
htp_participant_meta_data # 1,529 rows
htp_participant_meta_data |> count(ParticipantID) # 1,529 unique participants
htp_participant_meta_data |> count(Karyotype) # 1064 T21
htp_participant_meta_data |> count(Sex) # 783 Female
#
### 1.4.2 Read in HTP Visit/Event level meta data ----
htp_visit_meta_data <- htp_visit_meta_data_file |>
  read_tsv() |> 
  mutate(
    Sample_source_code = as_factor(Sample_source_code),
  )
#
htp_visit_meta_data # 2,338 rows
htp_visit_meta_data |> count(ParticipantID) |> arrange(-n) # 1,528 unique, participants; some have >1 visit
htp_visit_meta_data |> count(VisitID) |> arrange(-n) # 2,338 unique visits
#
### 1.4.3 Read in HTP Metab AQ data ----
htp_plasma_metab_AQ_data <- htp_plasma_metab_AQ_data_file |>
  read_tsv()
#
htp_plasma_metab_AQ_data |> distinct(ParticipantID) # 419 Participants
htp_plasma_metab_AQ_data |> distinct(VisitID) # 419 Visits/Plasma samples
htp_plasma_metab_AQ_data |> distinct(Compound_Name) # 2 metabolites
#


## 1.5 Calculate Kyn/Trp ratio ----
### 1.5.1 TOFA ----
tofa_KynTrp_ratio <- tofa_plasma_metab_AQ_data |> 
  filter(Compound_Name %in% c("kynurenine", "L-tryptophan")) |> 
  select(ParticipantID, VisitID, Compound_Name, Value) |> 
  pivot_wider(names_from = Compound_Name, values_from = Value) |> 
  mutate(
    Compound_Name = "KynTrp_ratio",
    Value = kynurenine / `L-tryptophan`,
    Units = "ratio",
    .keep = "unused"
  ) |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data)
#
#### Combine ratio with AQ data ----
tofa_plasma_metab_AQ_data_final <- bind_rows(
  tofa_plasma_metab_AQ_data,
  tofa_KynTrp_ratio
)
tofa_plasma_metab_AQ_data_final |> count(Compound_Name, Units)
#
### 1.5.2 HTP ----
htp_KynTrp_ratio <- htp_plasma_metab_AQ_data |> 
  filter(Compound_Name %in% c("Kynurenine", "Tryptophan")) |> # not really needed in this case
  select(ParticipantID, VisitID, Compound_Name, Value) |> 
  pivot_wider(names_from = Compound_Name, values_from = Value) |> 
  mutate(
    Compound_Name = "KynTrp_ratio",
    Value = Kynurenine / Tryptophan,
    Units = "ratio",
    .keep = "unused"
  )
#
#### Adjust Kyn/Trp ratio for Sex/Age/Sample source ----
htp_KynTrp_ratio_unadj <- htp_KynTrp_ratio |> 
  select(VisitID, Compound_Name, Value) |> 
  mutate(Value = log2(Value)) |> # need to log2 transform for batch correction
  pivot_wider(names_from = VisitID, values_from = Value) |> 
  column_to_rownames(var = "Compound_Name")
#
htp_KynTrp_ratio_SexAgeSource_adj <- htp_KynTrp_ratio_unadj |> 
  limma::removeBatchEffect(
    # get batch information, ensuring order matches unadj_data:
    batch = htp_KynTrp_ratio_unadj |> colnames() |> enframe(name = NULL, value = "VisitID") |> 
      inner_join(htp_visit_meta_data |> distinct(VisitID, Sample_source_code)) |> pull(Sample_source_code),
    # get sex information, ensuring order matches unadj_data:
    batch2 = htp_KynTrp_ratio_unadj |> colnames() |> enframe(name = NULL, value = "VisitID") |> 
      inner_join(htp_visit_meta_data) |> 
      inner_join(htp_participant_meta_data |> distinct(ParticipantID, Sex)) |> pull(Sex),
    # get age information, ensuring order matches unadj_data:
    covariates = htp_KynTrp_ratio_unadj |> colnames() |> enframe(name = NULL, value = "VisitID") |> 
      inner_join(htp_visit_meta_data |> distinct(VisitID, Age_years_at_visit)) |> pull(Age_years_at_visit),
    design = htp_KynTrp_ratio_unadj |> colnames() |> enframe(name = NULL, value = "VisitID") |> 
      # get karyotype information, ensuring order matches unadj_data:
      inner_join(htp_visit_meta_data) |> 
      inner_join(htp_participant_meta_data |> distinct(ParticipantID, Karyotype)) %>% 
      model.matrix(~ Karyotype, data = .)
  ) |> 
  as_tibble(rownames = "Compound_Name") |> # convert back to tibble
  pivot_longer(-Compound_Name, names_to = "VisitID", values_to = "Value_adj") |>  
  mutate(Value_adj = 2^Value_adj) # remove log2 transformation
#


# 2 Endpoint Stats ----
## 2.1 Check data distributions ----
### 2.1.1 check for outlier values ----
tofa_plasma_metab_AQ_data_final |> 
  mutate(
    extreme = rstatix::is_extreme(Value),
    outlier = rstatix::is_outlier(Value),
    .by = Compound_Name
  ) |> 
  count(Compound_Name, outlier, extreme)
#
tofa_plasma_metab_AQ_data_final |> 
  mutate(
    extreme = rstatix::is_extreme(Value),
    outlier = rstatix::is_outlier(Value),
    .by = Compound_Name
  ) |> 
  ggplot(aes("Value", Value)) +
  geom_sina(data = . %>% filter(outlier == FALSE & extreme == FALSE), maxwidth = 0.5, color = "grey") +
  geom_sina(data = . %>% filter(outlier == TRUE & extreme == FALSE), aes(color = "outlier"), maxwidth = 0.1) +
  geom_sina(data = . %>% filter(extreme == TRUE), aes(color = "extreme"), maxwidth = 0.1) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~ Compound_Name, scales = "free_y", nrow = 1) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    aspect.ratio = 1.3
  ) +
  labs(title = "Plasma Metabolite outliers")
#
### 2.1.2 Calc Differences from baseline ----
tofa_metab_diffs_2_weeks <- tofa_plasma_metab_AQ_data_final |> 
  filter(Event_Name %in% c("Baseline", "2 week")) |> 
  filter(!is.na(Value)) |>
  filter(is.finite(Value)) |> 
  add_count(Compound_Name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  select(ParticipantID, Event_Name, Compound_Name, Value) |> 
  pivot_wider(names_from = Event_Name, values_from = Value) |> 
  mutate(difference = `2 week` - Baseline) |> 
  select(ParticipantID, Compound_Name, difference) |> 
  mutate(
    extreme = rstatix::is_extreme(difference),
    outlier = rstatix::is_outlier(difference),
    .by = Compound_Name
  )
#
tofa_metab_diffs_8_weeks <- tofa_plasma_metab_AQ_data_final |> 
  filter(Event_Name %in% c("Baseline", "8 week")) |> 
  filter(!is.na(Value)) |>
  filter(is.finite(Value)) |> 
  add_count(Compound_Name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, Compound_Name) |> 
  select(ParticipantID, Event_Name, Compound_Name, Value) |> 
  pivot_wider(names_from = Event_Name, values_from = Value) |> 
  mutate(difference = `8 week` - Baseline) |> 
  select(ParticipantID, Compound_Name, difference) |> 
  mutate(
    extreme = rstatix::is_extreme(difference),
    outlier = rstatix::is_outlier(difference),
    .by = Compound_Name
  )
#
tofa_metab_diffs_16_weeks <- tofa_plasma_metab_AQ_data_final |> 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(Value)) |>
  filter(is.finite(Value)) |> 
  add_count(Compound_Name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, Compound_Name) |> 
  select(ParticipantID, Event_Name, Compound_Name, Value) |> 
  pivot_wider(names_from = Event_Name, values_from = Value) |> 
  mutate(difference = `16 week` - Baseline) |> 
  select(ParticipantID, Compound_Name, difference) |> 
  mutate(
    extreme = rstatix::is_extreme(difference),
    outlier = rstatix::is_outlier(difference),
    .by = Compound_Name
  )
#
tofa_metab_diffs_40_weeks <- tofa_plasma_metab_AQ_data_final |> 
  filter(Event_Name %in% c("Baseline", "40 week")) |> 
  filter(!is.na(Value)) |>
  filter(is.finite(Value)) |> 
  add_count(Compound_Name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, Compound_Name) |> 
  select(ParticipantID, Event_Name, Compound_Name, Value) |> 
  pivot_wider(names_from = Event_Name, values_from = Value) |> 
  mutate(difference = `40 week` - Baseline) |> 
  select(ParticipantID, Compound_Name, difference) |> 
  mutate(
    extreme = rstatix::is_extreme(difference),
    outlier = rstatix::is_outlier(difference),
    .by = Compound_Name
  )
#
### 2.1.3 Count extreme outliers in differences ----
tofa_metab_diffs_2_weeks |> 
  count(outlier, extreme)
tofa_metab_diffs_8_weeks |> 
  count(outlier, extreme)
tofa_metab_diffs_16_weeks |> 
  count(outlier, extreme)
tofa_metab_diffs_40_weeks |> 
  count(outlier, extreme)
#
### 2.1.4 Plot differences ---- 
tofa_metab_diffs_2_weeks |> 
  ggplot(aes(Compound_Name, difference, color = Compound_Name)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_sina(data = . %>% filter(extreme == FALSE), maxwidth = 0.5) +
  geom_sina(data = . %>% filter(extreme == TRUE), color = "red", maxwidth = 0.1) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~Compound_Name, scales = "free", nrow = 1) +
  theme(
    legend.title = element_blank(),
    aspect.ratio = 1.3
  ) +
  labs(title = "2 weeks: Distributions of differences vs. baseline", x = NULL)
#
tofa_metab_diffs_8_weeks |> 
  ggplot(aes(Compound_Name, difference, color = Compound_Name)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_sina(data = . %>% filter(extreme == FALSE), maxwidth = 0.5) +
  geom_sina(data = . %>% filter(extreme == TRUE), color = "red", maxwidth = 0.1) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~Compound_Name, scales = "free", nrow = 1) +
  theme(
    legend.title = element_blank(),
    aspect.ratio = 1.3
  ) +
  labs(title = "8 weeks: Distributions of differences vs. baseline", x = NULL)
#
tofa_metab_diffs_16_weeks |> 
  ggplot(aes(Compound_Name, difference, color = Compound_Name)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_sina(data = . %>% filter(extreme == FALSE), maxwidth = 0.5) +
  geom_sina(data = . %>% filter(extreme == TRUE), color = "red", maxwidth = 0.1) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~Compound_Name, scales = "free", nrow = 1) +
  theme(
    legend.title = element_blank(),
    aspect.ratio = 1.3
  ) +
  labs(title = "16 weeks: Distributions of differences vs. baseline", x = NULL)
#
tofa_metab_diffs_40_weeks |> 
  ggplot(aes(Compound_Name, difference, color = Compound_Name)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_sina(data = . %>% filter(extreme == FALSE), maxwidth = 0.5) +
  geom_sina(data = . %>% filter(extreme == TRUE), color = "red", maxwidth = 0.1) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~Compound_Name, scales = "free", nrow = 1) +
  theme(
    legend.title = element_blank(),
    aspect.ratio = 1.3
  ) +
  labs(title = "40 weeks: Distributions of differences vs. baseline", x = NULL)
#
### 2.1.5 Check if differences are normally distributed ----
tofa_metab_diffs_2_weeks |>
  ggpubr::ggqqplot("difference", title = "Q-Q plot: 2 week differences", facet.by = "Compound_Name",  scales = "free_y")
tofa_metab_diffs_8_weeks |>
  ggpubr::ggqqplot("difference", title = "Q-Q plot: 8 week differences", facet.by = "Compound_Name",  scales = "free_y")
tofa_metab_diffs_16_weeks |>
  ggpubr::ggqqplot("difference", title = "Q-Q plot: 16 week differences", facet.by = "Compound_Name",  scales = "free_y")
tofa_metab_diffs_40_weeks |>
  ggpubr::ggqqplot("difference", title = "Q-Q plot: 40 week differences", facet.by = "Compound_Name",  scales = "free_y")
#

## 2.2 Paired t tests  --------
# https://www.datanovia.com/en/lessons/t-test-in-r/#pstt
# Assumptions: the two groups are paired; No significant outliers in the difference between the two related groups; Normality. the difference of pairs follow a normal distribution
# Assessing equality of variances. Homogeneity of variances can be checked using
# the Levene’s test. Note that, by default, the t_test() function does not
# assume equal variances; instead of the standard Student’s t-test, it uses the
# Welch t-test by default, which is the considered the safer one. To use
# Student’s t-test, set var.equal = TRUE. The two methods give very similar
# results unless both the group sizes and the standard deviations are very
# different.
#

## 2.2.1 Run paired t tests 16 weeks ----
# Only 16 week timepoint is tested as a trial endpoint.
# Only kynurenine, quinolinic acid, and KynTrp_ratio are endpoints (tryptophan not tested)
tofa_plasma_metab_AQ_Ttest_res_16weeks <- tofa_plasma_metab_AQ_data_final |> 
  arrange(ParticipantID, Event_Name) |> # ensure correct sort order for rstatix::t_test()
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(Compound_Name != "L-tryptophan") |> 
  filter(!is.na(Value)) |>
  filter(is.finite(Value)) |> 
  add_count(Compound_Name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, Compound_Name) |> 
  group_by(Compound_Name) |> 
  rstatix::t_test(
    formula = Value ~ Event_Name,
    ref.group = "16 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    var.equal = TRUE, # set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) |> 
  mutate(BH_padj = p.adjust(p, method = "BH")) |>
  arrange(p)
#
### 2.2.2 Paired t test effect size ----
# The effect size for a paired-samples t-test can be calculated by dividing the
# mean difference by the standard deviation of the difference, as shown below.
# Cohen’s formula:
# d = mean(D)/sd(D), where D is the differences of the paired samples values.
tofa_plasma_metab_AQ_Ttest_effsize_16weeks <- tofa_plasma_metab_AQ_data_final |> 
  arrange(ParticipantID, Event_Name) |> # ensure correct sort order for rstatix::cohens_d()
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(Compound_Name != "L-tryptophan") |> 
  filter(!is.na(Value)) |>
  filter(is.finite(Value)) |> 
  add_count(Compound_Name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  group_by(Compound_Name) |> 
  rstatix::cohens_d(
    formula = Value ~ Event_Name,
    ref.group = "16 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    var.equal = TRUE, # set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
  )
#
### 2.2.3 Compile t test results ----
tofa_plasma_metab_AQ_Ttest_res_16weeks_full <- tofa_plasma_metab_AQ_Ttest_res_16weeks |> 
  inner_join(tofa_plasma_metab_AQ_Ttest_effsize_16weeks) |> 
  mutate(
    mean_diff = estimate,
  ) |> 
  inner_join(tofa_metab_diffs_16_weeks |> summarize(median_diff = median(difference), .by = Compound_Name)) |> 
  select(Compound_Name, mean_diff, median_diff, p, BH_padj, effsize, magnitude, group1, group2, n1, n2, everything())
#
tofa_plasma_metab_AQ_Ttest_res_16weeks_full
#
### 2.2.4 Export results ----
list(
  "ttest_results" = tofa_plasma_metab_AQ_Ttest_res_16weeks_full |> 
    select(
      Compound_Name,
      Timepoint = group1,
      n_pairs = n1,
      Mean_difference = mean_diff,
      Conf.low = conf.low,
      Conf.high = conf.high,
      Statistic = statistic,
      Effect_size = effsize,
      Magnitude = magnitude,
      pvalue = p,
      qvalue = BH_padj
    ) |> 
    mutate(
      Method = "Paired Student's t test, two-sided"
    )
) |> 
  export_excel(filename = "Ttest_results")
#


# 3 Sina plots individual metabolites -----
## 3.1 Baseline vs 16 weeks ---- 
tofa_plasma_metab_AQ_data_final |> 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(Value)) |> 
  filter(is.finite(Value)) |> 
  add_count(ParticipantID, name = "pair") |> 
  filter(pair >= 2) |>  # require at least 2 data points
  #
  group_by(Compound_Name) |> 
  group_split() |> 
  map( # using anonymous function: \(x)
    \(x) x |> 
      arrange(Compound_Name, Event_Name) |> 
      mutate( # add labels with n
        n = n_distinct(ParticipantID), 
        .by = Event_Name,
        label = paste0(str_extract(Event_Name, "^[[:alpha:]]|^\\d{1,2}"), "\n(n=", n, ")"),
        label = fct_inorder(label)
      ) |> 
      mutate(Event_Name = fct_drop(Event_Name)) |> # prevent empty x values
      complete(Compound_Name, ParticipantID, Event_Name) %>% # Prevent error from stat_signif()
      #
      {ggplot(data = ., aes(Event_Name, Value, color = Event_Name)) +
          geom_line(aes(group = ParticipantID), color = "grey90") +
          geom_sina(maxwidth = 0.5) +
          geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
          scale_color_manual(values = standard_colors) +
          facet_wrap(~ Compound_Name, scales = "free_y", nrow = 1) +
          theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            aspect.ratio = 1.5,
            panel.border = element_blank(),
            axis.line = element_line(colour = "black")
          ) +
          labs(
            x = NULL, 
            y = "Concentration (µM)"
          ) +
          # add sample numbers to labels
          scale_x_discrete(
            labels = . |>
              filter(!is.na(label)) |>
              distinct(Compound_Name, Event_Name, label) |>
              select(Event_Name, label) |>
              transmute(as.character(Event_Name), as.character(label)) |>
              deframe()
          )
      }
  ) |> 
  patchwork::wrap_plots() +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(
    title = "Endpoint Metabolites: Baseline vs. Treatment",
    subtitle = "nb: smaller N at 40 weeks"
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_individual_B_16_weeks", ".pdf")), device = cairo_pdf, width = 15, height = 5, units = "in")
#
## 3.2 Baseline vs 2/8/16/40 weeks ----
tofa_plasma_metab_AQ_data_final |> 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> 
  filter(!is.na(Value)) |> 
  filter(is.finite(Value)) |> 
  add_count(ParticipantID, name = "pair") |> 
  filter(pair >= 2) |>  # require at least 2 data points
  #
  group_by(Compound_Name) |> 
  group_split() |> 
  map( # using anonymous function: \(x)
    \(x) x |> 
      arrange(Compound_Name, Event_Name) |> 
      mutate( # add labels with n
        n = n_distinct(ParticipantID), 
        .by = Event_Name,
        label = paste0(str_extract(Event_Name, "^[[:alpha:]]|^\\d{1,2}"), "\n(n=", n, ")"),
        label = fct_inorder(label)
      ) |> 
      mutate(Event_Name = fct_drop(Event_Name)) |> # prevent empty x values
      complete(Compound_Name, ParticipantID, Event_Name) %>% # Prevent error from stat_signif()
      #
      {ggplot(data = ., aes(Event_Name, Value, color = Event_Name)) +
          geom_line(aes(group = ParticipantID), color = "grey90") +
          geom_sina(maxwidth = 0.5) +
          geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
          scale_color_manual(values = standard_colors) +
          facet_wrap(~ Compound_Name, scales = "free_y", nrow = 1) +
          theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            aspect.ratio = 1.5,
            panel.border = element_blank(),
            axis.line = element_line(colour = "black")
          ) +
          labs(
            x = NULL, 
            y = "Concentration (µM)"
          ) +
          # add sample numbers to labels
          scale_x_discrete(
            labels = . |>
              filter(!is.na(label)) |>
              distinct(Compound_Name, Event_Name, label) |>
              select(Event_Name, label) |>
              transmute(as.character(Event_Name), as.character(label)) |>
              deframe()
          )
      }
  ) |> 
  patchwork::wrap_plots() +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(
    title = "Endpoint Metabolites: Baseline vs. Treatment",
    subtitle = "nb: smaller N at 40 weeks"
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_individual_B_2_8_16_40_weeks", ".pdf")), device = cairo_pdf, width = 20, height = 5, units = "in")
#


# 4 Sina plots Kyn/Trp ratio -----
## 4.1 Baseline vs 16 ----
tofa_KynTrp_ratio |> 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(Value)) |> 
  filter(is.finite(Value)) |> 
  add_count(ParticipantID, name = "pair") |> 
  filter(pair >= 2) |>  # require at least 2 data points
  #
  group_by(Compound_Name) |> 
  group_split() |> 
  map( # using anonymous function: \(x)
    \(x) x |> 
      arrange(Compound_Name, Event_Name) |> 
      mutate( # add labels with n
        n = n_distinct(ParticipantID), 
        .by = Event_Name,
        label = paste0(str_extract(Event_Name, "^[[:alpha:]]|^\\d{1,2}"), "\n(n=", n, ")"),
        label = fct_inorder(label)
      ) |> 
      mutate(Event_Name = fct_drop(Event_Name)) |> # prevent empty x values
      complete(Compound_Name, ParticipantID, Event_Name) %>% # Prevent error from stat_signif()
      #
      {ggplot(data = ., aes(Event_Name, Value, color = Event_Name)) +
          geom_line(aes(group = ParticipantID), color = "grey90") +
          geom_sina(maxwidth = 0.5) +
          geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
          scale_color_manual(values = standard_colors) +
          facet_wrap(~ Compound_Name, scales = "free_y", nrow = 1) +
          theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            aspect.ratio = 1.5,
            panel.border = element_blank(),
            axis.line = element_line(colour = "black")
          ) +
          labs(
            x = NULL, 
            y = "Kyn:Trp ratio"
          ) +
          # add sample numbers to labels
          scale_x_discrete(
            labels = . |>
              filter(!is.na(label)) |>
              distinct(Compound_Name, Event_Name, label) |>
              select(Event_Name, label) |>
              transmute(as.character(Event_Name), as.character(label)) |>
              deframe()
          )
      }
  ) |> 
  patchwork::wrap_plots() +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(
    title = "Endpoint Metabolites: Baseline vs. Treatment",
    subtitle = "nb: smaller N at 40 weeks"
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_KynTrpRatio_B_16_weeks", ".pdf")), device = cairo_pdf, width = 15, height = 6, units = "in")
#
## 4.2 Baseline vs 2 weeks vs 8 weeks vs 16 vs 40 weeks ----
tofa_KynTrp_ratio |> 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> 
  filter(!is.na(Value)) |> 
  filter(is.finite(Value)) |> 
  add_count(ParticipantID, name = "pair") |> 
  filter(pair >= 2) |>  # require at least 2 data points
  #
  group_by(Compound_Name) |> 
  group_split() |> 
  map( # using anonymous function: \(x)
    \(x) x |> 
      arrange(Compound_Name, Event_Name) |> 
      mutate( # add labels with n
        n = n_distinct(ParticipantID), 
        .by = Event_Name,
        label = paste0(str_extract(Event_Name, "^[[:alpha:]]|^\\d{1,2}"), "\n(n=", n, ")"),
        label = fct_inorder(label)
      ) |> 
      mutate(Event_Name = fct_drop(Event_Name)) |> # prevent empty x values
      complete(Compound_Name, ParticipantID, Event_Name) %>% # Prevent error from stat_signif()
      #
      {ggplot(data = ., aes(Event_Name, Value, color = Event_Name)) +
          geom_line(aes(group = ParticipantID), color = "grey90") +
          geom_sina(maxwidth = 0.5) +
          geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
          scale_color_manual(values = standard_colors) +
          facet_wrap(~ Compound_Name, scales = "free_y", nrow = 1) +
          theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            aspect.ratio = 1.5,
            panel.border = element_blank(),
            axis.line = element_line(colour = "black")
          ) +
          labs(
            x = NULL, 
            y = "Kyn:Trp ratio"
          ) +
          # add sample numbers to labels
          scale_x_discrete(
            labels = . |>
              filter(!is.na(label)) |>
              distinct(Compound_Name, Event_Name, label) |>
              select(Event_Name, label) |>
              transmute(as.character(Event_Name), as.character(label)) |>
              deframe()
          )
      }
  ) |> 
  patchwork::wrap_plots() +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(
    title = "Endpoint Metabolites: Baseline vs. Treatment",
    subtitle = "nb: smaller N at 40 weeks"
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_KynTrpRatio_B_2_8_16_40_weeks", ".pdf")), device = cairo_pdf, width = 15, height = 5, units = "in")
#

## 4.3 Comparison to HTP -----
### 4.3.1 Sina plot: Ages 12-40 ----
# matching age range to clinical trial
htp_KynTrp_ratio_SexAgeSource_adj |> 
  inner_join(htp_visit_meta_data) |> 
  inner_join(htp_participant_meta_data) |> 
  filter(dplyr::between(Age_years_at_visit, 12, 40)) |> # keeps 294 rows
  mutate(extreme = rstatix::is_extreme(Value_adj), .by = Karyotype) |>
  filter(extreme == FALSE) |> # drops 3 rows
  ggplot(aes(Karyotype, Value_adj, color = Karyotype)) +
  scale_color_manual(values = c("Control" = "grey60", "T21" = "#009b4e")) +
  geom_sina(maxwidth = 0.5) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    aspect.ratio = 2.5,
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  labs(
    title = "HTP Kyn:Trp ratio",
    subtitle = "ages 12-40; extreme outliers removed",
    x = NULL
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_KynTrp_ratio_P4C_SexAgeSourceadj_12-40", ".pdf")), device = cairo_pdf, width = 10, height = 5, units = "in")
#
### 4.3.2 Unpaired Wilcox / Mann-Whitney ----
# matching age range to clinical trial
htp_KynTrp_ratio_SexAgeSource_adj |> 
  inner_join(htp_visit_meta_data) |> 
  inner_join(htp_participant_meta_data) |> 
  filter(dplyr::between(Age_years_at_visit, 12, 40)) |> # keeps 294 rows
  mutate(extreme = rstatix::is_extreme(Value_adj), .by = Karyotype) |>
  filter(extreme == FALSE) |> # drops 3 rows
  rstatix::wilcox_test(
    formula = Value_adj ~ Karyotype,
    ref.group = "T21",
    paired = FALSE,
    detailed = TRUE,
    p.adjust.method = "none"
  )
#


# 5 Sina Plots - Differences  --------
## 5.1 Baseline vs 16 weeks ----
tofa_metab_diffs_16_weeks |> 
  mutate(timepoint = "16 week") |> 
  filter(Compound_Name %in% c("kynurenine", "quinolinic acid", "KynTrp_ratio")) |> 
  mutate(Compound_Name = fct_relevel(Compound_Name, c("kynurenine", "quinolinic acid", "KynTrp_ratio"))) |> 
  ggplot(aes(timepoint, difference, color = timepoint)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_sina(maxwidth = 0.4) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  scale_color_manual(values = standard_colors) +
  facet_wrap(~Compound_Name, nrow = 1, scales = "free_y") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    aspect.ratio = 3
  ) +
  labs(
    title = "Plasma metab AQ: Differences vs Baseline"
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_diffs_16_weeks", ".pdf")), device = cairo_pdf, width = 20, height = 5, units = "in")
#
## 5.2 Baseline vs 2/8/16/40 weeks ----
bind_rows(
  tofa_metab_diffs_2_weeks |> mutate(timepoint = "2 week"),
  tofa_metab_diffs_8_weeks |> mutate(timepoint = "8 week"),
  tofa_metab_diffs_16_weeks |> mutate(timepoint = "16 week"),
  tofa_metab_diffs_40_weeks |> mutate(timepoint = "40 week")
) |> 
  mutate(timepoint = fct_relevel(timepoint, c("2 week", "8 week", "16 week", "40 week"))) |> 
  filter(Compound_Name %in% c("kynurenine", "quinolinic acid", "KynTrp_ratio")) |> 
  mutate(Compound_Name = fct_relevel(Compound_Name, c("kynurenine", "quinolinic acid", "KynTrp_ratio"))) |> 
  ggplot(aes(timepoint, difference, color = timepoint)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(aes(group = ParticipantID), color = "grey90") +
  geom_sina(maxwidth = 0.4) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  scale_color_manual(values = standard_colors) +
  facet_wrap(~Compound_Name, nrow = 1, scales = "free_y") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    aspect.ratio = 1.3
  ) +
  labs(
    title = "Plasma metab AQ:\nDistributions of differences vs Baseline",
    subtitle = "NB: smaller N at 40 weeks",
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_diffs_2_8_16_40_weeks", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#


# 6 AQ metabolites - Mixed effects linear regression models (non-stratified) ----
# Time as a categorical variable
#   No assumptions about trajectory
#   Directly compares each visit to baseline
# Random effects
#   Handles missing at random (MAR) data.
#   Accounts for within-subject correlation.
## 6.1 Set up models ----
### 6.1.1 Mixed effects LM: Event_name + 1|ParticipantID ----
tofa_plasma_metab_AQ_lm_mixedParticipantID <- tofa_plasma_metab_AQ_data_final |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Value)) |>
  filter(is.finite(Value)) |>
  select(Compound_Name, ParticipantID, Sex, Age_Baseline, Event_Name, Value) |> 
  nest(data = -Compound_Name) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Value ~ Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
### 6.1.2 Mixed effects LM: Sex + Event_name + 1|ParticipantID ----
tofa_plasma_metab_AQ_lm_fixedSex_mixedParticipantID <- tofa_plasma_metab_AQ_data_final |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Value)) |>
  filter(is.finite(Value)) |>
  select(Compound_Name, ParticipantID, Sex, Age_Baseline, Event_Name, Value) |> 
  nest(data = -Compound_Name) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Value ~ Sex + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
### 6.1.3 Mixed effects LM: Age + Event_name + 1|ParticipantID ----
tofa_plasma_metab_AQ_lm_fixedAge_mixedParticipantID <- tofa_plasma_metab_AQ_data_final |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Value)) |>
  filter(is.finite(Value)) |>
  select(Compound_Name, ParticipantID, Sex, Age_Baseline, Event_Name, Value) |> 
  nest(data = -Compound_Name) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Value ~ Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
### 6.1.4 Mixed effects LM: Sex + Age + Event_name + 1|ParticipantID ----
tofa_plasma_metab_AQ_lm_fixedSexAge_mixedParticipantID <- tofa_plasma_metab_AQ_data_final |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Value)) |>
  filter(is.finite(Value)) |>
  select(Compound_Name, ParticipantID, Sex, Age_Baseline, Event_Name, Value) |> 
  nest(data = -Compound_Name) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Value ~ Sex + Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#

## 6.2 Compare models ----
### 6.2.1 AIC/BIC ----
# Mixed models
tofa_plasma_metab_AQ_lm_mixedParticipantID |> unnest(glanced) |> select(Compound_Name, AIC, BIC)
tofa_plasma_metab_AQ_lm_fixedSex_mixedParticipantID |> unnest(glanced) |> select(Compound_Name, AIC, BIC)
tofa_plasma_metab_AQ_lm_fixedAge_mixedParticipantID |> unnest(glanced) |> select(Compound_Name, AIC, BIC)
tofa_plasma_metab_AQ_lm_fixedSexAge_mixedParticipantID |> unnest(glanced) |> select(Compound_Name, AIC, BIC)
#
### 6.2.2 Likelihood ratio tests ----
anova(
  tofa_plasma_metab_AQ_lm_mixedParticipantID |> pluck("fit", 1),
  tofa_plasma_metab_AQ_lm_fixedSex_mixedParticipantID |> pluck("fit", 1),
  test="LRT"
) |>  # refitting model(s) with ML (instead of REML)
  tidy() # p.value = 0.103
anova(
  tofa_plasma_metab_AQ_lm_mixedParticipantID |> pluck("fit", 1),
  tofa_plasma_metab_AQ_lm_fixedAge_mixedParticipantID |> pluck("fit", 1),
  test="LRT"
) |>  # refitting model(s) with ML (instead of REML)
  tidy() # p.value = 0.259
anova(
  tofa_plasma_metab_AQ_lm_mixedParticipantID |> pluck("fit", 1),
  tofa_plasma_metab_AQ_lm_fixedSexAge_mixedParticipantID |> pluck("fit", 1),
  test="LRT"
) |>  # refitting model(s) with ML (instead of REML)
  tidy() # p.value = 0.139
#

## 6.3 Model results ----
### 6.3.2 Mixed ParticipantID with fixed SexAge -----
tofa_plasma_metab_AQ_lm_fixedSexAge_mixedParticipantID |> 
  unnest(tidied) |> 
  select(Compound_Name, group, term, estimate, p.value)
# 
tofa_plasma_metab_AQ_lm_fixedSexAge_mixedParticipantID_results <- tofa_plasma_metab_AQ_lm_fixedSexAge_mixedParticipantID |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, "Event_Name")) |>
  select(Compound_Name, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = Compound_Name) |> # NB: by Compound_Name
  mutate(
    level = "",
    term = str_remove(term, "^Event_Name"),
    .after = Compound_Name,
  )
#


# 7 AQ metabolites -  Mixed effects linear regression models (stratified) ----
## 7.1 By Sex ----
lm_fixedAge_mixedParticipantID_by_Sex <- tofa_plasma_metab_AQ_data_final |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Value)) |>
  filter(is.finite(Value)) |>
  select(Compound_Name, ParticipantID, Sex, Age_Baseline, Event_Name, Value) |> 
  nest(data = -c(Compound_Name, Sex)) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Value ~ Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
lm_fixedAge_mixedParticipantID_by_Sex |> unnest(glanced) |> select(Compound_Name, Sex, AIC, BIC)
#
lm_fixedAge_mixedParticipantID_by_Sex_results <- lm_fixedAge_mixedParticipantID_by_Sex |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, "Event_Name")) |>
  select(Compound_Name, level = Sex, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = c(Compound_Name, level)) |> # NB: by Compound_Name + level
  mutate(term = str_remove(term, "^Event_Name")) |> 
  arrange(level)
#

## 7.2 By Age group ----
tofa_baseline_age_groups <- tofa_visit_meta_data |> 
  filter(Event_Name == "Baseline") |> 
  mutate(
    Age_group = case_when(
      Age_years_at_visit < 18 ~ "under18",
      Age_years_at_visit >= 18 ~ "18+"
    ),
    Age_group = fct_relevel(Age_group, c("under18", "18+"))
  ) |> 
  select(ParticipantID, Age_group)
tofa_baseline_age_groups |> count(Age_group)
#
lm_fixedSex_mixedParticipantID_by_Age_group <- tofa_plasma_metab_AQ_data_final |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  inner_join(tofa_baseline_age_groups) |> # add age groups info
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Value)) |>
  filter(is.finite(Value)) |>
  select(Compound_Name, ParticipantID, Sex, Age_group, Event_Name, Value) |> 
  nest(data = -c(Compound_Name, Age_group)) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Value ~ Sex + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
lm_fixedSex_mixedParticipantID_by_Age_group |> unnest(glanced) |> select(Compound_Name, Age_group, AIC, BIC) # most improved
#
lm_fixedSex_mixedParticipantID_by_Age_group_results <- lm_fixedSex_mixedParticipantID_by_Age_group |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, "Event_Name")) |>
  select(Compound_Name, level = Age_group, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = c(Compound_Name, level)) |> # NB: by Compound_Name + level
  mutate(term = str_remove(term, "^Event_Name")) |> 
  arrange(level)
#

## 7.3 By Obesity ----
tofa_baseline_obesity <- tofa_baseline_obesity_file |> 
  read_tsv()
tofa_baseline_obesity |> count(baseline_obesity_status)
#
lm_fixedSexAge_mixedParticipantID_by_Obesity <- tofa_plasma_metab_AQ_data_final |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  inner_join(tofa_baseline_obesity) |> # add Obesity info
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Value)) |>
  filter(is.finite(Value)) |>
  select(Compound_Name, ParticipantID, Sex, Age_Baseline, baseline_obesity_status, Event_Name, Value) |> 
  nest(data = -c(Compound_Name, baseline_obesity_status)) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Value ~ Sex + Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
lm_fixedSexAge_mixedParticipantID_by_Obesity |> unnest(glanced) |> select(Compound_Name, baseline_obesity_status, AIC, BIC) # most improved
#
lm_fixedSexAge_mixedParticipantID_by_Obesity_results <- lm_fixedSexAge_mixedParticipantID_by_Obesity |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, "Event_Name")) |>
  select(Compound_Name, level = baseline_obesity_status, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = c(Compound_Name, level)) |> # NB: by Compound_Name + level
  mutate(term = str_remove(term, "^Event_Name")) |> 
  arrange(level)
#

## 7.4 By COVID ----
# Slightly different models due to testing with different variable at 16 weeks vs 40 weeks
tofa_covid_history <- tofa_covid_history_file |> 
  read_tsv() |> 
  mutate(COVID_event_hx = fct_relevel(COVID_event_hx, c("no", "yes")))
tofa_covid_history |> filter(Event_Name == "2 week") |> 
  count(COVID_event_hx) # 0 participants had covid by 2 weeks - cannot test
tofa_covid_history |> filter(Event_Name == "8 week") |> 
  count(COVID_event_hx) # 2 participants had covid by 8 weeks - n too low to test
tofa_covid_history |> filter(Event_Name == "16 week") |> 
  count(COVID_event_hx) # 5 participants had covid by 2 weeks - small n but will test
tofa_covid_history |> filter(Event_Name == "40 week") |> 
  count(COVID_event_hx) # 10 participants had covid by 2 weeks - small n but will test
#
lm_fixedSexAge_mixedParticipantID_by_COVID <- bind_rows( # NEEDS 2 DIFFERENT MODELS FOR 16/40 WEEKS
  tofa_plasma_metab_AQ_data_final |> 
    inner_join(tofa_participant_meta_data) |> 
    inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
    inner_join(tofa_covid_history |> filter(Event_Name == "16 week") |> select(ParticipantID, COVID_event_hx)) |> # add covid info for 16 weeks
    filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
    filter(!is.na(Value)) |>
    filter(is.finite(Value)) |>
    select(Compound_Name, ParticipantID, Sex, Age_Baseline, COVID_event_hx, Event_Name, Value) |> 
    nest(data = -c(Compound_Name, COVID_event_hx)) |> 
    mutate(
      test = "16 week",
      fit = map(data, possibly(\(.x) lmerTest::lmer(Value ~ Sex + Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
      tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
      glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
      augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
    ),
  tofa_plasma_metab_AQ_data_final |> 
    inner_join(tofa_participant_meta_data) |> 
    inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
    inner_join(tofa_covid_history |> filter(Event_Name == "40 week") |> select(ParticipantID, COVID_event_hx)) |> # add covid info for 40 weeks
    filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
    filter(!is.na(Value)) |>
    filter(is.finite(Value)) |>
    select(Compound_Name, ParticipantID, Sex, Age_Baseline, COVID_event_hx, Event_Name, Value) |> 
    nest(data = -c(Compound_Name, COVID_event_hx)) |> 
    mutate(
      test = "40 week",
      fit = map(data, possibly(\(.x) lmerTest::lmer(Value ~ Sex + Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
      tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
      glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
      augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
    )
) |> 
  arrange(COVID_event_hx)
#
lm_fixedSexAge_mixedParticipantID_by_COVID |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, regex(test))) |> # keep only relevant tests
  select(Compound_Name, COVID_event_hx, term, AIC, BIC) # some improvements
#
lm_fixedSexAge_mixedParticipantID_by_COVID_results <- lm_fixedSexAge_mixedParticipantID_by_COVID |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, regex(test))) |> # keep only relevant tests
  filter(str_detect(term, "Event_Name")) |>
  select(Compound_Name, level = COVID_event_hx, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = c(Compound_Name, level)) |> # NB: by Compound_Name + level
  mutate(term = str_remove(term, "^Event_Name")) |> 
  arrange(level) |> 
  mutate(level = fct_recode(level, No = "no", Yes = "yes")) # relabel COVID levels
#


# 8 Export LM results ----
list(
  "LMM results" = list(
    "Overall" = tofa_plasma_metab_AQ_lm_fixedSexAge_mixedParticipantID_results,
    "Sex" = lm_fixedAge_mixedParticipantID_by_Sex_results,
    "Age_group" = lm_fixedSex_mixedParticipantID_by_Age_group_results,
    "Obesity" = lm_fixedSexAge_mixedParticipantID_by_Obesity_results,
    "COVID" = lm_fixedSexAge_mixedParticipantID_by_COVID_results
  ) |> 
    bind_rows(.id = "stratifier") |> 
    select(
      Stratifier = stratifier,
      Level = level,
      Compound_Name,
      Timepoint = term,
      n_observations = n_obs,
      Mean_difference = estimate,
      Conf.low = conf.low,
      Conf.high = conf.high,
      Statistic = statistic,
      pvalue = p.value,
      qvalue = BH_padj
    ) |> 
    mutate(
      Method = "Mixed effects linear model"
    )
) |> 
  export_excel(filename = "LMM_combined_results")
#


# 9 Forest plot(s) ----
## 9.2 Kynurenine ----
list(
  "Overall" = tofa_plasma_metab_AQ_lm_fixedSexAge_mixedParticipantID_results,
  "Sex" = lm_fixedAge_mixedParticipantID_by_Sex_results,
  "Age_group" = lm_fixedSex_mixedParticipantID_by_Age_group_results,
  "Obesity" = lm_fixedSexAge_mixedParticipantID_by_Obesity_results,
  "COVID" = lm_fixedSexAge_mixedParticipantID_by_COVID_results
) |> 
  bind_rows(.id = "stratifier") |> 
  mutate(stratifier = fct_relevel(stratifier, c("Overall", "Age_group", "Sex", "Obesity", "COVID"))) |> 
  arrange(stratifier) |> 
  filter(Compound_Name == "kynurenine") |> 
  mutate(comparison = paste0(stratifier, "|", level, "|", term)) |> 
  # Add blank separator rows:
  mutate(.row_id = row_number()) |> 
  group_by(stratifier) %>%
  group_modify(~ bind_rows(
    .x,
    tibble(
      .row_id   = min(.x$.row_id) - 0.5, # Subtract so blank rows are BEFORE each group when y-axis reversed
      comparison = paste0(unique(.y$stratifier), "                    "), # .y stores the grouping variable
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      BH_padj = NA_real_,
      color = NA_character_
    )
  )) |>
  ungroup() |>
  arrange(.row_id) |> 
  #
  mutate(comparison = fct_inorder(comparison)) |> 
  mutate(color = if_else(BH_padj < 0.1, "q < 0.1", "n.s.")) |> 
  ggplot(aes(estimate, comparison)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey 50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, height = 0.2)) +
  geom_point(shape = 15, aes(size = -log10(BH_padj), color = color)) +
  scale_color_manual(values = c("q < 0.1" = "red", "n.s." = "black")) +
  scale_size_continuous(range = c(1, 3), limits = c(0, 10)) + # set limits across mutliple endpoints
  scale_y_discrete(limits = rev) +
  theme(
    aspect.ratio = 3,
    legend.key = element_blank()
  ) +
  labs(
    title = "KYN: Treatment vs. Baseline",
    subtitle = "LMM",
    x = "Estimate",
    y = NULL
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "Forest_plot_LMMs_KYN", ".pdf")), device = cairo_pdf, width = 8, height = 5, units = "in")
#
## 9.2 Quinolinic acid ----
list(
  "Overall" = tofa_plasma_metab_AQ_lm_fixedSexAge_mixedParticipantID_results,
  "Sex" = lm_fixedAge_mixedParticipantID_by_Sex_results,
  "Age_group" = lm_fixedSex_mixedParticipantID_by_Age_group_results,
  "Obesity" = lm_fixedSexAge_mixedParticipantID_by_Obesity_results,
  "COVID" = lm_fixedSexAge_mixedParticipantID_by_COVID_results
) |> 
  bind_rows(.id = "stratifier") |> 
  mutate(stratifier = fct_relevel(stratifier, c("Overall", "Age_group", "Sex", "Obesity", "COVID"))) |> 
  arrange(stratifier) |> 
  filter(Compound_Name == "quinolinic acid") |> 
  mutate(comparison = paste0(stratifier, "|", level, "|", term)) |> 
  # Add blank separator rows:
  mutate(.row_id = row_number()) |> 
  group_by(stratifier) %>%
  group_modify(~ bind_rows(
    .x,
    tibble(
      .row_id   = min(.x$.row_id) - 0.5, # Subtract so blank rows are BEFORE each group when y-axis reversed
      comparison = paste0(unique(.y$stratifier), "                    "), # .y stores the grouping variable
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      BH_padj = NA_real_,
      color = NA_character_
    )
  )) |>
  ungroup() |>
  arrange(.row_id) |> 
  #
  mutate(comparison = fct_inorder(comparison)) |> 
  mutate(color = if_else(BH_padj < 0.1, "q < 0.1", "n.s.")) |> 
  ggplot(aes(estimate, comparison)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey 50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, height = 0.2)) +
  geom_point(shape = 15, aes(size = -log10(BH_padj), color = color)) +
  scale_color_manual(values = c("q < 0.1" = "red", "n.s." = "black")) +
  scale_size_continuous(range = c(1, 3), limits = c(0, 10)) + # set limits across mutliple endpoints
  scale_y_discrete(limits = rev) +
  theme(
    aspect.ratio = 3,
    legend.key = element_blank()
  ) +
  labs(
    title = "QA: Treatment vs. Baseline",
    subtitle = "LMM",
    x = "Estimate",
    y = NULL
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "Forest_plot_LMMs_QA", ".pdf")), device = cairo_pdf, width = 8, height = 5, units = "in")
#
## 9.3 KynTrp_ratio ----
list(
  "Overall" = tofa_plasma_metab_AQ_lm_fixedSexAge_mixedParticipantID_results,
  "Sex" = lm_fixedAge_mixedParticipantID_by_Sex_results,
  "Age_group" = lm_fixedSex_mixedParticipantID_by_Age_group_results,
  "Obesity" = lm_fixedSexAge_mixedParticipantID_by_Obesity_results,
  "COVID" = lm_fixedSexAge_mixedParticipantID_by_COVID_results
) |> 
  bind_rows(.id = "stratifier") |> 
  mutate(stratifier = fct_relevel(stratifier, c("Overall", "Age_group", "Sex", "Obesity", "COVID"))) |> 
  arrange(stratifier) |> 
  filter(Compound_Name == "KynTrp_ratio") |> 
  mutate(comparison = paste0(stratifier, "|", level, "|", term)) |> 
  # Add blank separator rows:
  mutate(.row_id = row_number()) |> 
  group_by(stratifier) %>%
  group_modify(~ bind_rows(
    .x,
    tibble(
      .row_id   = min(.x$.row_id) - 0.5, # Subtract so blank rows are BEFORE each group when y-axis reversed
      comparison = paste0(unique(.y$stratifier), "                    "), # .y stores the grouping variable
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      BH_padj = NA_real_,
      color = NA_character_
    )
  )) |>
  ungroup() |>
  arrange(.row_id) |> 
  #
  mutate(comparison = fct_inorder(comparison)) |> 
  mutate(color = if_else(BH_padj < 0.1, "q < 0.1", "n.s.")) |> 
  ggplot(aes(estimate, comparison)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey 50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, height = 0.2)) +
  geom_point(shape = 15, aes(size = -log10(BH_padj), color = color)) +
  scale_color_manual(values = c("q < 0.1" = "red", "n.s." = "black")) +
  scale_size_continuous(range = c(1, 3), limits = c(0, 10)) + # set limits across mutliple endpoints
  scale_y_discrete(limits = rev) +
  theme(
    aspect.ratio = 3,
    legend.key = element_blank()
  ) +
  labs(
    title = "Kyn/Trp: Treatment vs. Baseline",
    subtitle = "LMM",
    x = "Estimate",
    y = NULL
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "Forest_plot_LMMs_KynTrp_ratio", ".pdf")), device = cairo_pdf, width = 8, height = 5, units = "in")
#
## 9.3 Tryptophan ----
list(
  "Overall" = tofa_plasma_metab_AQ_lm_fixedSexAge_mixedParticipantID_results,
  "Sex" = lm_fixedAge_mixedParticipantID_by_Sex_results,
  "Age_group" = lm_fixedSex_mixedParticipantID_by_Age_group_results,
  "Obesity" = lm_fixedSexAge_mixedParticipantID_by_Obesity_results,
  "COVID" = lm_fixedSexAge_mixedParticipantID_by_COVID_results
) |> 
  bind_rows(.id = "stratifier") |> 
  mutate(stratifier = fct_relevel(stratifier, c("Overall", "Age_group", "Sex", "Obesity", "COVID"))) |> 
  arrange(stratifier) |> 
  filter(Compound_Name == "L-tryptophan") |> 
  mutate(comparison = paste0(stratifier, "|", level, "|", term)) |> 
  # Add blank separator rows:
  mutate(.row_id = row_number()) |> 
  group_by(stratifier) %>%
  group_modify(~ bind_rows(
    .x,
    tibble(
      .row_id   = min(.x$.row_id) - 0.5, # Subtract so blank rows are BEFORE each group when y-axis reversed
      comparison = paste0(unique(.y$stratifier), "                    "), # .y stores the grouping variable
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      BH_padj = NA_real_,
      color = NA_character_
    )
  )) |>
  ungroup() |>
  arrange(.row_id) |> 
  #
  mutate(comparison = fct_inorder(comparison)) |> 
  mutate(color = if_else(BH_padj < 0.1, "q < 0.1", "n.s.")) |> 
  ggplot(aes(estimate, comparison)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey 50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, height = 0.2)) +
  geom_point(shape = 15, aes(size = -log10(BH_padj), color = color)) +
  scale_color_manual(values = c("q < 0.1" = "red", "n.s." = "black")) +
  scale_size_continuous(range = c(1, 3), limits = c(0, 10)) + # set limits across mutliple endpoints
  scale_y_discrete(limits = rev) +
  theme(
    aspect.ratio = 3,
    legend.key = element_blank()
  ) +
  labs(
    title = "TRP: Treatment vs. Baseline",
    subtitle = "LMM",
    x = "Estimate",
    y = NULL
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "Forest_plot_LMMs_TRP", ".pdf")), device = cairo_pdf, width = 8, height = 5, units = "in")
#


################################################
# save workspace ----
save.image(file = here("rdata", paste0(out_file_prefix, ".RData")), compress = TRUE, safe = TRUE) # saves entire workspace (can be slow)
# To reload previously saved workspace:
# load(here("rdata", paste0(out_file_prefix, ".RData")))

# session_info ----
date()
sessionInfo()
################################################