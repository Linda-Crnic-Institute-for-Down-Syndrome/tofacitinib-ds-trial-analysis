################################################
# Title: Tofacitinib for Immune Skin Conditions in Down Syndrome: Analysis of Clinical Safety Labs
# Author(s):
#   - Matt Galbraith & Micah Donovan
# Affiliation(s):
#   - Linda Crnic Institute for Down Syndrome
#   - University of Colorado Anschutz
################################################

### Summary:  
# This workflow analyzes clinical safety labs acquired by non-fasting Compete Blood Count (CBC) 
# with diferential (diff), Comprehensive Metabolic Panel (CMP), Lipid panel, and select other values. 
# See README.md for more details.
#  

### Data type(s):
# Clinical trial (TOFA) datasets:
#    * Participant-level metadata; Available on request.
#    * Visit/Event-level metadata; Available on request.
#    * Baseline obesity status; Available on request.
#    * COVID-19 history; Available on request.
#    * Clinical safety labs; DOI: 10.5281/zenodo.20074415

# 0 General Setup -----
# RUN THIS FIRST TIME - Initialize and install packages with renv:
# renv::init(bioconductor = TRUE)
#
# To install the exact versions of all R packages base on renv.lock file (requires matching R version):
# renv::restore()

## 0.1 Load required libraries ----
library("readxl") # Used to read .xlsx files
library("openxlsx") # used for data export as Excel workbooks
library("tidyverse")
library("ggrepel") # required for labelling points
library("ggforce") # required for zooming and sina
library("skimr")
library("janitor")
library("rstatix")
library("coin")
library("conflicted")
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("count", "dplyr")
library("here") # generates path to current project directory
#

## 0.2 Set input files and other parameters ----
#
# Input data files
# Datasets used in this study can be obtained from the associated tofacitinib-ds-trial-analysis
# repositories (further details in README.md).
# Download each dataset to /data directory within this R project.
#
# Clinical trial datasets:
tofa_participant_meta_data_file <- here("data", "TOFA_Participant_metadata_zenodo_v1.txt") 
tofa_visit_meta_data_file <- here("data", "TOFA_Visit_metadata_zenodo_v1.txt") 
tofa_clinical_labs_data_file <- here("data", "TOFA_Clinical_Labs_Data_zenodo_v1.txt")
tofa_baseline_obesity_status_file <- here("data", "TOFA_Baseline_Obesity_Status_zenodo_v1.txt")
tofa_covid_event_history_file <- here("data", "TOFA_COVID_History_zenodo_v1.txt")
# 
# Other parameters:
out_file_prefix <- "TOFA_Clinical_Labs_" # should match this script title
standard_colors <- c("Baseline" = "#999999", "2 week" = "#c6dbef", "8 week" = "#9ecae1", "16 week" = "#6baed6", "40 week" = "#4292c6")
standard_colors2 <- c("2 weeks" = "#c6dbef", "8 weeks" = "#9ecae1", "16 weeks" = "#6baed6", "40 weeks" = "#4292c6")
# End required parameters ###
source(here("helper_functions.R")) # load helper functions
#

# 1.0 Read in data     ----
## 1.1 Participant level meta data    ----
tofa_participant_meta_data <- tofa_participant_meta_data_file %>% 
  read_tsv() %>% 
  mutate(
    Sex = fct_relevel(Sex, c("Female", "Male")) # set factor levels
  )
#
tofa_participant_meta_data # 47 rows
tofa_participant_meta_data %>% distinct(ParticipantID) # 47 participants in this table
tofa_participant_meta_data %>% count(Safety_eligible, Participant_notes) # 4 participants not eligible for Safety analyses
#

## 1.2 Visit level meta data ----
tofa_visit_meta_data <- tofa_visit_meta_data_file %>% 
  read_tsv() %>% 
  mutate(
    Event_Name = fct_relevel(Event_Name, c("Baseline", "2 week", "8 week", "16 week","40 week")) # set factor levels
  )
#

## 1.3 Participant baseline obesity status information   ----
# File listing participants as obese or non-obese at baseline
tofa_baseline_obesity_status <- tofa_baseline_obesity_status_file %>% 
  read_tsv() %>% 
  inner_join(tofa_participant_meta_data)
#
tofa_baseline_obesity_status # 43 rows
tofa_baseline_obesity_status %>% distinct(ParticipantID) # 43 participants in this table
#

## 1.4 Participant COVID-19 event history information   ----
# File listing whether a covid event was recorded for each between baseline and study time points
tofa_covid_event_history <- tofa_covid_event_history_file %>% 
  read_tsv() %>% 
  inner_join(tofa_participant_meta_data) %>% 
  filter(Safety_eligible == TRUE)
#
tofa_covid_event_history # 172 rows
tofa_baseline_obesity_status %>% distinct(ParticipantID) # 43 participants in this table
#

## 1.5 Join clinical safety labs to meta-data   ----
tofa_clinical_labs_data <- tofa_clinical_labs_data_file %>% 
  read_tsv() %>% 
  inner_join(tofa_visit_meta_data) %>% 
  inner_join(tofa_participant_meta_data)
#
tofa_clinical_labs_data %>% distinct(ParticipantID) # 43 participants in this table
tofa_clinical_labs_data %>% distinct(lab_name) # 46 laboratory names in this table
#

# 2 Stats ----
## 2.1 Differences from baseline ----
### 2.1.1 Calculate differences from baseline   ----
tofa_clinical_labs_diffs_2_weeks <- tofa_clinical_labs_data %>% 
  # NEW v0.9: keep TOFA0055
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "2 week")) %>% 
  filter(!is.na(lab_value)) %>% # drops 
  filter(is.finite(lab_value)) %>% 
  add_count(lab_name, ParticipantID, name = "pair") %>% 
  filter(pair == 2) %>%  # require pairs
  add_count(Event_Name, lab_name) %>% 
  group_by(lab_name) %>% 
  filter(!any(n < 2)) %>% # remove labs with n < X in any group
  select(ParticipantID, Event_Name, lab_name, lab_value) %>% 
  pivot_wider(names_from = Event_Name, values_from = lab_value) %>% 
  mutate(difference = `2 week` - Baseline) %>% 
  select(ParticipantID, lab_name, difference) %>% 
  group_by(lab_name) %>% 
  mutate(
    extreme = rstatix::is_extreme(difference),
    outlier = rstatix::is_outlier(difference)
  ) %>% 
  ungroup()
#
tofa_clinical_labs_diffs_8_weeks <- tofa_clinical_labs_data %>% 
  # NEW v0.9: keep TOFA0055
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "8 week")) %>% 
  filter(!is.na(lab_value)) %>% # drops 
  filter(is.finite(lab_value)) %>% 
  add_count(lab_name, ParticipantID, name = "pair") %>% 
  filter(pair == 2) %>%  # require pairs
  add_count(Event_Name, lab_name) %>% 
  group_by(lab_name) %>% 
  filter(!any(n < 2)) %>% # remove labs with n < X in any group
  select(ParticipantID, Event_Name, lab_name, lab_value) %>% 
  pivot_wider(names_from = Event_Name, values_from = lab_value) %>% 
  mutate(difference = `8 week` - Baseline) %>% 
  select(ParticipantID, lab_name, difference) %>% 
  group_by(lab_name) %>% 
  mutate(
    extreme = rstatix::is_extreme(difference),
    outlier = rstatix::is_outlier(difference)
  ) %>% 
  ungroup()
#
tofa_clinical_labs_diffs_16_weeks <- tofa_clinical_labs_data %>% 
  # NEW v0.9: keep TOFA0055
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "16 week")) %>% 
  filter(!is.na(lab_value)) %>% 
  filter(is.finite(lab_value)) %>% 
  add_count(lab_name, ParticipantID, name = "pair") %>% 
  filter(pair == 2) %>%  # require pairs
  add_count(Event_Name, lab_name) %>% 
  group_by(lab_name) %>% 
  filter(!any(n < 2)) %>% # remove labs with n < X in any group
  select(ParticipantID, Event_Name, lab_name, lab_value) %>% 
  pivot_wider(names_from = Event_Name, values_from = lab_value) %>% 
  mutate(difference = `16 week` - Baseline) %>% 
  select(ParticipantID, lab_name, difference) %>% 
  group_by(lab_name) %>% 
  mutate(
    extreme = rstatix::is_extreme(difference),
    outlier = rstatix::is_outlier(difference)
  ) %>% 
  ungroup()
#
tofa_clinical_labs_diffs_40_weeks <- tofa_clinical_labs_data %>% 
  # NEW v0.9: keep TOFA0055
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "40 week")) %>% 
  filter(!is.na(lab_value)) %>% 
  filter(is.finite(lab_value)) %>% 
  add_count(lab_name, ParticipantID, name = "pair") %>% 
  filter(pair == 2) %>%  # require pairs
  add_count(Event_Name, lab_name) %>% 
  group_by(lab_name) %>% 
  filter(!any(n < 2)) %>% # remove labs with n < X in any group
  select(ParticipantID, Event_Name, lab_name, lab_value) %>% 
  pivot_wider(names_from = Event_Name, values_from = lab_value) %>% 
  mutate(difference = `40 week` - Baseline) %>% 
  select(ParticipantID, lab_name, difference) %>% 
  group_by(lab_name) %>% 
  mutate(
    extreme = rstatix::is_extreme(difference),
    outlier = rstatix::is_outlier(difference)
  ) %>% 
  ungroup()
#

### 2.1.2 Count extreme outliers in differences ----
tofa_clinical_labs_diffs_2_weeks %>% 
  filter(extreme == TRUE) %>% 
  count(lab_name, name = "n_extreme") %>% 
  arrange(-n_extreme)
tofa_clinical_labs_diffs_8_weeks %>% 
  filter(extreme == TRUE) %>% 
  count(lab_name, name = "n_extreme") %>% 
  arrange(-n_extreme)
tofa_clinical_labs_diffs_16_weeks %>% 
  filter(extreme == TRUE) %>% 
  count(lab_name, name = "n_extreme") %>% 
  arrange(-n_extreme)
tofa_clinical_labs_diffs_40_weeks %>% 
  filter(extreme == TRUE) %>% 
  count(lab_name, name = "n_extreme") %>% 
  arrange(-n_extreme)
#

# 3.0 Paired t tests for Overall cohort  --------
## 3.1 Run paired t tests across timepoints   ----
### 3.1.1 2 weeks ----
tofa_clinical_labs_Ttest_res_2weeks <- tofa_clinical_labs_data %>% 
  # NEW v0.9: keep TOFA0055
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "2 week")) %>% 
  filter(!is.na(lab_value)) %>% 
  filter(is.finite(lab_value)) %>% 
  add_count(lab_name, ParticipantID, name = "pair") %>% 
  filter(pair == 2) %>%  # require pairs
  add_count(Event_Name, lab_name) %>% 
  group_by(lab_name) %>% 
  filter(!any(n < 2)) %>% # remove labs with n < X in any group
  rstatix::t_test(
    formula = lab_value ~ Event_Name,
    ref.group = "2 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    # var.equal = TRUE, set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) %>%
  left_join(tofa_clinical_labs_data %>% distinct(lab_name, test_name, lab_units, panel)) %>% 
  group_by(panel) %>% 
  mutate(BH_padj = p.adjust(p)) %>% 
  ungroup() %>% 
  arrange(p)
#

### 3.1.2 8 weeks ----
tofa_clinical_labs_Ttest_res_8weeks <- tofa_clinical_labs_data %>% 
  # NEW v0.9: keep TOFA0055
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "8 week")) %>% 
  filter(!is.na(lab_value)) %>% 
  filter(is.finite(lab_value)) %>% 
  add_count(lab_name, ParticipantID, name = "pair") %>% 
  filter(pair == 2) %>%  # require pairs
  add_count(Event_Name, lab_name) %>% 
  group_by(lab_name) %>% 
  filter(!any(n < 2)) %>% # remove labs with n < X in any group
  rstatix::t_test(
    formula = lab_value ~ Event_Name,
    ref.group = "8 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    # var.equal = TRUE, set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) %>%
  left_join(tofa_clinical_labs_data %>% distinct(lab_name, test_name, lab_units, panel)) %>% 
  group_by(panel) %>% 
  mutate(BH_padj = p.adjust(p)) %>% 
  ungroup() %>% 
  arrange(p)
#

### 3.1.3 16 weeks ----
tofa_clinical_labs_Ttest_res_16weeks <- tofa_clinical_labs_data %>% 
  # NEW v0.9: keep TOFA0055
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "16 week")) %>% 
  filter(!is.na(lab_value)) %>% 
  filter(is.finite(lab_value)) %>% 
  add_count(lab_name, ParticipantID, name = "pair") %>% 
  filter(pair == 2) %>%  # require pairs
  add_count(Event_Name, lab_name) %>% 
  group_by(lab_name) %>% 
  filter(!any(n < 2)) %>% # remove labs with n < X in any group
  rstatix::t_test(
    formula = lab_value ~ Event_Name,
    ref.group = "16 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    # var.equal = TRUE, set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) %>%
  left_join(tofa_clinical_labs_data %>% distinct(lab_name, test_name, lab_units, panel)) %>% 
  group_by(panel) %>% 
  mutate(BH_padj = p.adjust(p)) %>% 
  ungroup() %>% 
  arrange(p)
#

### 3.1.4 40 weeks ----
tofa_clinical_labs_Ttest_res_40weeks <- tofa_clinical_labs_data %>% 
  # NEW v0.9: keep TOFA0055
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "40 week")) %>% 
  filter(!is.na(lab_value)) %>% 
  filter(is.finite(lab_value)) %>% 
  add_count(lab_name, ParticipantID, name = "pair") %>% 
  filter(pair == 2) %>%  # require pairs
  add_count(Event_Name, lab_name) %>% 
  group_by(lab_name) %>% 
  filter(!any(n < 2)) %>% # remove labs with n < X in any group
  rstatix::t_test(
    formula = lab_value ~ Event_Name,
    ref.group = "40 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    # var.equal = TRUE, set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) %>%
  left_join(tofa_clinical_labs_data %>% distinct(lab_name, test_name, lab_units, panel)) %>% 
  group_by(panel) %>% 
  mutate(BH_padj = p.adjust(p)) %>% 
  ungroup() %>% 
  arrange(p)
#

## 3.2 Calculate paired t test effect sizes across time points    ----
# The effect size for a paired-samples t-test can be calculated by dividing the
# mean difference by the standard deviation of the difference, as shown below.
# Cohen’s formula:
# d = mean(D)/sd(D), where D is the differences of the paired samples values.
#
### 3.2.1 2 weeks ----
tofa_clinical_labs_Ttest_effsize_2weeks <- tofa_clinical_labs_data %>% 
  # NEW v0.9: keep TOFA0055
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "2 week")) %>% 
  filter(!is.na(lab_value)) %>% 
  filter(is.finite(lab_value)) %>% 
  add_count(lab_name, ParticipantID, name = "pair") %>% 
  filter(pair == 2) %>%  # require pairs
  add_count(Event_Name, lab_name) %>% 
  group_by(lab_name) %>% 
  filter(!any(n < 2)) %>% # remove labs with n < X in any group
  rstatix::cohens_d(
    formula = lab_value ~ Event_Name,
    ref.group = "2 week", # Seems to be the reverse of how we typically do this
    paired = TRUE
  )
#

### 3.2.2 8 weeks ----
tofa_clinical_labs_Ttest_effsize_8weeks <- tofa_clinical_labs_data %>% 
  # NEW v0.9: keep TOFA0055
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "8 week")) %>% 
  filter(!is.na(lab_value)) %>% 
  filter(is.finite(lab_value)) %>% 
  add_count(lab_name, ParticipantID, name = "pair") %>% 
  filter(pair == 2) %>%  # require pairs
  add_count(Event_Name, lab_name) %>% 
  group_by(lab_name) %>% 
  filter(!any(n < 2)) %>% # remove labs with n < X in any group
  rstatix::cohens_d(
    formula = lab_value ~ Event_Name,
    ref.group = "8 week", # Seems to be the reverse of how we typically do this
    paired = TRUE
  )
#

### 3.2.3 16 weeks ----
tofa_clinical_labs_Ttest_effsize_16weeks <- tofa_clinical_labs_data %>% 
  # NEW v0.9: keep TOFA0055
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "16 week")) %>% 
  filter(!is.na(lab_value)) %>% 
  filter(is.finite(lab_value)) %>% 
  add_count(lab_name, ParticipantID, name = "pair") %>% 
  filter(pair == 2) %>%  # require pairs
  add_count(Event_Name, lab_name) %>% 
  group_by(lab_name) %>% 
  filter(!any(n < 2)) %>% # remove labs with n < X in any group
  rstatix::cohens_d(
    formula = lab_value ~ Event_Name,
    ref.group = "16 week", # Seems to be the reverse of how we typically do this
    paired = TRUE
  )
#

### 3.2.4 40 weeks ----
tofa_clinical_labs_Ttest_effsize_40weeks <- tofa_clinical_labs_data %>% 
  # NEW v0.9: keep TOFA0055
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "40 week")) %>% 
  filter(!is.na(lab_value)) %>% 
  filter(is.finite(lab_value)) %>% 
  add_count(lab_name, ParticipantID, name = "pair") %>% 
  filter(pair == 2) %>%  # require pairs
  add_count(Event_Name, lab_name) %>% 
  group_by(lab_name) %>% 
  filter(!any(n < 2)) %>% # remove labs with n < X in any group
  rstatix::cohens_d(
    formula = lab_value ~ Event_Name,
    ref.group = "40 week", # Seems to be the reverse of how we typically do this
    paired = TRUE
  )
#

## 3.3 Compile t test results ----
# 2 weeks
tofa_clinical_labs_Ttest_res_2weeks_full <- tofa_clinical_labs_Ttest_res_2weeks %>% 
  inner_join(tofa_clinical_labs_Ttest_effsize_2weeks) %>% 
  mutate(
    mean_diff = estimate, 
  ) %>% 
  inner_join(tofa_clinical_labs_diffs_2_weeks %>% group_by(lab_name) %>% summarize(median_diff = median(difference))) %>% 
  select(panel, lab_name, test_name, lab_units, mean_diff, median_diff, p, BH_padj, effsize, magnitude, .y., group1, group2, everything())
#
tofa_clinical_labs_Ttest_res_2weeks_full
tofa_clinical_labs_Ttest_res_2weeks_full %>% filter(p < 0.05)
tofa_clinical_labs_Ttest_res_2weeks_full %>% filter(BH_padj < 0.1)
#
tofa_clinical_labs_Ttest_res_2weeks_full %>% filter(panel == "CBC")
tofa_clinical_labs_Ttest_res_2weeks_full %>% filter(panel == "Differential")
tofa_clinical_labs_Ttest_res_2weeks_full %>% filter(panel == "CMP")
tofa_clinical_labs_Ttest_res_2weeks_full %>% filter(panel == "Lipids (non-fasted)")
tofa_clinical_labs_Ttest_res_2weeks_full %>% filter(panel == "Other")
#

# 8 weeks
tofa_clinical_labs_Ttest_res_8weeks_full <- tofa_clinical_labs_Ttest_res_8weeks %>% 
  inner_join(tofa_clinical_labs_Ttest_effsize_8weeks) %>% 
  mutate(
    mean_diff = estimate, 
  ) %>% 
  inner_join(tofa_clinical_labs_diffs_8_weeks %>% group_by(lab_name) %>% summarize(median_diff = median(difference))) %>% 
  select(panel, lab_name, test_name, lab_units, mean_diff, median_diff, p, BH_padj, effsize, magnitude, .y., group1, group2, everything())
#
tofa_clinical_labs_Ttest_res_8weeks_full
tofa_clinical_labs_Ttest_res_8weeks_full %>% filter(p < 0.05)
tofa_clinical_labs_Ttest_res_8weeks_full %>% filter(BH_padj < 0.1)
#
tofa_clinical_labs_Ttest_res_8weeks_full %>% filter(panel == "CBC")
tofa_clinical_labs_Ttest_res_8weeks_full %>% filter(panel == "Differential")
tofa_clinical_labs_Ttest_res_8weeks_full %>% filter(panel == "CMP")
tofa_clinical_labs_Ttest_res_8weeks_full %>% filter(panel == "Lipids (non-fasted)")
tofa_clinical_labs_Ttest_res_8weeks_full %>% filter(panel == "Other")
#

# 16 weeks
tofa_clinical_labs_Ttest_res_16weeks_full <- tofa_clinical_labs_Ttest_res_16weeks %>% 
  inner_join(tofa_clinical_labs_Ttest_effsize_16weeks) %>% 
  mutate(
    mean_diff = estimate, 
  ) %>% 
  inner_join(tofa_clinical_labs_diffs_16_weeks %>% group_by(lab_name) %>% summarize(median_diff = median(difference))) %>% 
  select(panel, lab_name, test_name, lab_units, mean_diff, median_diff, p, BH_padj, effsize, magnitude, .y., group1, group2, everything())
#
tofa_clinical_labs_Ttest_res_16weeks_full
tofa_clinical_labs_Ttest_res_16weeks_full %>% filter(p < 0.05)
tofa_clinical_labs_Ttest_res_16weeks_full %>% filter(BH_padj < 0.1)
#
tofa_clinical_labs_Ttest_res_16weeks_full %>% filter(panel == "CBC")
tofa_clinical_labs_Ttest_res_16weeks_full %>% filter(panel == "Differential")
tofa_clinical_labs_Ttest_res_16weeks_full %>% filter(panel == "CMP")
tofa_clinical_labs_Ttest_res_16weeks_full %>% filter(panel == "Lipids (non-fasted)")
tofa_clinical_labs_Ttest_res_16weeks_full %>% filter(panel == "Other")
#

# 40 weeks
tofa_clinical_labs_Ttest_res_40weeks_full <- tofa_clinical_labs_Ttest_res_40weeks %>% 
  inner_join(tofa_clinical_labs_Ttest_effsize_40weeks) %>% 
  mutate(
    mean_diff = estimate, 
  ) %>% 
  inner_join(tofa_clinical_labs_diffs_40_weeks %>% group_by(lab_name) %>% summarize(median_diff = median(difference))) %>% 
  select(panel, lab_name, test_name, lab_units, mean_diff, median_diff, p, BH_padj, effsize, magnitude, .y., group1, group2, everything())
#
tofa_clinical_labs_Ttest_res_40weeks_full
tofa_clinical_labs_Ttest_res_40weeks_full %>% filter(p < 0.05)
tofa_clinical_labs_Ttest_res_40weeks_full %>% filter(BH_padj < 0.1)
#
tofa_clinical_labs_Ttest_res_40weeks_full %>% filter(panel == "CBC")
tofa_clinical_labs_Ttest_res_40weeks_full %>% filter(panel == "Differential")
tofa_clinical_labs_Ttest_res_40weeks_full %>% filter(panel == "CMP")
tofa_clinical_labs_Ttest_res_40weeks_full %>% filter(panel == "Lipids (non-fasted)")
tofa_clinical_labs_Ttest_res_40weeks_full %>% filter(panel == "Other")
#

## 3.4 Save results of overall T tests    ----
bind_rows(
  tofa_clinical_labs_Ttest_res_2weeks_full %>% mutate(timepoint = "2 weeks"),
  tofa_clinical_labs_Ttest_res_8weeks_full %>% mutate(timepoint = "8 weeks"),
  tofa_clinical_labs_Ttest_res_16weeks_full %>% mutate(timepoint = "16 weeks"),
  tofa_clinical_labs_Ttest_res_40weeks_full %>% mutate(timepoint = "40 weeks")
) %>% 
  select(timepoint, everything()) %>% 
  write_tsv(file = here("results", paste0(out_file_prefix, "tofa_clinical_labs_Ttest_results_overall.txt")))

# 4.0 Paired t tests stratified by Age group at baseline    ----
baseline_ages <- tofa_visit_meta_data %>% 
  filter(Event_Name == "Baseline") %>% 
  distinct(ParticipantID, Age_years_at_visit) %>% 
  rename(Baseline_Age = Age_years_at_visit) %>% 
  mutate(Baseline_Age_Group = case_when(
    Baseline_Age < 18 ~ "under18",
    Baseline_Age >= 18 ~ "18+"
  ))

## 4.1 Run paired t tests across timepoints   ----
### 4.1.1 2 weeks   ----
tofa_clinical_labs_Ttest_res_2weeks_age <- tofa_clinical_labs_data |> 
  inner_join(baseline_ages %>% select(ParticipantID, Baseline_Age_Group)) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "2 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, Baseline_Age_Group) |> 
  group_by(lab_name, Baseline_Age_Group) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::t_test(
    formula = lab_value ~ Event_Name,
    ref.group = "2 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    # var.equal = TRUE, set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) %>%
  left_join(tofa_clinical_labs_data %>% distinct(lab_name, test_name, lab_units, panel)) %>% 
  group_by(panel) |> 
  mutate(BH_padj = p.adjust(p)) |> 
  ungroup() |> 
  arrange(p)
#

### 4.1.2 8 weeks   ----
tofa_clinical_labs_Ttest_res_8weeks_age <- tofa_clinical_labs_data |> 
  inner_join(baseline_ages %>% select(ParticipantID, Baseline_Age_Group)) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "8 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, Baseline_Age_Group) |> 
  group_by(lab_name, Baseline_Age_Group) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::t_test(
    formula = lab_value ~ Event_Name,
    ref.group = "8 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    # var.equal = TRUE, set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) %>%
  left_join(tofa_clinical_labs_data %>% distinct(lab_name, test_name, lab_units, panel)) %>% 
  group_by(panel) |> 
  mutate(BH_padj = p.adjust(p)) |> 
  ungroup() |> 
  arrange(p)
#

### 4.1.3 16 weeks    ----
tofa_clinical_labs_Ttest_res_16weeks_age <- tofa_clinical_labs_data |> 
  inner_join(baseline_ages %>% select(ParticipantID, Baseline_Age_Group)) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, Baseline_Age_Group) |> 
  group_by(lab_name, Baseline_Age_Group) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::t_test(
    formula = lab_value ~ Event_Name,
    ref.group = "16 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    # var.equal = TRUE, set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) %>%
  left_join(tofa_clinical_labs_data %>% distinct(lab_name, test_name, lab_units, panel)) %>% 
  group_by(panel) |> 
  mutate(BH_padj = p.adjust(p)) |> 
  ungroup() |> 
  arrange(p)
#

### 4.1.4 40 weeks    ----
tofa_clinical_labs_Ttest_res_40weeks_age <- tofa_clinical_labs_data |> 
  inner_join(baseline_ages %>% select(ParticipantID, Baseline_Age_Group)) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "40 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, Baseline_Age_Group) |> 
  group_by(lab_name, Baseline_Age_Group) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::t_test(
    formula = lab_value ~ Event_Name,
    ref.group = "40 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    # var.equal = TRUE, set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) %>%
  left_join(tofa_clinical_labs_data %>% distinct(lab_name, test_name, lab_units, panel)) %>% 
  group_by(panel) |> 
  mutate(BH_padj = p.adjust(p)) |> 
  ungroup() |> 
  arrange(p)
#

## 4.2 Calculate paired t test effect sizes across time points    ----
### 4.2.1 2 weeks   ----
tofa_clinical_labs_Ttest_effsize_2weeks_age <- tofa_clinical_labs_data |> 
  inner_join(baseline_ages %>% select(ParticipantID, Baseline_Age_Group)) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "2 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, Baseline_Age_Group) |> 
  group_by(lab_name, Baseline_Age_Group) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::cohens_d(
    formula = lab_value ~ Event_Name,
    ref.group = "2 week", # Seems to be the reverse of how we typically do this
    paired = TRUE
  )
#

### 4.2.2 8 weeks   ----
tofa_clinical_labs_Ttest_effsize_8weeks_age <- tofa_clinical_labs_data |> 
  inner_join(baseline_ages %>% select(ParticipantID, Baseline_Age_Group)) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "8 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, Baseline_Age_Group) |> 
  group_by(lab_name, Baseline_Age_Group) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::cohens_d(
    formula = lab_value ~ Event_Name,
    ref.group = "8 week", # Seems to be the reverse of how we typically do this
    paired = TRUE
  )
#

### 4.2.3 16 weeks    ----
tofa_clinical_labs_Ttest_effsize_16weeks_age <- tofa_clinical_labs_data |> 
  inner_join(baseline_ages %>% select(ParticipantID, Baseline_Age_Group)) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, Baseline_Age_Group) |> 
  group_by(lab_name, Baseline_Age_Group) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::cohens_d(
    formula = lab_value ~ Event_Name,
    ref.group = "16 week", # Seems to be the reverse of how we typically do this
    paired = TRUE
  )
#

### 4.2.4 40 weeks    ----
tofa_clinical_labs_Ttest_effsize_40weeks_age <- tofa_clinical_labs_data |> 
  inner_join(baseline_ages %>% select(ParticipantID, Baseline_Age_Group)) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "40 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, Baseline_Age_Group) |> 
  group_by(lab_name, Baseline_Age_Group) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::cohens_d(
    formula = lab_value ~ Event_Name,
    ref.group = "40 week", # Seems to be the reverse of how we typically do this
    paired = TRUE
  )
#

## 4.3 Compile t test results   ----
# 2 weeks
tofa_clinical_labs_Ttest_res_2weeks_full_age <- tofa_clinical_labs_Ttest_res_2weeks_age |> 
  inner_join(tofa_clinical_labs_Ttest_effsize_2weeks_age) |> 
  mutate(
    mean_diff = estimate, 
  ) |> 
  inner_join(tofa_clinical_labs_diffs_2_weeks |> group_by(lab_name) |> summarize(median_diff = median(difference))) |> 
  select(panel, lab_name, test_name, lab_units, mean_diff, median_diff, p, BH_padj, effsize, magnitude, .y., group1, group2, everything())
#

# 8 weeks
tofa_clinical_labs_Ttest_res_8weeks_full_age <- tofa_clinical_labs_Ttest_res_8weeks_age |> 
  inner_join(tofa_clinical_labs_Ttest_effsize_8weeks_age) |> 
  mutate(
    mean_diff = estimate, 
  ) |> 
  inner_join(tofa_clinical_labs_diffs_8_weeks |> group_by(lab_name) |> summarize(median_diff = median(difference))) |> 
  select(panel, lab_name, test_name, lab_units, mean_diff, median_diff, p, BH_padj, effsize, magnitude, .y., group1, group2, everything())
#

# 16 weeks
tofa_clinical_labs_Ttest_res_16weeks_full_age <- tofa_clinical_labs_Ttest_res_16weeks_age |> 
  inner_join(tofa_clinical_labs_Ttest_effsize_16weeks_age) |> 
  mutate(
    mean_diff = estimate, 
  ) |> 
  inner_join(tofa_clinical_labs_diffs_16_weeks |> group_by(lab_name) |> summarize(median_diff = median(difference))) |> 
  select(panel, lab_name, test_name, lab_units, mean_diff, median_diff, p, BH_padj, effsize, magnitude, .y., group1, group2, everything())
#

# 40 weeks
tofa_clinical_labs_Ttest_res_40weeks_full_age <- tofa_clinical_labs_Ttest_res_40weeks_age |> 
  inner_join(tofa_clinical_labs_Ttest_effsize_40weeks_age) |> 
  mutate(
    mean_diff = estimate, 
  ) |> 
  inner_join(tofa_clinical_labs_diffs_40_weeks |> group_by(lab_name) |> summarize(median_diff = median(difference))) |> 
  select(panel, lab_name, test_name, lab_units, mean_diff, median_diff, p, BH_padj, effsize, magnitude, .y., group1, group2, everything())
#

## 4.4 Save results of age-stratified T tests    ----
bind_rows(
  tofa_clinical_labs_Ttest_res_2weeks_age %>% mutate(timepoint = "2 weeks"),
  tofa_clinical_labs_Ttest_res_8weeks_age %>% mutate(timepoint = "8 weeks"),
  tofa_clinical_labs_Ttest_res_16weeks_age %>% mutate(timepoint = "16 weeks"),
  tofa_clinical_labs_Ttest_res_40weeks_age %>% mutate(timepoint = "40 weeks")
) %>% 
  select(timepoint, everything()) %>% 
  write_tsv(file = here("results", paste0(out_file_prefix, "tofa_clinical_labs_Ttest_results_age_stratified.txt")))
#

# 5.0 Paired t tests stratified by Sex    ----
## 5.1 Run paired t tests across timepoints   ----
### 5.1.1 2 weeks   ----
tofa_clinical_labs_Ttest_res_2weeks_sex <- tofa_clinical_labs_data |> 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "2 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, Sex) |> 
  group_by(lab_name, Sex) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::t_test(
    formula = lab_value ~ Event_Name,
    ref.group = "2 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    # var.equal = TRUE, set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) %>%
  left_join(tofa_clinical_labs_data %>% distinct(lab_name, test_name, lab_units, panel)) %>% 
  group_by(panel) |> 
  mutate(BH_padj = p.adjust(p)) |> 
  ungroup() |> 
  arrange(p)
#

### 5.1.2 8 weeks   ----
tofa_clinical_labs_Ttest_res_8weeks_sex <- tofa_clinical_labs_data |> 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "8 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, Sex) |> 
  group_by(lab_name, Sex) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::t_test(
    formula = lab_value ~ Event_Name,
    ref.group = "8 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    # var.equal = TRUE, set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) %>%
  left_join(tofa_clinical_labs_data %>% distinct(lab_name, test_name, lab_units, panel)) %>% 
  group_by(panel) |> 
  mutate(BH_padj = p.adjust(p)) |> 
  ungroup() |> 
  arrange(p)
#

### 5.1.3 16 weeks    ----
tofa_clinical_labs_Ttest_res_16weeks_sex <- tofa_clinical_labs_data |> 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, Sex) |> 
  group_by(lab_name, Sex) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::t_test(
    formula = lab_value ~ Event_Name,
    ref.group = "16 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    # var.equal = TRUE, set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) %>%
  left_join(tofa_clinical_labs_data %>% distinct(lab_name, test_name, lab_units, panel)) %>% 
  group_by(panel) |> 
  mutate(BH_padj = p.adjust(p)) |> 
  ungroup() |> 
  arrange(p)
#

### 5.1.4 40 weeks    ----
tofa_clinical_labs_Ttest_res_40weeks_sex <- tofa_clinical_labs_data |> 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "40 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, Sex) |> 
  group_by(lab_name, Sex) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::t_test(
    formula = lab_value ~ Event_Name,
    ref.group = "40 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    # var.equal = TRUE, set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) %>%
  left_join(tofa_clinical_labs_data %>% distinct(lab_name, test_name, lab_units, panel)) %>% 
  group_by(panel) |> 
  mutate(BH_padj = p.adjust(p)) |> 
  ungroup() |> 
  arrange(p)
#

## 5.2 Calculate paired t test effect sizes across timepoints    ----
### 5.2.1 2 weeks   ----
tofa_clinical_labs_Ttest_effsize_2weeks_sex <- tofa_clinical_labs_data |> 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "2 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, Sex) |> 
  group_by(lab_name, Sex) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::cohens_d(
    formula = lab_value ~ Event_Name,
    ref.group = "2 week", # Seems to be the reverse of how we typically do this
    paired = TRUE
  )
#

### 5.2.2 8 weeks   ----
tofa_clinical_labs_Ttest_effsize_8weeks_sex <- tofa_clinical_labs_data |> 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "8 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, Sex) |> 
  group_by(lab_name, Sex) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::cohens_d(
    formula = lab_value ~ Event_Name,
    ref.group = "8 week", # Seems to be the reverse of how we typically do this
    paired = TRUE
  )
#

### 5.2.3 16 weeks    ----
tofa_clinical_labs_Ttest_effsize_16weeks_sex <- tofa_clinical_labs_data |> 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, Sex) |> 
  group_by(lab_name, Sex) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::cohens_d(
    formula = lab_value ~ Event_Name,
    ref.group = "16 week", # Seems to be the reverse of how we typically do this
    paired = TRUE
  )
#

### 5.2.4 40 weeks    ----
tofa_clinical_labs_Ttest_effsize_40weeks_sex <- tofa_clinical_labs_data |> 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "40 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, Sex) |> 
  group_by(lab_name, Sex) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::cohens_d(
    formula = lab_value ~ Event_Name,
    ref.group = "40 week", # Seems to be the reverse of how we typically do this
    paired = TRUE
  )
#

## 5.3 Compile t test results   ----
# 2 weeks
tofa_clinical_labs_Ttest_res_2weeks_full_sex <- tofa_clinical_labs_Ttest_res_2weeks_sex |> 
  inner_join(tofa_clinical_labs_Ttest_effsize_2weeks_sex) |> 
  mutate(
    mean_diff = estimate, 
  ) |> 
  inner_join(tofa_clinical_labs_diffs_2_weeks |> group_by(lab_name) |> summarize(median_diff = median(difference))) |> 
  select(panel, lab_name, test_name, lab_units, mean_diff, median_diff, p, BH_padj, effsize, magnitude, .y., group1, group2, everything())
#

# 8 weeks
tofa_clinical_labs_Ttest_res_8weeks_full_sex <- tofa_clinical_labs_Ttest_res_8weeks_sex |> 
  inner_join(tofa_clinical_labs_Ttest_effsize_8weeks_sex) |> 
  mutate(
    mean_diff = estimate, 
  ) |> 
  inner_join(tofa_clinical_labs_diffs_8_weeks |> group_by(lab_name) |> summarize(median_diff = median(difference))) |> 
  select(panel, lab_name, test_name, lab_units, mean_diff, median_diff, p, BH_padj, effsize, magnitude, .y., group1, group2, everything())
#

# 16 weeks
tofa_clinical_labs_Ttest_res_16weeks_full_sex <- tofa_clinical_labs_Ttest_res_16weeks_sex |> 
  inner_join(tofa_clinical_labs_Ttest_effsize_16weeks_sex) |> 
  mutate(
    mean_diff = estimate, 
  ) |> 
  inner_join(tofa_clinical_labs_diffs_16_weeks |> group_by(lab_name) |> summarize(median_diff = median(difference))) |> 
  select(panel, lab_name, test_name, lab_units, mean_diff, median_diff, p, BH_padj, effsize, magnitude, .y., group1, group2, everything())
#

# 40 weeks
tofa_clinical_labs_Ttest_res_40weeks_full_sex <- tofa_clinical_labs_Ttest_res_40weeks_sex |> 
  inner_join(tofa_clinical_labs_Ttest_effsize_40weeks_sex) |> 
  mutate(
    mean_diff = estimate, 
  ) |> 
  inner_join(tofa_clinical_labs_diffs_40_weeks |> group_by(lab_name) |> summarize(median_diff = median(difference))) |> 
  select(panel, lab_name, test_name, lab_units, mean_diff, median_diff, p, BH_padj, effsize, magnitude, .y., group1, group2, everything())
#

## 5.4 Save results of sex-stratified T tests    ----
bind_rows(
  tofa_clinical_labs_Ttest_res_2weeks_sex %>% mutate(timepoint = "2 weeks"),
  tofa_clinical_labs_Ttest_res_8weeks_sex %>% mutate(timepoint = "8 weeks"),
  tofa_clinical_labs_Ttest_res_16weeks_sex %>% mutate(timepoint = "16 weeks"),
  tofa_clinical_labs_Ttest_res_40weeks_sex %>% mutate(timepoint = "40 weeks")
) %>% 
  select(timepoint, everything()) %>% 
  write_tsv(file = here("results", paste0(out_file_prefix, "tofa_clinical_labs_Ttest_results_sex_stratified.txt")))
#

# 6.0 Paired t tests stratified by obesity status at baseline   ----
baseline_obesity_status <- tofa_baseline_obesity_status %>% 
  select(ParticipantID, baseline_obesity_status)

## 6.1 Run paired t tests across timepoints   ----
### 6.1.1 2 weeks   ----
tofa_clinical_labs_Ttest_res_2weeks_ob <- tofa_clinical_labs_data |> 
  inner_join(baseline_obesity_status %>% select(ParticipantID, baseline_obesity_status)) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "2 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, baseline_obesity_status) |> 
  group_by(lab_name, baseline_obesity_status) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::t_test(
    formula = lab_value ~ Event_Name,
    ref.group = "2 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    # var.equal = TRUE, set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) %>%
  left_join(tofa_clinical_labs_data %>% distinct(lab_name, test_name, lab_units, panel)) %>% 
  group_by(panel) |> 
  mutate(BH_padj = p.adjust(p)) |> 
  ungroup() |> 
  arrange(p)
#

### 6.1.2 8 weeks   ----
tofa_clinical_labs_Ttest_res_8weeks_ob <- tofa_clinical_labs_data |> 
  inner_join(baseline_obesity_status %>% select(ParticipantID, baseline_obesity_status)) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "8 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, baseline_obesity_status) |> 
  group_by(lab_name, baseline_obesity_status) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::t_test(
    formula = lab_value ~ Event_Name,
    ref.group = "8 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    # var.equal = TRUE, set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) %>%
  left_join(tofa_clinical_labs_data %>% distinct(lab_name, test_name, lab_units, panel)) %>% 
  group_by(panel) |> 
  mutate(BH_padj = p.adjust(p)) |> 
  ungroup() |> 
  arrange(p)
#

### 6.1.3 16 weeks    ----
tofa_clinical_labs_Ttest_res_16weeks_ob <- tofa_clinical_labs_data |> 
  inner_join(baseline_obesity_status %>% select(ParticipantID, baseline_obesity_status)) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, baseline_obesity_status) |> 
  group_by(lab_name, baseline_obesity_status) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::t_test(
    formula = lab_value ~ Event_Name,
    ref.group = "16 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    # var.equal = TRUE, set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) %>%
  left_join(tofa_clinical_labs_data %>% distinct(lab_name, test_name, lab_units, panel)) %>% 
  group_by(panel) |> 
  mutate(BH_padj = p.adjust(p)) |> 
  ungroup() |> 
  arrange(p)
#

### 6.1.4 40 weeks    ----
tofa_clinical_labs_Ttest_res_40weeks_ob <- tofa_clinical_labs_data |> 
  inner_join(baseline_obesity_status %>% select(ParticipantID, baseline_obesity_status)) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "40 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, baseline_obesity_status) |> 
  group_by(lab_name, baseline_obesity_status) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::t_test(
    formula = lab_value ~ Event_Name,
    ref.group = "40 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    # var.equal = TRUE, set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) %>%
  left_join(tofa_clinical_labs_data %>% distinct(lab_name, test_name, lab_units, panel)) %>% 
  group_by(panel) |> 
  mutate(BH_padj = p.adjust(p)) |> 
  ungroup() |> 
  arrange(p)
#

## 6.2 Calculate paired t test effect sizes across timepoints    ----
### 6.2.1 2 weeks   ----
tofa_clinical_labs_Ttest_effsize_2weeks_ob <- tofa_clinical_labs_data |> 
  inner_join(baseline_obesity_status %>% select(ParticipantID, baseline_obesity_status)) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "2 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, baseline_obesity_status) |> 
  group_by(lab_name, baseline_obesity_status) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::cohens_d(
    formula = lab_value ~ Event_Name,
    ref.group = "2 week", # Seems to be the reverse of how we typically do this
    paired = TRUE
  )
#

### 6.2.2 8 weeks   ----
tofa_clinical_labs_Ttest_effsize_8weeks_ob <- tofa_clinical_labs_data |> 
  inner_join(baseline_obesity_status %>% select(ParticipantID, baseline_obesity_status)) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "8 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, baseline_obesity_status) |> 
  group_by(lab_name, baseline_obesity_status) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::cohens_d(
    formula = lab_value ~ Event_Name,
    ref.group = "8 week", # Seems to be the reverse of how we typically do this
    paired = TRUE
  )
#

### 6.2.3 16 weeks    ----
tofa_clinical_labs_Ttest_effsize_16weeks_ob <- tofa_clinical_labs_data |> 
  inner_join(baseline_obesity_status %>% select(ParticipantID, baseline_obesity_status)) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, baseline_obesity_status) |> 
  group_by(lab_name, baseline_obesity_status) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::cohens_d(
    formula = lab_value ~ Event_Name,
    ref.group = "16 week", # Seems to be the reverse of how we typically do this
    paired = TRUE
  )
#

### 6.2.4 40 weeks    ----
tofa_clinical_labs_Ttest_effsize_40weeks_ob <- tofa_clinical_labs_data |> 
  inner_join(baseline_obesity_status %>% select(ParticipantID, baseline_obesity_status)) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "40 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, baseline_obesity_status) |> 
  group_by(lab_name, baseline_obesity_status) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::cohens_d(
    formula = lab_value ~ Event_Name,
    ref.group = "40 week", # Seems to be the reverse of how we typically do this
    paired = TRUE
  )
#

## 6.3 Compile t test results   ----
# 2 weeks
tofa_clinical_labs_Ttest_res_2weeks_full_ob <- tofa_clinical_labs_Ttest_res_2weeks_ob |> 
  inner_join(tofa_clinical_labs_Ttest_effsize_2weeks_ob) |> 
  mutate(
    mean_diff = estimate, 
  ) |> 
  inner_join(tofa_clinical_labs_diffs_2_weeks |> group_by(lab_name) |> summarize(median_diff = median(difference))) |> 
  select(panel, lab_name, test_name, lab_units, mean_diff, median_diff, p, BH_padj, effsize, magnitude, .y., group1, group2, everything())
#

# 8 weeks
tofa_clinical_labs_Ttest_res_8weeks_full_ob <- tofa_clinical_labs_Ttest_res_8weeks_ob |> 
  inner_join(tofa_clinical_labs_Ttest_effsize_8weeks_ob) |> 
  mutate(
    mean_diff = estimate, 
  ) |> 
  inner_join(tofa_clinical_labs_diffs_8_weeks |> group_by(lab_name) |> summarize(median_diff = median(difference))) |> 
  select(panel, lab_name, test_name, lab_units, mean_diff, median_diff, p, BH_padj, effsize, magnitude, .y., group1, group2, everything())
#

# 16 weeks
tofa_clinical_labs_Ttest_res_16weeks_full_ob <- tofa_clinical_labs_Ttest_res_16weeks_ob |> 
  inner_join(tofa_clinical_labs_Ttest_effsize_16weeks_ob) |> 
  mutate(
    mean_diff = estimate, 
  ) |> 
  inner_join(tofa_clinical_labs_diffs_16_weeks |> group_by(lab_name) |> summarize(median_diff = median(difference))) |> 
  select(panel, lab_name, test_name, lab_units, mean_diff, median_diff, p, BH_padj, effsize, magnitude, .y., group1, group2, everything())
#

# 40 weeks
tofa_clinical_labs_Ttest_res_40weeks_full_ob <- tofa_clinical_labs_Ttest_res_40weeks_ob |> 
  inner_join(tofa_clinical_labs_Ttest_effsize_40weeks_ob) |> 
  mutate(
    mean_diff = estimate, 
  ) |> 
  inner_join(tofa_clinical_labs_diffs_40_weeks |> group_by(lab_name) |> summarize(median_diff = median(difference))) |> 
  select(panel, lab_name, test_name, lab_units, mean_diff, median_diff, p, BH_padj, effsize, magnitude, .y., group1, group2, everything())
#

## 6.4 Save results of sex-stratified T tests    ----
bind_rows(
  tofa_clinical_labs_Ttest_res_2weeks_ob %>% mutate(timepoint = "2 weeks"),
  tofa_clinical_labs_Ttest_res_8weeks_ob %>% mutate(timepoint = "8 weeks"),
  tofa_clinical_labs_Ttest_res_16weeks_ob %>% mutate(timepoint = "16 weeks"),
  tofa_clinical_labs_Ttest_res_40weeks_ob %>% mutate(timepoint = "40 weeks")
) %>% 
  select(timepoint, everything()) %>% 
  write_tsv(file = here("results", paste0(out_file_prefix, "tofa_clinical_labs_Ttest_results_obesity_stratified.txt")))
#

# 7.0 Paired t tests stratified by history of covid event between baseline and timepoint    ----
## 7.1 Run paired t tests across timepoints   ----
covid_status_16w <- tofa_covid_event_history %>% 
  filter(Event_Name == "16 week") %>% 
  select(ParticipantID, COVID_event_hx)
#
covid_status_40w <- tofa_covid_event_history %>% 
  filter(Event_Name == "40 week") %>% 
  select(ParticipantID, COVID_event_hx)
#

### 7.1.1 16 weeks   ----
tofa_clinical_labs_Ttest_res_16weeks_covid <- tofa_clinical_labs_data |> 
  filter(Safety_eligible == TRUE) %>% 
  inner_join(covid_status_16w) %>% 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, COVID_event_hx) |> 
  group_by(lab_name, COVID_event_hx) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::t_test(
    formula = lab_value ~ Event_Name,
    ref.group = "16 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    # var.equal = TRUE, set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) %>%
  left_join(tofa_clinical_labs_data %>% distinct(lab_name, test_name, lab_units, panel)) %>% 
  group_by(panel) |> 
  mutate(BH_padj = p.adjust(p)) |> 
  ungroup() |> 
  arrange(p)
#

### 7.1.2 40 weeks    ----
tofa_clinical_labs_Ttest_res_40weeks_covid <- tofa_clinical_labs_data |> 
  filter(Safety_eligible == TRUE) %>% 
  inner_join(covid_status_40w) %>% 
  filter(Event_Name %in% c("Baseline", "40 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, COVID_event_hx) |> 
  group_by(lab_name, COVID_event_hx) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::t_test(
    formula = lab_value ~ Event_Name,
    ref.group = "40 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    # var.equal = TRUE, set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) %>%
  left_join(tofa_clinical_labs_data %>% distinct(lab_name, test_name, lab_units, panel)) %>% 
  group_by(panel) |> 
  mutate(BH_padj = p.adjust(p)) |> 
  ungroup() |> 
  arrange(p)
#

## 7.2 Calculate paired t test effect sizes across timepoints    ----
### 7.2.1 16 weeks   ----
tofa_clinical_labs_Ttest_effsize_16weeks_covid <- tofa_clinical_labs_data |> 
  inner_join(covid_status_16w) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, COVID_event_hx) |> 
  group_by(lab_name, COVID_event_hx) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::cohens_d(
    formula = lab_value ~ Event_Name,
    ref.group = "16 week", # Seems to be the reverse of how we typically do this
    paired = TRUE
  )
#

### 7.2.2 40 weeks   ----
tofa_clinical_labs_Ttest_effsize_40weeks_covid <- tofa_clinical_labs_data |> 
  inner_join(covid_status_40w) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name %in% c("Baseline", "40 week")) |> 
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, lab_name, COVID_event_hx) |> 
  group_by(lab_name, COVID_event_hx) |> 
  filter(!any(n < 2)) |> # remove labs with n < X in any group
  rstatix::cohens_d(
    formula = lab_value ~ Event_Name,
    ref.group = "40 week", # Seems to be the reverse of how we typically do this
    paired = TRUE
  )
#

## 7.3 Compile t test results   ----
# 16 weeks
tofa_clinical_labs_Ttest_res_16weeks_full_covid <- tofa_clinical_labs_Ttest_res_16weeks_covid |> 
  inner_join(tofa_clinical_labs_Ttest_effsize_16weeks_covid) |> 
  mutate(
    mean_diff = estimate, 
  ) |> 
  inner_join(tofa_clinical_labs_diffs_16_weeks |> group_by(lab_name) |> summarize(median_diff = median(difference))) |> 
  select(panel, lab_name, test_name, lab_units, mean_diff, median_diff, p, BH_padj, effsize, magnitude, .y., group1, group2, everything())
#

# 40 weeks
tofa_clinical_labs_Ttest_res_40weeks_full_covid <- tofa_clinical_labs_Ttest_res_40weeks_covid |> 
  inner_join(tofa_clinical_labs_Ttest_effsize_40weeks_covid) |> 
  mutate(
    mean_diff = estimate, 
  ) |> 
  inner_join(tofa_clinical_labs_diffs_40_weeks |> group_by(lab_name) |> summarize(median_diff = median(difference))) |> 
  select(panel, lab_name, test_name, lab_units, mean_diff, median_diff, p, BH_padj, effsize, magnitude, .y., group1, group2, everything())
#

## 7.4 Save results of sex-stratified T tests    ----
bind_rows(
  tofa_clinical_labs_Ttest_res_16weeks_covid %>% mutate(timepoint = "16 weeks"),
  tofa_clinical_labs_Ttest_res_40weeks_covid %>% mutate(timepoint = "40 weeks")
) %>% 
  select(timepoint, everything()) %>% 
  write_tsv(file = here("results", paste0(out_file_prefix, "tofa_clinical_labs_Ttest_results_covid_stratified.txt")))
#

# 8.0 Plots   ----
## 8.1 Sina plots: Clinical lab values    ----
labs_of_interest <- c(
  "WBC",
  "NEUTROPHIL",
  "Hgb",
  "AST",
  "ALT",
  "HDL-C",
  "Urea_Nitrogen_BUN"
)
#
tofa_clinical_labs_data |> 
  filter(lab_name  %in% labs_of_interest) %>% 
  mutate(lab_name = paste0(lab_name, "[", lab_units, "]")) %>% 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> 
  #
  filter(!is.na(lab_value)) |> 
  filter(is.finite(lab_value)) |> 
  add_count(lab_name, ParticipantID, name = "pair") |> 
  filter(pair >= 2) |>  # require at least 2 data points
  #
  group_by(lab_name) |> 
  group_split() |> 
  map( # using anonymous function: \(x)
    \(x) x |> 
      arrange(lab_name, Event_Name) |> 
      mutate( # add labels with n
        n = n_distinct(ParticipantID), 
        .by = c(lab_name, Event_Name),
        label = paste0(str_extract(Event_Name, "^[[:alpha:]]|^\\d{1,2}"), "\n(", n, ")"),
        label = fct_inorder(label)
      ) |> 
      mutate(Event_Name = fct_drop(Event_Name)) |> # prevent empty x values
      complete(lab_name, ParticipantID, Event_Name) %>% # Prevent error from stat_signif()
      #
      {ggplot(data = ., aes(Event_Name, lab_value, color = Event_Name)) +
          scale_color_manual(values = standard_colors) +
          geom_line(aes(group = ParticipantID), color = "gray50", alpha = 0.3, linewidth = 0.75) +
          geom_sina(maxwidth = 0.3) +
          geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
          facet_wrap(~ lab_name, scales = "free_y", nrow = 1) +
          theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            aspect.ratio = 1.3,
            panel.border = element_blank(),
            axis.line = element_line(colour = "black")
          ) +
          labs(
            x = NULL,
            y = "lab_units"
          ) +
          scale_x_discrete(
            labels = . |>
              filter(!is.na(label)) |>
              distinct(lab_name, Event_Name, label) |>
              select(Event_Name, label) |>
              transmute(as.character(Event_Name), as.character(label)) |>
              deframe()
          )
      }
  ) |> 
  patchwork::wrap_plots() +
  patchwork::plot_layout(
    # guides = "collect", 
    nrow = 1
  ) +
  patchwork::plot_annotation(
    title = "Clinical labs: Baseline vs. Treatment"
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_plots_all", ".pdf")), device = cairo_pdf, width = 120, height = 6.5, units = "in", limitsize =  FALSE)
#

## 8.2 Sina plot: Stratified clinical lab differences from baseline    ----
# Example: AST in Obese and Non-Obese
diffs_combined_ob <- bind_rows(
  tofa_clinical_labs_diffs_2_weeks |> mutate(timepoint = "2 weeks"),
  tofa_clinical_labs_diffs_8_weeks |> mutate(timepoint = "8 weeks"),
  tofa_clinical_labs_diffs_16_weeks |> mutate(timepoint = "16 weeks"),
  tofa_clinical_labs_diffs_40_weeks |> mutate(timepoint = "40 weeks")
) %>% 
  mutate(timepoint = fct_relevel(timepoint, c("2 weeks", "8 weeks", "16 weeks", "40 weeks"))) %>%
  inner_join(baseline_obesity_status) %>% 
  mutate(baseline_obesity_status = fct_relevel(baseline_obesity_status, c("Non-obese", "Obese")))
#
diffs_combined_ob %>% 
  filter(lab_name == "AST") %>% 
  filter(!is.na(difference)) |> 
  filter(is.finite(difference)) |> 
  add_count(baseline_obesity_status, ParticipantID, name = "pair") |> 
  filter(pair >= 2) |>  # require at least 2 data points
  group_by(baseline_obesity_status) |> 
  group_split() |> 
  map( # using anonymous function: \(x)
    \(x) x |> 
      arrange(baseline_obesity_status, timepoint) |> 
      mutate( # add labels with n
        n = n_distinct(ParticipantID), 
        .by = c(baseline_obesity_status, timepoint),
        label = paste0(str_extract(timepoint, "^[[:alpha:]]|^\\d{1,2}"), "\n(n=", n, ")"),
        label = fct_inorder(label)
      ) |> 
      mutate(timepoint = fct_drop(timepoint)) |> # prevent empty x values
      mutate(baseline_obesity_status = fct_drop(baseline_obesity_status)) |> # prevent empty x values
      complete(baseline_obesity_status, ParticipantID, timepoint) %>% # Prevent error from stat_signif()
      #
      {ggplot(data = ., aes(timepoint, difference, color = timepoint)) +
          scale_color_manual(values = standard_colors2) +
          geom_sina(maxwidth = 0.3) +
          geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
          geom_hline(aes(yintercept = 0), linetype = 2) +
          geom_line(aes(group = ParticipantID), color = "gray50", alpha = 0.3, linewidth = 0.75) +
          facet_wrap(~ baseline_obesity_status, nrow = 1) +
          scale_y_continuous(limits = c(-8, 21)) +
          theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            aspect.ratio = 1.3,
            panel.border = element_blank(),
            axis.line = element_line(colour = "black")
          ) +
          labs(
            x = NULL,
            y = "difference"
          ) +
          scale_x_discrete(
            labels = . |>
              filter(!is.na(label)) |>
              distinct(baseline_obesity_status, timepoint, label) |>
              select(timepoint, label) |>
              transmute(as.character(timepoint), as.character(label)) |>
              deframe()
          )
      }
  ) |> 
  patchwork::wrap_plots() +
  patchwork::plot_layout(
    nrow = 1
  ) +
  patchwork::plot_annotation(
    title = "AST: Differences Treatment vs. Baseline",
    subtitle = "BMI stratified"
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_plots_AST_diffs_by_BMI", ".pdf")), device = cairo_pdf, width = 25, height = 6.5, units = "in")
#

## 8.3 Forest plot: Statistics from T test comparing timepoints vs. baseline    ----
# Example: select Clinical Labs from CBC + Differential
# Set plot parameters
row_unit_height <- 0.18
base_height <- 1.5
plot_width <- 6.5
gap_size <- 5   # 4 timepoints + 1 spacer

# Select labs of interest
cbc_diff_labs_select <- c(
  "WBC", 
  "NEUTROPHIL",
  "LYMPHOCYTE",
  "RBC",
  "percent_Hct",
  "Hgb"
)
#
# Render plot
cbc_fp_dat <- bind_rows(
  tofa_clinical_labs_Ttest_res_2weeks_full,
  tofa_clinical_labs_Ttest_res_8weeks_full,
  tofa_clinical_labs_Ttest_res_16weeks_full,
  tofa_clinical_labs_Ttest_res_40weeks_full
) %>% 
  rename(timepoint = group1) %>% 
  filter(panel %in% c("Differential", "CBC")) |> 
  filter(lab_name %in% cbc_diff_labs_select) |> 
  mutate(
    lab_name = factor(lab_name, levels = cbc_diff_labs_select),
    timepoint = fct_relevel(timepoint, "2 week", "8 week", "16 week", "40 week")
  ) |> 
  arrange(lab_name, timepoint) |> 
  group_by(lab_name) |> 
  mutate(tp_index = row_number()) |> 
  ungroup() |> 
  mutate(
    lab_index = as.integer(lab_name),
    y_pos = lab_index * gap_size + tp_index,
    row_label = paste(lab_name, timepoint, sep = " | "),
    color = if_else(p < 0.05, "p < 0.05", "n.s.")
  ) 
#
cbc_fp_dat %>% 
  ggplot(aes(estimate, y_pos)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point(shape = 15, aes(size = -log10(p), color = color)) +
  scale_y_reverse(
    breaks = cbc_fp_dat$y_pos,
    labels = cbc_fp_dat$row_label
  ) +
  scale_color_manual(values = c("p < 0.05" = "red", "n.s." = "black")) +
  theme(
    legend.key = element_blank(),
    aspect.ratio = 1.2
  ) +
  labs(
    title = "CBC + Differential: Treatment vs. Baseline",
    subtitle = "T test",
    x = "Estimate",
    y = NULL
  )
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