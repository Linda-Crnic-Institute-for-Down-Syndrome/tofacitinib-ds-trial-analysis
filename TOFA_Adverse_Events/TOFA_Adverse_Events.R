################################################
# Title: Tofacitinib for Immune Skin Conditions in Down Syndrome: Analysis of Adverse Events
# Author(s):
#   - Micah Donovan
# Affiliation(s):
#   - Linda Crnic Institute for Down Syndrome
#   - University of Colorado Anschutz
################################################

### Summary:  
# This workflow analyzes Adverse Event data based on physician assessment of grade, 
# severity, and relatedness to tofacitinib.
# See README.md for more details.
#  

### Data type(s):
# Clinical trial (TOFA) datasets:
#    * Participant-level metadata; Available on request.
#    * Visit/Event-level metadata; Available on request.
#    * Baseline obesity status; Available on request.
#    * COVID-19 history; Available on request.
#    * Adverse Events from baseline to 16 weeks; Available on request.
#    * Adverse Events from baseline to 40 weeks; Available on request.
#


# 0 General Setup -----
# RUN THIS FIRST TIME - Initialize and install packages with renv:
# renv::init(bioconductor = TRUE)
#
# To install the exact versions of all R packages base on renv.lock file (requires matching R version):
# renv::restore()

## 0.1 Load required libraries ----
library("knitr")
library("DT") # Used for sortable, searchable datatables in reports
library("readxl") # Used to read .xlsx files
library("openxlsx") # used for data export as Excel workbooks
library("tidyverse")
library("magrittr")
library("ggrepel") # required for labelling genes
library("ggforce") # required for zooming and sina
library("plotly") # required for interactive plots/IFNA1
library("tictoc") # timer
library("skimr")
library("janitor")
library("ggplot2")
library("here")
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
tofa_ae_16w_data_file <- here("data", "TOFA_Adverse_Events_Data_16weeks_zenodo_v1.txt") 
tofa_ae_40w_data_file <- here("data", "TOFA_Adverse_Events_Data_40weeks_zenodo_v1.txt") 
tofa_baseline_obesity_status_file <- here("data", "TOFA_Baseline_Obesity_Status_zenodo_v1.txt")
tofa_covid_event_history_file <- here("data", "TOFA_COVID_History_zenodo_v1.txt")
# 
# Other parameters:
out_file_prefix <- "TOFA_Adverse_Events_" # should match this script title
#
source(here("helper_functions.R")) # load helper functions
#

# 1.0 Read in data     ----
## 1.1 Participant level meta data    ----
tofa_participant_meta_data <- tofa_participant_meta_data_file |> 
  read_tsv() |> 
  mutate(
    Sex = fct_relevel(Sex, c("Female", "Male")) # set factor levels
  )
#
tofa_participant_meta_data # 47 rows
tofa_participant_meta_data |> distinct(ParticipantID) # 47 participants in this table
tofa_participant_meta_data |> count(Safety_eligible, Participant_notes) # 4 participants not eligible for Safety analyses
#

## 1.2 Visit level meta data ----
tofa_visit_meta_data <- tofa_visit_meta_data_file |> 
  read_tsv() |> 
  mutate(
    Event_Name = fct_relevel(Event_Name, c("Baseline", "2 week", "8 week", "16 week","40 week")) # set factor levels
  )
#
tofa_visit_meta_data # 202 rows
tofa_visit_meta_data |> distinct(ParticipantID) # 43 Participants with samples
#

## 1.3 Participant baseline obesity status information   ----
# File listing participants as obese or non-obese at baseline
tofa_baseline_obesity_status <- tofa_baseline_obesity_status_file %>% 
  read_tsv() %>% 
  inner_join(tofa_participant_meta_data)
#
tofa_baseline_obesity_status # 43 rows
tofa_baseline_obesity_status |> distinct(ParticipantID) # 43 participants in this table
#

## 1.4 Participant COVID-19 event history information   ----
# File listing whether a covid event was recorded for each between baseline and study time points
tofa_covid_event_history <- tofa_covid_event_history_file %>% 
  read_tsv() %>% 
  inner_join(tofa_participant_meta_data) %>% 
  filter(Safety_eligible == TRUE)
#
tofa_covid_event_history # 172 rows
tofa_baseline_obesity_status |> distinct(ParticipantID) # 43 participants in this table
#

## 1.5 Tofacitinib Adverse Events data    ----
### 1.5.1 16 week Tofacitinib Adverse Event data    ----
tofa_ae_16w_data <- tofa_ae_16w_data_file %>% 
  read_tsv() %>% 
  inner_join(tofa_participant_meta_data) %>% 
  filter(Safety_eligible == TRUE)
#
tofa_ae_16w_data # 2,247 rows
tofa_ae_16w_data %>% skimr::skim()
tofa_ae_16w_data %>% distinct(ParticipantID) # 43 participants in this table
tofa_ae_16w_data %>% distinct(AE) # 39 adverse event types in this table
#

### 1.5.2 40 week Tofacitinib Adverse Event data    ----
tofa_ae_40w_data <- tofa_ae_40w_data_file %>% 
  read_tsv() %>% 
  inner_join(tofa_participant_meta_data) %>% 
  filter(Safety_eligible == TRUE)
#
tofa_ae_40w_data # 2,582 rows
tofa_ae_40w_data %>% skimr::skim()
tofa_ae_40w_data %>% distinct(ParticipantID) # 43 participants in this table
tofa_ae_40w_data %>% distinct(AE) # 42 adverse event types in this table
#

# 2.0 Calculate Incidence Rates for both time points    ----
## 2.1 16 week Incidence Rates    ----
total_n_16w <- length(unique(tofa_ae_16w_data %>% 
                           pull(ParticipantID)))
#
irs_16w <- tofa_ae_16w_data %>%
  filter(AE_occurred == 1) %>%
  mutate(Grade = if_else(is.na(Grade), 4, Grade)) %>% # Set placeholder Grade 4 for NAs
  group_by(AE_Event_Name, ParticipantID) %>%
  slice_max(order_by = Grade, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(AE) %>%
  summarise(n_participants = n_distinct(ParticipantID), .groups = 'drop') %>%
  mutate(incidence = round(n_participants / total_n_16w * 100, digits = 1)) %>%
  arrange(-incidence)
irs_16w
#
irs_16w %>% write_tsv(file = here("results", paste0(out_file_prefix, "irs_16w.txt")))
#

## 2.2 40 week Incidence Rates    ----
total_n_40w <- length(unique(tofa_ae_40w_data %>% 
                               pull(ParticipantID)))
#
irs_40w <- tofa_ae_40w_data %>%
  filter(AE_occurred == 1) %>%
  mutate(Grade = if_else(is.na(Grade), 4, Grade)) %>% # Set placeholder Grade 4 for NAs
  group_by(AE_Event_Name, ParticipantID) %>%
  slice_max(order_by = Grade, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(AE) %>%
  summarise(n_participants = n_distinct(ParticipantID), .groups = 'drop') %>%
  mutate(incidence = round(n_participants / total_n_40w * 100, digits = 1)) %>%
  arrange(-incidence)
irs_40w
#
irs_40w %>% write_tsv(file = here("results", paste0(out_file_prefix, "irs_40w.txt")))
#

# 3.0 Calculate Event Rates for both time points    ----
## 3.1 16 week Event Rates    ----
ers_16w <- tofa_ae_16w_data %>%
  filter(AE_occurred == 1) %>%
  group_by(AE) %>%
  summarise(
    n_events = n(),  # Each row is one AE occurrence
    .groups = "drop"
  ) %>%
  mutate(
    event_rate_per_100 = round(n_events / total_n_16w * 100, digits = 1)
  ) %>%
  arrange(desc(event_rate_per_100))
ers_16w
#
ers_16w %>% write_tsv(file = here("results", paste0(out_file_prefix, "ers_16w.txt")))
#

## 3.2 40 week Event Rates    ----
ers_40w <- tofa_ae_40w_data %>%
  filter(AE_occurred == 1) %>%
  group_by(AE) %>%
  summarise(
    n_events = n(),  # Each row is one AE occurrence
    .groups = "drop"
  ) %>%
  mutate(
    event_rate_per_100 = round(n_events / total_n_40w * 100, digits = 1)
  ) %>%
  arrange(desc(event_rate_per_100))
ers_40w
#
ers_40w %>% write_tsv(file = here("results", paste0(out_file_prefix, "ers_40w.txt")))
#

# 4.0 Calculate Exposure-Adjusted Incidence Rates (EAIR) for both time points   ----
## 4.1 16 week EAIR   ----
eairs_16w <- tofa_ae_16w_data %>%
  group_by(ParticipantID, AE) %>% 
  slice_min(AE_exposure_number, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% 
  group_by(AE) %>%
  summarise(
    n_participants = sum(AE_occurred),
    total_patient_years = sum(AE_exposure_years, na.rm = TRUE),
    EAIR_per_100py = round((n_participants/total_patient_years) * 100, digits = 1)
  ) %>% 
  arrange(-EAIR_per_100py)
eairs_16w
#
eairs_16w %>% write_tsv(file = here("results", paste0(out_file_prefix, "eairs_16w.txt")))
#

## 4.2 40 week EAIR   ----
eairs_40w <- tofa_ae_40w_data %>%
  group_by(ParticipantID, AE) %>% 
  slice_min(AE_exposure_number, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% 
  group_by(AE) %>%
  summarise(
    n_participants = sum(AE_occurred),
    total_patient_years = sum(AE_exposure_years, na.rm = TRUE),
    EAIR_per_100py = round((n_participants/total_patient_years) * 100, digits = 1)
  ) %>% 
  arrange(-EAIR_per_100py)
eairs_40w
#
eairs_40w %>% write_tsv(file = here("results", paste0(out_file_prefix, "eairs_40w.txt")))
#

# 5.0 Calculate Exposure-Adjusted Event Rates (EAER) for both time points   ----
## 5.1 16 week EAER   ----
eaers_16w <- tofa_ae_16w_data %>%
  group_by(ParticipantID, AE) %>%
  summarise(
    n_events = sum(AE_occurred),
    .groups = "drop"
  ) %>%
  left_join(tofa_ae_16w_data %>% 
              distinct(ParticipantID, participant_years_constrained)) %>% 
  group_by(AE) %>%
  summarise(
    total_participant_years = sum(participant_years_constrained),
    n_events = sum(n_events),
    EAER_per_100py = round((n_events / total_participant_years) * 100, digits = 1)
  ) %>%
  arrange(desc(n_events))
eaers_16w
#
eaers_16w %>% write_tsv(file = here("results", paste0(out_file_prefix, "eaers_16w.txt")))
#

## 5.2 40 week EAER   ----
eaers_40w <- tofa_ae_40w_data %>%
  group_by(ParticipantID, AE) %>%
  summarise(
    n_events = sum(AE_occurred),
    .groups = "drop"
  ) %>%
  left_join(tofa_ae_40w_data %>% 
              distinct(ParticipantID, participant_years_constrained)) %>% 
  group_by(AE) %>%
  summarise(
    total_participant_years = sum(participant_years_constrained),
    n_events = sum(n_events),
    EAER_per_100py = round((n_events / total_participant_years) * 100, digits = 1)
  ) %>%
  arrange(desc(n_events))
eaers_40w
#
eaers_40w %>% write_tsv(file = here("results", paste0(out_file_prefix, "eaers_40w.txt")))
#

# 6.0 Stratified analyses of Incidence Rates by age group, sex, obesity status, and COVID-19 events   ----
# Create data object with stratification groups per participant
participant_strata <- tofa_visit_meta_data %>% 
  inner_join(tofa_participant_meta_data %>% 
               select(ParticipantID, Safety_eligible, Sex)) %>% 
  filter(Safety_eligible == TRUE) %>% 
  filter(Event_Name == "Baseline") %>% 
  distinct(ParticipantID, Age_years_at_visit, Sex) %>% 
  rename(Baseline_Age = Age_years_at_visit) %>% 
  mutate(Baseline_Age_Group = case_when(
    Baseline_Age < 18 ~ "Pediatric",
    Baseline_Age >= 18 ~ "Adult"
  )) %>% 
  inner_join(tofa_baseline_obesity_status %>% distinct(ParticipantID, baseline_obesity_status)) %>% 
  left_join(tofa_covid_event_history %>%
              filter(Event_Name %in% c("16 week", "40 week")) %>% 
              pivot_wider(names_from = Event_Name, values_from = COVID_event_hx) %>% 
              rename(covid_event_by_16w = '16 week') %>% 
              rename(covid_event_by_40w = '40 week')) %>% 
  select(ParticipantID, Baseline_Age_Group, Sex, baseline_obesity_status, covid_event_by_16w, covid_event_by_40w)
participant_strata
#

# Reduce AE data to 1 row per (ParticipantID, AE) with AE_occurred in {0,1}
## 16 weeks
ae_incidence_16w <- tofa_ae_16w_data %>%
  mutate(
    Grade = if_else(is.na(Grade), 4, Grade) # placeholder choice
  ) %>%
  group_by(AE_Event_Name, ParticipantID) %>%
  slice_max(order_by = Grade, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(ParticipantID, AE = AE, AE_occurred = AE_occurred) %>%
  group_by(ParticipantID, AE) %>%
  summarise(AE_occurred = max(AE_occurred), .groups = "drop")
ae_incidence_16w
#

## 40 weeks
ae_incidence_40w <- tofa_ae_40w_data %>%
  mutate(
    Grade = if_else(is.na(Grade), 4, Grade) # placeholder choice
  ) %>%
  group_by(AE_Event_Name, ParticipantID) %>%
  slice_max(order_by = Grade, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(ParticipantID, AE = AE, AE_occurred = AE_occurred) %>%
  group_by(ParticipantID, AE) %>%
  summarise(AE_occurred = max(AE_occurred), .groups = "drop")
ae_incidence_40w
#

## 6.1 Fishers exact test for incidence by age groups   ----
### 6.1.1 16 weeks    ----
fisher_16w_age <- participant_strata %>%
  select(ParticipantID, Baseline_Age_Group) %>% 
  left_join(ae_incidence_16w, by = "ParticipantID") %>%
  group_by(AE) %>%
  complete(ParticipantID, fill = list(AE_occurred = 0L)) %>%
  ungroup() %>%
  select(ParticipantID, Baseline_Age_Group, AE, AE_occurred) %>%
  group_by(AE) %>%
  summarise(
    adult_yes   = sum(Baseline_Age_Group == "Adult"   & AE_occurred == 1),
    adult_no    = sum(Baseline_Age_Group == "Adult"   & AE_occurred == 0),
    pediatric_yes = sum(Baseline_Age_Group == "Pediatric" & AE_occurred == 1),
    pediatric_no  = sum(Baseline_Age_Group == "Pediatric" & AE_occurred == 0),
    .groups = "drop"
  ) %>%
  mutate(
    fisher = pmap(
      list(adult_yes, adult_no, pediatric_yes, pediatric_no), # OR = Adult vs Pediatric
      ~ fisher.test(matrix(c(..1, ..2, ..3, ..4), nrow = 2, byrow = TRUE))
    ),
    p_value    = map_dbl(fisher, ~ .x$p.value),
    odds_ratio = map_dbl(fisher, ~ unname(.x$estimate)),
    conf_low   = map_dbl(fisher, ~ .x$conf.int[1]),
    conf_high  = map_dbl(fisher, ~ .x$conf.int[2])
  ) %>%
  select(AE, adult_yes, adult_no, pediatric_yes, pediatric_no, odds_ratio, conf_low, conf_high, p_value) %>%
  arrange(p_value)
fisher_16w_age
#
fisher_16w_age %>% write_tsv(file = here("results", paste0(out_file_prefix, "fisher_16w_age.txt")))
#

### 6.1.2 40 weeks    ----
fisher_40w_age <- participant_strata %>%
  select(ParticipantID, Baseline_Age_Group) %>% 
  left_join(ae_incidence_40w, by = "ParticipantID") %>%
  group_by(AE) %>%
  complete(ParticipantID, fill = list(AE_occurred = 0L)) %>%
  ungroup() %>%
  select(ParticipantID, Baseline_Age_Group, AE, AE_occurred) %>%
  group_by(AE) %>%
  summarise(
    adult_yes   = sum(Baseline_Age_Group == "Adult"   & AE_occurred == 1),
    adult_no    = sum(Baseline_Age_Group == "Adult"   & AE_occurred == 0),
    pediatric_yes = sum(Baseline_Age_Group == "Pediatric" & AE_occurred == 1),
    pediatric_no  = sum(Baseline_Age_Group == "Pediatric" & AE_occurred == 0),
    .groups = "drop"
  ) %>%
  mutate(
    fisher = pmap(
      list(adult_yes, adult_no, pediatric_yes, pediatric_no), # OR = Adult vs Pediatric
      ~ fisher.test(matrix(c(..1, ..2, ..3, ..4), nrow = 2, byrow = TRUE))
    ),
    p_value    = map_dbl(fisher, ~ .x$p.value),
    odds_ratio = map_dbl(fisher, ~ unname(.x$estimate)),
    conf_low   = map_dbl(fisher, ~ .x$conf.int[1]),
    conf_high  = map_dbl(fisher, ~ .x$conf.int[2])
  ) %>%
  select(AE, adult_yes, adult_no, pediatric_yes, pediatric_no, odds_ratio, conf_low, conf_high, p_value) %>%
  arrange(p_value)
fisher_40w_age
#
fisher_40w_age %>% write_tsv(file = here("results", paste0(out_file_prefix, "fisher_40w_age.txt")))
#

## 6.2 Fishers exact test for incidence by sex    ----
### 6.2.1 16 weeks    ----
fisher_16w_sex <- participant_strata %>%
  select(ParticipantID, Sex) %>% 
  left_join(ae_incidence_16w, by = "ParticipantID") %>%
  group_by(AE) %>%
  complete(ParticipantID, fill = list(AE_occurred = 0L)) %>%
  ungroup() %>%
  select(ParticipantID, Sex, AE, AE_occurred) %>%
  group_by(AE) %>%
  summarise(
    male_yes   = sum(Sex == "Male"   & AE_occurred == 1),
    male_no    = sum(Sex == "Male"   & AE_occurred == 0),
    female_yes = sum(Sex == "Female" & AE_occurred == 1),
    female_no  = sum(Sex == "Female" & AE_occurred == 0),
    .groups = "drop"
  ) %>%
  mutate(
    fisher = pmap(
      list(male_yes, male_no, female_yes, female_no), # OR = Males vs Females
      ~ fisher.test(matrix(c(..1, ..2, ..3, ..4), nrow = 2, byrow = TRUE))
    ),
    p_value    = map_dbl(fisher, ~ .x$p.value),
    odds_ratio = map_dbl(fisher, ~ unname(.x$estimate)),
    conf_low   = map_dbl(fisher, ~ .x$conf.int[1]),
    conf_high  = map_dbl(fisher, ~ .x$conf.int[2])
  ) %>%
  select(AE, male_yes, male_no, female_yes, female_no, odds_ratio, conf_low, conf_high, p_value) %>%
  arrange(p_value)
fisher_16w_sex
#
fisher_16w_sex %>% write_tsv(file = here("results", paste0(out_file_prefix, "fisher_16w_sex.txt")))
#

### 6.2.2 40 weeks    ----
fisher_40w_sex <- participant_strata %>%
  select(ParticipantID, Sex) %>% 
  left_join(ae_incidence_40w, by = "ParticipantID") %>%
  group_by(AE) %>%
  complete(ParticipantID, fill = list(AE_occurred = 0L)) %>%
  ungroup() %>%
  select(ParticipantID, Sex, AE, AE_occurred) %>%
  group_by(AE) %>%
  summarise(
    male_yes   = sum(Sex == "Male"   & AE_occurred == 1),
    male_no    = sum(Sex == "Male"   & AE_occurred == 0),
    female_yes = sum(Sex == "Female" & AE_occurred == 1),
    female_no  = sum(Sex == "Female" & AE_occurred == 0),
    .groups = "drop"
  ) %>%
  mutate(
    fisher = pmap(
      list(male_yes, male_no, female_yes, female_no), # OR = Males vs Females
      ~ fisher.test(matrix(c(..1, ..2, ..3, ..4), nrow = 2, byrow = TRUE))
    ),
    p_value    = map_dbl(fisher, ~ .x$p.value),
    odds_ratio = map_dbl(fisher, ~ unname(.x$estimate)),
    conf_low   = map_dbl(fisher, ~ .x$conf.int[1]),
    conf_high  = map_dbl(fisher, ~ .x$conf.int[2])
  ) %>%
  select(AE, male_yes, male_no, female_yes, female_no, odds_ratio, conf_low, conf_high, p_value) %>%
  arrange(p_value)
fisher_40w_sex
#
fisher_40w_sex %>% write_tsv(file = here("results", paste0(out_file_prefix, "fisher_40w_sex.txt")))
#

## 6.3 Fishers exact test for incidence by obesity status at baseline    ----
### 6.3.1 16 weeks    ----
fisher_16w_obesity <- participant_strata %>%
  select(ParticipantID, baseline_obesity_status) %>% 
  left_join(ae_incidence_16w, by = "ParticipantID") %>%
  group_by(AE) %>%
  complete(ParticipantID, fill = list(AE_occurred = 0L)) %>%
  ungroup() %>%
  select(ParticipantID, baseline_obesity_status, AE, AE_occurred) %>%
  group_by(AE) %>%
  summarise(
    obese_yes   = sum(baseline_obesity_status == "Obese"   & AE_occurred == 1),
    obese_no    = sum(baseline_obesity_status == "Obese"   & AE_occurred == 0),
    nonobese_yes = sum(baseline_obesity_status == "Non-obese" & AE_occurred == 1),
    nonobese_no  = sum(baseline_obesity_status == "Non-obese" & AE_occurred == 0),
    .groups = "drop"
  ) %>%
  mutate(
    fisher = pmap(
      list(obese_yes, obese_no, nonobese_yes, nonobese_no), # OR = Males vs Females
      ~ fisher.test(matrix(c(..1, ..2, ..3, ..4), nrow = 2, byrow = TRUE))
    ),
    p_value    = map_dbl(fisher, ~ .x$p.value),
    odds_ratio = map_dbl(fisher, ~ unname(.x$estimate)),
    conf_low   = map_dbl(fisher, ~ .x$conf.int[1]),
    conf_high  = map_dbl(fisher, ~ .x$conf.int[2])
  ) %>%
  select(AE, obese_yes, obese_no, nonobese_yes, nonobese_no, odds_ratio, conf_low, conf_high, p_value) %>%
  arrange(p_value)
fisher_16w_obesity
#
fisher_16w_obesity %>% write_tsv(file = here("results", paste0(out_file_prefix, "fisher_16w_obesity.txt")))
#

### 6.3.2 40 weeks    ----
fisher_40w_obesity <- participant_strata %>%
  select(ParticipantID, baseline_obesity_status) %>% 
  left_join(ae_incidence_40w, by = "ParticipantID") %>%
  group_by(AE) %>%
  complete(ParticipantID, fill = list(AE_occurred = 0L)) %>%
  ungroup() %>%
  select(ParticipantID, baseline_obesity_status, AE, AE_occurred) %>%
  group_by(AE) %>%
  summarise(
    obese_yes   = sum(baseline_obesity_status == "Obese"   & AE_occurred == 1),
    obese_no    = sum(baseline_obesity_status == "Obese"   & AE_occurred == 0),
    nonobese_yes = sum(baseline_obesity_status == "Non-obese" & AE_occurred == 1),
    nonobese_no  = sum(baseline_obesity_status == "Non-obese" & AE_occurred == 0),
    .groups = "drop"
  ) %>%
  mutate(
    fisher = pmap(
      list(obese_yes, obese_no, nonobese_yes, nonobese_no), # OR = Males vs Females
      ~ fisher.test(matrix(c(..1, ..2, ..3, ..4), nrow = 2, byrow = TRUE))
    ),
    p_value    = map_dbl(fisher, ~ .x$p.value),
    odds_ratio = map_dbl(fisher, ~ unname(.x$estimate)),
    conf_low   = map_dbl(fisher, ~ .x$conf.int[1]),
    conf_high  = map_dbl(fisher, ~ .x$conf.int[2])
  ) %>%
  select(AE, obese_yes, obese_no, nonobese_yes, nonobese_no, odds_ratio, conf_low, conf_high, p_value) %>%
  arrange(p_value)
fisher_40w_obesity
#
fisher_40w_obesity %>% write_tsv(file = here("results", paste0(out_file_prefix, "fisher_40w_obesity.txt")))
#

## 6.4 Fishers exact test for incidence by history of covid event between baseline and timepoint    ----
### 6.4.1 16 weeks    ----
fisher_16w_covid <- participant_strata %>%
  select(ParticipantID, covid_event_by_16w) %>% 
  left_join(ae_incidence_16w, by = "ParticipantID") %>%
  filter(!is.na(covid_event_by_16w)) %>% 
  group_by(AE) %>%
  summarise(
    covid_yes   = sum(covid_event_by_16w == "yes"   & AE_occurred == 1),
    covid_no    = sum(covid_event_by_16w == "yes"   & AE_occurred == 0),
    noncovid_yes = sum(covid_event_by_16w == "no" & AE_occurred == 1),
    noncovid_no  = sum(covid_event_by_16w == "no" & AE_occurred == 0),
    .groups = "drop"
  ) %>%
  mutate(
    fisher = pmap(
      list(covid_yes, covid_no, noncovid_yes, noncovid_no), # OR = covid vs non-covid
      ~ fisher.test(matrix(c(..1, ..2, ..3, ..4), nrow = 2, byrow = TRUE))
    ),
    p_value    = map_dbl(fisher, ~ .x$p.value),
    odds_ratio = map_dbl(fisher, ~ unname(.x$estimate)),
    conf_low   = map_dbl(fisher, ~ .x$conf.int[1]),
    conf_high  = map_dbl(fisher, ~ .x$conf.int[2])
  ) %>%
  select(AE, covid_yes, covid_no, noncovid_yes, noncovid_no, odds_ratio, conf_low, conf_high, p_value) %>%
  arrange(p_value)
fisher_16w_covid
#
fisher_16w_covid %>% write_tsv(file = here("results", paste0(out_file_prefix, "fisher_16w_covid.txt")))
#

### 6.4.2 40 weeks    ----
fisher_40w_covid <- participant_strata %>%
  select(ParticipantID, covid_event_by_40w) %>% 
  left_join(ae_incidence_40w, by = "ParticipantID") %>%
  filter(!is.na(covid_event_by_40w)) %>% 
  group_by(AE) %>%
  summarise(
    covid_yes   = sum(covid_event_by_40w == "yes"   & AE_occurred == 1),
    covid_no    = sum(covid_event_by_40w == "yes"   & AE_occurred == 0),
    noncovid_yes = sum(covid_event_by_40w == "no" & AE_occurred == 1),
    noncovid_no  = sum(covid_event_by_40w == "no" & AE_occurred == 0),
    .groups = "drop"
  ) %>%
  mutate(
    fisher = pmap(
      list(covid_yes, covid_no, noncovid_yes, noncovid_no), # OR = covid vs non-covid
      ~ fisher.test(matrix(c(..1, ..2, ..3, ..4), nrow = 2, byrow = TRUE))
    ),
    p_value    = map_dbl(fisher, ~ .x$p.value),
    odds_ratio = map_dbl(fisher, ~ unname(.x$estimate)),
    conf_low   = map_dbl(fisher, ~ .x$conf.int[1]),
    conf_high  = map_dbl(fisher, ~ .x$conf.int[2])
  ) %>%
  select(AE, covid_yes, covid_no, noncovid_yes, noncovid_no, odds_ratio, conf_low, conf_high, p_value) %>%
  arrange(p_value)
fisher_40w_covid
#
fisher_40w_covid %>% write_tsv(file = here("results", paste0(out_file_prefix, "fisher_40w_covid.txt")))
#

# 7.0 Select plots    ----
## 7.1 Stacked bar chart: Most common adverse events by grade (16 weeks)   ----
tofa_ae_16w_data %>% 
  filter(AE %in% (
    irs_16w %>% filter(n_participants >= 2) %>%  pull(AE)
  )) %>% 
  filter(AE_occurred == 1) %>% 
  group_by(AE, ParticipantID) %>% 
  slice_max(order_by = Grade, n = 1, with_ties = FALSE) %>% 
  ungroup() %>% 
  select(AE, Grade) %>%
  group_by(AE) %>% 
  count(Grade) %>%
  ungroup() %>% 
  pivot_wider(names_from = Grade, 
              values_from = n,
              names_prefix = "Grade_",
              values_fill = 0) %>% 
  mutate(
    no_ae = total_n_16w - (Grade_1 + Grade_2 + Grade_3),  # Calculate "no_ae"
    severity_score = Grade_2 * 2 + Grade_1  # Assign weights to grades (Grade_2 > Grade_1)
  ) %>%  
  pivot_longer(cols = c(-AE, -severity_score),
               names_to = "Grade",
               values_to = "count") %>% 
  mutate(Grade = fct_relevel(Grade, rev(c("no_ae", "Grade_1", "Grade_2", "Grade_3"))))  %>%
  group_by(AE) %>%
  mutate(no_ae_count = count[Grade == "no_ae"]) %>%
  ungroup() %>% 
  mutate(AE = fct_relevel(AE, c(
    irs_16w %>% filter(n_participants >= 2) %>%  pull(AE)
  ))) %>% 
  ggplot(aes(x = AE, y = count, fill = Grade)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5, colour = "black") +
  labs(title = "AE Grade Distribution by Adverse Event",
       subtitle = "16 weeks; possibly related",
       x = "Adverse Event", 
       y = "Count") +
  scale_fill_manual(values = c("no_ae" = "#fff7ec","Grade_1" = "#fdd49e", "Grade_2" = "#fc8d59", "Grade_3" = "#d7301f")) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    aspect.ratio = 0.4,
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Subtitle font size and style
    axis.line = element_line(colour = "black")
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "stacked_bar_chart_IncidenceCounts_16weeks", ".pdf")), device = cairo_pdf, width = 15, height = 6.5, units = "in")
#

## 7.2 Bar chart: Exposure adjusted event rates (16 weeks)   ----
eaers_16w %>% 
  filter(n_events >= 2) %>% 
  select(AE, EAER_per_100py) %>% 
  mutate(Event_Name = "16 week") %>% 
  arrange(-EAER_per_100py) %>% 
  mutate(AE = fct_inorder(AE)) %>% 
  ggplot(aes(x = AE, y = EAER_per_100py, fill = Event_Name)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(
    title = "Exposure Adjusted Event Rates",
    subtitle = "16 weeks; AEs with incidence in >= 2",
  ) +
  scale_fill_manual(values = c("16 week" = "#6baed6"
                               # "40 week" = "#4292c6"
  )) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    aspect.ratio = 0.4,
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),  # Subtitle font size and style
    axis.line = element_line(colour = "black")
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "bar_chart_EAER_16weeks_filtered_by_n2", ".pdf")), device = cairo_pdf, width = 15, height = 6.5, units = "in")
#

## 7.3 Pie chart: Proportion of adverse events by grade (16 weeks)   ----
tofa_ae_16w_data %>% 
  filter(AE_occurred == 1) %>% 
  count(Grade) %>% 
  mutate(total = sum(n),
         proportion = n/total,
         Grade = paste0("Grade ", Grade) %>% as.character(),
         Grade = fct_relevel(Grade, c("Grade 1", "Grade 2", "Grade 3")))  %>% 
  ggplot(aes(x = "", y = proportion, fill = Grade)) +
  geom_bar(width = 1, stat = "identity", colour = "black", alpha = 0.8) +
  coord_polar("y") +
  theme_void() +  # Remove unnecessary axes
  labs(title = "AE Grade Distributions") +
  scale_fill_manual(values = c("Grade 1" = "#fdd49e", "Grade 2" = "#fc8d59", "Grade 3" = "#d7301f")) +
  geom_text(aes(label = sprintf("%.3f", proportion)), 
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme(
    strip.text = element_text(size = 3)  # Control facet label (AE) size
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "pie_chart_AEs_by_grade_16w", ".pdf")), device = cairo_pdf, width = 4, height = 3.5, units = "in")
#

## 7.4 Heatmap: Incidence 16 vs. 40 weeks    ----
### 7.4.1 Combine and complete data object for Incidence Rates at 16 and 40 weeks   ----
irs_dat <- bind_rows(
  irs_16w %>% mutate(timepoint = "16_weeks"),
  irs_40w %>% mutate(timepoint = "40_weeks")
) %>%
  complete(AE, timepoint = c("16_weeks", "40_weeks")) %>% 
  mutate(n_participants = if_else(is.na(n_participants), 0, n_participants),
         incidence = if_else(is.na(incidence), 0, incidence))
irs_dat
#

### 7.4.2 Heatmap for Incidence at 16 and 40 weeks    ----
irs_dat %>% 
  mutate(timepoint = fct_relevel(timepoint, c("16_weeks", "40_weeks"))) %>%
  mutate(AE = fct_relevel(AE, c(
    irs_dat %>% 
      filter(timepoint == "16_weeks") %>% 
      arrange(-incidence) %>% 
      pull(AE)
  ))) %>% 
  tidy_hm(
    row = timepoint,
    col = AE,
    value = incidence,
    rotate_text = TRUE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    row_size = 5,
    col_size = 5,
    scale_type = "continuous",
    palette = "Reds",
    palette_rev = FALSE,
    col_dend = FALSE
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "heatmap_incidence_16_vs_40", ".pdf")), device = cairo_pdf, width = 25, height = 5, units = "in")
#

## 7.5 Heatmap: EAER 16 vs. 40 weeks   ----
### 7.5.1 Combine and complete data object for EAERs at 16 and 40 weeks   ----
eaer_dat <- bind_rows(
  eaers_16w %>% mutate(timepoint = "16_weeks"),
  eaers_40w %>% mutate(timepoint = "40_weeks")
) %>%
  complete(AE, timepoint = c("16_weeks", "40_weeks"))
eaer_dat
#

### 7.5.2 Heatmap for EAER at 16 and 40 weeks    ----
eaer_dat %>% 
  mutate(AE = fct_relevel(AE, c(
    irs_dat %>% # Factor to match order of Incidence heatmap (see 7.4) 
      filter(timepoint == "16_weeks") %>% 
      arrange(-incidence) %>% 
      pull(AE)
  ))) %>% 
  tidy_hm(
    row = timepoint,
    col = AE,
    value = EAER_per_100py,
    rotate_text = TRUE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    row_size = 5,
    col_size = 5,
    scale_type = "continuous",
    palette = "Blues",
    palette_rev = FALSE,
    col_dend = FALSE
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "heatmap_eaer_16_vs_40", ".pdf")), device = cairo_pdf, width = 25, height = 5, units = "in")
#

## 7.6 Forest plot: AEs in participants with COVID Events vs. those without (40 weeks)   ----
fisher_40w_covid %>% 
  arrange(p_value) %>% 
  filter(!AE %in% c(
    "Infections and infestations - Other, specify", # Removed because mapped directly to COVID
    "CPK increased", # Removed because inf OR and p > 0.05
    "Chills" # Removed because inf OR and p > 0.05
  )) %>% 
  mutate(
    filter_out = case_when(
      AE == "Productive cough" ~ "No",
      !is.finite(odds_ratio) ~ "Yes",
      TRUE ~ "No"
    )
  ) %>% 
  filter(filter_out == "No") %>% 
  head(8) %>% 
  arrange(odds_ratio, p_value) %>%
  mutate(
    AE = factor(AE, levels = unique(AE)),
    color = if_else(p_value < 0.05, "p < 0.05", "n.s.")
  ) %>% 
  ggplot(aes(odds_ratio, AE)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
  geom_errorbar(aes(xmin = conf_low, xmax = conf_high), height = 0.2) +
  geom_point(shape = 15, aes(size = -log10(p_value), color = color)) +
  scale_color_manual(values = c("p < 0.05" = "red", "n.s." = "black")) +
  scale_x_continuous(
    trans = "log",
    labels = \(x) sprintf("%.1f", x)
  ) +
  scale_size_continuous(range = c(2, 4)) +
  theme(
    aspect.ratio = 1.2,
    legend.key = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "COVID stratified AEs: 40 weeks",
    subtitle = "Fisher's exact test",
    x = "Odds ratio",
    y = NULL
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "ForestPlot_COVID_AEs_40weeks", ".pdf")), device = cairo_pdf, width = 10, height = 4, units = "in")
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