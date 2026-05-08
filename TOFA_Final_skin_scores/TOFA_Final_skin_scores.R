################################################
# Title: Tofacitinib for Immune Skin Conditions in Down Syndrome: Analysis of Dermatological/skin scores 
# Author(s):
#   - Matthew Galbraith
# Affiliation(s):
#   - Linda Crnic Institute for Down Syndrome
#   - University of Colorado Anschutz
################################################

### Summary:  
# Dermatological scores for skin conditions, assigned by board-certified dermatologists or participants and/or caregivers (DLQI).
# See README.md for more details.
#  

### Data type(s):
# Clinical trial (TOFA) datasets:
#    * Participant-level metadata; Available on request.
#    * Visit/Event-level metadata; Available on request.
#    * Baseline obesity status; Available on request.
#    * COVID-19 history; Available on request.
#    * Dermatological scores; DOI: 10.5281/zenodo.20077742
#


# 0.1 General Setup -----
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
tofa_participant_meta_data_file <- here("data/TOFA_Participant_metadata_zenodo_v1.txt") # Source: Available on request
tofa_visit_meta_data_file <- here("data/TOFA_Visit_metadata_zenodo_v1.txt") # Source: Available on request
tofa_baseline_obesity_file <- here("data/TOFA_Baseline_Obesity_Status_zenodo_v1.txt") # Source: Available on request
tofa_covid_history_file <- here("data/TOFA_COVID_History_zenodo_v1.txt") # Source: Available on request
tofa_skin_scores_data_file <- here("data/TOFA_Skin_scores_data_zenodo_v1.txt.gz") # Source: 10.5281/zenodo.20077742
#
# Other parameters:
standard_colors <- c("Baseline" = "#999999", "2 week" = "#c6dbef", "8 week" = "#9ecae1", "16 week" = "#6baed6", "40 week" = "#4292c6")
#
out_file_prefix <- "TOFA_Final_skin_scores_" # should match this script title
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

## 1.2 Read in TOFA Skin scores data ----
tofa_skin_scores_data <- tofa_skin_scores_data_file |>
  read_tsv() |> 
  mutate(
    Event_Name = fct_relevel(Event_Name, c("Baseline", "8 week", "16 week","40 week")) # set factor levels
  )
#
tofa_skin_scores_data # 486 rows
tofa_skin_scores_data |> distinct(ParticipantID) # 42 Participants
tofa_skin_scores_data |> distinct(VisitID) # 157 VisitIDs
tofa_skin_scores_data |> distinct(Score_name) # 5 scores
#

## 1.3 Join with meta data, filter eligible etc, filter to endpoint cytokines ----
tofa_skin_scores_data <- tofa_skin_scores_data |> 
  inner_join(tofa_visit_meta_data) |> # returns 486 of 486 rows
  inner_join(tofa_participant_meta_data) |> # returns 486 of 486 rows
  filter(Endpoint_eligible == TRUE) # returns 486 of 486 rows
#
### 1.3.1 ParticipantID vs Event_Name summary ---- 
tofa_skin_scores_data |> 
  distinct(ParticipantID, VisitID, Event_Name) |> 
  arrange(Event_Name) |> 
  pivot_wider(names_from = Event_Name, values_from = VisitID) |> 
  print(n = Inf)
# Note: Not all scores were collected at all time points.
# Note: Not all Participants were assessed for all scores and not all completed extension to 40 weeks.
# No DLQI scores at 8 weeks.


# 2 Calculate 'Percent Baseline' scores ----
tofa_skin_scores_perc_baseline <- tofa_skin_scores_data |> 
  filter(Event_Name %in% c("Baseline")) |> 
  select(ParticipantID, Score_name, Baseline_value = Score_value, Qualifying_Skin_Condition) |> 
  inner_join(
    tofa_skin_scores_data |> 
      filter(Event_Name %in% c("Baseline", "8 week", "16 week", "40 week")) |> # keep baseline to allow paired tests
      select(ParticipantID, VisitID, Event_Name, Score_name, Score_value)
  ) |> 
  mutate(
    perc_baseline = Score_value / Baseline_value * 100,
    # Fix NaN values when Baseline_value and Score_value are both 0:
    perc_baseline = if_else(is.nan(perc_baseline), 100, perc_baseline),
    # Combine disease specific qualifying scores for endpoint:
    Score_combined = case_when(
      Score_name == "IGA" ~ "IGA_perc_baseline",
      Score_name == "DLQI" ~ "DLQI_perc_baseline",
      .default = "Combined_disease_scores_perc_baseline"
    ),
    Score_combined = fct_relevel(Score_combined, c("IGA_perc_baseline", "DLQI_perc_baseline", "Combined_disease_scores_perc_baseline"))
  )
#


# 3 Endpoint Stats - 'Percent baseline' ----
## 3.1 Check data distributions ----
### 3.1.1 check for outlier values ----
tofa_skin_scores_perc_baseline |> 
  mutate(
    extreme = rstatix::is_extreme(perc_baseline),
    outlier = rstatix::is_outlier(perc_baseline),
    .by = Score_combined
  ) |> 
  count(Score_combined, outlier, extreme)
#
tofa_skin_scores_perc_baseline |> 
  mutate(
    extreme = rstatix::is_extreme(perc_baseline),
    outlier = rstatix::is_outlier(perc_baseline),
    .by = Score_combined
  ) |> 
  ggplot(aes("perc_baseline", perc_baseline)) +
  geom_sina(data = . %>% filter(outlier == FALSE & extreme == FALSE), maxwidth = 0.5, color = "grey") +
  geom_sina(data = . %>% filter(outlier == TRUE & extreme == FALSE), aes(color = "outlier"), maxwidth = 0.1) +
  geom_sina(data = . %>% filter(extreme == TRUE), aes(color = "extreme"), maxwidth = 0.1) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~ Score_combined, scales = "free_y", nrow = 1) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    aspect.ratio = 1.3
  ) +
  labs(title = "Skin score perc baseline outliers")
#
### 3.1.2 Calc Differences from baseline ----
tofa_skin_perc_baseline_diffs_16_weeks <- tofa_skin_scores_perc_baseline |> 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(perc_baseline)) |> 
  filter(is.finite(perc_baseline)) |> 
  add_count(Score_combined, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  select(ParticipantID, Event_Name, Score_combined, perc_baseline) |> 
  pivot_wider(names_from = Event_Name, values_from = perc_baseline) |> 
  mutate(difference = `16 week` - Baseline) |> 
  select(ParticipantID, Score_combined, difference) |> 
  mutate(
    extreme = rstatix::is_extreme(difference),
    outlier = rstatix::is_outlier(difference),
    .by = Score_combined
  )
#
### 3.1.3 Count extreme outliers in differences ----
tofa_skin_perc_baseline_diffs_16_weeks |> 
  count(Score_combined, outlier, extreme)
#
### 3.1.4 Plot differences ----
tofa_skin_perc_baseline_diffs_16_weeks |> 
  ggplot(aes(Score_combined, difference, color = Score_combined)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_sina(data = . %>% filter(outlier == FALSE & extreme == FALSE), maxwidth = 0.5, color = "grey") +
  geom_sina(data = . %>% filter(outlier == TRUE & extreme == FALSE), aes(color = "outlier"), maxwidth = 0.1) +
  geom_sina(data = . %>% filter(extreme == TRUE), aes(color = "extreme"), maxwidth = 0.1) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~ Score_combined, scales = "free", nrow = 1) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    aspect.ratio = 1.3
  ) +
  labs(title = "16 weeks: Distributions of differences vs. baseline", x = NULL)
#
### 3.1.5 Check if differences are normally distributed ----
tofa_skin_perc_baseline_diffs_16_weeks |>
  ggpubr::ggqqplot("difference", title = "Q-Q plot: 16 week differences", facet.by = "Score_combined",  scales = "free_y")
#

## 3.2  Unpaired t test (NOT USED) ----
tofa_skin_scores_perc_baseline |> 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(perc_baseline)) |> 
  filter(is.finite(perc_baseline)) |> 
  add_count(Score_combined, ParticipantID, name = "pair") |> 
  filter(pair >= 2) |> # require pairs
  filter(Event_Name %in% c("16 week")) |> 
  group_by(Score_combined) |> 
  rstatix::t_test(
    formula = perc_baseline ~ 1,
    paired = FALSE,
    mu = 100, # Null hypothesis - no change from baseline
    var.equal = TRUE, # set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) |> 
  # need to 'reverse' estimate:
  mutate(estimate = estimate - 100)
# Same p-values as paired t test below
#

## 3.3 Paired t tests  --------
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
### 3.3.1 Run paired t tests 16 weeks ----
# Only 16 week timepoint is tested as a trial endpoint.
tofa_skin_perc_baseline_Ttest_res_16weeks <- tofa_skin_scores_perc_baseline |> 
  arrange(ParticipantID, Event_Name) |> # ensure correct sort order for rstatix::t_test()
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(perc_baseline)) |> 
  filter(is.finite(perc_baseline)) |> 
  add_count(Score_combined, ParticipantID, name = "pair") |> 
  filter(pair >= 2) |>  # require pairs
  group_by(Score_combined) |> 
  rstatix::t_test(
    formula = perc_baseline ~ Event_Name,
    ref.group = "16 week",
    paired = TRUE,
    var.equal = TRUE, # set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  )
# Same p-values as unpaired t test above
#

### 3.3.2 Paired t test effect sizes ----
# The effect size for a paired-samples t-test can be calculated by dividing the
# mean difference by the standard deviation of the difference, as shown below.
# Cohen’s formula:
# d = mean(D)/sd(D), where D is the differences of the paired samples values.
tofa_skin_perc_baseline_Ttest_effsize_16weeks <- tofa_skin_scores_perc_baseline |> 
  arrange(ParticipantID, Event_Name) |> # ensure correct sort order for rstatix::cohens_d()
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(perc_baseline)) |> 
  filter(is.finite(perc_baseline)) |> 
  add_count(Score_combined, ParticipantID, name = "pair") |> 
  filter(pair >= 2) |>  # require pairs
  group_by(Score_combined) |> 
  rstatix::cohens_d(
    formula = perc_baseline ~ Event_Name,
    ref.group = "16 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    var.equal = TRUE, # set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
  )
#

### 3.3.3 Compile t test results ----
# 16 weeks
tofa_skin_perc_baseline_Ttest_res_16weeks_full <- tofa_skin_perc_baseline_Ttest_res_16weeks |> 
  inner_join(tofa_skin_perc_baseline_Ttest_effsize_16weeks) |> 
  mutate(
    mean_diff = estimate,
  ) |> 
  select(Score_combined, .y., group1, group2, mean_diff, p, effsize, magnitude, everything())
#
tofa_skin_perc_baseline_Ttest_res_16weeks_full
#

### 3.3.4 Export results ----
list(
  "ttests_perc_baseline" = tofa_skin_perc_baseline_Ttest_res_16weeks_full |> 
    select(
      Score_name = Score_combined,
      Timepoint = group1,
      n_pairs = n1,
      Mean_difference = mean_diff,
      Conf.low = conf.low,
      Conf.high = conf.high,
      Statistic = statistic,
      Effect_size = effsize,
      Magnitude = magnitude,
      pvalue = p
    ) |> 
    mutate(
      qvalue = "", # placeholder, not used here
      Method = "Paired Student's t test, two-sided"
    )
) |> 
  export_excel(filename = "Ttest_results_perc_baseline")
#


# 4 Stats - Individual score values ----
## 4.1 Check data distributions ----
### 4.1.1 check for outlier values ----
tofa_skin_scores_data |> 
  mutate(
    extreme = rstatix::is_extreme(Score_value),
    outlier = rstatix::is_outlier(Score_value),
    .by = Score_name
  ) |> 
  count(Score_name, outlier, extreme)
#
tofa_skin_scores_data |> 
  mutate(
    extreme = rstatix::is_extreme(Score_value),
    outlier = rstatix::is_outlier(Score_value),
    .by = Score_value
  ) |> 
  ggplot(aes("Score_value", Score_value)) +
  geom_sina(data = . %>% filter(outlier == FALSE & extreme == FALSE), maxwidth = 0.5, color = "grey") +
  geom_sina(data = . %>% filter(outlier == TRUE & extreme == FALSE), aes(color = "outlier"), maxwidth = 0.1) +
  geom_sina(data = . %>% filter(extreme == TRUE), aes(color = "extreme"), maxwidth = 0.1) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~ Score_name, scales = "free_y", nrow = 1) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    aspect.ratio = 1.3
  ) +
  labs(title = "Skin score outliers")
#
### 4.1.2 Calc Differences from baseline ----
tofa_skin_diffs_16_weeks <- tofa_skin_scores_data |> 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(Score_value)) |> 
  filter(is.finite(Score_value)) |> 
  add_count(Score_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  select(ParticipantID, Event_Name, Score_name, Score_value) |> 
  pivot_wider(names_from = Event_Name, values_from = Score_value) |> 
  mutate(difference = `16 week` - Baseline) |> 
  select(ParticipantID, Score_name, difference) |> 
  mutate(
    extreme = rstatix::is_extreme(difference),
    outlier = rstatix::is_outlier(difference),
    .by = Score_name
  )
#
### 4.1.3 Count extreme outliers in differences ----
tofa_skin_diffs_16_weeks |> 
  count(Score_name, outlier, extreme)
#
### 4.1.4 Plot differences ----
tofa_skin_diffs_16_weeks |> 
  ggplot(aes(Score_name, difference, color = Score_name)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_sina(data = . %>% filter(outlier == FALSE & extreme == FALSE), maxwidth = 0.5, color = "grey") +
  geom_sina(data = . %>% filter(outlier == TRUE & extreme == FALSE), aes(color = "outlier"), maxwidth = 0.1) +
  geom_sina(data = . %>% filter(extreme == TRUE), aes(color = "extreme"), maxwidth = 0.1) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~ Score_name, scales = "free", nrow = 1) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    aspect.ratio = 1.3
  ) +
  labs(title = "16 weeks: Distributions of differences vs. baseline", x = NULL)
#
### 4.1.5 Check if differences are normally distributed ----
tofa_skin_diffs_16_weeks |>
  ggpubr::ggqqplot("difference", title = "Q-Q plot: 16 week differences", facet.by = "Score_name",  scales = "free_y")
#

## 4.2 Paired t tests  --------
### 4.2.1 Run paired t tests 16 weeks ----
tofa_skin_scores_Ttest_res_16weeks <- tofa_skin_scores_data |> 
  arrange(ParticipantID, Event_Name) |> # ensure correct sort order for rstatix::t_test()
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(Score_value)) |> 
  filter(is.finite(Score_value)) |> 
  add_count(Score_name, ParticipantID, name = "pair") |> 
  filter(pair >= 2) |>  # require pairs
  group_by(Score_name) |> 
  rstatix::t_test(
    formula = Score_value ~ Event_Name,
    ref.group = "16 week",
    paired = TRUE,
    var.equal = TRUE, # set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) |> 
  mutate(BH_padj = p.adjust(p, method = "BH"))
#

### 4.2.2 Paired t test effect sizes ----
tofa_skin_scores_Ttest_effsize_16weeks <- tofa_skin_scores_data |> 
  arrange(ParticipantID, Event_Name) |> # ensure correct sort order for rstatix::cohens_d()
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(Score_value)) |> 
  filter(is.finite(Score_value)) |> 
  add_count(Score_name, ParticipantID, name = "pair") |> 
  filter(pair >= 2) |>  # require pairs
  group_by(Score_name) |> 
  rstatix::cohens_d(
    formula = Score_value ~ Event_Name,
    ref.group = "16 week", # Seems to be the reverse of how we typically do this
    paired = TRUE,
    var.equal = TRUE, # set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
  )
#

### 4.2.3 Compile t test results ----
# 16 weeks
tofa_skin_scores_Ttest_res_16weeks_full <- tofa_skin_scores_Ttest_res_16weeks |> 
  inner_join(tofa_skin_scores_Ttest_effsize_16weeks) |> 
  mutate(
    mean_diff = estimate,
  ) |> 
  select(Score_name, .y., group1, group2, mean_diff, p, BH_padj, effsize, magnitude, everything()) |> 
  mutate()
#
tofa_skin_scores_Ttest_res_16weeks_full
#

### 4.2.4 Export results ----
list(
  "ttests_ind_scores" = tofa_skin_scores_Ttest_res_16weeks_full |> 
    select(
      Score_name,
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
  export_excel(filename = "Ttest_results_individual_scores")
#



# 5 Paired plots  --------
## 5.1 Baseline vs 16 weeks ----
tofa_skin_scores_data |> 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(Score_value)) |> 
  filter(is.finite(Score_value)) |> 
  add_count(Score_name, ParticipantID, name = "pair") |> 
  filter(pair >= 2) |>  # require at least 2 data points
  #
  group_by(Score_name) |> 
  group_split() |> 
  map( # using anonymous function: \(x)
    \(x) x |> 
      arrange(Score_name, Event_Name) |> 
      mutate( # add labels with n
        n = n_distinct(ParticipantID), 
        .by = c(Score_name, Event_Name),
        label = paste0(str_extract(Event_Name, "^[[:alpha:]]|^\\d{1,2}"), "\n(n=", n, ")"),
        label = fct_inorder(label)
      ) |> 
      mutate(Event_Name = fct_drop(Event_Name)) |> # prevent empty x values
      mutate(Score_name = fct_drop(Score_name)) |> # prevent empty x values
      complete(Score_name, ParticipantID, Event_Name) %>% # Prevent error from stat_signif()
      #
      {ggplot(data = ., aes(Event_Name, Score_value, color = Event_Name)) +
          scale_color_manual(values = standard_colors) +
          geom_line(aes(group = ParticipantID), color = "grey90") +
          geom_point() +
          stat_summary(fun = "median", fun.min = "median", fun.max= "median", geom = "crossbar", size= 0.3, width = 0.4, color = "black", show.legend = FALSE) +
          facet_wrap(~ Score_name, scales = "free_y", nrow = 1) +
          theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            aspect.ratio = 1.5,
            panel.border = element_blank(),
            axis.line = element_line(colour = "black")
          ) +
          labs(
            x = NULL
          ) +
          scale_x_discrete(
            labels = . |>
              filter(!is.na(label)) |>
              distinct(Score_name, Event_Name, label) |>
              select(Event_Name, label) |>
              transmute(as.character(Event_Name), as.character(label)) |>
              deframe()
          )
      }
  ) |> 
  patchwork::wrap_plots() +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(
    title = "Skin scores: Baseline vs. Treatment"
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "paired_plots_B_16_weeks", ".pdf")), device = cairo_pdf, width = 20, height = 6, units = "in")
#



# 6 Sina plots - score values  --------
## 6.1 Baseline vs 16 weeks ----
tofa_skin_scores_data |> 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(Score_value)) |> 
  filter(is.finite(Score_value)) |> 
  add_count(Score_name, ParticipantID, name = "pair") |> 
  filter(pair >= 2) |>  # require at least 2 data points
  #
  group_by(Score_name) |> 
  group_split() |> 
  map( # using anonymous function: \(x)
    \(x) x |> 
      arrange(Score_name, Event_Name) |> 
      mutate( # add labels with n
        n = n_distinct(ParticipantID), 
        .by = c(Score_name, Event_Name),
        label = paste0(str_extract(Event_Name, "^[[:alpha:]]|^\\d{1,2}"), "\n(n=", n, ")"),
        label = fct_inorder(label)
      ) |> 
      mutate(Event_Name = fct_drop(Event_Name)) |> # prevent empty x values
      mutate(Score_name = fct_drop(Score_name)) |> # prevent empty x values
      complete(Score_name, ParticipantID, Event_Name) %>% # Prevent error from stat_signif()
      #
      {ggplot(data = ., aes(Event_Name, Score_value, color = Event_Name)) +
          scale_color_manual(values = standard_colors) +
          geom_line(aes(group = ParticipantID), color = "grey90") +
          geom_sina(maxwidth = 0.3) +
          geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
          # stat_summary(fun = "median", fun.min = "median", fun.max= "median", geom = "crossbar", size= 0.3, width = 0.4, color = "black", show.legend = FALSE) +
          facet_wrap(~ Score_name, scales = "free_y", nrow = 1) +
          theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            aspect.ratio = 1.5,
            panel.border = element_blank(),
            axis.line = element_line(colour = "black")
          ) +
          labs(
            x = NULL
          ) +
          
          # add sample numbers to labels
          scale_x_discrete(
            labels = . |>
              filter(!is.na(label)) |>
              distinct(Score_name, Event_Name, label) |>
              select(Event_Name, label) |>
              transmute(as.character(Event_Name), as.character(label)) |>
              deframe()
          )
      }
  ) |> 
  patchwork::wrap_plots() +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(
    title = "Skin scores: Baseline vs. Treatment"
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_plots_B_16_weeks", ".pdf")), device = cairo_pdf, width = 20, height = 6, units = "in")
#
## 6.2 Baseline vs 8 weeks vs 16 vs 40 weeks ----
tofa_skin_scores_data |> 
  filter(Event_Name %in% c("Baseline", "8 week", "16 week", "40 week")) |> 
  filter(!is.na(Score_value)) |> 
  filter(is.finite(Score_value)) |> 
  add_count(Score_name, ParticipantID, name = "pair") |> 
  filter(pair >= 2) |>  # require at least 2 data points
  #
  group_by(Score_name) |> 
  group_split() |> 
  map( # using anonymous function: \(x)
    \(x) x |> 
      arrange(Score_name, Event_Name) |> 
      mutate( # add labels with n
        n = n_distinct(ParticipantID), 
        .by = c(Score_name, Event_Name),
        label = paste0(str_extract(Event_Name, "^[[:alpha:]]|^\\d{1,2}"), "\n(n=", n, ")"),
        label = fct_inorder(label)
      ) |> 
      mutate(Event_Name = fct_drop(Event_Name)) |> # prevent empty x values
      mutate(Score_name = fct_drop(Score_name)) |> # prevent empty x values
      complete(Score_name, ParticipantID, Event_Name) %>% # Prevent error from stat_signif()
      #
      {ggplot(data = ., aes(Event_Name, Score_value, color = Event_Name)) +
          scale_color_manual(values = standard_colors) +
          geom_line(aes(group = ParticipantID), color = "grey90") +
          geom_sina(maxwidth = 0.5) +
          geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
          facet_wrap(~ Score_name, scales = "free_y", nrow = 1) +
          theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            aspect.ratio = 1.5,
            panel.border = element_blank(),
            axis.line = element_line(colour = "black")
          ) +
          labs(
            x = NULL
          ) +
          
          # add sample numbers to labels
          scale_x_discrete(
            labels = . |>
              filter(!is.na(label)) |>
              distinct(Score_name, Event_Name, label) |>
              select(Event_Name, label) |>
              transmute(as.character(Event_Name), as.character(label)) |>
              deframe()
          )
      }
  ) |> 
  patchwork::wrap_plots() +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(
    title = "Skin scores: Baseline vs. Treatment"
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_plots_B_8_16_40_weeks", ".pdf")), device = cairo_pdf, width = 20, height = 5, units = "in")
#


# 7 Sina plots - 'Percent Baseline'
## 7.1 Sina plots 'Percent Baseline' combined MSS/SALT ----
# Have to combine separate filters to get correct numbers for combined disease specific groups
bind_rows(
  tofa_skin_scores_perc_baseline |> 
    filter(Event_Name %in% c("Baseline")) |> 
    filter(!is.na(perc_baseline)) |>
    filter(is.finite(perc_baseline)),
  tofa_skin_scores_perc_baseline |> 
    filter(Event_Name %in% c("Baseline", "8 week")) |> 
    filter(!is.na(perc_baseline)) |>
    filter(is.finite(perc_baseline)) |> 
    add_count(Score_combined, ParticipantID, name = "pair") |> 
    filter(pair >= 2) |>  # require pairs
    add_count(Event_Name, Score_combined) |> 
    group_by(Score_combined) |> 
    filter(!any(n < 2)) |> 
    filter(Event_Name != "Baseline"),
  tofa_skin_scores_perc_baseline |> 
    filter(Event_Name %in% c("Baseline", "16 week")) |> 
    filter(!is.na(perc_baseline)) |>
    filter(is.finite(perc_baseline)) |> 
    add_count(Score_combined, ParticipantID, name = "pair") |> 
    filter(pair >= 2) |>  # require pairs
    add_count(Event_Name, Score_combined) |> 
    group_by(Score_combined) |> 
    filter(!any(n < 2)) |> 
    filter(Event_Name != "Baseline"),
  tofa_skin_scores_perc_baseline |> 
    filter(Event_Name %in% c("Baseline", "40 week")) |> 
    filter(!is.na(perc_baseline)) |> 
    filter(is.finite(perc_baseline)) |> 
    # # CUSTOM - need to remove TOFA0072v2 as no matching v10E causes t_test error for Disease_specific
    # filter(TOFA_LabID != "TOFA0072v2") |> 
    add_count(Score_combined, ParticipantID, name = "pair") |> 
    filter(pair >= 2) |>  # require pairs
    add_count(Event_Name, Score_combined) |> 
    group_by(Score_combined) |> 
    filter(!any(n < 2)) |> 
    filter(Event_Name != "Baseline")
) |> 
  group_by(Score_combined) |> 
  group_split() |> 
  map( # using anonymous function: \(x)
    \(x) x |> 
      arrange(Score_combined, Event_Name) |> 
      mutate( # add labels with n
        n = n_distinct(ParticipantID, Score_name), # NOTE the difference here
        .by = c(Score_combined, Event_Name),
        label = paste0(str_extract(Event_Name, "^[[:alpha:]]|^\\d{1,2}"), "\n(n=", n, ")"),
        label = fct_inorder(label)
      ) |> 
      mutate(Event_Name = fct_drop(Event_Name)) |> # prevent empty x values
      mutate(Score_combined = fct_drop(Score_combined)) |> # prevent empty x values
      complete(Score_combined, ParticipantID, Event_Name) %>% # Prevent error from stat_signif()
      #
      {ggplot(data = ., aes(Event_Name, perc_baseline, color = Event_Name)) +
          geom_hline(yintercept = 100, linetype = 2) +
          geom_line(aes(group = ParticipantID), color = "grey90") +
          geom_sina(maxwidth = 0.5) +
          scale_color_manual(values = standard_colors) +
          geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
          facet_wrap(~ Score_combined, scales = "free_y", nrow = 1) +
          theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            aspect.ratio = 1.5,
            panel.border = element_blank(),
            axis.line = element_line(colour = "black")
          ) +
          labs(x = NULL) +
          # add sample numbers to labels
          scale_x_discrete(
            labels = . |>
              filter(!is.na(label)) |>
              distinct(Score_combined, Event_Name, label) |>
              select(Event_Name, label) |>
              transmute(as.character(Event_Name), as.character(label)) |>
              deframe()
          )
      }
  ) |> 
  patchwork::wrap_plots() +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(
    title = "Skin scores: Percent Baseline",
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_plots_combined_PERC_BASELINE_8_16_40_weeks", ".pdf")), device = cairo_pdf, width = 20, height = 6, units = "in")
#


# 8 Improvement threshold proportions -----
## 8.1 IGA 0/1 ----
# % with IGA 0 or 1 at each time point 
IGA_0_1_long <- tofa_skin_scores_data |>
  filter(Event_Name %in% c("Baseline", "8 week", "16 week", "40 week")) |>
  mutate(Event_Name = fct_drop(Event_Name)) |> 
  select(Event_Name, Score_name, Score_value) |>
  filter(Score_name == "IGA", !is.na(Score_value)) |>
  count(Event_Name, response = Score_value <= 1) |>
  complete(Event_Name, response, fill = list(n = 0)) |>
  mutate(
    perc = n / sum(n) * 100,
    .by = Event_Name
  ) |>
  filter(response == TRUE) |> 
  select(-response) |> 
  mutate(measure = "IGA_0_1", .before = everything())
#
## 8.2 DLQI 0/1 ----
# % with DLQI 0 or 1 at each time point 
DLQI_0_1_long <- tofa_skin_scores_data |>
  filter(Event_Name %in% c("Baseline", "16 week", "40 week")) |> # not collected at week 8
  mutate(Event_Name = fct_drop(Event_Name)) |> 
  select(Event_Name, Score_name, Score_value) |>
  filter(Score_name == "DLQI", !is.na(Score_value)) |>
  count(Event_Name, response = Score_value <= 1) |>
  complete(Event_Name, response, fill = list(n = 0)) |>
  mutate(
    perc = n / sum(n) * 100,
    .by = Event_Name
  ) |>
  filter(response == TRUE) |> 
  select(-response) |> 
  mutate(measure = "DLQI_0_1", .before = everything())
#
## 8.3 MSS 50 and 90 ----
# % with with at least 50 or 90% improvement at each time point 
MSS_50_long <- tofa_skin_scores_perc_baseline |>
  filter(Event_Name %in% c("8 week", "16 week", "40 week")) |> # no Baseline for % change
  mutate(Event_Name = fct_drop(Event_Name)) |> 
  select(Event_Name, Score_name, perc_baseline) |>
  filter(Score_name == "MSS", !is.na(perc_baseline)) |>
  count(Event_Name, response = perc_baseline <= 50) |>
  complete(Event_Name, response, fill = list(n = 0)) |>
  mutate(
    perc = n / sum(n) * 100,
    .by = Event_Name
  ) |>
  filter(response == TRUE) |> 
  select(-response) |> 
  mutate(measure = "MSS_50", .before = everything())
#
MSS_90_long <- tofa_skin_scores_perc_baseline |>
  filter(Event_Name %in% c("8 week", "16 week", "40 week")) |> # no Baseline for % change
  mutate(Event_Name = fct_drop(Event_Name)) |> 
  select(Event_Name, Score_name, perc_baseline) |>
  filter(Score_name == "MSS", !is.na(perc_baseline)) |>
  count(Event_Name, response = perc_baseline <= 10) |>
  complete(Event_Name, response, fill = list(n = 0)) |>
  mutate(
    perc = n / sum(n) * 100,
    .by = Event_Name
  ) |>
  filter(response == TRUE) |> 
  select(-response) |> 
  mutate(measure = "MSS_90", .before = everything())
#
## 8.4 SALT 50 and 90 ----
# % with with at least 50 or 90% improvement at each time point 
SALT_50_long <- tofa_skin_scores_perc_baseline |>
  filter(Event_Name %in% c("8 week", "16 week", "40 week")) |> # no Baseline for % change
  mutate(Event_Name = fct_drop(Event_Name)) |> 
  select(Event_Name, Score_name, perc_baseline) |>
  filter(Score_name == "SALT", !is.na(perc_baseline)) |>
  count(Event_Name, response = perc_baseline <= 50) |>
  complete(Event_Name, response, fill = list(n = 0)) |>
  mutate(
    perc = n / sum(n) * 100,
    .by = Event_Name
  ) |>
  filter(response == TRUE) |> 
  select(-response) |> 
  mutate(measure = "SALT_50", .before = everything())
#
SALT_90_long <- tofa_skin_scores_perc_baseline |>
  filter(Event_Name %in% c("8 week", "16 week", "40 week")) |> # no Baseline for % change
  mutate(Event_Name = fct_drop(Event_Name)) |> 
  select(Event_Name, Score_name, perc_baseline) |>
  filter(Score_name == "SALT", !is.na(perc_baseline)) |>
  count(Event_Name, response = perc_baseline <= 10) |>
  complete(Event_Name, response, fill = list(n = 0)) |>
  mutate(
    perc = n / sum(n) * 100,
    .by = Event_Name
  ) |>
  filter(response == TRUE) |> 
  select(-response) |> 
  mutate(measure = "SALT_90", .before = everything())
#

## 8.5 Plot together as heatmap ----
bind_rows(
  IGA_0_1_long,
  DLQI_0_1_long,
  MSS_50_long,
  MSS_90_long,
  SALT_50_long,
  SALT_90_long
) |> 
  mutate(measure = fct_relevel(measure, rev(c("IGA_0_1", "DLQI_0_1", "MSS_50", "MSS_90", "SALT_50", "SALT_90")))) |> 
  ggplot(aes(Event_Name, measure, fill = perc)) +
  geom_tile() +
  geom_text(
    aes(label = paste0(round(perc), "%")),
    size = 3
  ) +
  scale_fill_distiller(palette = "Blues", direction = 0) +
  theme(aspect.ratio = 1) +
  labs(
    title = "",
    x = NULL, y = NULL
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "heatmap_improvement_threshold_perc", ".pdf")), device = cairo_pdf, width = 8, height = 3, units = "in")
#


# 9 Mixed effects linear regression models (non-stratified) ----
# Time as a categorical variable
#   No assumptions about trajectory
#   Directly compares each visit to baseline
# Random effects
#   Handles missing at random (MAR) data.
#   Accounts for within-subject correlation.
## 9.1 Set up models ----
### 9.1.1 Mixed effects LM: Event_name + 1|ParticipantID ----
tofa_skin_scores_lm_mixedParticipantID <- tofa_skin_scores_data |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Score_value)) |>
  filter(is.finite(Score_value)) |>
  select(Score_name, ParticipantID, Sex, Age_Baseline, Event_Name, Score_value) |> 
  nest(data = -Score_name) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Score_value ~ Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
### 9.1.2 Mixed effects LM: Sex + Event_name + 1|ParticipantID ----
tofa_skin_scores_lm_fixedSex_mixedParticipantID <- tofa_skin_scores_data |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Score_value)) |>
  filter(is.finite(Score_value)) |>
  select(Score_name, ParticipantID, Sex, Age_Baseline, Event_Name, Score_value) |> 
  nest(data = -Score_name) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Score_value ~ Sex + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
### 9.1.3 Mixed effects LM: Age + Event_name + 1|ParticipantID ----
tofa_skin_scores_lm_fixedAge_mixedParticipantID <- tofa_skin_scores_data |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Score_value)) |>
  filter(is.finite(Score_value)) |>
  select(Score_name, ParticipantID, Sex, Age_Baseline, Event_Name, Score_value) |> 
  nest(data = -Score_name) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Score_value ~ Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  ) # 1 warning
#
### 9.1.4 Mixed effects LM: Sex + Age + Event_name + 1|ParticipantID (PREFERRED) ----
tofa_skin_scores_lm_fixedSexAge_mixedParticipantID <- tofa_skin_scores_data |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Score_value)) |>
  filter(is.finite(Score_value)) |>
  select(Score_name, ParticipantID, Sex, Age_Baseline, Event_Name, Score_value) |> 
  nest(data = -Score_name) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Score_value ~ Sex + Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#

## 9.2 Compare models ----
### 9.2.1 AIC/BIC ----
# Mixed models
tofa_skin_scores_lm_mixedParticipantID |> unnest(glanced) |> select(Score_name, AIC, BIC)
tofa_skin_scores_lm_fixedSex_mixedParticipantID |> unnest(glanced) |> select(Score_name, AIC, BIC)
tofa_skin_scores_lm_fixedAge_mixedParticipantID |> unnest(glanced) |> select(Score_name, AIC, BIC)
tofa_skin_scores_lm_fixedSexAge_mixedParticipantID |> unnest(glanced) |> select(Score_name, AIC, BIC)
#
### 9.2.2 Likelihood ratio tests ----
anova(
  tofa_skin_scores_lm_mixedParticipantID |> pluck("fit", 1),
  tofa_skin_scores_lm_fixedSex_mixedParticipantID |> pluck("fit", 1),
  test="LRT"
) |>  # refitting model(s) with ML (instead of REML)
  tidy() # p.value = 0.813
anova(
  tofa_skin_scores_lm_mixedParticipantID |> pluck("fit", 1),
  tofa_skin_scores_lm_fixedAge_mixedParticipantID |> pluck("fit", 1),
  test="LRT"
) |>  # refitting model(s) with ML (instead of REML)
  tidy() # p.value = 0.515
anova(
  tofa_skin_scores_lm_mixedParticipantID |> pluck("fit", 1),
  tofa_skin_scores_lm_fixedSexAge_mixedParticipantID |> pluck("fit", 1),
  test="LRT"
) |>  # refitting model(s) with ML (instead of REML)
  tidy() # p.value = 0.791
#

## 9.3 Model results ----
### 9.3.1 Mixed ParticipantID with fixed SexAge -----
tofa_skin_scores_lm_fixedSexAge_mixedParticipantID |> 
  unnest(tidied) |>
  select(Score_name, group, term, estimate, p.value)
#
tofa_skin_scores_lm_fixedSexAge_mixedParticipantID_results <- tofa_skin_scores_lm_fixedSexAge_mixedParticipantID |> 
  unnest(c(tidied, glanced)) |>
  filter(str_detect(term, "Event_Name")) |>
  select(Score_name, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = Score_name) |> # NB: by Score_name
  mutate(
    level = "",
    term = str_remove(term, "^Event_Name"),
    .after = Score_name,
  )
#


# 10 AQ metabolites -  Mixed effects linear regression models (stratified) ----
## 10.1 By Sex ----
lm_fixedAge_mixedParticipantID_by_Sex <- tofa_skin_scores_data |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Score_value)) |>
  filter(is.finite(Score_value)) |>
  select(Score_name, ParticipantID, Sex, Age_Baseline, Event_Name, Score_value) |> 
  nest(data = -c(Score_name, Sex)) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Score_value ~ Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  ) # 1 warning(s)
#
lm_fixedAge_mixedParticipantID_by_Sex |> unnest(glanced) |> select(Score_name, Sex, AIC, BIC)
#
lm_fixedAge_mixedParticipantID_by_Sex_results <- lm_fixedAge_mixedParticipantID_by_Sex |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, "Event_Name")) |>
  select(Score_name, level = Sex, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = c(Score_name, level)) |> # NB: by Score_name + level
  mutate(term = str_remove(term, "^Event_Name")) |> 
  arrange(level)
#

## 10.2 By Age group ----
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
lm_fixedSex_mixedParticipantID_by_Age_group <- tofa_skin_scores_data |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  inner_join(tofa_baseline_age_groups) |> # add age groups info
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Score_value)) |>
  filter(is.finite(Score_value)) |>
  select(Score_name, ParticipantID, Sex, Age_group, Event_Name, Score_value) |> 
  nest(data = -c(Score_name, Age_group)) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Score_value ~ Sex + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
lm_fixedSex_mixedParticipantID_by_Age_group |> unnest(glanced) |> select(Score_name, Age_group, AIC, BIC)
#
lm_fixedSex_mixedParticipantID_by_Age_group_results <- lm_fixedSex_mixedParticipantID_by_Age_group |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, "Event_Name")) |>
  select(Score_name, level = Age_group, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = c(Score_name, level)) |> # NB: by Score_name + level
  mutate(term = str_remove(term, "^Event_Name")) |> 
  arrange(level)
#

## 10.3 By Obesity ----
tofa_baseline_obesity <- tofa_baseline_obesity_file |> 
  read_tsv()
tofa_baseline_obesity |> count(baseline_obesity_status)
#
lm_fixedSexAge_mixedParticipantID_by_Obesity <- tofa_skin_scores_data |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  inner_join(tofa_baseline_obesity) |> # add Obesity info
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Score_value)) |>
  filter(is.finite(Score_value)) |>
  select(Score_name, ParticipantID, Sex, Age_Baseline, baseline_obesity_status, Event_Name, Score_value) |> 
  nest(data = -c(Score_name, baseline_obesity_status)) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Score_value ~ Sex + Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
lm_fixedSexAge_mixedParticipantID_by_Obesity |> unnest(glanced) |> select(Score_name, baseline_obesity_status, AIC, BIC)
#
lm_fixedSexAge_mixedParticipantID_by_Obesity_results <- lm_fixedSexAge_mixedParticipantID_by_Obesity |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, "Event_Name")) |>
  select(Score_name, level = baseline_obesity_status, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = c(Score_name, level)) |> # NB: by Score_name + level
  mutate(term = str_remove(term, "^Event_Name")) |> 
  arrange(level)
#

## 10.4 By COVID ----
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
  tofa_skin_scores_data |> 
    inner_join(tofa_participant_meta_data) |> 
    inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
    inner_join(tofa_covid_history |> filter(Event_Name == "16 week") |> select(ParticipantID, COVID_event_hx)) |> # add covid info for 16 weeks
    filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
    filter(!is.na(Score_value)) |>
    filter(is.finite(Score_value)) |>
    select(Score_name, ParticipantID, Sex, Age_Baseline, COVID_event_hx, Event_Name, Score_value) |> 
    nest(data = -c(Score_name, COVID_event_hx)) |> 
    mutate(
      test = "16 week",
      fit = map(data, possibly(\(.x) lmerTest::lmer(Score_value ~ Sex + Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
      tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
      glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
      augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
    ),
  tofa_skin_scores_data |> 
    inner_join(tofa_participant_meta_data) |> 
    inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
    inner_join(tofa_covid_history |> filter(Event_Name == "40 week") |> select(ParticipantID, COVID_event_hx)) |> # add covid info for 16 weeks
    filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
    filter(!is.na(Score_value)) |>
    filter(is.finite(Score_value)) |>
    select(Score_name, ParticipantID, Sex, Age_Baseline, COVID_event_hx, Event_Name, Score_value) |> 
    nest(data = -c(Score_name, COVID_event_hx)) |> 
    mutate(
      test = "40 week",
      fit = map(data, possibly(\(.x) lmerTest::lmer(Score_value ~ Sex + Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
      tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
      glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
      augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
    ) # 2 warning(s)
) |> 
  arrange(COVID_event_hx)
#
lm_fixedSexAge_mixedParticipantID_by_COVID |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, regex(test))) |> # keep only relevant tests
  select(Score_name, COVID_event_hx, term, AIC, BIC) # some improvements
#
lm_fixedSexAge_mixedParticipantID_by_COVID_results <- lm_fixedSexAge_mixedParticipantID_by_COVID |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, regex(test))) |> # keep only relevant tests
  filter(str_detect(term, "Event_Name")) |>
  select(Score_name, level = COVID_event_hx, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = c(Score_name, level)) |> # NB: by Score_name + level
  mutate(term = str_remove(term, "^Event_Name")) |> 
  arrange(level) |> 
  mutate(level = fct_recode(level, No = "no", Yes = "yes")) # relabel COVID levels
#


# 11 Export LM results ----
list(
  "LMM results" = list(
    "Overall" = tofa_skin_scores_lm_fixedSexAge_mixedParticipantID_results,
    "Sex" = lm_fixedAge_mixedParticipantID_by_Sex_results,
    "Age_group" = lm_fixedSex_mixedParticipantID_by_Age_group_results,
    "Obesity" = lm_fixedSexAge_mixedParticipantID_by_Obesity_results,
    "COVID" = lm_fixedSexAge_mixedParticipantID_by_COVID_results
  ) |> 
    bind_rows(.id = "stratifier") |> 
    select(
      Stratifier = stratifier,
      Level = level,
      Score_name = Score_name,
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


# 12 Forest plot(s) ----
# 12.1 IGA ----
list(
  "Overall" = tofa_skin_scores_lm_fixedSexAge_mixedParticipantID_results,
  "Sex" = lm_fixedAge_mixedParticipantID_by_Sex_results,
  "Age_group" = lm_fixedSex_mixedParticipantID_by_Age_group_results,
  "Obesity" = lm_fixedSexAge_mixedParticipantID_by_Obesity_results,
  "COVID" = lm_fixedSexAge_mixedParticipantID_by_COVID_results
) |> 
  bind_rows(.id = "stratifier") |> 
  mutate(stratifier = fct_relevel(stratifier, c("Overall", "Age_group", "Sex", "Obesity", "COVID"))) |> 
  arrange(stratifier) |> 
  filter(Score_name == "IGA") |> 
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
  scale_size_continuous(range = c(1, 3), limits = c(0, 10)) + # set limits across mutliple scores
  scale_y_discrete(limits = rev) +
  theme(
    aspect.ratio = 3,
    legend.key = element_blank()
  ) +
  labs(
    title = "IGA: Treatment vs. Baseline",
    subtitle = "LMM",
    x = "Estimate",
    y = NULL
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "Forest_plot_LMMs_IGA", ".pdf")), device = cairo_pdf, width = 8, height = 5, units = "in")
#
# 12.2 DLQI ----
list(
  "Overall" = tofa_skin_scores_lm_fixedSexAge_mixedParticipantID_results,
  "Sex" = lm_fixedAge_mixedParticipantID_by_Sex_results,
  "Age_group" = lm_fixedSex_mixedParticipantID_by_Age_group_results,
  "Obesity" = lm_fixedSexAge_mixedParticipantID_by_Obesity_results,
  "COVID" = lm_fixedSexAge_mixedParticipantID_by_COVID_results
) |> 
  bind_rows(.id = "stratifier") |> 
  mutate(stratifier = fct_relevel(stratifier, c("Overall", "Age_group", "Sex", "Obesity", "COVID"))) |> 
  arrange(stratifier) |> 
  filter(Score_name == "DLQI") |> 
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
  scale_size_continuous(range = c(1, 3), limits = c(0, 10)) + # set limits across mutliple scores
  scale_y_discrete(limits = rev) +
  theme(
    aspect.ratio = 3,
    legend.key = element_blank()
  ) +
  labs(
    title = "DLQI: Treatment vs. Baseline",
    subtitle = "LMM",
    x = "Estimate",
    y = NULL
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "Forest_plot_LMMs_DLQI", ".pdf")), device = cairo_pdf, width = 8, height = 5, units = "in")
#
# 12.3 MSS ----
list(
  "Overall" = tofa_skin_scores_lm_fixedSexAge_mixedParticipantID_results,
  "Sex" = lm_fixedAge_mixedParticipantID_by_Sex_results,
  "Age_group" = lm_fixedSex_mixedParticipantID_by_Age_group_results,
  "Obesity" = lm_fixedSexAge_mixedParticipantID_by_Obesity_results,
  "COVID" = lm_fixedSexAge_mixedParticipantID_by_COVID_results
) |> 
  bind_rows(.id = "stratifier") |> 
  mutate(stratifier = fct_relevel(stratifier, c("Overall", "Age_group", "Sex", "Obesity", "COVID"))) |> 
  arrange(stratifier) |> 
  filter(Score_name == "MSS") |> 
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
  scale_size_continuous(range = c(1, 3), limits = c(0, 10)) + # set limits across mutliple scores
  scale_y_discrete(limits = rev) +
  theme(
    aspect.ratio = 3,
    legend.key = element_blank()
  ) +
  labs(
    title = "MSS: Treatment vs. Baseline",
    subtitle = "LMM",
    x = "Estimate",
    y = NULL
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "Forest_plot_LMMs_MSS", ".pdf")), device = cairo_pdf, width = 8, height = 5, units = "in")
#
# 12.4 SALT ----
list(
  "Overall" = tofa_skin_scores_lm_fixedSexAge_mixedParticipantID_results,
  "Sex" = lm_fixedAge_mixedParticipantID_by_Sex_results,
  "Age_group" = lm_fixedSex_mixedParticipantID_by_Age_group_results,
  "Obesity" = lm_fixedSexAge_mixedParticipantID_by_Obesity_results,
  "COVID" = lm_fixedSexAge_mixedParticipantID_by_COVID_results
) |> 
  bind_rows(.id = "stratifier") |> 
  mutate(stratifier = fct_relevel(stratifier, c("Overall", "Age_group", "Sex", "Obesity", "COVID"))) |> 
  arrange(stratifier) |> 
  filter(Score_name == "SALT") |> 
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
  scale_size_continuous(range = c(1, 3), limits = c(0, 10)) + # set limits across mutliple scores
  scale_y_discrete(limits = rev) +
  theme(
    aspect.ratio = 3,
    legend.key = element_blank()
  ) +
  labs(
    title = "SALT: Treatment vs. Baseline",
    subtitle = "LMM",
    x = "Estimate",
    y = NULL
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "Forest_plot_LMMs_SALT", ".pdf")), device = cairo_pdf, width = 8, height = 5, units = "in")
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
