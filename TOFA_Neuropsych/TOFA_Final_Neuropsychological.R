################################################
# Title: Tofacitinib for Immune Skin Conditions in Down Syndrome: Analysis of Neuropsychological Assessment scores
# Author(s):
#   - Matthew Galbraith
# Affiliation(s):
#   - Linda Crnic Institute for Down Syndrome
#   - University of Colorado Anschutz
################################################

### Summary:  
# Analysis of Neuropsychological Assessment scores.
# See README.md for more details.
#  

### Data type(s):
# Clinical trial (TOFA) datasets:
#    * Participant-level metadata; Available on request.
#    * Visit/Event-level metadata; Available on request.
#    * Baseline obesity status; Available on request.
#    * COVID-19 history; Available on request.
#    * Neuropsychological Assessment scores; DOI: 10.5281/zenodo.20080323
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
tofa_participant_meta_data_file <- here("data/TOFA_Participant_metadata_zenodo_v1.txt") # Source: Available on request
tofa_visit_meta_data_file <- here("data/TOFA_Visit_metadata_zenodo_v1.txt") # Source: Available on request
tofa_baseline_obesity_file <- here("data/TOFA_Baseline_Obesity_Status_zenodo_v1.txt") # Source: Available on request
tofa_covid_history_file <- here("data/TOFA_COVID_History_zenodo_v1.txt") # Source: Available on request
tofa_neuropsychological_scores_file <- here("data/TOFA_Neuropsychological_Endpoint_scores_data_zenodo_v1.txt.gz") # Source: 10.5281/zenodo.20080323
#
# Other parameters:
standard_colors <- c("Baseline" = "#999999", "2 week" = "#c6dbef", "8 week" = "#9ecae1", "16 week" = "#6baed6", "40 week" = "#4292c6")
#
out_file_prefix <- "TOFA_Final_neuropsychological_" # should match this script title
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

## 1.2 Read in TOFA Neuropsych data ----
tofa_neuropsychological_scores <- tofa_neuropsychological_scores_file |>
  read_tsv() |> 
  mutate(
    Event_Name = fct_relevel(Event_Name, c("Baseline", "16 week","40 week")) # set factor levels
  )
#
tofa_neuropsychological_scores # 2,185 rows
tofa_neuropsychological_scores |> distinct(ParticipantID) # 42 Participants
tofa_neuropsychological_scores |> distinct(VisitID) # 115 VisitIDs
tofa_neuropsychological_scores |> distinct(Score_name) # 19 scores
#

## 1.3 Join with meta data, filter eligible etc, filter to endpoint cytokines ----
tofa_neuropsychological_scores <- tofa_neuropsychological_scores |> 
  inner_join(tofa_visit_meta_data) |> # returns 2,185 of 2,185 rows
  inner_join(tofa_participant_meta_data) |> # returns 2,185 of 2,185 rows
  filter(Endpoint_eligible == TRUE) # returns 2,185 of 2,185 rows
#
### 1.3.1 ParticipantID vs Event_Name summary ---- 
tofa_neuropsychological_scores |> 
  distinct(ParticipantID, VisitID, Event_Name) |> 
  arrange(Event_Name) |> 
  pivot_wider(names_from = Event_Name, values_from = VisitID) |> 
  print(n = Inf)
# Note: Some Participants have missing samples and not all completed extension to 40 weeks.


# 2 Endpoint Stats ----
## 2.1 Check data distributions ----
### 2.1.1 check for outlier values ----
tofa_neuropsychological_scores |> 
  mutate(
    extreme = rstatix::is_extreme(Score_value),
    outlier = rstatix::is_outlier(Score_value),
    .by = Score_name
  ) |> 
  ggplot(aes(Score_name, Score_value, color = Score_name)) +
  geom_sina(data = . %>% filter(outlier == FALSE & extreme == FALSE), maxwidth = 0.5, color = "grey") +
  geom_sina(data = . %>% filter(outlier == TRUE & extreme == FALSE), aes(color = "outlier"), maxwidth = 0.1) +
  geom_sina(data = . %>% filter(extreme == TRUE), aes(color = "extreme"), maxwidth = 0.1) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~ Score_name, scales = "free") +
  coord_cartesian(clip = 'off') + # prevent label clipping
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    aspect.ratio = 1.3
  ) +
  labs(title = "Distributions of Score Values + outliers")
#
### 2.1.2 Calc Differences from baseline ----
tofa_neuropsychological_16week_diffs <- tofa_neuropsychological_scores |> 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(Score_value)) |> 
  filter(is.finite(Score_value)) |> 
  add_count(Score_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, Score_name) |> 
  group_by(Score_name) |> 
  filter(!any(n < 2)) |> # remove scores with n < X in any group
  select(ParticipantID, Event_Name, Assessment, Score_name, Score_value) |> 
  pivot_wider(names_from = Event_Name, values_from = Score_value) |> 
  mutate(difference = `16 week` - Baseline) |> 
  select(Assessment, ParticipantID, Score_name, difference) |> 
  mutate(
    extreme = rstatix::is_extreme(difference),
    outlier = rstatix::is_outlier(difference)
  ) |> 
  ungroup()
#
tofa_neuropsychological_40week_diffs <- tofa_neuropsychological_scores |> 
  filter(Event_Name %in% c("Baseline", "40 week")) |> 
  filter(!is.na(Score_value)) |> 
  filter(is.finite(Score_value)) |> 
  add_count(Score_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, Score_name) |> 
  group_by(Score_name) |> 
  filter(!any(n < 2)) |> # remove scores with n < X in any group
  select(ParticipantID, Event_Name, Assessment, Score_name, Score_value) |> 
  pivot_wider(names_from = Event_Name, values_from = Score_value) |> 
  mutate(difference = `40 week` - Baseline) |> 
  select(Assessment, ParticipantID, Score_name, difference) |> 
  mutate(
    extreme = rstatix::is_extreme(difference),
    outlier = rstatix::is_outlier(difference)
  ) |> 
  ungroup()
#
### 2.1.3 Count extreme outliers in differences ----
tofa_neuropsychological_16week_diffs |> 
  count(outlier, extreme)
tofa_neuropsychological_40week_diffs |> 
  count(outlier, extreme)
#
### 2.1.4 Plot differences ---- 
tofa_neuropsychological_16week_diffs |> 
  ggplot(aes(Score_name, difference, color = Score_name)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_sina(data = . %>% filter(outlier == FALSE & extreme == FALSE), maxwidth = 0.5, color = "grey") +
  geom_sina(data = . %>% filter(outlier == TRUE & extreme == FALSE), aes(color = "outlier"), maxwidth = 0.1) +
  geom_sina(data = . %>% filter(extreme == TRUE), aes(color = "extreme"), maxwidth = 0.1) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~ Score_name, scales = "free") +
  coord_cartesian(clip = 'off') + # prevent label clipping
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    aspect.ratio = 1.3,
  ) +
  labs(title = "Distributions of 16 week vs Baseline differences")
#
tofa_neuropsychological_40week_diffs |> 
  ggplot(aes(Score_name, difference, color = Score_name)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_sina(data = . %>% filter(outlier == FALSE & extreme == FALSE), maxwidth = 0.5, color = "grey") +
  geom_sina(data = . %>% filter(outlier == TRUE & extreme == FALSE), aes(color = "outlier"), maxwidth = 0.1) +
  geom_sina(data = . %>% filter(extreme == TRUE), aes(color = "extreme"), maxwidth = 0.1) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~ Score_name, scales = "free") +
  coord_cartesian(clip = 'off') + # prevent label clipping
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    aspect.ratio = 1.3,
  ) +
  labs(title = "Distributions of 40 week vs Baseline differences")
#
### 2.1.5 Check if differences are normally distributed ----
tofa_neuropsychological_16week_diffs |>
  ggpubr::ggqqplot("difference", title = "Q-Q plot: 16 week differences", facet.by = "Score_name",  scales = "free_y")
tofa_neuropsychological_40week_diffs |>
  ggpubr::ggqqplot("difference", title = "Q-Q plot: 40 week differences", facet.by = "Score_name",  scales = "free_y")
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

### 2.2.1 Run paired t tests 16 weeks ----
# Only 16 week timepoint is tested as a trial endpoint.
tofa_neuropsychological_Ttest_res_16weeks <- tofa_neuropsychological_scores |> 
  arrange(ParticipantID, Event_Name) |> # ensure correct sort order for rstatix::t_test()
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(Score_value)) |> 
  filter(is.finite(Score_value)) |> 
  add_count(Score_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, Score_name) |> 
  group_by(Assessment, Score_name) |> 
  filter(!any(n < 2)) |> # remove scores with n < X in any group
  rstatix::t_test(
    formula = Score_value ~ Event_Name,
    ref.group = "16 week",
    paired = TRUE,
    var.equal = TRUE, # set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  ) |> 
  mutate(
    BH_padj = p.adjust(p, method = "BH"),
    .by = Assessment # adjusts within each Assessment group
    ) |> 
  arrange(p)
#
### 2.2.2 Paired t test effect size ----
# The effect size for a paired-samples t-test can be calculated by dividing the
# mean difference by the standard deviation of the difference, as shown below.
# Cohen’s formula:
# d = mean(D)/sd(D), where D is the differences of the paired samples values.
tofa_neuropsychological_Ttest_effsize_16weeks <- tofa_neuropsychological_scores |> 
  arrange(ParticipantID, Event_Name) |> # ensure correct sort order for rstatix::cohens_d()
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(Score_value)) |> 
  filter(is.finite(Score_value)) |> 
  add_count(Score_name, ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  add_count(Event_Name, Score_name) |> 
  group_by(Assessment, Score_name) |> 
  filter(!any(n < 2)) |> # remove scores with n < X in any group
  rstatix::cohens_d(
    formula = Score_value ~ Event_Name,
    ref.group = "16 week",
    paired = TRUE,
    var.equal = TRUE # set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
  )
#
### 2.2.3 Compile t test results ----
tofa_neuropsychological_Ttest_res_16weeks_full <- tofa_neuropsychological_Ttest_res_16weeks |> 
  inner_join(tofa_neuropsychological_Ttest_effsize_16weeks) |> 
  mutate(
    mean_diff = estimate,
  ) |> 
  inner_join(tofa_neuropsychological_16week_diffs |> summarize(median_diff = median(difference), .by = Score_name)) |> 
  select(Assessment, Score_name, mean_diff, median_diff, p, BH_padj, effsize, magnitude, group1, group2, n1, n2, everything()) |> 
  mutate(Assessment = fct_relevel(Assessment, c("leiter", "cantab", "kbit", "ppvt", "sobc", "nepsyii", "promis"))) |> 
  arrange(Assessment, Score_name)
#
tofa_neuropsychological_Ttest_res_16weeks_full
#
### 2.2.4 Export results ----
list(
  "ttest_results" = tofa_neuropsychological_Ttest_res_16weeks_full |> 
    select(
      Assessment,
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
  export_excel(filename = "Ttest_results")
#


# 3 Paired plots - score values --------
#
## 3.1 Baseline vs. 16 weeks - All ----
tofa_neuropsychological_scores |> 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  #
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
  patchwork::plot_layout(nrow = 2, guides = "collect") +
  patchwork::plot_annotation(
    title = "Baseline vs. Treatment"
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "paired_plots_B_16_weeks", ".pdf")), device = cairo_pdf, width = 35, limitsize = FALSE, height = 6, units = "in")
#
## 3.2 Baseline vs. 16 weeks - Leiter ----
tofa_neuropsychological_scores |> 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(Assessment == "leiter") |> 
  #
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
  patchwork::plot_layout(
    # guides = "collect", 
    nrow = 1
  ) +
  patchwork::plot_annotation(
    title = "Leiter: Baseline vs. Treatment",
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "paired_plots_Leiter_B_16_weeks", ".pdf")), device = cairo_pdf, width = 20, limitsize = FALSE, height = 6, units = "in")
#


# 4 Sina plots - score values --------
## 4.1 Baseline vs. 16 weeks - Leiter ----
tofa_neuropsychological_scores |> 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(Assessment == "leiter") |> 
  #
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
  patchwork::plot_layout(
    # guides = "collect", 
    nrow = 1
  ) +
  patchwork::plot_annotation(
    title = "Leiter: Baseline vs. Treatment"
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_plots_Leiter_B_16_weeks", ".pdf")), device = cairo_pdf, width = 20, limitsize = FALSE, height = 6, units = "in")
#


# 5 Sina plots - Differences --------
## 5.1 Baseline vs. 16 weeks - ALL  ----
bind_rows(
  tofa_neuropsychological_16week_diffs |> mutate(Event_Name = "16 week"),
  tofa_neuropsychological_40week_diffs |> mutate(Event_Name = "40 week")
) |> 
  mutate(Event_Name = fct_relevel(Event_Name, c("16 week", "40 week"))) |> 
  #
  group_by(Score_name) |> 
  group_split() %>%
  # .[1] %>% # for testing
  map( # using anonymous function: \(x)
    \(x) { x2 <- x |> 
      arrange(Score_name, Event_Name) |> 
      mutate( # add labels with n
        n = n_distinct(ParticipantID), 
        .by = c(Score_name, Event_Name),
        label = paste0(str_extract(Event_Name, "^[[:alpha:]]|^\\d{1,2}"), "\n(n=", n, ")"),
        label = fct_inorder(label)
      ) |> 
      mutate(Event_Name = fct_drop(Event_Name)) |> # prevent empty x values
      complete(Score_name, ParticipantID, Event_Name) # Prevent error from stat_signif()
      #
    x2 %>%
      {ggplot(data = ., aes(Event_Name, difference, color = Event_Name)) +
          geom_hline(yintercept = 0, linetype = 2) +
          geom_line(aes(group = ParticipantID), color = "grey90") +
          geom_sina(maxwidth = 0.3) +
          geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
          scale_color_manual(values = standard_colors) +
          geom_text_repel(
            data = x2 |> mutate(TOFA_label =  if_else(ParticipantID %in% c("TOFA0001", "TOFA0011"), ParticipantID, "")),
            aes(label = TOFA_label),
            size = 2,
            nudge_x = 0.25,
            min.segment.length = 0,
            color = "black"
          ) +
          facet_wrap(~ Score_name, scales = "free_y", nrow = 1) +
          coord_cartesian(clip = 'off') + # prevent label clipping
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
  }) |> 
  patchwork::wrap_plots() +
  patchwork::plot_layout(
    # guides = "collect", 
    nrow = 2
  ) +
  patchwork::plot_annotation(
    title = "Endpoints scores DIFFERENCES: Baseline vs. Treatment"
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_plots_DIFFS_ENDPOINTS_16_40_weeks_DSRD_labels", ".pdf")), device = cairo_pdf, width = 70, limitsize = FALSE, height = 6, units = "in")
#
## 5.2 Baseline vs. 16 weeks - leiter ----
bind_rows(
  tofa_neuropsychological_16week_diffs |> mutate(Event_Name = "16 week"),
  tofa_neuropsychological_40week_diffs |> mutate(Event_Name = "40 week")
) |> 
  mutate(Event_Name = fct_relevel(Event_Name, c("16 week", "40 week"))) |> 
  filter(Assessment == "leiter") |> 
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
      complete(Score_name, ParticipantID, Event_Name) %>% # Prevent error from stat_signif()
      #
      {ggplot(data = ., aes(Event_Name, difference, color = Event_Name)) +
          geom_hline(yintercept = 0, linetype = 2) +
          geom_line(aes(group = ParticipantID), color = "grey90") +
          geom_sina(maxwidth = 0.3) +
          geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
          scale_color_manual(values = standard_colors) +
          facet_wrap(~ Score_name, scales = "free_y", nrow = 1) +
          coord_cartesian(clip = 'off') + # prevent label clipping
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
  patchwork::plot_layout(
    # guides = "collect", 
    nrow = 1
  ) +
  patchwork::plot_annotation(
    title = "leiter scores DIFFERENCES: Baseline vs. Treatment",
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_plots_leiter_DIFFS_ENDPOINTS_16_40_weeks", ".pdf")), device = cairo_pdf, width = 20, limitsize = FALSE, height = 6, units = "in")
#


# 6 Dumbbell plots ----
### 6.1. 16 weeks ----
tofa_neuropsychological_scores |> 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  group_by(Score_name) |> 
  group_split() %>%
  # .[1] |> # for testing
  map( # using anonymous function: \(x)
    \(x) {x2 <- x |> 
      filter(!is.na(Score_value)) |>
      filter(is.finite(Score_value)) |>
      add_count(Score_name, ParticipantID, name = "pair") |>
      filter(pair == 2) |> # require pairs
      filter(Event_Name == "Baseline") |>
      select(ParticipantID, Score_name, Baseline = Score_value) |>
      inner_join(
        x |>
          filter(Event_Name %in% c("Baseline", "16 week")) |>
          filter(!is.na(Score_value)) |>
          filter(is.finite(Score_value)) |>
          add_count(Score_name, ParticipantID, name = "pair") |>
          filter(pair == 2) |> # require pairs
          filter(Event_Name == "16 week") |>
          select(ParticipantID, Score_name, Week16 = Score_value)
      ) |>
      mutate(
        diff = Week16 - Baseline,
        direction = case_when(
          diff > 0 ~ "Increased",
          diff < 0 ~ "Decreased",
          diff == 0 ~ "Unchanged"
        ),
        direction = fct_relevel(direction, c("Decreased", "Unchanged", "Increased"))
      ) |>
      arrange(direction, Baseline, diff) |>
      mutate(ParticipantID = fct_inorder(ParticipantID))
      x2 |> ggplot() +
      ggarchery::geom_arrowsegment(
        aes(y = ParticipantID, yend = ParticipantID, x = Baseline, xend = Week16, col = direction, linewidth = abs(diff), alpha = abs(diff)),
        arrows = arrow(angle = 20, length = unit(0.03, "npc"), type = "open"),
        arrow_positions = 0.6
      ) +
      scale_color_manual(values = c("Decreased" = "#0571b0", "Unchanged" = "#cccccc", "Increased" = "#ca0020")) +
      scale_linewidth(range = c(0.25, 0.75)) +
      scale_alpha(range = c(0.2, 0.75)) +
      ggnewscale::new_scale_color() +
      geom_point(aes(Baseline, ParticipantID, color = "Baseline")) +
      scale_color_manual(values = standard_colors) +
      geom_point(aes(Week16, ParticipantID, color = "16 week")) +
        geom_text(
          data = x2 |> mutate(label =  if_else(ParticipantID %in% c("TOFA0001", "TOFA0011"), ParticipantID, "")),
          aes(Week16, ParticipantID, label = label),,
          size = 2,
          nudge_x = 2
        ) +
      facet_wrap(~ Score_name, scales = "free_x", nrow = 1) +
      coord_cartesian(clip = 'off') + # prevent label clipping
      theme(
        aspect.ratio = 1.5,
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      ) +
      labs(
        x = "Units",
        y = NULL
      )
      }) |> 
  patchwork::wrap_plots() +
  patchwork::plot_layout(
    guides = "collect",
    nrow = 2
  ) +
  patchwork::plot_annotation(
    title = "Change in score values: 16 weeks vs Baseline",
    subtitle = "Endpoint measures; arranged by: direction, Baseline, diff",
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "dumbbell_16_weeks_ENDPOINTS_DSRD_labels", ".pdf")), device = cairo_pdf, width = 35, limitsize = FALSE, height = 10, units = "in")
#




# 9 Differences Heatmap(s) with mean improvement score ----
#
### 9.1 Endpoints 16 weeks vs Baseline ----
hm_dat_diffs <- bind_rows(
  tofa_neuropsychological_16week_diffs |> mutate(Event_Name = "16 week"),
  tofa_neuropsychological_40week_diffs |> mutate(Event_Name = "40 week")
) |> 
  inner_join(tofa_participant_meta_data |> select(ParticipantID, Age_consent_years, Sex)) |> 
  mutate(Event_Name = fct_relevel(Event_Name, c("16 week", "40 week"))) |> 
  mutate(scaled_diff = difference / sd(difference), .by = c(Score_name)) |> 
  mutate(ParticipantID = paste0(ParticipantID, "_", Event_Name)) |> 
  mutate(
    improvement = case_when(
      Score_name == "ne_total_err" ~  scaled_diff * -1, # lower is better
      Score_name == "ne_total_time" ~  scaled_diff * -1, # lower is better
      Score_name == "le_as_raw" ~  scaled_diff, # higher is better
      Score_name == "le_rm_gsv" ~  scaled_diff, # higher is better
      Score_name == "le_so_gsv" ~  scaled_diff, # higher is better
      Score_name == "le_ns_effect_raw" ~  scaled_diff, # higher is better
      Score_name == "le_fm_gsv" ~  scaled_diff, # higher is better
      Score_name == "pr_pos_affect_tscore" ~  scaled_diff, # higher is better
      Score_name == "pr_anxiety_tscore" ~  scaled_diff * -1, # lower is better
      Score_name == "pr_gen_concern_tscore" ~  scaled_diff * -1, # lower is better
      Score_name == "pr_depression_tscore" ~  scaled_diff * -1, # lower is better
      Score_name == "ca_sspfsl" ~  scaled_diff, # higher is better
      Score_name == "ca_palnpr212" ~  scaled_diff, # higher is better
      Score_name == "ca_rtismrt" ~  scaled_diff * -1, # lower is better
      Score_name == "so_median_rt" ~  scaled_diff * -1, # lower is better
      Score_name == "kbit_iq_composite_standard" ~  scaled_diff, # higher is better
      Score_name == "kbit2_raw_verbal" ~  scaled_diff, # higher is better
      Score_name == "pp_raw" ~  scaled_diff, # higher is better
      .default = scaled_diff
    )
  ) |> 
  mutate(mean_improvement = mean(improvement), .by = c(ParticipantID)) |> 
  arrange(Event_Name, mean_improvement) |> 
  mutate(ParticipantID = fct_inorder(ParticipantID))
#
hm_diffs_lim <- hm_dat_diffs |>
  pull(scaled_diff) |>
  abs() |>
  max() |>
  round(2)
breaks_diffs <- seq(-hm_diffs_lim, hm_diffs_lim, length.out = 11)
hm_diffs_palette <- circlize::colorRamp2(
  breaks_diffs,
  RColorBrewer::brewer.pal(11, "RdBu") |> rev()
)
#
mean_improvement_lim <- hm_dat_diffs |> 
  pull(mean_improvement) |> 
  abs() |>
  max() |>
  round(2)
mean_improvement_palette <- circlize::colorRamp2(
  seq(-mean_improvement_lim, mean_improvement_lim, length.out = 11),
  RColorBrewer::brewer.pal(11, "PRGn") |> rev()
)
#
age_max <- hm_dat_diffs |> 
  pull(Age_consent_years) |> 
  max() |>
  round(2)
age_min <- hm_dat_diffs |> 
  pull(Age_consent_years) |> 
  min() |>
  round(2)
age_palette <- circlize::colorRamp2(
  seq(age_min, age_max, length.out = 9),
  RColorBrewer::brewer.pal(9, "Blues")
)
#
# plot heatmap
hm_dat_diffs |>
  group_by(Event_Name, Assessment) |> 
    tidyHeatmap::heatmap(
      Score_name,
      ParticipantID,
      scaled_diff,
      palette_value = hm_diffs_palette,
      heatmap_legend_param = list(color_bar = "continuous", at = seq(-hm_diffs_lim, hm_diffs_lim, length.out = 5)),
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      row_title = NULL,
      show_column_names = TRUE,
      column_title = NULL,
      border = TRUE
    ) |>
  tidyHeatmap::annotation_tile(mean_improvement, palette = mean_improvement_palette, border = TRUE) |> 
  tidyHeatmap::annotation_tile(Age_consent_years, palette = age_palette, border = TRUE) |> 
  tidyHeatmap::wrap_heatmap()
ggsave(filename = here("plots", paste0(out_file_prefix, "heatmap_improvements_ENDPOINTS_PRGn_withAge", ".pdf")), device = cairo_pdf, width = 10, height = 5, units = "in")
#


# 10 Baseline scores vs Differences correlations ----
## 10.1 Spearman correlations -----
tofa_neuropsychological_endpoints_spearman <- tofa_neuropsychological_scores |> 
  filter(Event_Name == "Baseline") %>% 
  rename(Baseline = Score_value) %>% 
  select(ParticipantID, Score_name, Baseline) |> 
  inner_join(
    tofa_neuropsychological_16week_diffs |> select(ParticipantID, Assessment, Score_name, diff_16weeks = difference)
  ) |> 
  left_join(
    tofa_neuropsychological_40week_diffs |> select(ParticipantID, Assessment, Score_name, diff_40weeks = difference)
  ) |> 
  nest(data = -c(Assessment, Score_name)) |> 
  mutate(
    cor_16weeks = map(data, \(x) rstatix::cor_test(x, Baseline, diff_16weeks, method = "spearman")),
    cor_40weeks = map(data, \(x) rstatix::cor_test(x, Baseline, diff_40weeks, method = "spearman"))
  )
#
tofa_neuropsychological_endpoints_spearman |> 
  unnest(cor_16weeks) |> 
  mutate(BH_padj = p.adjust(p, method = "BH")) |>
  arrange(p) |> 
  filter(BH_padj < 0.1)
tofa_neuropsychological_endpoints_spearman |> 
  unnest(cor_40weeks) |> 
  mutate(BH_padj = p.adjust(p, method = "BH")) |>
  arrange(p) |> 
  filter(BH_padj < 0.1)
#

## 10.2 Heatmap -----
cor_hm_dat <- bind_rows(
  tofa_neuropsychological_endpoints_spearman |> 
    unnest(cor_16weeks) |> 
    select(-c(data, cor_40weeks)) |> 
    mutate(Event_Name = "16 weeks") |> 
    mutate(BH_padj = p.adjust(p, method = "BH")),
  tofa_neuropsychological_endpoints_spearman |> 
    unnest(cor_40weeks) |> 
    select(-c(data, cor_16weeks)) |> 
    mutate(Event_Name = "40 weeks") |> 
    mutate(BH_padj = p.adjust(p, method = "BH"))
) |> 
  rename(rho = cor)
#
cor_hm_lim <- cor_hm_dat |>
  pull(rho) |>
  abs() |>
  max() |>
  round(2)
cor_hm_breaks <- seq(-cor_hm_lim, cor_hm_lim, length.out = 11)
cor_hm_palette <- circlize::colorRamp2(
  cor_hm_breaks,
  RColorBrewer::brewer.pal(11, "RdBu") |> rev()
)
# plot heatmap
cor_hm_dat |>
  arrange(rho) |> 
  # mutate(Score_name = fct_inorder(Score_name)) |> 
  group_by(Assessment) |> 
  tidyHeatmap::heatmap(
    Score_name,
    Event_Name,
    rho,
    palette_value = cor_hm_palette,
    heatmap_legend_param = list(color_bar = "continuous", at = seq(-cor_hm_lim, cor_hm_lim, length.out = 5)),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    row_title = NULL,
    show_column_names = TRUE,
    column_title = NULL,
    border = TRUE
  ) |> 
  tidyHeatmap::layer_asterisk(BH_padj < 0.1) |>
  tidyHeatmap::wrap_heatmap() +
labs(
  title = "neuropsychological Endpoints: Correlation between Baseline and Differences",
)
ggsave(filename = here("plots", paste0(out_file_prefix, "heatmap_spearman_Baseline_vs_Diffs_", ".pdf")), device = cairo_pdf, width = 4.5, height = 5, units = "in")
#

## 10.3 Scatter plots with rho + p -----
tofa_neuropsychological_scores |> 
  filter(Event_Name == "Baseline") %>% 
  rename(Baseline = Score_value) %>% 
  select(ParticipantID, Score_name, Baseline) |> 
  inner_join(
    tofa_neuropsychological_16week_diffs |> select(ParticipantID, Assessment, Score_name, diff_16weeks = difference)
  ) |> 
  left_join(
    tofa_neuropsychological_40week_diffs |> select(ParticipantID, Assessment, Score_name, diff_40weeks = difference)
  ) |> 
  ggplot(aes(Baseline, diff_16weeks)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Score_name, scales = "free", nrow = 2) +
  coord_cartesian(clip = 'off') + # prevent label clipping
  theme(aspect.ratio = 1) +
  ggpubr::stat_cor(method = "spearman", cor.coef.name = "rho", alternative = "two.sided", p.accuracy = 0.05)
ggsave(filename = here("plots", paste0(out_file_prefix, "scatter_16_weeks_Endpoint_diffs", ".pdf")), device = cairo_pdf, width = 35, limitsize = FALSE, height = 10, units = "in")
#


# 11 Mixed effects linear regression models (non-stratified) ----
# Time as a categorical variable
#   No assumptions about trajectory
#   Directly compares each visit to baseline
# Random effects
#   Handles missing at random (MAR) data.
#   Accounts for within-subject correlation.
## 11.1 Set up models ----
### 11.1.1 Mixed effects LM: Event_name + 1|ParticipantID ----
tofa_neuropsychological_lm_mixedParticipantID <- tofa_neuropsychological_scores |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Score_value)) |>
  filter(is.finite(Score_value)) |>
  select(Assessment, Score_name, ParticipantID, Sex, Age_Baseline, Event_Name, Score_value) |> 
  nest(data = -c(Assessment, Score_name)) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Score_value ~ Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
### 11.1.2 Mixed effects LM: Sex + Event_name + 1|ParticipantID ----
tofa_neuropsychological_lm_fixedSex_mixedParticipantID <- tofa_neuropsychological_scores |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Score_value)) |>
  filter(is.finite(Score_value)) |>
  select(Assessment, Score_name, ParticipantID, Sex, Age_Baseline, Event_Name, Score_value) |> 
  nest(data = -c(Assessment, Score_name)) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Score_value ~ Sex + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  ) 
#
### 11.1.3 Mixed effects LM: Age + Event_name + 1|ParticipantID ----
tofa_neuropsychological_lm_fixedAge_mixedParticipantID <- tofa_neuropsychological_scores |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Score_value)) |>
  filter(is.finite(Score_value)) |>
  select(Assessment, Score_name, ParticipantID, Sex, Age_Baseline, Event_Name, Score_value) |> 
  nest(data = -c(Assessment, Score_name)) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Score_value ~ Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  ) 
#
### 11.1.4 Mixed effects LM: Sex + Age + Event_name + 1|ParticipantID ----
tofa_neuropsychological_lm_fixedSexAge_mixedParticipantID <- tofa_neuropsychological_scores |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Score_value)) |>
  filter(is.finite(Score_value)) |>
  select(Assessment, Score_name, ParticipantID, Sex, Age_Baseline, Event_Name, Score_value) |> 
  nest(data = -c(Assessment, Score_name)) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Score_value ~ Sex + Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
## 11.2 Compare models ----
### 11.2.1 AIC/BIC ----
tofa_neuropsychological_lm_mixedParticipantID |> unnest(glanced) |> select(Score_name, AIC, BIC)
tofa_neuropsychological_lm_fixedSex_mixedParticipantID |> unnest(glanced) |> select(Score_name, AIC, BIC)
tofa_neuropsychological_lm_fixedAge_mixedParticipantID |> unnest(glanced) |> select(Score_name, AIC, BIC)
tofa_neuropsychological_lm_fixedSexAge_mixedParticipantID |> unnest(glanced) |> select(Score_name, AIC, BIC)
#
### 11.2.2 Likelihood ratio tests ----
anova(
  tofa_neuropsychological_lm_mixedParticipantID |> pluck("fit", 1),
  tofa_neuropsychological_lm_fixedSex_mixedParticipantID |> pluck("fit", 1),
  test="LRT"
) |>  # refitting model(s) with ML (instead of REML)
  tidy() # p.value = 0.959
anova(
  tofa_neuropsychological_lm_mixedParticipantID |> pluck("fit", 1),
  tofa_neuropsychological_lm_fixedAge_mixedParticipantID |> pluck("fit", 1),
  test="LRT"
) |>  # refitting model(s) with ML (instead of REML)
  tidy() # p.value = 0.0184
anova(
  tofa_neuropsychological_lm_mixedParticipantID |> pluck("fit", 1),
  tofa_neuropsychological_lm_fixedSexAge_mixedParticipantID |> pluck("fit", 1),
  test="LRT"
) |>  # refitting model(s) with ML (instead of REML)
  tidy() # p.value = 0.0621
#

## 11.3 Model results ----
### 11.3.1 Mixed ParticipantID with fixed SexAge -----
tofa_neuropsychological_lm_fixedSexAge_mixedParticipantID |> 
  unnest(tidied) |> 
  select(Score_name, group, term, estimate, p.value)
# 
tofa_neuropsychological_lm_fixedSexAge_mixedParticipantID_results <- tofa_neuropsychological_lm_fixedSexAge_mixedParticipantID |> 
  unnest(c(tidied, glanced)) |>
  filter(str_detect(term, "Event_Name")) |>
  select(Assessment, Score_name, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = Score_name) |> # NB: by Score_name
  mutate(
    level = "",
    term = str_remove(term, "^Event_Name"),
    .after = Score_name,
  )
#
tofa_neuropsychological_lm_fixedSexAge_mixedParticipantID_results 
#


# 12 Mixed effects linear regression models (stratified) ----
## 12.1 By Sex ----
lm_fixedAge_mixedParticipantID_by_Sex <- tofa_neuropsychological_scores |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Score_value)) |>
  filter(is.finite(Score_value)) |>
  select(Assessment, Score_name, ParticipantID, Sex, Age_Baseline, Event_Name, Score_value) |> 
  nest(data = -c(Assessment, Score_name, Sex)) |> 
  mutate(
    fit = map(data, possibly(\(.x) lmerTest::lmer(Score_value ~ Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
lm_fixedAge_mixedParticipantID_by_Sex |> unnest(glanced) |> select(Score_name, Sex, AIC, BIC)
#
lm_fixedAge_mixedParticipantID_by_Sex_results <- lm_fixedAge_mixedParticipantID_by_Sex |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, "Event_Name")) |>
  select(Assessment, Score_name, level = Sex, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = c(Score_name, level)) |> # NB: by Score_name + level
  mutate(term = str_remove(term, "^Event_Name")) |> 
  arrange(level)
#

## 12.2 By Age group ----
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
lm_fixedSex_mixedParticipantID_by_Age_group <- tofa_neuropsychological_scores |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  inner_join(tofa_baseline_age_groups) |> # add age groups info
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Score_value)) |>
  filter(is.finite(Score_value)) |>
  select(Assessment, Score_name, ParticipantID, Sex, Age_group, Event_Name, Score_value) |> 
  nest(data = -c(Assessment, Score_name, Age_group)) |> 
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
  select(Assessment, Score_name, level = Age_group, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = c(Score_name, level)) |> # NB: by Score_name + level
  mutate(term = str_remove(term, "^Event_Name")) |> 
  arrange(level)
#

## 12.3 By Obesity ----
tofa_baseline_obesity <- tofa_baseline_obesity_file |> 
  read_tsv()
tofa_baseline_obesity |> count(baseline_obesity_status)
#
lm_fixedSexAge_mixedParticipantID_by_Obesity <- tofa_neuropsychological_scores |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  inner_join(tofa_baseline_obesity) |> # add Obesity info
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
  filter(!is.na(Score_value)) |>
  filter(is.finite(Score_value)) |>
  select(Assessment, Score_name, ParticipantID, Sex, Age_Baseline, baseline_obesity_status, Event_Name, Score_value) |> 
  nest(data = -c(Assessment, Score_name, baseline_obesity_status)) |> 
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
  select(Assessment, Score_name, level = baseline_obesity_status, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = c(Score_name, level)) |> # NB: by Score_name + level
  mutate(term = str_remove(term, "^Event_Name")) |> 
  arrange(level)
#

## 12.4 By COVID ----
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
  tofa_neuropsychological_scores |> 
    inner_join(tofa_participant_meta_data) |> 
    inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
    inner_join(tofa_covid_history |> filter(Event_Name == "16 week") |> select(ParticipantID, COVID_event_hx)) |> # add covid info for 16 weeks
    filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
    filter(!is.na(Score_value)) |>
    filter(is.finite(Score_value)) |>
    select(Assessment, Score_name, ParticipantID, Sex, Age_Baseline, COVID_event_hx, Event_Name, Score_value) |> 
    nest(data = -c(Assessment, Score_name, COVID_event_hx)) |> 
    mutate(
      test = "16 week",
      fit = map(data, possibly(\(.x) lmerTest::lmer(Score_value ~ Sex + Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
      tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
      glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
      augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
    ),
  tofa_neuropsychological_scores |> 
    inner_join(tofa_participant_meta_data) |> 
    inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
    inner_join(tofa_covid_history |> filter(Event_Name == "40 week") |> select(ParticipantID, COVID_event_hx)) |> # add covid info for 40 weeks
    filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> # keep only timepoints of interest
    filter(!is.na(Score_value)) |>
    filter(is.finite(Score_value)) |>
    select(Assessment, Score_name, ParticipantID, Sex, Age_Baseline, COVID_event_hx, Event_Name, Score_value) |> 
    nest(data = -c(Assessment, Score_name, COVID_event_hx)) |> 
    mutate(
      test = "40 week",
      fit = map(data, possibly(\(.x) lmerTest::lmer(Score_value ~ Sex + Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE))),
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
  select(Score_name, COVID_event_hx, term, AIC, BIC)
#
lm_fixedSexAge_mixedParticipantID_by_COVID_results <- lm_fixedSexAge_mixedParticipantID_by_COVID |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, regex(test))) |> # keep only relevant tests
  filter(str_detect(term, "Event_Name")) |>
  select(Assessment, Score_name, level = COVID_event_hx, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = c(Score_name, level)) |> # NB: by Score_name + level
  mutate(term = str_remove(term, "^Event_Name")) |> 
  arrange(level) |> 
  mutate(level = fct_recode(level, No = "no", Yes = "yes")) # relabel COVID levels
#


# 13 Export LM results ----
list(
  "LMM results" = list(
    "Overall" = tofa_neuropsychological_lm_fixedSexAge_mixedParticipantID_results,
    "Sex" = lm_fixedAge_mixedParticipantID_by_Sex_results,
    "Age_group" = lm_fixedSex_mixedParticipantID_by_Age_group_results,
    "Obesity" = lm_fixedSexAge_mixedParticipantID_by_Obesity_results,
    "COVID" = lm_fixedSexAge_mixedParticipantID_by_COVID_results
  ) |> 
    bind_rows(.id = "stratifier") |> 
    mutate(Assessment = fct_relevel(Assessment, c("leiter", "cantab", "kbit", "ppvt", "sobc", "nepsyii", "promis"))) |> 
    arrange(Assessment) |> 
    select(
      Stratifier = stratifier,
      Level = level,
      Assessment,
      Score_name,
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
  export_excel(filename = "LMM_combined_results_ALL")
#

# 14 Forest plot(s) ----
## 14.1 kb_iq_composite_standard ----
list(
  "Overall" = tofa_neuropsychological_lm_fixedSexAge_mixedParticipantID_results,
  "Sex" = lm_fixedAge_mixedParticipantID_by_Sex_results,
  "Age_group" = lm_fixedSex_mixedParticipantID_by_Age_group_results,
  "Obesity" = lm_fixedSexAge_mixedParticipantID_by_Obesity_results,
  "COVID" = lm_fixedSexAge_mixedParticipantID_by_COVID_results
) |>
  bind_rows(.id = "stratifier") |>
  mutate(stratifier = fct_relevel(stratifier, c("Overall", "Age_group", "Sex", "Obesity", "COVID"))) |>
  arrange(stratifier) |>
  filter(Score_name == "kb_iq_composite_standard") |>
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
  scale_size_continuous(range = c(2, 4)) +
  scale_y_discrete(limits = rev) +
  theme(
    aspect.ratio = 1.5,
    legend.key = element_blank()
  ) +
  labs(
    title = "kb_iq_composite_standard: Treatment vs. Baseline",
    subtitle = "LMM",
    x = "Estimate",
    y = NULL
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "Forest_plot_LMMs_kb_iq_composite", ".pdf")), device = cairo_pdf, width = 8, height = 8, units = "in")
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
