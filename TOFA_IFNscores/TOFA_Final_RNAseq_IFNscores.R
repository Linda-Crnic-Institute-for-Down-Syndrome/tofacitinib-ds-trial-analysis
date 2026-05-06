################################################
# Title: Tofacitinib for Immune Skin Conditions in Down Syndrome: Analysis of Interferon (IFN) scores 
# Author(s):
#   - Matthew Galbraith
# Affiliation(s):
#   - Linda Crnic Institute for Down Syndrome
#   - University of Colorado Anschutz
################################################

### Summary:  
# IFN signaling scores based on gene expression, measured by RNAseq of whole blood.
# Paired-end strand-specific globin-depleted polyA+ libraries generated and sequenced by Novogene.
# See README.md for more details.
#  

### Data type(s):
# Clinical trial (TOFA) datasets:
#    * Participant-level metadata; Available on request.
#    * Visit/Event-level metadata; Available on request.
#    * Baseline obesity status; Available on request.
#    * COVID-19 history; Available on request.
#    * PAXgene whole blood RNAseq data (RPKMs); DOI: 10.5281/zenodo.19954573
#
# Human Trisome Project (HTP) datasets:
#    * Participant-level metadata; DOI: 10.5281/zenodo.19962380
#    * Visit/Event-level metadata; DOI: 10.5281/zenodo.19962380
#    * PAXgene whole blood RNAseq data (RPKMs); DOI: 10.5281/zenodo.20044079
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
tofa_rpkms_file <- here("data/TOFA_PAXgene_RNAseq_RPKMs_zenodo_v1.txt.gz") # Source: 10.5281/zenodo.19954573
tofa_participant_meta_data_file <- here("data/TOFA_Participant_metadata_zenodo_v1.txt") # Source: Available on request
tofa_visit_meta_data_file <- here("data/TOFA_Visit_metadata_zenodo_v1.txt") # Source: Available on request
tofa_baseline_obesity_file <- here("data/TOFA_Baseline_Obesity_Status_zenodo_v1.txt") # Source: Available on request
tofa_covid_history_file <- here("data/TOFA_COVID_History_zenodo_v1.txt") # Source: Available on request
# Human Trisome Project datasets:
htp_participant_meta_data_file <- here("data/HTP_Participant_metadata_zenodo_v1.txt") # Source: 10.5281/zenodo.19962380
htp_visit_meta_data_file <- here("data/HTP_Visit_metadata_zenodo_v1.txt") # Source: 10.5281/zenodo.19962380
htp_rpkms_data_file <- here("data/HTP_PAXgene_RNAseq_RPKMs_long_zenodo_v1.txt.gz") # Source: 10.5281/zenodo.20044079
#
# Other parameters:
standard_colors <- c("Baseline" = "#999999", "2 week" = "#c6dbef", "8 week" = "#9ecae1", "16 week" = "#6baed6", "40 week" = "#4292c6")
#
out_file_prefix <- "TOFA_RNAseq_IFNscores_" # should match this script title
#
source(here("helper_functions.R")) # load helper functions
#

# 1 Read in data ----
## 1.1 Read in TOFA RPKMs data ----
tofa_rpkms <- tofa_rpkms_file |> 
  read_tsv() |> 
  rename(RPKM = Value)
#
tofa_rpkms # 11,707,766 rows
tofa_rpkms |> distinct(ParticipantID) # 42 Participants
tofa_rpkms |> distinct(VisitID) # 193 VisitIDs = PAXgene Whole Blood RNA samples
tofa_rpkms |> distinct(EnsemblID) # 60,662 EnsemblIDs
#


## 1.2. Read in TOFA meta data ----
### 1.2.2 Participant level meta data ----
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
### 1.2.3 Event/Visit level meta data ----
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

## 1.3 Join TOFA meta data with RPKMs ----
tofa_rpkms <- tofa_rpkms |> 
  inner_join(tofa_visit_meta_data) |> # returns 11,707,766 of 11,707,766 rows
  inner_join(tofa_participant_meta_data) |> # returns 11,707,766 of 11,707,766 rows
  filter(Endpoint_eligible == TRUE) # ensure only endpoint eligible participants
#
### 1.3.1 ParticipantID vs Event_Name summary ---- 
tofa_rpkms |> 
  distinct(ParticipantID, VisitID, Event_Name) |> 
  arrange(Event_Name) |> 
  pivot_wider(names_from = Event_Name, values_from = VisitID) |> 
  print(n = Inf)
# Note: Not all Participants completed extension to 40 weeks and some have missing samples.

## 1.4 Read in HTP data ----
### 1.4.1 Read in HTP RPKMs data ----
htp_rpkms_data <- htp_rpkms_data_file |> 
  read_tsv() |> 
  rename(RPKM = Value)
#
### 1.4.2 Read in HTP Participant level meta data ----
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
### 1.4.3 Read in HTP Visit/Event level meta data ----
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


# 2 Data exploration ----
## 2.1 Plot individual genes -----
### 2.1.1 Sina plot(s): Baseline vs 16 weeks ----
tofa_rpkms |> 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(RPKM)) |> 
  filter(is.finite(RPKM)) |> 
  add_count(VisitID, name = "pair") |> 
  filter(pair >= 2) |>  # require at least 2 data points
  # Select genes:
  filter(Gene_name %in% c("RSAD2", "IFI44L", "ISG15", "BPGM", "GMPR", "IFI27")) |> 
  #
  group_by(Gene_name) |> 
  group_split() |> 
  map( # using anonymous function: \(x)
    \(x) x |> 
      arrange(Gene_name, Event_Name) |> 
      mutate( # add labels with n
        n = n_distinct(VisitID), 
        .by = Event_Name,
        label = paste0(str_extract(Event_Name, "^[[:alpha:]]|^\\d{1,2}"), "\n(n=", n, ")"),
        label = fct_inorder(label)
      ) |> 
      mutate(Event_Name = fct_drop(Event_Name)) |> # prevent empty x values
      complete(Gene_name, VisitID, Event_Name) %>% # Prevent error from stat_signif()
      #
      {ggplot(data = ., aes(Event_Name, log2(RPKM), color = Event_Name)) +
          geom_line(aes(group = VisitID), color = "grey90") +
          geom_sina(maxwidth = 0.5) +
          geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
          scale_color_manual(values = standard_colors) +
          facet_wrap(~ Gene_name, scales = "free_y", nrow = 1) +
          theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            aspect.ratio = 2.5,
            panel.border = element_blank(),
            axis.line = element_line(colour = "black")
          ) +
          labs(x = NULL) +
          scale_x_discrete(
            labels = . |>
              filter(!is.na(label)) |>
              distinct(Gene_name, Event_Name, label) |>
              select(Event_Name, label) |>
              transmute(as.character(Event_Name), as.character(label)) |>
              deframe()
          )
      }
  ) |> 
  patchwork::wrap_plots() +
  patchwork::plot_layout(nrow = 1) +
  patchwork::plot_annotation(
    title = "Individual genes: Baseline vs. Treatment"
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_individual_B_16_weeks", ".pdf")), device = cairo_pdf, width = 25, height = 4, units = "in")
#


# 3 Calculate RPKM Z-scores ----
## 3.1 TOFA -----
tofa_zscores_unadj <- tofa_rpkms |> 
  mutate(
    mean = RPKM %>% log2() %>% na_if(-Inf) %>% mean(na.rm = TRUE), # need to deal with 0 counts
    sd = RPKM %>% log2() %>% na_if(-Inf) %>% sd(na.rm = TRUE),  # need to deal with 0 counts
    .by = EnsemblID
  ) |> 
  mutate(
    zscore = (log2(RPKM) - mean) / sd
  ) %>%
  arrange(EnsemblID)
#
## 3.2 HTP -----
# Here using RPKMs adjusted for Sex, Age, and Sample source using limma::removeBatchEffect()
htp_zscores <- htp_rpkms_data |> 
  inner_join(htp_visit_meta_data) |> 
  select(ParticipantID, VisitID, EnsemblID, Gene_name, RPKM = Value_SexAgeSource_adj) |> 
  mutate(
    mean = RPKM %>% log2() %>% na_if(-Inf) %>% mean(na.rm = TRUE), # need to deal with 0 counts
    sd = RPKM %>% log2() %>% na_if(-Inf) %>% sd(na.rm = TRUE),  # need to deal with 0 counts
    .by = EnsemblID
  ) |> 
  mutate(
    zscore = (log2(RPKM) - mean) / sd
  ) %>%
  arrange(EnsemblID)
#


# 4 Calculate IFN scores ----
# Interferon scores are calculated from 16 interferon-stimulated genes (ISGs) that are elevated in people with trisomy 21
isgs_16 <- c("CCL5", "IFI27", "GMPR", "FCGR1A", "USP18", "BPGM", "GZMA", "CMKLR1", "PLSCR1", "IFI44", "CD274", "IFI44L", "IRF7", "RSAD2", "ISG15", "IFITM1")
#
## 4.1 TOFA -----
tofa_IFN_score_unadj_16genes <- tofa_zscores_unadj |> 
  filter(Gene_name %in% isgs_16) |>
  # need to remove these to prevent -Inf sum scores:
  filter(is.finite(zscore)) |> # drops 2 rows: IFI27 and USP18
  #
  summarize(
    IFN_score = sum(zscore),
    n_genes = n_distinct(EnsemblID),
    .by = c(ParticipantID, VisitID, Event_Name)
  ) |> 
  arrange(VisitID, Event_Name)
#
### Export
tofa_IFN_score_unadj_16genes |> 
  write_tsv(file = here("results", paste0(out_file_prefix, "IFN_score_unadj_16genes", ".txt.gz")))
#
## 4.2 HTP -----
htp_IFN_score_16genes <- htp_zscores |> 
  filter(Gene_name %in% isgs_16) |>
  # need to remove these to prevent -Inf sum scores:
  # drops 0 rows:
  filter(is.finite(zscore)) |> 
  #
  summarize(
    IFN_score = sum(zscore),
    n_genes = n_distinct(EnsemblID),
    .by = c(ParticipantID, VisitID)
  )
# 
#### Sina plot: Ages 12-40 ----
# matching age range to clinical trial
htp_IFN_score_16genes |> 
  inner_join(htp_visit_meta_data) |> 
  inner_join(htp_participant_meta_data) |> 
  filter(dplyr::between(Age_years_at_visit, 12, 40)) |> # keeps 282 rows
  mutate(extreme = rstatix::is_extreme(IFN_score), .by = Karyotype) |>
  filter(extreme == FALSE) |> # remove extreme outliers: drops 7 rows
  ggplot(aes(Karyotype, IFN_score, color = Karyotype)) +
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
    title = "HTP IFN scores",
    subtitle = "ages 12-40; extreme outliers removed",
    x = NULL
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_IFN_scores_HTP_12-40", ".pdf")), device = cairo_pdf, width = 10, height = 5, units = "in")
#
#### Unpaired Wilcox / Mann-Whitney ----
# matching age range to clinical trial
htp_IFN_score_16genes |> 
  inner_join(htp_visit_meta_data) |> 
  inner_join(htp_participant_meta_data) |> 
  filter(dplyr::between(Age_years_at_visit, 12, 40)) |> # keeps 282 rows
  mutate(extreme = rstatix::is_extreme(IFN_score), .by = Karyotype) |>
  filter(extreme == FALSE) |> # remove extreme outliers: drops 7 rows
  rstatix::wilcox_test(
    formula = IFN_score ~ Karyotype,
    ref.group = "T21",
    paired = FALSE,
    detailed = TRUE,
    p.adjust.method = "none"
  )
#


# 5 Stats ----
## 5.1 Check data distributions ----
### 5.1.1 check for outlier values ----
tofa_IFN_score_unadj_16genes |> 
  mutate(
    extreme = rstatix::is_extreme(IFN_score),
    outlier = rstatix::is_outlier(IFN_score)
  ) |> 
  count(outlier, extreme)
#
tofa_IFN_score_unadj_16genes |> 
  mutate(
    extreme = rstatix::is_extreme(IFN_score),
    outlier = rstatix::is_outlier(IFN_score)
  ) |> 
  ggplot(aes("IFN_score", IFN_score)) +
  geom_sina(data = . %>% filter(outlier == FALSE & extreme == FALSE), maxwidth = 0.5, color = "grey") +
  geom_sina(data = . %>% filter(outlier == TRUE & extreme == FALSE), aes(color = "outlier"), maxwidth = 0.1) +
  geom_sina(data = . %>% filter(extreme == TRUE), aes(color = "extreme"), maxwidth = 0.1) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    aspect.ratio = 1.3
  ) +
  labs(title = "IFN_score outliers")
#
### 5.1.2 Calc Differences from baseline ----
tofa_IFN_scores_diffs_2_weeks <- tofa_IFN_score_unadj_16genes |> 
  filter(Event_Name %in% c("Baseline", "2 week")) |> 
  filter(!is.na(IFN_score)) |>
  filter(is.finite(IFN_score)) |> 
  add_count(ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  select(ParticipantID, Event_Name, IFN_score) |> 
  pivot_wider(names_from = Event_Name, values_from = IFN_score) |> 
  mutate(difference = `2 week` - Baseline) |> 
  select(ParticipantID, difference) |> 
  mutate(
    extreme = rstatix::is_extreme(difference),
    outlier = rstatix::is_outlier(difference)
  )
#
tofa_IFN_scores_diffs_8_weeks <- tofa_IFN_score_unadj_16genes |> 
  filter(Event_Name %in% c("Baseline", "8 week")) |> 
  filter(!is.na(IFN_score)) |>
  filter(is.finite(IFN_score)) |> 
  add_count(ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  select(ParticipantID, Event_Name, IFN_score) |> 
  pivot_wider(names_from = Event_Name, values_from = IFN_score) |> 
  mutate(difference = `8 week` - Baseline) |> 
  select(ParticipantID, difference) |> 
  mutate(
    extreme = rstatix::is_extreme(difference),
    outlier = rstatix::is_outlier(difference)
  )
#
tofa_IFN_scores_diffs_16_weeks <- tofa_IFN_score_unadj_16genes |> 
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(IFN_score)) |>
  filter(is.finite(IFN_score)) |> 
  add_count(ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  select(ParticipantID, Event_Name, IFN_score) |> 
  pivot_wider(names_from = Event_Name, values_from = IFN_score) |> 
  mutate(difference = `16 week` - Baseline) |> 
  select(ParticipantID, difference) |> 
  mutate(
    extreme = rstatix::is_extreme(difference),
    outlier = rstatix::is_outlier(difference)
  )
#
tofa_IFN_scores_diffs_40_weeks <- tofa_IFN_score_unadj_16genes |> 
  filter(Event_Name %in% c("Baseline", "40 week")) |> 
  filter(!is.na(IFN_score)) |>
  filter(is.finite(IFN_score)) |> 
  add_count(ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  select(ParticipantID, Event_Name, IFN_score) |> 
  pivot_wider(names_from = Event_Name, values_from = IFN_score) |> 
  mutate(difference = `40 week` - Baseline) |> 
  select(ParticipantID, difference) |> 
  mutate(
    extreme = rstatix::is_extreme(difference),
    outlier = rstatix::is_outlier(difference)
  )
#
### 5.1.3 Check extreme outliers in differences ----
tofa_IFN_scores_diffs_2_weeks |> 
  count(outlier, extreme)
tofa_IFN_scores_diffs_8_weeks |> 
  count(outlier, extreme)
tofa_IFN_scores_diffs_16_weeks |> 
  count(outlier, extreme)
tofa_IFN_scores_diffs_40_weeks |> 
  count(outlier, extreme)
#
### 5.1.4 Plot differences from Baseline ---- 
sina_diffs_2wks <- tofa_IFN_scores_diffs_2_weeks |> 
  ggplot(aes(x = "IFN score", difference)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_sina(data = . %>% filter(outlier == FALSE & extreme == FALSE), maxwidth = 0.5, color = "grey") +
  geom_sina(data = . %>% filter(outlier == TRUE), aes(color = "outlier"), maxwidth = 0.1) +
  geom_sina(data = . %>% filter(extreme == TRUE), aes(color = "extreme"), maxwidth = 0.1) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  theme(
    legend.title = element_blank(),
    aspect.ratio = 1.3
  ) +
  labs(title = "2 weeks", x = NULL)
#
sina_diffs_8wks <- tofa_IFN_scores_diffs_8_weeks |> 
  ggplot(aes(x = "IFN score", difference)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_sina(data = . %>% filter(outlier == FALSE & extreme == FALSE), maxwidth = 0.5, color = "grey") +
  geom_sina(data = . %>% filter(outlier == TRUE), aes(color = "outlier"), maxwidth = 0.1) +
  geom_sina(data = . %>% filter(extreme == TRUE), aes(color = "extreme"), maxwidth = 0.1) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  theme(
    legend.title = element_blank(),
    aspect.ratio = 1.3
  ) +
  labs(title = "8 weeks", x = NULL)
#
sina_diffs_16wks <- tofa_IFN_scores_diffs_16_weeks |> 
  ggplot(aes(x = "IFN score", difference)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_sina(data = . %>% filter(outlier == FALSE & extreme == FALSE), maxwidth = 0.5, color = "grey") +
  geom_sina(data = . %>% filter(outlier == TRUE), aes(color = "outlier"), maxwidth = 0.1) +
  geom_sina(data = . %>% filter(extreme == TRUE), aes(color = "extreme"), maxwidth = 0.1) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  theme(
    legend.title = element_blank(),
    aspect.ratio = 1.3
  ) +
  labs(title = "16 weeks", x = NULL)
#
sina_diffs_40wks <- tofa_IFN_scores_diffs_40_weeks |> 
  ggplot(aes(x = "IFN score", difference)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_sina(data = . %>% filter(outlier == FALSE & extreme == FALSE), maxwidth = 0.5, color = "grey") +
  geom_sina(data = . %>% filter(outlier == TRUE), aes(color = "outlier"), maxwidth = 0.1) +
  geom_sina(data = . %>% filter(extreme == TRUE), aes(color = "extreme"), maxwidth = 0.1) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  theme(
    legend.title = element_blank(),
    aspect.ratio = 1.3
  ) +
  labs(title = "40 weeks", x = NULL)
#
sina_diffs_2wks + expand_limits(y = c(-50, 30)) +
  sina_diffs_8wks + expand_limits(y = c(-50, 30)) +
  sina_diffs_16wks + expand_limits(y = c(-50, 30)) +
  sina_diffs_40wks + expand_limits(y = c(-50, 30)) +
  plot_layout(guides = "collect", nrow = 1) +
  plot_annotation(title = "Distributions of differences vs. baseline")
#
### 5.1.5 Check if differences are normally distributed ----
tofa_IFN_scores_diffs_2_weeks |>
  ggpubr::ggqqplot("difference", title = "Q-Q plot: 2 week differences")
tofa_IFN_scores_diffs_8_weeks |>
  ggpubr::ggqqplot("difference", title = "Q-Q plot: 8 week differences")
tofa_IFN_scores_diffs_16_weeks |>
  ggpubr::ggqqplot("difference", title = "Q-Q plot: 16 week differences")
tofa_IFN_scores_diffs_40_weeks |>
  ggpubr::ggqqplot("difference", title = "Q-Q plot: 40 week differences")
#

# 5.2 Paired t test for 16 week endpoint  --------
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
### 5.2.1 Run paired t tests 16 weeks ----
# Only 16 week timepoint is tested as a trial endpoint.
tofa_IFN_scores_Ttest_res_16weeks <- tofa_IFN_score_unadj_16genes |> 
  arrange(ParticipantID, Event_Name) |> # ensure correct sort order for rstatix::t_test
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(IFN_score)) |>
  filter(is.finite(IFN_score)) |> 
  add_count(ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  mutate(score_name = "IFN score") |> 
  group_by(score_name) |> 
  rstatix::t_test(
    formula = IFN_score ~ Event_Name,
    ref.group = "16 week",
    paired = TRUE,
    var.equal = TRUE, # set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
    detailed = TRUE,
    p.adjust.method = "none"
  )
#
### 5.2.2 Calculate t test effect size ----
# The effect size for a paired-samples t-test can be calculated by dividing the
# mean difference by the standard deviation of the difference, as shown below.
# Cohen’s formula:
# d = mean(D)/sd(D), where D is the differences of the paired samples values.
tofa_IFN_scores_Ttest_effsize_16weeks <- tofa_IFN_score_unadj_16genes |> 
  arrange(ParticipantID, Event_Name) |> # ensure correct sort order for rstatix::cohens_d
  filter(Event_Name %in% c("Baseline", "16 week")) |> 
  filter(!is.na(IFN_score)) |>
  filter(is.finite(IFN_score)) |> 
  add_count(ParticipantID, name = "pair") |> 
  filter(pair == 2) |>  # require pairs
  mutate(score_name = "IFN score") |> 
  group_by(score_name) |> 
  rstatix::cohens_d(
    formula = IFN_score ~ Event_Name,
    ref.group = "16 week",
    paired = TRUE,
    var.equal = TRUE, # set this to use Student's vs Welch t-test; The two methods give very similar results unless both the group sizes and the standard deviations are very different.
  )
#
### 5.2.3 Compile t test results ----
tofa_IFN_scores_Ttest_res_16weeks_full <- tofa_IFN_scores_Ttest_res_16weeks |> 
  inner_join(tofa_IFN_scores_Ttest_effsize_16weeks) |> 
  mutate(mean_diff = estimate) |> 
  inner_join(tofa_IFN_scores_diffs_16_weeks |> mutate(score_name = "IFN score") |> group_by(score_name) |> summarize(median_diff = median(difference))) |> 
  select(score_name, mean_diff, median_diff, p, effsize, magnitude, .y., group1, group2, everything())
#
tofa_IFN_scores_Ttest_res_16weeks_full
#
### 5.2.4 Export t test results ----
list(
  "ttest_results" = tofa_IFN_scores_Ttest_res_16weeks_full |> 
    select(
      Score_name = score_name,
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
  export_excel(filename = "Ttest_results")
#


# # 6 Paired plots (NOT USED)  --------
# ## 6.1 Baseline vs 2 weeks vs 8 weeks vs 16 vs 40 weeks ----
# tofa_IFN_score_unadj_16genes |> 
#   filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> 
#   filter(!is.na(IFN_score)) |> 
#   filter(is.finite(IFN_score)) |> 
#   add_count(ParticipantID, name = "pair") |> 
#   filter(pair >= 2) |>  # require at least 2 data points
#   ggplot(aes(Event_Name, IFN_score, color = Event_Name)) +
#   scale_color_manual(values = standard_colors) +
#   geom_line(aes(group = ParticipantID), color = "grey90") +
#   geom_point() +
#   theme(
#     axis.text.x = element_blank(),
#     legend.title = element_blank(),
#     legend.position = "bottom",
#     aspect.ratio = 1.3,
#     panel.border = element_blank(),
#     axis.line = element_line(colour = "black")
#   ) +
#   labs(
#     title = "IFN scores: Baseline vs. TOFA",
#     subtitle = "NB: smaller N at 40 weeks",
#     x = NULL
#   )
# ggsave(filename = here("plots", paste0(out_file_prefix, "paired_plots_B_2_8_16_40_weeks", ".pdf")), device = cairo_pdf, width = 25, height = 3, units = "in")
# #

# 7 Sina plots - IFN scores  --------
## 7.2 Baseline vs 2/8/16/40 weeks  ----
tofa_IFN_score_unadj_16genes |> 
  filter(Event_Name %in% c("Baseline", "2 week", "8 week", "16 week", "40 week")) |> 
  filter(!is.na(IFN_score)) |> 
  filter(is.finite(IFN_score)) |> 
  add_count(ParticipantID, name = "pair") |> 
  filter(pair >= 2) |>  # require at least 2 data points
  #
  mutate(Score = "IFN score") |> 
  group_by(Score) |>
  group_split() |>
  map( # using anonymous function: \(x)
    \(x) x |> 
      arrange(Score, Event_Name) |>
      arrange(Event_Name) |> 
      mutate( # add labels with n
        n = n_distinct(ParticipantID),
        .by = Event_Name,
        label = paste0(str_extract(Event_Name, "^[[:alpha:]]|^\\d{1,2}"), "\n(n=", n, ")"),
        # label = paste0(str_extract(Event_Name, "^[[:alpha:]]|^\\d{1,2}"), " (n=", n, ")"),
        label = fct_inorder(label)
      ) |>
      mutate(Event_Name = fct_drop(Event_Name)) |> # prevent empty x values
      complete(Score, ParticipantID, Event_Name) %>% # Prevent error from stat_signif()
      #
      {ggplot(data = ., aes(Event_Name, IFN_score, color = Event_Name)) +
          scale_color_manual(values = standard_colors) +
          geom_line(aes(group = ParticipantID), color = "grey90") +
          geom_sina(maxwidth = 0.5) +
          geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
          theme(
            legend.title = element_blank(),
            legend.position = "bottom",
            aspect.ratio = 1.5,
            panel.border = element_blank(),
            axis.line = element_line(colour = "black")
          ) +
          labs(
            title = "IFN scores: Baseline vs. TOFA",
            subtitle = "nb: smaller N at 40 weeks",
            x = NULL
          ) +
          scale_x_discrete(
            labels = . |>
              filter(!is.na(label)) |>
              distinct(Score, Event_Name, label) |>
              select(Event_Name, label) |>
              transmute(as.character(Event_Name), as.character(label)) |>
              deframe()
          )
      }
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_plot_signif_B_2_8_16_40_weeks", ".pdf")), device = cairo_pdf, width = 15, height = 5, units = "in")
#


# 8 Sina Plots - Differences  --------
## 8.1 Baseline vs 16 weeks ----
tofa_IFN_scores_diffs_16_weeks |> 
  mutate(
    timepoint = "16 week",
    score_name = "IFN score"
         ) |> 
  ggplot(aes(timepoint, difference, color = timepoint)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_sina(maxwidth = 0.4) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  scale_color_manual(values = standard_colors) +
  facet_wrap(~score_name, nrow = 1, scales = "free_y") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    aspect.ratio = 3
  ) +
  labs(
    title = "IFN score: Difference from Baseline",
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_diffs_16_weeks", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#
## 8.2 Baseline vs 2/8/16/40 weeks ----
bind_rows(
  tofa_IFN_scores_diffs_2_weeks |> mutate(timepoint = "2 week"),
  tofa_IFN_scores_diffs_8_weeks |> mutate(timepoint = "8 week"),
  tofa_IFN_scores_diffs_16_weeks |> mutate(timepoint = "16 week"),
  tofa_IFN_scores_diffs_40_weeks |> mutate(timepoint = "40 week")
) |> 
  mutate(timepoint = fct_relevel(timepoint, c("2 week", "8 week", "16 week", "40 week"))) |> 
  ggplot(aes(timepoint, difference, color = timepoint)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(aes(group = ParticipantID), color = "grey90") +
  geom_sina(maxwidth = 0.4) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  scale_color_manual(values = standard_colors) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    aspect.ratio = 1.3
  ) +
  labs(
    title = "IFN scores: Differences vs Baseline",
    subtitle = "NB: smaller N at 40 weeks",
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_diffs_2_8_16_40_weeks", ".pdf")), device = cairo_pdf, width = 10, height = 6, units = "in")
#


# 12 Mixed effects linear regression models (non-stratified) ----
# Time as a categorical variable
#   No assumptions about trajectory
#   Directly compares each visit to baseline
# Random effects
#   Handles missing at random (MAR) data.
#   Accounts for within-subject correlation.
## 12.2 Set up models ----
### 12.2.1 Mixed effects LM: Event_name + 1|ParticipantID ----
tofa_IFN_scores_lm_mixedParticipantID <- tofa_IFN_score_unadj_16genes |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(!is.na(IFN_score)) |>
  filter(is.finite(IFN_score)) |> 
  mutate(score_name = "IFN score") |> 
  select(score_name, ParticipantID, Sex, Age_Baseline, Event_Name, Value = IFN_score) |> 
  nest(data = -score_name) |> 
  mutate(
    fit = map(data, ~ lmerTest::lmer(Value ~ Event_Name + (1 | ParticipantID), data = .x, REML = TRUE)),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
### 12.2.2 Mixed effects LM: Sex + Event_name + 1|ParticipantID ----
tofa_IFN_scores_lm_fixedSex_mixedParticipantID <- tofa_IFN_score_unadj_16genes |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(!is.na(IFN_score)) |>
  filter(is.finite(IFN_score)) |> 
  mutate(score_name = "IFN score") |> 
  select(score_name, ParticipantID, Sex, Age_Baseline, Event_Name, Value = IFN_score) |> 
  nest(data = -score_name) |> 
  mutate(
    fit = map(data, ~ lmerTest::lmer(Value ~ Sex + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE)),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
### 12.2.3 Mixed effects LM: Age + Event_name + 1|ParticipantID ----
tofa_IFN_scores_lm_fixedAge_mixedParticipantID <- tofa_IFN_score_unadj_16genes |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(!is.na(IFN_score)) |>
  filter(is.finite(IFN_score)) |> 
  mutate(score_name = "IFN score") |> 
  select(score_name, ParticipantID, Sex, Age_Baseline, Event_Name, Value = IFN_score) |> 
  nest(data = -score_name) |> 
  mutate(
    fit = map(data, ~ lmerTest::lmer(Value ~ Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE)),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
### 12.2.4 Mixed effects LM: Sex + Age + Event_name + 1|ParticipantID (PREFERRED) ----
tofa_IFN_scores_lm_fixedSexAge_mixedParticipantID <- tofa_IFN_score_unadj_16genes |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(!is.na(IFN_score)) |>
  filter(is.finite(IFN_score)) |> 
  mutate(score_name = "IFN score") |> 
  select(score_name, ParticipantID, Sex, Age_Baseline, Event_Name, Value = IFN_score) |> 
  nest(data = -score_name) |> 
  mutate(
    fit = map(data, ~ lmerTest::lmer(Value ~ Sex + Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE)),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#

## 12.3 Compare models ----
### 12.3.1 AIC/BIC ----
# Mixed models
tofa_IFN_scores_lm_mixedParticipantID |> unnest(glanced) |> select(AIC, BIC) # 1292
tofa_IFN_scores_lm_fixedSex_mixedParticipantID |> unnest(glanced) |> select(AIC, BIC) # 1287; decent improvement
tofa_IFN_scores_lm_fixedAge_mixedParticipantID |> unnest(glanced) |> select(AIC, BIC) # 1296; no improvement
tofa_IFN_scores_lm_fixedSexAge_mixedParticipantID |> unnest(glanced) |> select(AIC, BIC) # 1291; minor improvement; PREFERRED
#
### 12.3.2 Likelihood ratio tests ----
anova(
  tofa_IFN_scores_lm_mixedParticipantID |> pluck("fit", 1),
  tofa_IFN_scores_lm_fixedSex_mixedParticipantID |> pluck("fit", 1),
  test="LRT"
) |>  # refitting model(s) with ML (instead of REML)
  tidy() # p.value = 0.0376
anova(
  tofa_IFN_scores_lm_mixedParticipantID |> pluck("fit", 1),
  tofa_IFN_scores_lm_fixedAge_mixedParticipantID |> pluck("fit", 1),
  test="LRT"
) |>  # refitting model(s) with ML (instead of REML)
  tidy() # p.value = 0.540
anova(
  tofa_IFN_scores_lm_mixedParticipantID |> pluck("fit", 1),
  tofa_IFN_scores_lm_fixedSexAge_mixedParticipantID |> pluck("fit", 1),
  test="LRT"
) |>  # refitting model(s) with ML (instead of REML)
  tidy() # p.value = 0.0987
#


## 12.4 Model results ----
### 12.4.2 Mixed ParticipantID with fixed SexAge -----
tofa_IFN_scores_lm_fixedSexAge_mixedParticipantID |> 
  unnest(tidied) |> 
  select(score_name, group, term, estimate, p.value)
#
tofa_IFN_scores_lm_fixedSexAge_mixedParticipantID_results <- tofa_IFN_scores_lm_fixedSexAge_mixedParticipantID |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, "Event_Name")) |>
  select(score_name, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = score_name) |> # NB: by score_name
  mutate(
    level = "",
    term = str_remove(term, "^Event_Name"),
    .after = score_name,
  )
#
tofa_IFN_scores_lm_fixedSexAge_mixedParticipantID_results
#


# 13 Mixed effects linear regression models (stratified) ----
## 13.1 By Sex ----
lm_fixedAge_mixedParticipantID_by_Sex <- tofa_IFN_score_unadj_16genes |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  filter(!is.na(IFN_score)) |>
  filter(is.finite(IFN_score)) |> 
  mutate(score_name = "IFN score") |> 
  select(score_name, ParticipantID, Sex, Age_Baseline, Event_Name, Value = IFN_score) |> 
  nest(data = -c(score_name, Sex)) |> 
  mutate(
    fit = map(data, ~ lmerTest::lmer(Value ~ Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE)),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
lm_fixedAge_mixedParticipantID_by_Sex |> unnest(glanced) |> select(Sex, AIC, BIC) # 577 / 699; no improvement
#
lm_fixedAge_mixedParticipantID_by_Sex_results <- lm_fixedAge_mixedParticipantID_by_Sex |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, "Event_Name")) |>
  select(score_name, level = Sex, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = c(score_name, level)) |> # NB: by score_name + level
  mutate(term = str_remove(term, "^Event_Name")) |> 
  arrange(level)
#

## 13.2 By Age group ----
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
lm_fixedSex_mixedParticipantID_by_Age_group <- tofa_IFN_score_unadj_16genes |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  inner_join(tofa_baseline_age_groups) |> # add age groups info
  filter(!is.na(IFN_score)) |>
  filter(is.finite(IFN_score)) |> 
  mutate(score_name = "IFN score") |> 
  select(score_name, ParticipantID, Sex, Age_group, Event_Name, Value = IFN_score) |> 
  nest(data = -c(score_name, Age_group)) |> 
  mutate(
    fit = map(data, ~ lmerTest::lmer(Value ~ Sex + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE)),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
lm_fixedSex_mixedParticipantID_by_Age_group |> unnest(glanced) |> select(Age_group, AIC, BIC) # 1011 / 257; minor improvement
#
lm_fixedSex_mixedParticipantID_by_Age_group_results <- lm_fixedSex_mixedParticipantID_by_Age_group |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, "Event_Name")) |>
  select(score_name, level = Age_group, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = c(score_name, level)) |> # NB: by score_name + level
  mutate(term = str_remove(term, "^Event_Name")) |> 
  arrange(level)
#

## 13.3 By Obesity ----
tofa_baseline_obesity <- tofa_baseline_obesity_file |> 
  read_tsv()
tofa_baseline_obesity |> count(baseline_obesity_status)
#
lm_fixedSexAge_mixedParticipantID_by_Obesity <- tofa_IFN_score_unadj_16genes |> 
  inner_join(tofa_participant_meta_data) |> 
  inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
  inner_join(tofa_baseline_obesity) |> # add Obesity info
  filter(!is.na(IFN_score)) |>
  filter(is.finite(IFN_score)) |> 
  mutate(score_name = "IFN score") |> 
  select(score_name, ParticipantID, Sex, Age_Baseline, baseline_obesity_status, Event_Name, Value = IFN_score) |> 
  nest(data = -c(score_name, baseline_obesity_status)) |> 
  mutate(
    fit = map(data, ~ lmerTest::lmer(Value ~ Sex + Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE)),
    tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
    glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
    augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
  )
#
lm_fixedSexAge_mixedParticipantID_by_Obesity |> unnest(glanced) |> select(baseline_obesity_status, AIC, BIC) # 550 / 714; minor improvement
#
lm_fixedSexAge_mixedParticipantID_by_Obesity_results <- lm_fixedSexAge_mixedParticipantID_by_Obesity |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, "Event_Name")) |>
  select(score_name, level = baseline_obesity_status, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = c(score_name, level)) |> # NB: by score_name + level
  mutate(term = str_remove(term, "^Event_Name")) |> 
  arrange(level)
#

## 13.4 By COVID ----
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
  tofa_IFN_score_unadj_16genes |> 
    inner_join(tofa_participant_meta_data) |> 
    inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
    inner_join(tofa_covid_history |> filter(Event_Name == "16 week") |> select(ParticipantID, COVID_event_hx)) |> # add covid info for 16 weeks
    filter(Event_Name %in% c("Baseline", "16 week", "40 week")) |> # keep all participants for mixed effects
    filter(!is.na(IFN_score)) |>
    filter(is.finite(IFN_score)) |> 
    mutate(score_name = "IFN score") |> 
    select(score_name, ParticipantID, Sex, Age_Baseline, COVID_event_hx, Event_Name, Value = IFN_score) |> 
    nest(data = -c(score_name, COVID_event_hx)) |> 
    mutate(
      test = "16 week",
      fit = map(data, ~ lmerTest::lmer(Value ~ Sex + Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE)),
      tidied = map(fit, \(.x) broom.mixed::tidy(.x, conf.int = TRUE)), # see ?broom.mixed:::tidy.lme
      glanced = map(fit, broom.mixed::glance), # see ?broom.mixed:::glance.lme
      augmented = map(fit, broom.mixed::augment), # see ?broom.mixed:::augment.lme
    ),
  tofa_IFN_score_unadj_16genes |> 
    inner_join(tofa_participant_meta_data) |> 
    inner_join(tofa_visit_meta_data |> filter(Event_Name == "Baseline") |> select(ParticipantID, Age_Baseline = Age_years_at_visit)) |> 
    inner_join(tofa_covid_history |> filter(Event_Name == "40 week") |> select(ParticipantID, COVID_event_hx)) |> # add covid info for 40 weeks
    filter(Event_Name %in% c("Baseline", "16 week", "40 week")) |> # keep all participants for mixed effects
    filter(!is.na(IFN_score)) |>
    filter(is.finite(IFN_score)) |> 
    mutate(score_name = "IFN score") |> 
    select(score_name, ParticipantID, Sex, Age_Baseline, COVID_event_hx, Event_Name, Value = IFN_score) |> 
    nest(data = -c(score_name, COVID_event_hx)) |> 
    mutate(
      test = "40 week",
      fit = map(data, ~ lmerTest::lmer(Value ~ Sex + Age_Baseline + Event_Name + (1 | ParticipantID), data = .x, REML = TRUE)),
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
  select(COVID_event_hx, term, AIC, BIC)
#
lm_fixedSexAge_mixedParticipantID_by_COVID_results <- lm_fixedSexAge_mixedParticipantID_by_COVID |> 
  unnest(c(tidied, glanced)) |> 
  filter(str_detect(term, regex(test))) |> # keep only relevant tests
  filter(str_detect(term, "Event_Name")) |>
  select(score_name, level = COVID_event_hx, term, n_obs = nobs, estimate, conf.low, conf.high, statistic, p.value) |> 
  mutate(BH_padj = p.adjust(p.value, method = "BH"), .by = c(score_name, level)) |> # NB: by score_name + level
  mutate(term = str_remove(term, "^Event_Name")) |> 
  arrange(level) |> 
  mutate(level = fct_recode(level, No = "no", Yes = "yes")) # relabel COVID levels
#


# 14 Export LM results ----
list(
  "LMM results" = list(
    "Overall" = tofa_IFN_scores_lm_fixedSexAge_mixedParticipantID_results,
    "Sex" = lm_fixedAge_mixedParticipantID_by_Sex_results,
    "Age_group" = lm_fixedSex_mixedParticipantID_by_Age_group_results,
    "Obesity" = lm_fixedSexAge_mixedParticipantID_by_Obesity_results,
    "COVID" = lm_fixedSexAge_mixedParticipantID_by_COVID_results
  ) |> 
    bind_rows(.id = "stratifier") |> 
    select(
      Stratifier = stratifier,
      Level = level,
      Score_name = score_name,
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

# 15 LM Forest plot(s) ----
list(
  "Overall" = tofa_IFN_scores_lm_fixedSexAge_mixedParticipantID_results,
  "Sex" = lm_fixedAge_mixedParticipantID_by_Sex_results,
  "Age_group" = lm_fixedSex_mixedParticipantID_by_Age_group_results,
  "Obesity" = lm_fixedSexAge_mixedParticipantID_by_Obesity_results,
  "COVID" = lm_fixedSexAge_mixedParticipantID_by_COVID_results
) |> 
  bind_rows(.id = "stratifier") |> 
  mutate(stratifier = fct_relevel(stratifier, c("Overall", "Age_group", "Sex", "Obesity", "COVID"))) |> 
  arrange(stratifier) |> 
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
    title = "IFN score: Treatment vs. Baseline",
    subtitle = "LMM",
    x = "Estimate",
    y = NULL
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "Forest_plot_LMMs", ".pdf")), device = cairo_pdf, width = 8, height = 5, units = "in")
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


