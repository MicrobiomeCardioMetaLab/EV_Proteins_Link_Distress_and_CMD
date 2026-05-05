# Table S1: Baseline characteristics of the discovery and replication cohorts.
# This script summarizes demographic characteristics, lifestyle factors, emotional distress scores,
# and cardiometabolic disease status by anxiety and depression groups.
# Outputs include cohort-specific Table S1 files and binary SAS/SDS group counts.

rm(list = ls())
library(dplyr)
library(readr)

project_dir <- "/path/to/project"
data_dir <- file.path(project_dir, "data")
out_dir <- file.path(project_dir, "results", "TableS1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cohort_files <- list(
  Discovery  = file.path(data_dir, "discovery_score_covariates_cmd.csv"),
  Replication = file.path(data_dir, "replication_score_covariates_cmd.csv")
)

case_code <- 2
control_code <- 1

# Format continuous variables as mean (SD).
format_mean_sd <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
# Format binary variables as n (%).
format_n_pct <- function(x, denom) sprintf("%d (%.1f%%)", sum(x == 1, na.rm = TRUE), 100 * sum(x == 1, na.rm = TRUE) / denom)

# Harmonize demographic, lifestyle, emotional distress, and CMD variables.
prepare_baseline_data <- function(dat) {
  dat %>%
    mutate(
      Female = case_when(Gender == 1 ~ 1, Gender == 2 ~ 0, TRUE ~ NA_real_),
      Male = case_when(Gender == 2 ~ 1, Gender == 1 ~ 0, TRUE ~ NA_real_),
      Current_smoker = case_when(Smoking_status_now %in% c(1, 2) ~ 1, Smoking_status_now %in% c(3, 4) ~ 0, TRUE ~ NA_real_),
      Current_alcohol_drinker = case_when(Alcohol_intake_status == 1 ~ 1, Alcohol_intake_status %in% c(2, 3) ~ 0, TRUE ~ NA_real_),
      High_education = case_when(Education_level >= 5 ~ 1, Education_level %in% c(1, 2, 3, 4) ~ 0, TRUE ~ NA_real_),
      T2D_case = case_when(T2D == case_code ~ 1, T2D == control_code ~ 0, TRUE ~ NA_real_),
      Hypertension_case = case_when(Hypertension == case_code ~ 1, Hypertension == control_code ~ 0, TRUE ~ NA_real_),
      Dyslipidemia_case = case_when(Dyslipidemia == case_code ~ 1, Dyslipidemia == control_code ~ 0, TRUE ~ NA_real_),
      MetS_case = case_when(MetS == case_code ~ 1, MetS == control_code ~ 0, TRUE ~ NA_real_),
      Anxiety = case_when(SAS_score >= 50 ~ 1, SAS_score < 50 ~ 0, TRUE ~ NA_real_),
      Depression = case_when(SDS_score >= 50 ~ 1, SDS_score < 50 ~ 0, TRUE ~ NA_real_)
    )
}

# Summarize baseline characteristics for one cohort or subgroup.
summarise_baseline <- function(dat) {
  n <- nrow(dat)
  tibble(
    Characteristics = c(
      "N",
      "SAS score, mean (s.d.)",
      "SDS score, mean (s.d.)",
      "Age, mean (s.d.), years",
      "Female sex, n (%)",
      "Male sex, n (%)",
      "BMI, mean (s.d.), kg m-2",
      "High education level, n (%)",
      "Current smoker, n (%)",
      "Current alcohol drinker, n (%)",
      "Type 2 diabetes, n (%)",
      "Hypertension, n (%)",
      "Dyslipidemia, n (%)",
      "Metabolic syndrome, n (%)"
    ),
    Value = c(
      as.character(n),
      format_mean_sd(dat$SAS_score),
      format_mean_sd(dat$SDS_score),
      format_mean_sd(dat$Age),
      format_n_pct(dat$Female, n),
      format_n_pct(dat$Male, n),
      format_mean_sd(dat$BMI),
      format_n_pct(dat$High_education, n),
      format_n_pct(dat$Current_smoker, n),
      format_n_pct(dat$Current_alcohol_drinker, n),
      format_n_pct(dat$T2D_case, n),
      format_n_pct(dat$Hypertension_case, n),
      format_n_pct(dat$Dyslipidemia_case, n),
      format_n_pct(dat$MetS_case, n)
    )
  )
}

# Generate Table S1 and SAS/SDS group counts for one cohort.
make_table_s1 <- function(dat, cohort_name) {
  dat <- prepare_baseline_data(dat)
  groups <- list(
    All = dat,
    Anxiety_Healthy = dat %>% filter(Anxiety == 0),
    Anxiety_Case = dat %>% filter(Anxiety == 1),
    Depression_Healthy = dat %>% filter(Depression == 0),
    Depression_Case = dat %>% filter(Depression == 1)
  )
  out <- Reduce(
    function(x, y) left_join(x, y, by = "Characteristics"),
    lapply(names(groups), function(g) {
      summarise_baseline(groups[[g]]) %>% rename(!!g := Value)
    })
  )
  counts <- tibble(
    Cohort = cohort_name,
    Group = c("Anxiety Healthy", "Anxiety Case", "Depression Healthy", "Depression Case"),
    N = c(sum(dat$Anxiety == 0, na.rm = TRUE), sum(dat$Anxiety == 1, na.rm = TRUE), sum(dat$Depression == 0, na.rm = TRUE), sum(dat$Depression == 1, na.rm = TRUE)),
    Percent = round(100 * N / nrow(dat), 1)
  )
  list(table = out, counts = counts)
}

all_tables <- list()
all_counts <- list()

for (cohort_name in names(cohort_files)) {
  dat <- read_csv(cohort_files[[cohort_name]], show_col_types = FALSE)
  res <- make_table_s1(dat, cohort_name)
  write_csv(res$table, file.path(out_dir, paste0("TableS1_", cohort_name, "_baseline_characteristics.csv")))
  write_csv(res$counts, file.path(out_dir, paste0("TableS1_", cohort_name, "_SAS_SDS_group_counts.csv")))
  all_tables[[cohort_name]] <- res$table %>% mutate(Cohort = cohort_name, .before = 1)
  all_counts[[cohort_name]] <- res$counts
}

write_csv(bind_rows(all_tables), file.path(out_dir, "TableS1_baseline_characteristics_all_cohorts.csv"))
write_csv(bind_rows(all_counts), file.path(out_dir, "TableS1_SAS_SDS_binary_group_counts_all_cohorts.csv"))
