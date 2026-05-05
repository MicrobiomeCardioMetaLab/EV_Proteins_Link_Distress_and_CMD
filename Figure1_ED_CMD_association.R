# Figure 1: Association between emotional distress and cardiometabolic disorders.
# This script fits multivariable logistic regression models for anxiety/depression status
# against each cardiometabolic disorder in discovery, replication, and pooled cohorts.
# Outputs include odds ratios, 95% confidence intervals, P values, and FDR-adjusted results.

rm(list = ls())
library(dplyr)
library(readr)

project_dir <- "/path/to/project"
data_dir <- file.path(project_dir, "data")
out_dir <- file.path(project_dir, "results", "Figure1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

score_files <- list(
  Discovery = file.path(data_dir, "discovery_score_covariates.csv"),
  Replication = file.path(data_dir, "replication_score_covariates.csv")
)

cmd_files <- list(
  Discovery = file.path(data_dir, "discovery_cardiometabolic_disease.csv"),
  Replication = file.path(data_dir, "replication_cardiometabolic_disease.csv")
)

cmd_vars <- c("T2D", "MetS", "Dyslipidemia", "Hypertension")
covars <- c("Age_scaled", "Gender", "BMI_scaled", "Smoking_status_now", "Alcohol_intake_status", "Education_level")
case_code <- 2
control_code <- 1

# Merge score and CMD data and prepare model-ready variables.
prepare_model_data <- function(score_df, cmd_df) {
  cmd_keep <- intersect(cmd_vars, names(cmd_df))
  score_df %>%
    inner_join(cmd_df %>% select(sample_ID, all_of(cmd_keep)), by = "sample_ID") %>%
    mutate(
      Anxiety = case_when(SAS_score >= 50 ~ 1, SAS_score < 50 ~ 0, TRUE ~ NA_real_),
      Depression = case_when(SDS_score >= 50 ~ 1, SDS_score < 50 ~ 0, TRUE ~ NA_real_),
      Age_scaled = as.numeric(scale(Age)),
      BMI_scaled = as.numeric(scale(BMI)),
      Gender = factor(Gender),
      Smoking_status_now = factor(Smoking_status_now),
      Alcohol_intake_status = factor(Alcohol_intake_status),
      Education_level = factor(Education_level),
      across(all_of(cmd_keep), ~ case_when(.x == case_code ~ 1, .x == control_code ~ 0, TRUE ~ NA_real_))
    )
}

# Fit one adjusted logistic regression model for one ED-CMD pair.
run_one_logistic <- function(dat, ed_var, cmd_var, cohort_name) {
  use_vars <- c(cmd_var, ed_var, covars)
  df <- dat[, use_vars, drop = FALSE]
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 50 || length(unique(df[[cmd_var]])) < 2 || length(unique(df[[ed_var]])) < 2) {
    return(tibble())
  }
  names(df)[names(df) == cmd_var] <- "Y"
  names(df)[names(df) == ed_var] <- "X"
  fit <- tryCatch(glm(Y ~ X + Age_scaled + Gender + BMI_scaled + Smoking_status_now + Alcohol_intake_status + Education_level, data = df, family = binomial), error = function(e) NULL)
  if (is.null(fit) || !("X" %in% names(coef(fit)))) return(tibble())
  beta <- unname(coef(fit)["X"])
  se <- summary(fit)$coefficients["X", "Std. Error"]
  tibble(
    Cohort = cohort_name,
    ED = ed_var,
    CMD = cmd_var,
    N = nrow(df),
    Cases = sum(df$Y == 1, na.rm = TRUE),
    Controls = sum(df$Y == 0, na.rm = TRUE),
    Beta = beta,
    SE = se,
    OR = exp(beta),
    CI_lower = exp(beta - 1.96 * se),
    CI_upper = exp(beta + 1.96 * se),
    P_value = summary(fit)$coefficients["X", "Pr(>|z|)"]
  )
}

# Run all ED-CMD logistic regression models within one cohort.
run_logistic_set <- function(dat, cohort_name) {
  bind_rows(lapply(c("Anxiety", "Depression"), function(ed) {
    bind_rows(lapply(intersect(cmd_vars, names(dat)), function(cmd) run_one_logistic(dat, ed, cmd, cohort_name)))
  }))
}

cohort_data <- list()
for (cohort_name in names(score_files)) {
  score_df <- read_csv(score_files[[cohort_name]], show_col_types = FALSE)
  cmd_df <- read_csv(cmd_files[[cohort_name]], show_col_types = FALSE)
  cohort_data[[cohort_name]] <- prepare_model_data(score_df, cmd_df)
}

overall_data <- bind_rows(cohort_data, .id = "Source_cohort")
all_results <- bind_rows(
  bind_rows(lapply(names(cohort_data), function(cohort_name) run_logistic_set(cohort_data[[cohort_name]], cohort_name))),
  run_logistic_set(overall_data, "Overall")
)

all_results <- all_results %>%
  group_by(Cohort, ED) %>%
  mutate(FDR_by_ED = p.adjust(P_value, method = "BH")) %>%
  ungroup() %>%
  mutate(FDR_global = p.adjust(P_value, method = "BH")) %>%
  arrange(ED, CMD, factor(Cohort, levels = c("Discovery", "Replication", "Overall")))

write_csv(all_results, file.path(out_dir, "Figure1_ED_CMD_logistic_results.csv"))
write_csv(all_results %>% filter(Cohort == "Discovery", FDR_by_ED < 0.05), file.path(out_dir, "Figure1_ED_CMD_logistic_results_discovery_FDR005.csv"))
