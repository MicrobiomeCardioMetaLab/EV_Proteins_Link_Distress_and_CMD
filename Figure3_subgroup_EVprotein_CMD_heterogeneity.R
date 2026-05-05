# Figure 3: Subgroup heterogeneity of CMD-associated EV proteins by emotional distress status.
# This script evaluates whether associations between CMD-related EV proteins and CMD phenotypes
# differ between participants with and without anxiety/depression.
# Outputs include subgroup-specific correlations and heterogeneity test results.

rm(list = ls())
library(dplyr)
library(readr)
library(tidyr)
project_dir <- "/path/to/project"
data_dir <- file.path(project_dir, "data")
out_dir <- file.path(project_dir, "results", "Figure3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cmd_ev_file <- file.path(project_dir, "results", "Figure2", "Figure2_EVprotein_CMD_discovery_FDR005_merged.csv")

score_files <- list(
  Discovery = file.path(data_dir, "discovery_score_covariates.csv"),
  Replication = file.path(data_dir, "replication_score_covariates.csv")
)

cmd_files <- list(
  Discovery = file.path(data_dir, "discovery_cardiometabolic_disease.csv"),
  Replication = file.path(data_dir, "replication_cardiometabolic_disease.csv")
)

ev_protein_files <- list(
  Discovery = file.path(data_dir, "discovery_EV_protein_log_RSN.csv"),
  Replication = file.path(data_dir, "replication_EV_protein_log_RSN.csv")
)

covars <- c("Age", "Gender", "BMI", "Smoking_status_now", "Alcohol_intake_status", "Education_level")
cmd_vars <- c("T2D", "MetS", "Dyslipidemia", "Hypertension")
case_code <- 2
control_code <- 1
min_total_n <- 10
min_cases <- 5
min_controls <- 5

# Merge score, CMD, and EV protein data and prepare subgroup variables.
prepare_data <- function(score_df, cmd_df, protein_df) {
  cmd_keep <- intersect(cmd_vars, names(cmd_df))
  score_df %>%
    inner_join(cmd_df %>% select(sample_ID, all_of(cmd_keep)), by = "sample_ID") %>%
    inner_join(protein_df, by = "sample_ID") %>%
    mutate(
      Anxiety = case_when(SAS_score >= 50 ~ 1, SAS_score < 50 ~ 0, TRUE ~ NA_real_),
      Depression = case_when(SDS_score >= 50 ~ 1, SDS_score < 50 ~ 0, TRUE ~ NA_real_),
      across(all_of(cmd_keep), ~ case_when(.x == case_code ~ 1, .x == control_code ~ 0, TRUE ~ NA_real_)),
      Gender = factor(Gender),
      Smoking_status_now = factor(Smoking_status_now),
      Alcohol_intake_status = factor(Alcohol_intake_status),
      Education_level = factor(Education_level)
    )
}

# Remove covariate effects from one EV protein feature.
residualize_feature <- function(data, feature, covariates) {
  keep <- c(feature, covariates)
  df <- data[, keep, drop = FALSE]
  idx <- complete.cases(df)
  res <- rep(NA_real_, nrow(data))
  if (sum(idx) < min_total_n) return(res)
  model_df <- cbind(y = df[[feature]][idx], df[idx, covariates, drop = FALSE])
  fit <- tryCatch(lm(y ~ ., data = model_df), error = function(e) NULL)
  if (is.null(fit)) return(res)
  res[idx] <- resid(fit)
  res
}

# Estimate CMD-protein Spearman correlation within one ED subgroup.
run_stratified_spearman <- function(dat, outcome, protein, strat_var, strat_value) {
  if (!(outcome %in% names(dat)) || !(protein %in% names(dat)) || !(strat_var %in% names(dat))) return(tibble())
  df0 <- dat %>% filter(.data[[strat_var]] == strat_value)
  protein_resid <- residualize_feature(df0, protein, covars)
  df <- tibble(outcome_value = df0[[outcome]], protein_resid = protein_resid) %>% filter(complete.cases(.))
  cases <- sum(df$outcome_value == 1, na.rm = TRUE)
  controls <- sum(df$outcome_value == 0, na.rm = TRUE)
  if (nrow(df) < min_total_n || cases < min_cases || controls < min_controls || length(unique(df$outcome_value)) < 2) {
    return(tibble(outcome = outcome, protein = protein, strat_var = strat_var, strat_value = strat_value, n = nrow(df), cases = cases, controls = controls, rho = NA_real_, p = NA_real_, fisher_z = NA_real_))
  }
  ct <- tryCatch(cor.test(df$protein_resid, df$outcome_value, method = "spearman", exact = FALSE), error = function(e) NULL)
  if (is.null(ct)) {
    return(tibble(outcome = outcome, protein = protein, strat_var = strat_var, strat_value = strat_value, n = nrow(df), cases = cases, controls = controls, rho = NA_real_, p = NA_real_, fisher_z = NA_real_))
  }
  rho <- unname(ct$estimate)
  rho_clip <- pmin(pmax(rho, -0.999999), 0.999999)
  tibble(outcome = outcome, protein = protein, strat_var = strat_var, strat_value = strat_value, n = nrow(df), cases = cases, controls = controls, rho = rho, p = ct$p.value, fisher_z = atanh(rho_clip))
}

# Test heterogeneity between two independent subgroup correlations.
compare_independent_correlations <- function(r_yes, n_yes, r_no, n_no) {
  if (any(is.na(c(r_yes, n_yes, r_no, n_no))) || n_yes <= 3 || n_no <= 3) {
    return(tibble(z_diff = NA_real_, Q = NA_real_, df = 1, p_heterogeneity = NA_real_, I2 = NA_real_))
  }
  r_yes <- pmin(pmax(r_yes, -0.999999), 0.999999)
  r_no <- pmin(pmax(r_no, -0.999999), 0.999999)
  z_diff <- (atanh(r_yes) - atanh(r_no)) / sqrt(1 / (n_yes - 3) + 1 / (n_no - 3))
  Q <- z_diff^2
  tibble(z_diff = z_diff, Q = Q, df = 1, p_heterogeneity = 2 * pnorm(abs(z_diff), lower.tail = FALSE), I2 = ifelse(Q > 1, 100 * (Q - 1) / Q, 0))
}

# Standardize input CMD-protein pair tables for downstream analysis.
standardize_cmd_protein_table <- function(x) {
  outcome_col <- intersect(c("Trait", "CMD", "CMDs", "outcome"), names(x))[1]
  protein_col <- intersect(c("protein", "Protein"), names(x))[1]
  x %>% transmute(outcome = .data[[outcome_col]], protein = .data[[protein_col]]) %>% distinct()
}

# Run subgroup heterogeneity analysis for all CMD-protein pairs in one cohort.
run_subgroup_analysis_one_cohort <- function(cohort_name, cmd_protein_pairs) {
  dat <- prepare_data(
    read_csv(score_files[[cohort_name]], show_col_types = FALSE),
    read_csv(cmd_files[[cohort_name]], show_col_types = FALSE),
    read_csv(ev_protein_files[[cohort_name]], show_col_types = FALSE)
  )
  stratum_results <- bind_rows(lapply(seq_len(nrow(cmd_protein_pairs)), function(i) {
    outcome_i <- cmd_protein_pairs$outcome[i]
    protein_i <- cmd_protein_pairs$protein[i]
    bind_rows(
      run_stratified_spearman(dat, outcome_i, protein_i, "Anxiety", 1),
      run_stratified_spearman(dat, outcome_i, protein_i, "Anxiety", 0),
      run_stratified_spearman(dat, outcome_i, protein_i, "Depression", 1),
      run_stratified_spearman(dat, outcome_i, protein_i, "Depression", 0)
    )
  })) %>% mutate(Cohort = cohort_name, .before = 1)
  heterogeneity_results <- stratum_results %>%
    mutate(stratum = ifelse(strat_value == 1, "case", "control")) %>%
    select(Cohort, outcome, protein, strat_var, stratum, n, cases, controls, rho, p, fisher_z) %>%
    pivot_wider(names_from = stratum, values_from = c(n, cases, controls, rho, p, fisher_z), values_fn = ~ .x[1]) %>%
    rowwise() %>%
    mutate(tmp = list(compare_independent_correlations(rho_case, n_case, rho_control, n_control))) %>%
    unnest_wider(tmp) %>%
    ungroup() %>%
    group_by(Cohort, strat_var) %>%
    mutate(FDR_heterogeneity = p.adjust(p_heterogeneity, method = "BH")) %>%
    ungroup()
  list(stratum = stratum_results, heterogeneity = heterogeneity_results)
}

cmd_protein_pairs <- read_csv(cmd_ev_file, show_col_types = FALSE) %>% standardize_cmd_protein_table()

all_stratum <- list()
all_heterogeneity <- list()

for (cohort_name in names(score_files)) {
  res <- run_subgroup_analysis_one_cohort(cohort_name, cmd_protein_pairs)
  all_stratum[[cohort_name]] <- res$stratum
  all_heterogeneity[[cohort_name]] <- res$heterogeneity
}

stratum_df <- bind_rows(all_stratum)
heterogeneity_df <- bind_rows(all_heterogeneity)

write_csv(stratum_df, file.path(out_dir, "Figure3_CMD_EVprotein_psych_stratified_residSpearman_results.csv"))
write_csv(heterogeneity_df, file.path(out_dir, "Figure3_CMD_EVprotein_psych_stratified_heterogeneity_results.csv"))
write_csv(heterogeneity_df %>% filter(Cohort == "Discovery", FDR_heterogeneity < 0.05), file.path(out_dir, "Figure3_discovery_stratified_heterogeneity_FDR005.csv"))
