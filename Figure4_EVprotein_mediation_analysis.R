# Figure 4: Mediation analysis of EV proteins linking emotional distress to cardiometabolic disorders.
# This script constructs ED-EV protein-CMD candidate triplets and performs covariate-adjusted
# mediation analysis in discovery and replication cohorts.
# Outputs include ACME, ADE, total effects, mediation proportions, FDR-adjusted mediation results,
# and replicated mediator triplets.

rm(list = ls())
library(dplyr)
library(readr)
library(mediation)
set.seed(123)
project_dir <- "/path/to/project"
data_dir <- file.path(project_dir, "data")
out_dir <- file.path(project_dir, "results", "Figure4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

xy_file <- file.path(project_dir, "results", "Figure1", "Figure1_ED_CMD_logistic_results.csv")
xm_file <- file.path(project_dir, "results", "Figure2", "Figure2_EVprotein_ED_discovery_FDR005_merged.csv")
my_file <- file.path(project_dir, "results", "Figure2", "Figure2_EVprotein_CMD_discovery_FDR005_merged.csv")

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

covariates <- c("Age", "Gender", "BMI", "Smoking_status_now", "Alcohol_intake_status", "Education_level")
case_code <- 2
control_code <- 1
sims_n <- 999

# Standardize ED-protein association results for candidate triplet construction.
standardize_xm <- function(x) {
  trait_col <- intersect(c("Trait", "ED", "group"), names(x))[1]
  protein_col <- intersect(c("protein", "Protein"), names(x))[1]
  rho_col <- intersect(c("spearman_rho_discovery", "spearman_rho", "rho_discovery"), names(x))[1]
  p_col <- intersect(c("p_spearman_discovery", "p_spearman", "spearman_p_discovery"), names(x))[1]
  fdr_col <- intersect(c("FDR_by_trait_discovery", "FDR_by_trait", "adj_p_by_protein_disc"), names(x))[1]
  x %>%
    transmute(X = .data[[trait_col]], M = .data[[protein_col]], rho_xm = .data[[rho_col]], p_xm = .data[[p_col]], fdr_xm = .data[[fdr_col]]) %>%
    mutate(X = case_when(X == "SAS" ~ "Anxiety", X == "SDS" ~ "Depression", TRUE ~ X)) %>%
    distinct()
}

# Standardize protein-CMD association results for candidate triplet construction.
standardize_my <- function(x) {
  trait_col <- intersect(c("Trait", "CMD", "CMDs"), names(x))[1]
  protein_col <- intersect(c("protein", "Protein"), names(x))[1]
  rho_col <- intersect(c("spearman_rho_discovery", "spearman_rho", "rho_discovery"), names(x))[1]
  p_col <- intersect(c("p_spearman_discovery", "p_spearman", "spearman_p_discovery"), names(x))[1]
  fdr_col <- intersect(c("FDR_by_trait_discovery", "FDR_by_trait", "adj_p_by_protein_disc"), names(x))[1]
  x %>%
    transmute(Y = .data[[trait_col]], M = .data[[protein_col]], rho_my = .data[[rho_col]], p_my = .data[[p_col]], fdr_my = .data[[fdr_col]]) %>%
    distinct()
}

# Build ED-EV protein-CMD triplets supported by discovery associations.
build_candidate_triplets <- function() {
  xy <- read_csv(xy_file, show_col_types = FALSE)
  xm <- read_csv(xm_file, show_col_types = FALSE) %>% standardize_xm()
  my <- read_csv(my_file, show_col_types = FALSE) %>% standardize_my()
  xy_sig <- xy %>%
    filter(Cohort == "Discovery", FDR_by_ED < 0.05) %>%
    transmute(X = ED, Y = CMD, OR_xy = OR, P_xy = P_value, FDR_xy = FDR_by_ED)
  xm_sig <- xm %>% filter(fdr_xm < 0.05)
  my_sig <- my %>% filter(fdr_my < 0.05)
  xy_sig %>%
    inner_join(xm_sig, by = "X") %>%
    inner_join(my_sig, by = c("Y", "M")) %>%
    distinct(X, M, Y, .keep_all = TRUE)
}

# Merge phenotype, covariate, CMD, and EV protein data for mediation models.
prepare_mediation_data <- function(score_df, cmd_df, protein_df, triplets) {
  x_need <- unique(triplets$X)
  y_need <- unique(triplets$Y)
  m_need <- intersect(unique(triplets$M), names(protein_df))
  score_df %>%
    mutate(
      Anxiety = case_when(SAS_score >= 50 ~ 1, SAS_score < 50 ~ 0, TRUE ~ NA_real_),
      Depression = case_when(SDS_score >= 50 ~ 1, SDS_score < 50 ~ 0, TRUE ~ NA_real_),
      Gender = factor(Gender),
      Smoking_status_now = factor(Smoking_status_now),
      Alcohol_intake_status = factor(Alcohol_intake_status),
      Education_level = factor(Education_level)
    ) %>%
    select(sample_ID, all_of(covariates), any_of(x_need)) %>%
    inner_join(cmd_df %>% mutate(across(any_of(y_need), ~ case_when(.x == case_code ~ 1, .x == control_code ~ 0, TRUE ~ NA_real_))) %>% select(sample_ID, any_of(y_need)), by = "sample_ID") %>%
    inner_join(protein_df %>% select(sample_ID, any_of(m_need)), by = "sample_ID")
}

# Run mediation analysis and return NULL if model fitting fails.
safe_mediate_summary <- function(model_m, model_y, treat, mediator, control.value = NULL, treat.value = NULL) {
  tryCatch({
    fit <- if (is.null(control.value) || is.null(treat.value)) {
      mediate(model_m, model_y, treat = treat, mediator = mediator, boot = TRUE, sims = sims_n)
    } else {
      mediate(model_m, model_y, treat = treat, mediator = mediator, boot = TRUE, sims = sims_n, control.value = control.value, treat.value = treat.value)
    }
    summary(fit)
  }, error = function(e) NULL)
}

# Safely extract one mediation summary statistic.
safe_extract <- function(x, nm) {
  if (is.null(x) || !(nm %in% names(x))) return(NA_real_)
  unname(x[[nm]])
}

# Create a placeholder row when a mediation model cannot be fitted.
empty_mediation_row <- function(x, m, y, n, note) {
  tibble(
    X = x, M = m, Y = y, n = n,
    acme_forward = NA_real_, acme_p_forward = NA_real_, ade_forward = NA_real_, ade_p_forward = NA_real_, total_forward = NA_real_, total_p_forward = NA_real_, prop_med_forward = NA_real_, prop_med_p_forward = NA_real_,
    acme_inverse = NA_real_, acme_p_inverse = NA_real_, ade_inverse = NA_real_, ade_p_inverse = NA_real_, total_inverse = NA_real_, total_p_inverse = NA_real_, prop_med_inverse = NA_real_, prop_med_p_inverse = NA_real_,
    note = note
  )
}

# Fit one mediator model and one outcome model for a single triplet.
run_one_mediation <- function(all_df, x, m, y) {
  cols_needed <- unique(c(x, m, y, covariates))
  cols_needed <- cols_needed[cols_needed %in% names(all_df)]
  tmp <- all_df %>% select(all_of(cols_needed))
  names(tmp)[names(tmp) == x] <- "X"
  names(tmp)[names(tmp) == m] <- "M"
  names(tmp)[names(tmp) == y] <- "Y"
  tmp <- tmp %>% filter(complete.cases(.))
  if (nrow(tmp) < 50) return(empty_mediation_row(x, m, y, nrow(tmp), "skip: n < 50"))
  if (length(unique(tmp$X)) < 2) return(empty_mediation_row(x, m, y, nrow(tmp), "skip: X has < 2 levels"))
  if (length(unique(tmp$Y)) < 2) return(empty_mediation_row(x, m, y, nrow(tmp), "skip: Y has < 2 levels"))
  tmp$M <- as.numeric(scale(tmp$M))
  cov_in_df <- covariates[covariates %in% names(tmp)]
  cov_formula <- if (length(cov_in_df) > 0) paste("+", paste(cov_in_df, collapse = " + ")) else ""
  model_m_forward <- tryCatch(lm(as.formula(paste("M ~ X", cov_formula)), data = tmp), error = function(e) NULL)
  model_y_forward <- tryCatch(glm(as.formula(paste("Y ~ X + M", cov_formula)), data = tmp, family = binomial), error = function(e) NULL)
  sum_forward <- if (!is.null(model_m_forward) && !is.null(model_y_forward)) safe_mediate_summary(model_m_forward, model_y_forward, "X", "M") else NULL
  model_m_inverse <- tryCatch(glm(as.formula(paste("X ~ M", cov_formula)), data = tmp, family = binomial), error = function(e) NULL)
  model_y_inverse <- tryCatch(glm(as.formula(paste("Y ~ M + X", cov_formula)), data = tmp, family = binomial), error = function(e) NULL)
  sum_inverse <- if (!is.null(model_m_inverse) && !is.null(model_y_inverse)) safe_mediate_summary(model_m_inverse, model_y_inverse, "M", "X", control.value = 0, treat.value = 1) else NULL
  tibble(
    X = x, M = m, Y = y, n = nrow(tmp),
    acme_forward = safe_extract(sum_forward, "d.avg"), acme_p_forward = safe_extract(sum_forward, "d.avg.p"),
    ade_forward = safe_extract(sum_forward, "z.avg"), ade_p_forward = safe_extract(sum_forward, "z.avg.p"),
    total_forward = safe_extract(sum_forward, "tau.coef"), total_p_forward = safe_extract(sum_forward, "tau.p"),
    prop_med_forward = safe_extract(sum_forward, "n.avg"), prop_med_p_forward = safe_extract(sum_forward, "n.avg.p"),
    acme_inverse = safe_extract(sum_inverse, "d.avg"), acme_p_inverse = safe_extract(sum_inverse, "d.avg.p"),
    ade_inverse = safe_extract(sum_inverse, "z.avg"), ade_p_inverse = safe_extract(sum_inverse, "z.avg.p"),
    total_inverse = safe_extract(sum_inverse, "tau.coef"), total_p_inverse = safe_extract(sum_inverse, "tau.p"),
    prop_med_inverse = safe_extract(sum_inverse, "n.avg"), prop_med_p_inverse = safe_extract(sum_inverse, "n.avg.p"),
    note = ifelse(is.null(sum_forward) | is.null(sum_inverse), "mediate error in forward or inverse model", "")
  )
}

# Run mediation analyses for all candidate triplets in one cohort.
run_mediation_for_cohort <- function(cohort_name, triplets) {
  all_df <- prepare_mediation_data(
    read_csv(score_files[[cohort_name]], show_col_types = FALSE),
    read_csv(cmd_files[[cohort_name]], show_col_types = FALSE),
    read_csv(ev_protein_files[[cohort_name]], show_col_types = FALSE),
    triplets
  )
  bind_rows(lapply(seq_len(nrow(triplets)), function(i) run_one_mediation(all_df, triplets$X[i], triplets$M[i], triplets$Y[i]))) %>%
    mutate(Cohort = cohort_name, .before = 1) %>%
    mutate(FDR_acme_forward = p.adjust(acme_p_forward, method = "BH"), FDR_acme_inverse = p.adjust(acme_p_inverse, method = "BH")) %>%
    left_join(triplets, by = c("X", "M", "Y"))
}

triplets <- build_candidate_triplets()
write_csv(triplets, file.path(out_dir, "Figure4_candidate_X_M_Y_triplets.csv"))

mediation_results <- bind_rows(lapply(names(score_files), function(cohort_name) run_mediation_for_cohort(cohort_name, triplets)))
write_csv(mediation_results, file.path(out_dir, "Figure4_EVprotein_mediation_results.csv"))

mediation_merged <- mediation_results %>%
  select(Cohort, X, M, Y, n, acme_forward, acme_p_forward, FDR_acme_forward, acme_inverse, acme_p_inverse, FDR_acme_inverse, note) %>%
  pivot_wider(names_from = Cohort, values_from = c(n, acme_forward, acme_p_forward, FDR_acme_forward, acme_inverse, acme_p_inverse, FDR_acme_inverse, note), names_sep = "_")
write_csv(mediation_merged, file.path(out_dir, "Figure4_EVprotein_mediation_discovery_replication_merged.csv"))

mediator_proteins <- mediation_results %>%
  filter(Cohort == "Discovery", FDR_acme_forward < 0.05, is.na(acme_p_inverse) | acme_p_inverse > 0.05) %>%
  distinct(M) %>%
  rename(protein = M)
write_csv(mediator_proteins, file.path(out_dir, "Figure4_significant_forward_mediator_EVproteins.csv"))
