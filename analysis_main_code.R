#### Main analysis codes for the role of EV proteins in emotional distress and cardiometabolic disorders study
#### For questions, please contact Shuai Yang (shuai_yang8@163.com)




#----- 1. Association between emotional distress and cardiometabolic disorders -----
# This script fits multivariable logistic regression models for anxiety/depression status
# against each cardiometabolic disorder in discovery and replication cohorts.
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

write_csv(all_results, file.path(out_dir, "ED_CMD_logistic_results.csv"))
write_csv(all_results %>% filter(Cohort == "Discovery", FDR_by_ED < 0.05), file.path(out_dir, "ED_CMD_logistic_results_discovery_FDR005.csv"))




#----- 2. EV protein associations with emotional distress and cardiometabolic disorders -----
# This script tests covariate-adjusted EV protein associations with anxiety, depression,
# and cardiometabolic disorders, merges discovery and replication results, annotates proteins,
# and compares disease associations between EV proteins and corresponding plasma proteins.
# Outputs include Wilcoxon/Spearman association tables, validated results, enrichment-ready protein sets,
# and EV-versus-plasma effect heterogeneity results.

rm(list = ls())
library(dplyr)
library(readr)
library(tidyr)
library(purrr)

project_dir <- "/path/to/project"
data_dir <- file.path(project_dir, "data")
out_dir <- file.path(project_dir, "results", "Figure2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

plasma_protein_files <- list(
  Discovery = file.path(data_dir, "discovery_plasma_protein_log_RSN.csv"),
  Replication = file.path(data_dir, "replication_plasma_protein_log_RSN.csv")
)

annotation_file <- file.path(data_dir, "protein_annotation.csv")

cmd_vars <- c("T2D", "MetS", "Dyslipidemia", "Hypertension")
covars <- c("Age", "Gender", "BMI", "Smoking_status_now", "Alcohol_intake_status", "Education_level")
case_code <- 2
control_code <- 1

# Convert covariates to analysis-ready factor variables.
prepare_covariates <- function(dat) {
  dat %>%
    mutate(
      Gender = factor(Gender),
      Smoking_status_now = factor(Smoking_status_now),
      Alcohol_intake_status = factor(Alcohol_intake_status),
      Education_level = factor(Education_level)
    )
}

# Remove covariate effects from one protein feature using linear regression.
residualize_one_feature <- function(data, feature, covariates) {
  keep <- c(feature, covariates)
  df <- data[, keep, drop = FALSE]
  idx <- complete.cases(df)
  res <- rep(NA_real_, nrow(data))
  if (sum(idx) < 10) return(res)
  model_df <- cbind(y = df[[feature]][idx], df[idx, covariates, drop = FALSE])
  fit <- tryCatch(lm(y ~ ., data = model_df), error = function(e) NULL)
  if (is.null(fit)) return(res)
  res[idx] <- resid(fit)
  res
}

# Test one binary trait against all residualized protein features.
run_binary_feature_association <- function(data, trait, trait_label, feature_names, covariates) {
  bind_rows(lapply(feature_names, function(feature) {
    feature_resid <- residualize_one_feature(data, feature, covariates)
    df <- tibble(trait = data[[trait]], feature_resid = feature_resid) %>% filter(!is.na(trait), !is.na(feature_resid))
    if (nrow(df) < 10 || length(unique(df$trait)) < 2) return(tibble())
    n_control <- sum(df$trait == 0)
    n_case <- sum(df$trait == 1)
    if (n_control < 3 || n_case < 3) return(tibble())
    wt <- tryCatch(wilcox.test(feature_resid ~ trait, data = df, exact = FALSE), error = function(e) NULL)
    ct <- tryCatch(cor.test(df$trait, df$feature_resid, method = "spearman", exact = FALSE), error = function(e) NULL)
    if (is.null(wt) || is.null(ct)) return(tibble())
    tibble(
      Trait = trait_label,
      protein = feature,
      n_control = n_control,
      n_case = n_case,
      median_control = median(df$feature_resid[df$trait == 0], na.rm = TRUE),
      median_case = median(df$feature_resid[df$trait == 1], na.rm = TRUE),
      mean_control = mean(df$feature_resid[df$trait == 0], na.rm = TRUE),
      mean_case = mean(df$feature_resid[df$trait == 1], na.rm = TRUE),
      diff_median_case_minus_control = median_case - median_control,
      diff_mean_case_minus_control = mean_case - mean_control,
      wilcox_statistic = unname(wt$statistic),
      p_wilcox = wt$p.value,
      spearman_rho = unname(ct$estimate),
      p_spearman = ct$p.value,
      sample_size = nrow(df),
      direction = case_when(spearman_rho > 0 ~ "Higher_in_case", spearman_rho < 0 ~ "Lower_in_case", TRUE ~ "Equal")
    )
  }))
}

# Run EV/plasma protein association analyses for anxiety and depression.
run_ed_protein_associations <- function(score_df, protein_df, cohort_name, protein_type) {
  protein_names <- setdiff(names(protein_df), "sample_ID")
  dat <- score_df %>%
    inner_join(protein_df, by = "sample_ID") %>%
    prepare_covariates() %>%
    mutate(
      Anxiety = case_when(SAS_score >= 50 ~ 1, SAS_score < 50 ~ 0, TRUE ~ NA_real_),
      Depression = case_when(SDS_score >= 50 ~ 1, SDS_score < 50 ~ 0, TRUE ~ NA_real_)
    )
  bind_rows(
    run_binary_feature_association(dat, "Anxiety", "Anxiety", protein_names, covars),
    run_binary_feature_association(dat, "Depression", "Depression", protein_names, covars)
  ) %>%
    mutate(Cohort = cohort_name, Protein_type = protein_type, .before = 1) %>%
    group_by(Cohort, Protein_type, Trait) %>%
    mutate(FDR_by_trait = p.adjust(p_spearman, method = "BH")) %>%
    ungroup() %>%
    mutate(FDR_global = p.adjust(p_spearman, method = "BH"))
}

# Run EV/plasma protein association analyses for CMD phenotypes.
run_cmd_protein_associations <- function(score_df, cmd_df, protein_df, cohort_name, protein_type) {
  protein_names <- setdiff(names(protein_df), "sample_ID")
  cmd_keep <- intersect(cmd_vars, names(cmd_df))
  dat <- score_df %>%
    inner_join(cmd_df %>% select(sample_ID, all_of(cmd_keep)), by = "sample_ID") %>%
    inner_join(protein_df, by = "sample_ID") %>%
    prepare_covariates() %>%
    mutate(across(all_of(cmd_keep), ~ case_when(.x == case_code ~ 1, .x == control_code ~ 0, TRUE ~ NA_real_)))
  bind_rows(lapply(cmd_keep, function(cmd) run_binary_feature_association(dat, cmd, cmd, protein_names, covars))) %>%
    mutate(Cohort = cohort_name, Protein_type = protein_type, .before = 1) %>%
    group_by(Cohort, Protein_type, Trait) %>%
    mutate(FDR_by_trait = p.adjust(p_spearman, method = "BH")) %>%
    ungroup() %>%
    mutate(FDR_global = p.adjust(p_spearman, method = "BH"))
}

# Merge discovery and replication association results by key variables.
merge_discovery_replication <- function(res, key_cols) {
  disc <- res %>% filter(Cohort == "Discovery") %>% rename_with(~ paste0(.x, "_discovery"), .cols = -all_of(key_cols))
  rep <- res %>% filter(Cohort == "Replication") %>% rename_with(~ paste0(.x, "_replication"), .cols = -all_of(key_cols))
  inner_join(disc, rep, by = key_cols) %>%
    mutate(
      replicated_nominal = sign(spearman_rho_discovery) == sign(spearman_rho_replication) & p_spearman_replication < 0.05,
      replicated_wilcox = sign(spearman_rho_discovery) == sign(spearman_rho_replication) & p_wilcox_replication < 0.05
    )
}

# Add protein annotation when an annotation table is available.
add_annotation <- function(res, annotation_file) {
  if (!file.exists(annotation_file)) return(res)
  anno <- read_csv(annotation_file, show_col_types = FALSE)
  if ("GENE_name" %in% names(anno)) {
    anno <- anno %>% distinct(GENE_name, .keep_all = TRUE)
    return(left_join(res, anno, by = c("protein" = "GENE_name")))
  }
  if ("protein" %in% names(anno)) {
    anno <- anno %>% distinct(protein, .keep_all = TRUE)
    return(left_join(res, anno, by = "protein"))
  }
  res
}

ev_ed_results <- bind_rows(lapply(names(score_files), function(cohort_name) {
  run_ed_protein_associations(read_csv(score_files[[cohort_name]], show_col_types = FALSE), read_csv(ev_protein_files[[cohort_name]], show_col_types = FALSE), cohort_name, "EV")
}))

ev_cmd_results <- bind_rows(lapply(names(score_files), function(cohort_name) {
  run_cmd_protein_associations(read_csv(score_files[[cohort_name]], show_col_types = FALSE), read_csv(cmd_files[[cohort_name]], show_col_types = FALSE), read_csv(ev_protein_files[[cohort_name]], show_col_types = FALSE), cohort_name, "EV")
}))

plasma_ed_results <- bind_rows(lapply(names(score_files), function(cohort_name) {
  run_ed_protein_associations(read_csv(score_files[[cohort_name]], show_col_types = FALSE), read_csv(plasma_protein_files[[cohort_name]], show_col_types = FALSE), cohort_name, "Plasma")
}))

plasma_cmd_results <- bind_rows(lapply(names(score_files), function(cohort_name) {
  run_cmd_protein_associations(read_csv(score_files[[cohort_name]], show_col_types = FALSE), read_csv(cmd_files[[cohort_name]], show_col_types = FALSE), read_csv(plasma_protein_files[[cohort_name]], show_col_types = FALSE), cohort_name, "Plasma")
}))

ev_ed_results <- add_annotation(ev_ed_results, annotation_file)
ev_cmd_results <- add_annotation(ev_cmd_results, annotation_file)
plasma_ed_results <- add_annotation(plasma_ed_results, annotation_file)
plasma_cmd_results <- add_annotation(plasma_cmd_results, annotation_file)

write_csv(ev_ed_results, file.path(out_dir, "EVprotein_ED_associations.csv"))
write_csv(ev_cmd_results, file.path(out_dir, "EVprotein_CMD_associations.csv"))
write_csv(plasma_ed_results, file.path(out_dir, "plasma_protein_ED_associations.csv"))
write_csv(plasma_cmd_results, file.path(out_dir, "plasma_protein_CMD_associations.csv"))

ed_merge <- merge_discovery_replication(ev_ed_results, c("Protein_type", "Trait", "protein"))
cmd_merge <- merge_discovery_replication(ev_cmd_results, c("Protein_type", "Trait", "protein"))

write_csv(ed_merge, file.path(out_dir, "EVprotein_ED_discovery_replication_merged.csv"))
write_csv(cmd_merge, file.path(out_dir, "EVprotein_CMD_discovery_replication_merged.csv"))
write_csv(ed_merge %>% filter(FDR_by_trait_discovery < 0.05), file.path(out_dir, "EVprotein_ED_discovery_FDR005_merged.csv"))
write_csv(cmd_merge %>% filter(FDR_by_trait_discovery < 0.05), file.path(out_dir, "EVprotein_CMD_discovery_FDR005_merged.csv"))

shared_ev_proteins <- ed_merge %>%
  filter(FDR_by_trait_discovery < 0.05) %>%
  transmute(ED = Trait, protein) %>%
  inner_join(cmd_merge %>% filter(FDR_by_trait_discovery < 0.05) %>% transmute(CMD = Trait, protein), by = "protein") %>%
  distinct()
write_csv(shared_ev_proteins, file.path(out_dir, "ED_CMD_shared_significant_EVproteins.csv"))

if (requireNamespace("clusterProfiler", quietly = TRUE) && requireNamespace("ReactomePA", quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  gene_symbol <- unique(shared_ev_proteins$protein)
  gene_entrez <- clusterProfiler::bitr(gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  reactome_res <- ReactomePA::enrichPathway(gene = unique(gene_entrez$ENTREZID), organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.20, readable = TRUE)
  write_csv(as.data.frame(reactome_res), file.path(out_dir, "shared_EVproteins_Reactome_enrichment.csv"))
}

# Compute Spearman correlation with missing-value and sample-size checks.
safe_spearman <- function(x, y) {
  out <- tryCatch(suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)), error = function(e) NULL)
  if (is.null(out)) return(list(r = NA_real_, p = NA_real_))
  list(r = unname(out$estimate), p = out$p.value)
}

# Compare dependent overlapping correlations for EV and plasma counterparts.
run_cocor_overlap <- function(r_ev, r_plasma, r_ev_plasma, n) {
  out <- tibble(diff_rho_EV_minus_plasma = r_ev - r_plasma, p_heterogeneity = NA_real_, hittner_z = NA_real_, zou_ci_low = NA_real_, zou_ci_high = NA_real_)
  if (!requireNamespace("cocor", quietly = TRUE)) return(out)
  if (!all(is.finite(c(r_ev, r_plasma, r_ev_plasma, n))) || n < 4 || any(abs(c(r_ev, r_plasma, r_ev_plasma)) >= 1)) return(out)
  cr <- tryCatch(cocor::cocor.dep.groups.overlap(r.jk = r_ev, r.jh = r_plasma, r.kh = r_ev_plasma, n = n, alternative = "two.sided", test = c("hittner2003", "zou2007")), error = function(e) NULL)
  if (is.null(cr)) return(out)
  ht <- tryCatch(cocor::as.htest(cr), error = function(e) NULL)
  if (is.null(ht)) return(out)
  if (!is.null(ht$hittner2003)) {
    out$hittner_z <- unname(ht$hittner2003$statistic)
    out$p_heterogeneity <- ht$hittner2003$p.value
  }
  if (!is.null(ht$zou2007) && !is.null(ht$zou2007$conf.int)) {
    out$zou_ci_low <- unname(ht$zou2007$conf.int[1])
    out$zou_ci_high <- unname(ht$zou2007$conf.int[2])
  }
  out
}

# Compare trait associations between EV proteins and matched plasma proteins.
compare_ev_plasma_effects <- function(score_df, cmd_df, ev_df, plasma_df, traits, trait_type) {
  common_proteins <- intersect(setdiff(names(ev_df), "sample_ID"), setdiff(names(plasma_df), "sample_ID"))
  dat <- score_df %>%
    inner_join(cmd_df, by = "sample_ID") %>%
    inner_join(ev_df %>% select(sample_ID, all_of(common_proteins)), by = "sample_ID") %>%
    inner_join(plasma_df %>% select(sample_ID, all_of(common_proteins)), by = "sample_ID", suffix = c("_EV", "_plasma")) %>%
    prepare_covariates() %>%
    mutate(
      Anxiety = case_when(SAS_score >= 50 ~ 1, SAS_score < 50 ~ 0, TRUE ~ NA_real_),
      Depression = case_when(SDS_score >= 50 ~ 1, SDS_score < 50 ~ 0, TRUE ~ NA_real_),
      across(any_of(cmd_vars), ~ case_when(.x == case_code ~ 1, .x == control_code ~ 0, TRUE ~ NA_real_))
    )
  bind_rows(lapply(traits, function(trait) {
    bind_rows(lapply(common_proteins, function(protein) {
      ev_name <- paste0(protein, "_EV")
      plasma_name <- paste0(protein, "_plasma")
      if (!(ev_name %in% names(dat)) || !(plasma_name %in% names(dat))) return(tibble())
      ev_resid <- residualize_one_feature(dat, ev_name, covars)
      plasma_resid <- residualize_one_feature(dat, plasma_name, covars)
      df <- tibble(trait_value = dat[[trait]], EV = ev_resid, Plasma = plasma_resid) %>% filter(complete.cases(.))
      if (nrow(df) < 10 || length(unique(df$trait_value)) < 2) return(tibble())
      cor_ev <- safe_spearman(df$trait_value, df$EV)
      cor_plasma <- safe_spearman(df$trait_value, df$Plasma)
      cor_ep <- safe_spearman(df$EV, df$Plasma)
      het <- run_cocor_overlap(cor_ev$r, cor_plasma$r, cor_ep$r, nrow(df))
      tibble(Trait_type = trait_type, Trait = trait, protein = protein, n = nrow(df), rho_EV = cor_ev$r, p_EV = cor_ev$p, rho_plasma = cor_plasma$r, p_plasma = cor_plasma$p, rho_EV_plasma = cor_ep$r, p_EV_plasma = cor_ep$p) %>% bind_cols(het)
    }))
  })) %>% mutate(FDR_heterogeneity = p.adjust(p_heterogeneity, method = "BH"))
}

score_disc <- read_csv(score_files$Discovery, show_col_types = FALSE)
cmd_disc <- read_csv(cmd_files$Discovery, show_col_types = FALSE)
ev_disc <- read_csv(ev_protein_files$Discovery, show_col_types = FALSE)
plasma_disc <- read_csv(plasma_protein_files$Discovery, show_col_types = FALSE)

ed_heterogeneity <- compare_ev_plasma_effects(score_disc, cmd_disc, ev_disc, plasma_disc, c("Anxiety", "Depression"), "ED")
cmd_heterogeneity <- compare_ev_plasma_effects(score_disc, cmd_disc, ev_disc, plasma_disc, intersect(cmd_vars, names(cmd_disc)), "CMD")

write_csv(ed_heterogeneity, file.path(out_dir, "EV_vs_plasma_ED_effect_heterogeneity.csv"))
write_csv(cmd_heterogeneity, file.path(out_dir, "EV_vs_plasma_CMD_effect_heterogeneity.csv"))




#----- 3. Subgroup heterogeneity of CMD-associated EV proteins by emotional distress status -----
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

cmd_ev_file <- file.path(project_dir, "results", "Figure2", "EVprotein_CMD_discovery_FDR005_merged.csv")

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

write_csv(stratum_df, file.path(out_dir, "CMD_EVprotein_psych_stratified_residSpearman_results.csv"))
write_csv(heterogeneity_df, file.path(out_dir, "CMD_EVprotein_psych_stratified_heterogeneity_results.csv"))
write_csv(heterogeneity_df %>% filter(Cohort == "Discovery", FDR_heterogeneity < 0.05), file.path(out_dir, "discovery_stratified_heterogeneity_FDR005.csv"))




#----- 4. Mediation analysis of EV proteins linking emotional distress to cardiometabolic disorders -----
# This script constructs ED-EV protein-CMD candidate triplets and performs
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

xy_file <- file.path(project_dir, "results", "Figure1", "ED_CMD_logistic_results.csv")
xm_file <- file.path(project_dir, "results", "Figure2", "EVprotein_ED_discovery_FDR005_merged.csv")
my_file <- file.path(project_dir, "results", "Figure2", "EVprotein_CMD_discovery_FDR005_merged.csv")

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
write_csv(triplets, file.path(out_dir, "candidate_X_M_Y_triplets.csv"))

mediation_results <- bind_rows(lapply(names(score_files), function(cohort_name) run_mediation_for_cohort(cohort_name, triplets)))
write_csv(mediation_results, file.path(out_dir, "EVprotein_mediation_results.csv"))

mediation_merged <- mediation_results %>%
  select(Cohort, X, M, Y, n, acme_forward, acme_p_forward, FDR_acme_forward, acme_inverse, acme_p_inverse, FDR_acme_inverse, note) %>%
  pivot_wider(names_from = Cohort, values_from = c(n, acme_forward, acme_p_forward, FDR_acme_forward, acme_inverse, acme_p_inverse, FDR_acme_inverse, note), names_sep = "_")
write_csv(mediation_merged, file.path(out_dir, "EVprotein_mediation_discovery_replication_merged.csv"))

mediator_proteins <- mediation_results %>%
  filter(Cohort == "Discovery", FDR_acme_forward < 0.05, is.na(acme_p_inverse) | acme_p_inverse > 0.05) %>%
  distinct(M) %>%
  rename(protein = M)
write_csv(mediator_proteins, file.path(out_dir, "significant_forward_mediator_EVproteins.csv"))




#----- 5. Association between gut microbiome features and mediator EV proteins -----
# This script residualizes microbiome and mediator EV protein features against covariates,
# tests microbe-protein associations, and evaluates discovery-to-replication consistency.
# Outputs include species/pathway associations with mediator EV proteins and replicated pairs.

rm(list = ls())
library(dplyr)
library(readr)
library(tidyr)
library(purrr)

project_dir <- "/path/to/project"
data_dir <- file.path(project_dir, "data")
out_dir <- file.path(project_dir, "results", "Figure5")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mediator_file <- file.path(project_dir, "results", "Figure4", "significant_forward_mediator_EVproteins.csv")

covariate_files <- list(
  Discovery = file.path(data_dir, "discovery_score_covariates.csv"),
  Replication = file.path(data_dir, "replication_score_covariates.csv")
)

ev_protein_files <- list(
  Discovery = file.path(data_dir, "discovery_EV_protein_log_RSN.csv"),
  Replication = file.path(data_dir, "replication_EV_protein_log_RSN.csv")
)

species_files <- list(
  Discovery = file.path(data_dir, "discovery_species_relative_abundance_clr.csv"),
  Replication = file.path(data_dir, "replication_species_relative_abundance_clr.csv")
)

pathway_files <- list(
  Discovery = file.path(data_dir, "discovery_pathway_relative_abundance_clr.csv"),
  Replication = file.path(data_dir, "replication_pathway_relative_abundance_clr.csv")
)

species_annotation_file <- file.path(data_dir, "species_annotation.csv")
pathway_annotation_file <- file.path(data_dir, "pathway_annotation.csv")

covars <- c("Age", "Gender", "BMI", "Smoking_status_now", "Alcohol_intake_status", "Education_level")
min_n <- 10

# Convert covariates to analysis-ready factor variables.
prepare_covariates <- function(dat) {
  dat %>%
    mutate(
      Gender = factor(Gender),
      Smoking_status_now = factor(Smoking_status_now),
      Alcohol_intake_status = factor(Alcohol_intake_status),
      Education_level = factor(Education_level)
    )
}

# Remove covariate effects from one microbiome or protein feature.
residualize_one <- function(data, feature, covariates) {
  keep <- c(feature, covariates)
  df <- data[, keep, drop = FALSE]
  idx <- complete.cases(df)
  res <- rep(NA_real_, nrow(data))
  if (sum(idx) < min_n) return(res)
  model_df <- cbind(y = df[[feature]][idx], df[idx, covariates, drop = FALSE])
  fit <- tryCatch(lm(y ~ ., data = model_df), error = function(e) NULL)
  if (is.null(fit)) return(res)
  res[idx] <- resid(fit)
  res
}

# Residualize a set of features and export the residual matrix.
residualize_matrix <- function(feature_df, covar_df, feature_names, out_prefix) {
  dat <- covar_df %>% prepare_covariates() %>% inner_join(feature_df %>% select(sample_ID, all_of(feature_names)), by = "sample_ID")
  out <- tibble(sample_ID = dat$sample_ID)
  for (feature in feature_names) out[[feature]] <- residualize_one(dat, feature, covars)
  write_csv(out, file.path(out_dir, paste0(out_prefix, "_adjusted_residuals.csv")))
  out
}

# Test Spearman associations between microbiome features and mediator EV proteins.
run_microbe_protein_association <- function(microbe_resid, protein_resid, microbe_type, cohort_name) {
  microbe_names <- setdiff(names(microbe_resid), "sample_ID")
  protein_names <- setdiff(names(protein_resid), "sample_ID")
  assoc_data <- inner_join(microbe_resid, protein_resid, by = "sample_ID", suffix = c("_microbe", "_protein"))
  bind_rows(lapply(microbe_names, function(microbe) {
    bind_rows(lapply(protein_names, function(protein) {
      df <- assoc_data[, c(microbe, protein), drop = FALSE]
      df <- df[complete.cases(df), , drop = FALSE]
      if (nrow(df) < min_n || sd(df[[microbe]], na.rm = TRUE) == 0 || sd(df[[protein]], na.rm = TRUE) == 0) return(tibble())
      ct <- tryCatch(cor.test(df[[microbe]], df[[protein]], method = "spearman", exact = FALSE), error = function(e) NULL)
      if (is.null(ct)) return(tibble())
      tibble(Cohort = cohort_name, Microbe_type = microbe_type, microbe = microbe, protein = protein, sample_size = nrow(df), R = unname(ct$estimate), P = ct$p.value)
    }))
  })) %>%
    mutate(FDR_global = p.adjust(P, method = "BH")) %>%
    group_by(Cohort, Microbe_type, microbe) %>%
    mutate(FDR_by_microbe = p.adjust(P, method = "BH")) %>%
    ungroup() %>%
    group_by(Cohort, Microbe_type, protein) %>%
    mutate(FDR_by_protein = p.adjust(P, method = "BH")) %>%
    ungroup()
}

# Select mediator EV proteins from the mediation result table.
select_mediator_proteins <- function(protein_df) {
  mediator_tbl <- read_csv(mediator_file, show_col_types = FALSE)
  mediator_proteins <- intersect(mediator_tbl$protein, names(protein_df))
  protein_df %>% select(sample_ID, all_of(mediator_proteins))
}

# Run residualization and microbe-protein association analysis for one cohort.
run_one_cohort <- function(cohort_name, microbe_type, microbe_file) {
  covar_df <- read_csv(covariate_files[[cohort_name]], show_col_types = FALSE)
  protein_df <- read_csv(ev_protein_files[[cohort_name]], show_col_types = FALSE) %>% select_mediator_proteins()
  microbe_df <- read_csv(microbe_file, show_col_types = FALSE)
  microbe_names <- setdiff(names(microbe_df), "sample_ID")
  protein_names <- setdiff(names(protein_df), "sample_ID")
  microbe_resid <- residualize_matrix(microbe_df, covar_df, microbe_names, paste0("Figure5_", cohort_name, "_", microbe_type))
  protein_resid <- residualize_matrix(protein_df, covar_df, protein_names, paste0("Figure5_", cohort_name, "_mediator_EVproteins"))
  run_microbe_protein_association(microbe_resid, protein_resid, microbe_type, cohort_name)
}

species_results <- bind_rows(lapply(names(species_files), function(cohort_name) run_one_cohort(cohort_name, "Species", species_files[[cohort_name]])))
pathway_results <- bind_rows(lapply(names(pathway_files), function(cohort_name) run_one_cohort(cohort_name, "Pathway", pathway_files[[cohort_name]])))

write_csv(species_results, file.path(out_dir, "species_mediator_EVprotein_associations.csv"))
write_csv(pathway_results, file.path(out_dir, "pathway_mediator_EVprotein_associations.csv"))

# Merge discovery and replication microbe-protein association results.
merge_discovery_replication <- function(res) {
  disc <- res %>% filter(Cohort == "Discovery") %>% rename_with(~ paste0(.x, "_discovery"), .cols = -c(Microbe_type, microbe, protein))
  rep <- res %>% filter(Cohort == "Replication") %>% rename_with(~ paste0(.x, "_replication"), .cols = -c(Microbe_type, microbe, protein))
  inner_join(disc, rep, by = c("Microbe_type", "microbe", "protein")) %>%
    mutate(replicated_nominal = sign(R_discovery) == sign(R_replication) & P_replication < 0.05)
}

species_merged <- merge_discovery_replication(species_results)
pathway_merged <- merge_discovery_replication(pathway_results)

write_csv(species_merged, file.path(out_dir, "species_mediator_EVprotein_discovery_replication_merged.csv"))
write_csv(pathway_merged, file.path(out_dir, "pathway_mediator_EVprotein_discovery_replication_merged.csv"))
write_csv(species_merged %>% filter(FDR_by_microbe_discovery < 0.05, replicated_nominal), file.path(out_dir, "species_mediator_EVprotein_replicated_hits.csv"))
write_csv(pathway_merged %>% filter(FDR_by_microbe_discovery < 0.05, replicated_nominal), file.path(out_dir, "pathway_mediator_EVprotein_replicated_hits.csv"))




#----- 6. Genetic regulation of mediator EV proteins -----
# Analysis content: this script prepares phenotype residuals for 56 mediator EV proteins,
# generates PLINK pQTL GWAS job scripts for discovery and validation cohorts, merges
# discovery-significant SNPs with validation results, annotates independent clumped lead pQTLs,
# summarizes cis/trans pQTL patterns, and tests downstream SNP-protein and optional
# microbe-SNP-protein interaction models.

rm(list = ls())
library(data.table)
library(dplyr)
library(readr)
library(tidyr)
project_dir <- "/path/to/project"
data_dir <- file.path(project_dir, "data")
out_dir <- file.path(project_dir, "results", "Figure6")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The functions below define each analysis step. Run the required steps manually in sequence.

n_mediator_proteins <- 56
genomewide_p <- 5e-8
bonferroni_p <- genomewide_p / n_mediator_proteins
replication_p <- 0.05
manhattan_extract_p <- 1e-3
clump_p1 <- 5e-8
clump_p2 <- 1e-5
clump_r2 <- 0.05
clump_kb <- 500
min_n <- 50

plink_work_dir <- file.path(project_dir, "plink_56_mediator_EVprotein_pQTL")
plink_script_dir <- file.path(plink_work_dir, "slurm_jobs")
plink_log_dir <- file.path(plink_work_dir, "log")
plink_bfile <- "/path/to/genotype/plink_binary_prefix"
conda_activate_plink <- "source /path/to/miniconda3/bin/activate gwas_env"
conda_activate_r <- "source /path/to/miniconda3/bin/activate r_env"
slurm_partition <- "compute_partition"

plink_cohorts <- tibble::tribble(
  ~Cohort,       ~pheno_file,                                                                  ~result_dir,           ~sig5e8_dir,                     ~sigbonf_dir,                      ~clump_dir,
  "Discovery",  file.path(plink_work_dir, "discovery_56mediator_EVprotein_residuals.txt"),   "result_discovery",   "result_sig_5x10-8_discovery",   "result_sig_bonf_discovery",      "gwas_discovery_clumped",
  "Validation", file.path(plink_work_dir, "validation_56mediator_EVprotein_residuals.txt"),  "result_validation",  "result_sig_5x10-8_validation",  "result_sig_bonf_validation",     "gwas_validation_clumped"
)

pqtl_summary_file <- file.path(data_dir, "mediator_EVprotein_pQTL_summary.csv")
protein_gene_position_file <- file.path(data_dir, "protein_gene_positions_GRCh38.csv")
genomic_lambda_file <- file.path(data_dir, "genomic_inflation_lambda_all_pQTL.csv")
lead_snp_file <- file.path(data_dir, "lead_pQTLs_annotated.tsv")
snp_annotation_file <- file.path(data_dir, "lead_pQTLs_gene_annotated.tsv")

covariate_files <- list(
  Discovery = file.path(data_dir, "discovery_score_covariates.csv"),
  Validation = file.path(data_dir, "validation_score_covariates.csv")
)

ev_protein_files <- list(
  Discovery = file.path(data_dir, "discovery_EV_protein_log_RSN.csv"),
  Validation = file.path(data_dir, "validation_EV_protein_log_RSN.csv")
)

genotype_files <- list(
  Discovery = file.path(data_dir, "discovery_lead_pQTL_genotypes_additive.csv"),
  Validation = file.path(data_dir, "validation_lead_pQTL_genotypes_additive.csv")
)

species_residual_files <- list(
  Discovery = file.path(project_dir, "results", "Figure5", "Discovery_Species_adjusted_residuals.csv"),
  Validation = file.path(project_dir, "results", "Figure5", "Validation_Species_adjusted_residuals.csv")
)

pathway_residual_files <- list(
  Discovery = file.path(project_dir, "results", "Figure5", "Figure5_Discovery_Pathway_adjusted_residuals.csv"),
  Validation = file.path(project_dir, "results", "Figure5", "Figure5_Validation_Pathway_adjusted_residuals.csv")
)

protein_residual_files <- list(
  Discovery = file.path(project_dir, "results", "Figure6", "Figure6_Discovery_mediator_EVprotein_residuals_for_pQTL.csv"),
  Validation = file.path(project_dir, "results", "Figure6", "Figure6_Validation_mediator_EVprotein_residuals_for_pQTL.csv")
)

species_protein_pair_file <- file.path(project_dir, "results", "Figure5", "Figure5_species_mediator_EVprotein_replicated_hits.csv")
pathway_protein_pair_file <- file.path(project_dir, "results", "Figure5", "Figure5_pathway_mediator_EVprotein_replicated_hits.csv")
mediator_file <- file.path(project_dir, "results", "Figure4", "Figure4_significant_forward_mediator_EVproteins.csv")

covars <- c("Age", "Gender", "BMI", "Smoking_status_now", "Alcohol_intake_status", "Education_level")

# Standardize protein identifiers for matching across annotation tables.
std_key <- function(x) toupper(gsub("[.]", "-", x))

# Quote a shell argument safely for generated SLURM scripts.
shell_quote <- function(x) shQuote(x, type = "sh")

# Make file-safe names for job and output prefixes.
safe_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

# Return the first available column from a candidate list.
first_existing_col <- function(dat, candidates, required = TRUE) {
  hit <- intersect(candidates, names(dat))[1]
  if (is.na(hit) && required) stop("Required column not found: ", paste(candidates, collapse = ", "))
  hit
}

# Convert covariates to analysis-ready factor variables.
prepare_covariates <- function(dat) {
  dat %>%
    mutate(
      Gender = factor(Gender),
      Smoking_status_now = factor(Smoking_status_now),
      Alcohol_intake_status = factor(Alcohol_intake_status),
      Education_level = factor(Education_level)
    )
}

# Remove covariate effects from one mediator EV protein.
residualize_one <- function(data, feature, covariates) {
  keep <- c(feature, covariates)
  df <- data[, keep, drop = FALSE]
  idx <- complete.cases(df)
  res <- rep(NA_real_, nrow(data))
  if (sum(idx) < min_n) return(res)
  model_df <- cbind(y = df[[feature]][idx], df[idx, covariates, drop = FALSE])
  fit <- tryCatch(lm(y ~ ., data = model_df), error = function(e) NULL)
  if (is.null(fit)) return(res)
  res[idx] <- resid(fit)
  res
}

# Generate residualized mediator EV protein matrices for pQTL analysis.
make_protein_residuals <- function(cohort_name) {
  mediator_tbl <- read_csv(mediator_file, show_col_types = FALSE)
  protein_df <- read_csv(ev_protein_files[[cohort_name]], show_col_types = FALSE)
  target_col <- first_existing_col(mediator_tbl, c("protein", "Protein", "mediator", "GENE_name"))
  target_proteins <- intersect(unique(mediator_tbl[[target_col]]), names(protein_df))
  if (length(target_proteins) == 0) stop("No mediator EV proteins were found in the protein matrix for ", cohort_name, ".")
  covar_df <- read_csv(covariate_files[[cohort_name]], show_col_types = FALSE) %>% prepare_covariates()
  dat <- covar_df %>% inner_join(protein_df %>% select(sample_ID, all_of(target_proteins)), by = "sample_ID")
  out <- tibble(sample_ID = dat$sample_ID)
  for (protein in target_proteins) out[[protein]] <- residualize_one(dat, protein, covars)
  write_csv(out, protein_residual_files[[cohort_name]])
  out
}

# Write a PLINK phenotype file with FID/IID plus residualized mediator EV proteins.
write_plink_pheno_file <- function(cohort_name, pheno_file) {
  resid_df <- if (file.exists(protein_residual_files[[cohort_name]])) {
    read_csv(protein_residual_files[[cohort_name]], show_col_types = FALSE)
  } else {
    make_protein_residuals(cohort_name)
  }
  pheno <- resid_df %>% mutate(FID = sample_ID, IID = sample_ID) %>% select(FID, IID, -sample_ID, everything())
  dir.create(dirname(pheno_file), recursive = TRUE, showWarnings = FALSE)
  write_tsv(pheno, pheno_file, na = "NA")
  pheno_file
}

# Extract mediator protein names from a PLINK phenotype file, assuming columns 1-2 are FID/IID.
get_plink_pheno_names <- function(pheno_file) {
  if (!file.exists(pheno_file)) stop("Phenotype file does not exist: ", pheno_file)
  header <- names(fread(pheno_file, nrows = 0))
  if (length(header) <= 2) stop("Phenotype file must contain FID, IID, and at least one protein column.")
  header[-c(1, 2)]
}

# Generate one SLURM file per mediator protein for PLINK linear pQTL analysis and clumping.
write_plink_slurm_scripts <- function(cohort_name, pheno_file, result_dir, sig5e8_dir, sigbonf_dir, clump_dir,
                                      submit = FALSE) {
  proteins <- get_plink_pheno_names(pheno_file)
  dir.create(plink_script_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plink_log_dir, recursive = TRUE, showWarnings = FALSE)
  full_result_dir <- file.path(plink_work_dir, result_dir)
  full_sig5e8_dir <- file.path(plink_work_dir, sig5e8_dir)
  full_sigbonf_dir <- file.path(plink_work_dir, sigbonf_dir)
  full_clump_dir <- file.path(plink_work_dir, clump_dir)
  dir.create(full_result_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(full_sig5e8_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(full_sigbonf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(full_clump_dir, recursive = TRUE, showWarnings = FALSE)
  
  job_files <- character()
  for (protein in proteins) {
    prefix <- safe_name(protein)
    job_file <- file.path(plink_script_dir, paste0("run_", prefix, "_", tolower(cohort_name), ".slurm"))
    out_prefix <- file.path(full_result_dir, paste0(prefix, "_pQTL"))
    sig5e8_file <- file.path(full_sig5e8_dir, paste0(prefix, "_pQTL.assoc.linear.sig"))
    sigbonf_file <- file.path(full_sigbonf_dir, paste0(prefix, "_pQTL.assoc.linear.sig"))
    clump_prefix <- file.path(full_clump_dir, paste0(prefix, "_pQTL_clumped"))
    log_out <- file.path(plink_log_dir, paste0(prefix, "_", tolower(cohort_name), ".out"))
    log_err <- file.path(plink_log_dir, paste0(prefix, "_", tolower(cohort_name), ".err"))
    
    script <- c(
      "#!/bin/bash",
      paste0("#SBATCH -J pQTL_", prefix),
      paste0("#SBATCH -p ", slurm_partition),
      "#SBATCH -N 1",
      "#SBATCH -n 1",
      paste0("#SBATCH -o ", shell_quote(log_out)),
      paste0("#SBATCH -e ", shell_quote(log_err)),
      "set -euo pipefail",
      conda_activate_plink,
      paste0("cd ", shell_quote(plink_work_dir)),
      paste0("mkdir -p ", shell_quote(full_result_dir), " ", shell_quote(full_sig5e8_dir), " ", shell_quote(full_sigbonf_dir), " ", shell_quote(full_clump_dir)),
      paste0("echo '[PLINK] processing ", protein, " in ", cohort_name, " cohort'"),
      paste(
        "plink",
        "--bfile", shell_quote(plink_bfile),
        "--pheno", shell_quote(pheno_file),
        "--pheno-name", shell_quote(protein),
        "--linear --allow-no-sex --threads 1",
        "--out", shell_quote(out_prefix)
      ),
      paste0("awk 'NR==1 || $9 < ", format(genomewide_p, scientific = TRUE), "' ", shell_quote(paste0(out_prefix, ".assoc.linear")), " > ", shell_quote(sig5e8_file)),
      paste0("awk 'NR==1 || $9 < ", format(bonferroni_p, scientific = TRUE), "' ", shell_quote(paste0(out_prefix, ".assoc.linear")), " > ", shell_quote(sigbonf_file)),
      paste(
        "plink",
        "--bfile", shell_quote(plink_bfile),
        "--clump", shell_quote(paste0(out_prefix, ".assoc.linear")),
        "--clump-p1", clump_p1,
        "--clump-p2", clump_p2,
        "--clump-r2", clump_r2,
        "--clump-kb", clump_kb,
        "--threads 1",
        "--out", shell_quote(clump_prefix)
      )
    )
    writeLines(script, job_file)
    job_files <- c(job_files, job_file)
    if (isTRUE(submit)) system2("sbatch", job_file)
  }
  
  submit_file <- file.path(plink_script_dir, paste0("submit_all_", tolower(cohort_name), "_pQTL_jobs.sh"))
  writeLines(c("#!/bin/bash", paste("sbatch", shell_quote(job_files))), submit_file)
  Sys.chmod(c(job_files, submit_file), mode = "0755")
  tibble(Cohort = cohort_name, protein = proteins, job_file = job_files)
}

# Merge discovery-significant pQTLs with the corresponding validation-cohort PLINK output.
merge_discovery_validation_pqtl <- function(discovery_sig_dir, validation_result_dir, out_dir,
                                            p_thr = genomewide_p, validation_p_thr = replication_p) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(discovery_sig_dir, pattern = "_pQTL\\.assoc\\.linear\\.sig$", full.names = TRUE)
  if (length(files) == 0) return(data.table())
  
  merged_list <- list()
  for (f in files) {
    protein_name <- sub("_pQTL\\.assoc\\.linear\\.sig$", "", basename(f))
    val_file <- file.path(validation_result_dir, paste0(protein_name, "_pQTL.assoc.linear"))
    if (!file.exists(val_file)) next
    
    dis <- fread(f)
    if (!all(c("SNP", "P") %in% names(dis))) next
    sig_dis <- dis[!is.na(P) & P <= p_thr]
    if (nrow(sig_dis) == 0) next
    
    val <- fread(val_file)
    if (!("SNP" %in% names(val))) next
    merged <- merge(sig_dis, val, by = "SNP", all.x = TRUE, suffixes = c("_Discovery", "_Validation"))
    merged[, protein := protein_name]
    
    beta_d <- first_existing_col(merged, c("BETA_Discovery", "BETA.x", "BETA"), required = FALSE)
    beta_v <- first_existing_col(merged, c("BETA_Validation", "BETA.y"), required = FALSE)
    p_v <- first_existing_col(merged, c("P_Validation", "P.y"), required = FALSE)
    if (!is.na(beta_d) && !is.na(beta_v) && !is.na(p_v)) {
      merged[, replicated_nominal := !is.na(get(p_v)) & get(p_v) < validation_p_thr & sign(get(beta_d)) == sign(get(beta_v))]
    } else {
      merged[, replicated_nominal := NA]
    }
    
    out_file <- file.path(out_dir, paste0(protein_name, "_sigSNP_validation.txt"))
    fwrite(merged, out_file, sep = "\t")
    merged_list[[protein_name]] <- merged
  }
  
  all_merged <- rbindlist(merged_list, use.names = TRUE, fill = TRUE)
  if (nrow(all_merged) > 0) {
    fwrite(all_merged, file.path(out_dir, "merged_sig_5x10-8_snp_dis_val_all_proteins.txt"), sep = "\t")
    fwrite(all_merged[replicated_nominal == TRUE], file.path(out_dir, "merged_sig_5x10-8_snp_dis_val_all_proteins_replicated.txt"), sep = "\t")
  }
  all_merged
}

# Merge per-protein discovery-validation pQTL files into one table.
combine_merged_pqtl_files <- function(indir, out_file) {
  files <- list.files(indir, pattern = "_sigSNP_validation\\.txt$", full.names = TRUE)
  if (length(files) == 0) return(data.table())
  all_list <- lapply(files, function(f) {
    protein_name <- sub("_sigSNP_validation\\.txt$", "", basename(f))
    dt <- fread(f)
    dt[, protein := protein_name]
    dt
  })
  all_merged <- rbindlist(all_list, use.names = TRUE, fill = TRUE)
  fwrite(all_merged, out_file, sep = "\t")
  all_merged
}

# Annotate independent lead pQTLs from PLINK clumped output files.
annotate_clumped_lead_snps <- function(merged_file, clump_dir, out_file, lead_out_file) {
  if (!file.exists(merged_file)) stop("Merged pQTL file does not exist: ", merged_file)
  merged <- fread(merged_file)
  if (!all(c("SNP", "protein") %in% names(merged))) stop("Merged pQTL file must contain SNP and protein columns.")
  
  clump_files <- list.files(clump_dir, pattern = "_pQTL_clumped\\.clumped$", full.names = TRUE)
  merged[, is_lead := FALSE]
  for (f in clump_files) {
    protein_name <- sub("_pQTL_clumped\\.clumped$", "", basename(f))
    clump <- tryCatch(fread(f, fill = TRUE), error = function(e) data.table())
    if (nrow(clump) == 0 || !("SNP" %in% names(clump))) next
    lead_snps <- unique(clump$SNP[!is.na(clump$SNP)])
    merged[protein == protein_name & SNP %in% lead_snps, is_lead := TRUE]
  }
  fwrite(merged, out_file, sep = "\t")
  fwrite(merged[is_lead == TRUE], lead_out_file, sep = "\t")
  merged
}

# Extract p < 1e-3 pQTL records from all discovery GWAS files for downstream Manhattan-plot input.
extract_pqtl_records_for_manhattan <- function(result_dir, out_file, p_thr = manhattan_extract_p) {
  files <- list.files(result_dir, pattern = "_pQTL\\.assoc\\.linear$", full.names = TRUE)
  if (length(files) == 0) return(data.table())
  all_list <- list()
  for (f in files) {
    dat <- fread(f)
    required <- c("SNP", "CHR", "BP", "P")
    if (!all(required %in% names(dat))) next
    prefix <- sub("_pQTL\\.assoc\\.linear$", "", basename(f))
    man_dat <- dat[, .(SNP, CHR = as.numeric(CHR), BP = as.numeric(BP), P = as.numeric(P))]
    man_dat <- man_dat[!is.na(CHR) & !is.na(BP) & !is.na(P) & P < p_thr]
    if (nrow(man_dat) > 0) {
      man_dat[, Protein := prefix]
      all_list[[prefix]] <- man_dat
    }
  }
  all_sig_dat <- rbindlist(all_list, use.names = TRUE, fill = TRUE)
  if (nrow(all_sig_dat) > 0) fwrite(all_sig_dat, out_file, sep = "\t")
  all_sig_dat
}

# Estimate genomic inflation lambda from PLINK linear-association P values.
calculate_lambda_from_plink <- function(result_dir, out_file) {
  files <- list.files(result_dir, pattern = "_pQTL\\.assoc\\.linear$", full.names = TRUE)
  if (length(files) == 0) return(data.table())
  lambda_list <- lapply(files, function(f) {
    protein_name <- sub("_pQTL\\.assoc\\.linear$", "", basename(f))
    dat <- fread(f, select = "P")
    p <- dat$P[!is.na(dat$P) & dat$P > 0 & dat$P <= 1]
    if (length(p) == 0) return(data.table())
    lambda <- median(qchisq(1 - p, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
    data.table(protein = protein_name, n_snps = length(p), lambda = lambda)
  })
  lambda_df <- rbindlist(lambda_list, use.names = TRUE, fill = TRUE)
  fwrite(lambda_df, out_file, sep = "\t")
  lambda_df
}

# Annotate each pQTL as cis or trans using genomic distance to the protein-coding gene.
annotate_cis_trans <- function(pqtl_df, gene_pos_df) {
  gene_name_col <- first_existing_col(gene_pos_df, c("Genes", "GENE_name", "protein", "Protein"))
  pqtl_df2 <- pqtl_df %>% mutate(match_key = std_key(protein))
  gene_pos2 <- gene_pos_df %>%
    mutate(match_key = std_key(.data[[gene_name_col]])) %>%
    rename_with(~ gsub(" ", "_", .x)) %>%
    distinct(match_key, .keep_all = TRUE)
  joined <- pqtl_df2 %>% left_join(gene_pos2, by = "match_key") %>% select(-match_key)
  chr_col <- first_existing_col(joined, c("protein_CHR", "Gene_CHR", "CHR"))
  start_col <- first_existing_col(joined, c("Gene_start__(GRCh38)", "Gene_start_(GRCh38)", "Gene_start", "protein_Gene_start__(GRCh38)", "protein_Gene_start_(GRCh38)"))
  end_col <- first_existing_col(joined, c("Gene_end_(GRCh38)", "Gene_end", "protein_Gene_end_(GRCh38)"))
  pqtl_chr_col <- first_existing_col(joined, c("pQTL_CHR", "CHR_pQTL", "CHR_Discovery", "CHR"))
  pqtl_bp_col <- first_existing_col(joined, c("pQTL_BP", "BP_pQTL", "BP_Discovery", "BP"))
  joined %>%
    mutate(
      protein_center = (.data[[start_col]] + .data[[end_col]]) / 2,
      is_cis = case_when(
        as.character(.data[[chr_col]]) == as.character(.data[[pqtl_chr_col]]) & abs(protein_center - .data[[pqtl_bp_col]]) <= 1e6 ~ "cis-pQTL",
        !is.na(.data[[chr_col]]) & !is.na(.data[[pqtl_chr_col]]) ~ "trans-pQTL",
        TRUE ~ NA_character_
      )
    )
}

# Summarize genomic inflation lambda values across pQTL scans.
summarise_lambda <- function(lambda_df) {
  lambda_col <- first_existing_col(lambda_df, c("lambda", "genomic_lambda"))
  tibble(
    n_tests = sum(!is.na(lambda_df[[lambda_col]])),
    median_lambda = median(lambda_df[[lambda_col]], na.rm = TRUE),
    mean_lambda = mean(lambda_df[[lambda_col]], na.rm = TRUE),
    sd_lambda = sd(lambda_df[[lambda_col]], na.rm = TRUE),
    min_lambda = min(lambda_df[[lambda_col]], na.rm = TRUE),
    max_lambda = max(lambda_df[[lambda_col]], na.rm = TRUE)
  )
}

# Summarize how many proteins are associated with each SNP and vice versa.
summarise_pqtl_multiplicity <- function(lead_df) {
  snp_col <- first_existing_col(lead_df, c("SNP", "ID", "pQTL", "variant"))
  protein_col <- first_existing_col(lead_df, c("protein", "Protein"))
  dat <- lead_df %>% distinct(.data[[snp_col]], .data[[protein_col]], .keep_all = TRUE)
  snp_per_protein <- dat %>% group_by(SNP = .data[[snp_col]]) %>% summarise(n_proteins = n_distinct(.data[[protein_col]]), .groups = "drop")
  protein_per_snp <- dat %>% group_by(protein = .data[[protein_col]]) %>% summarise(n_snps = n_distinct(.data[[snp_col]]), .groups = "drop")
  list(
    snp_per_protein = snp_per_protein,
    protein_per_snp = protein_per_snp,
    snp_frequency = snp_per_protein %>% count(n_proteins, name = "frequency"),
    protein_frequency = protein_per_snp %>% mutate(n_snps_bin = ifelse(n_snps > 10, ">10", as.character(n_snps))) %>% count(n_snps_bin, name = "frequency")
  )
}

# Locate the genotype dosage column corresponding to one lead SNP.
find_snp_column <- function(genotype_df, snp_id) {
  hit <- names(genotype_df)[names(genotype_df) == snp_id | startsWith(names(genotype_df), paste0(snp_id, "_"))]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# Test lead SNP associations with residualized mediator EV proteins.
run_snp_protein_association <- function(cohort_name, lead_df) {
  genotype_df <- read_csv(genotype_files[[cohort_name]], show_col_types = FALSE)
  protein_resid <- if (file.exists(protein_residual_files[[cohort_name]])) read_csv(protein_residual_files[[cohort_name]], show_col_types = FALSE) else make_protein_residuals(cohort_name)
  snp_col <- first_existing_col(lead_df, c("SNP", "ID", "pQTL", "variant"))
  protein_col <- first_existing_col(lead_df, c("protein", "Protein"))
  pairs <- lead_df %>% transmute(SNP = .data[[snp_col]], protein = .data[[protein_col]]) %>% distinct()
  common_id <- intersect(genotype_df$sample_ID, protein_resid$sample_ID)
  genotype_df <- genotype_df[match(common_id, genotype_df$sample_ID), ]
  protein_resid <- protein_resid[match(common_id, protein_resid$sample_ID), ]
  bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
    snp <- pairs$SNP[i]
    protein <- pairs$protein[i]
    snp_name <- find_snp_column(genotype_df, snp)
    if (is.na(snp_name) || !(protein %in% names(protein_resid))) return(tibble())
    df <- tibble(y = protein_resid[[protein]], g = genotype_df[[snp_name]]) %>% filter(complete.cases(.))
    if (nrow(df) < min_n || length(unique(df$g)) < 2) return(tibble())
    fit <- tryCatch(lm(y ~ g, data = df), error = function(e) NULL)
    if (is.null(fit) || !("g" %in% rownames(summary(fit)$coefficients))) return(tibble())
    tab <- summary(fit)$coefficients
    tibble(Cohort = cohort_name, protein = protein, SNP = snp, N = nrow(df), Beta = tab["g", "Estimate"], SE = tab["g", "Std. Error"], T = tab["g", "t value"], P = tab["g", "Pr(>|t|)"])
  })) %>% mutate(FDR = p.adjust(P, method = "BH"))
}

# Fit optional microbe-SNP interaction models for mediator EV proteins.
run_microbe_snp_protein_interaction <- function(cohort_name, pair_file, microbe_residual_file, output_prefix) {
  if (!file.exists(pair_file) || !file.exists(microbe_residual_file) || !file.exists(protein_residual_files[[cohort_name]])) return(tibble())
  microbe_pairs <- read_csv(pair_file, show_col_types = FALSE)
  lead_df <- read_tsv(lead_snp_file, show_col_types = FALSE)
  snp_col <- first_existing_col(lead_df, c("SNP", "ID", "pQTL", "variant"))
  protein_col <- first_existing_col(lead_df, c("protein", "Protein"))
  snp_pairs <- lead_df %>% transmute(protein = .data[[protein_col]], SNP = .data[[snp_col]]) %>% distinct()
  protein_col_microbe <- first_existing_col(microbe_pairs, c("protein", "protein_discovery"))
  microbe_col <- first_existing_col(microbe_pairs, c("microbe", "species_discovery", "pathway_discovery"))
  pairs <- microbe_pairs %>% transmute(protein = .data[[protein_col_microbe]], microbe = .data[[microbe_col]]) %>% distinct() %>% inner_join(snp_pairs, by = "protein")
  microbe_df <- read_csv(microbe_residual_file, show_col_types = FALSE)
  protein_df <- read_csv(protein_residual_files[[cohort_name]], show_col_types = FALSE)
  genotype_df <- read_csv(genotype_files[[cohort_name]], show_col_types = FALSE)
  common_id <- Reduce(intersect, list(microbe_df$sample_ID, protein_df$sample_ID, genotype_df$sample_ID))
  microbe_df <- microbe_df[match(common_id, microbe_df$sample_ID), ]
  protein_df <- protein_df[match(common_id, protein_df$sample_ID), ]
  genotype_df <- genotype_df[match(common_id, genotype_df$sample_ID), ]
  res <- bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
    protein <- pairs$protein[i]
    microbe <- pairs$microbe[i]
    snp <- pairs$SNP[i]
    snp_name <- find_snp_column(genotype_df, snp)
    if (is.na(snp_name) || !(protein %in% names(protein_df)) || !(microbe %in% names(microbe_df))) return(tibble())
    df <- tibble(y = as.numeric(scale(protein_df[[protein]])), m = as.numeric(scale(microbe_df[[microbe]])), g = genotype_df[[snp_name]]) %>% filter(complete.cases(.))
    if (nrow(df) < min_n || length(unique(df$g)) < 2) return(tibble())
    fit <- tryCatch(lm(y ~ m * g, data = df), error = function(e) NULL)
    if (is.null(fit) || !("m:g" %in% rownames(summary(fit)$coefficients))) return(tibble())
    tab <- summary(fit)$coefficients
    beta <- tab["m:g", "Estimate"]
    se <- tab["m:g", "Std. Error"]
    tibble(Cohort = cohort_name, protein = protein, microbe = microbe, SNP = snp, N = nrow(df), Beta_interaction = beta, SE_interaction = se, CI95_low = beta - 1.96 * se, CI95_high = beta + 1.96 * se, T_interaction = tab["m:g", "t value"], P_interaction = tab["m:g", "Pr(>|t|)"])
  })) %>% mutate(FDR_interaction = p.adjust(P_interaction, method = "BH"))
  write_csv(res, file.path(out_dir, paste0(output_prefix, "_", cohort_name, "_microbe_SNP_protein_interaction.csv")))
  res
}