# Figure 2: EV protein associations with emotional distress and cardiometabolic disorders.
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

write_csv(ev_ed_results, file.path(out_dir, "Figure2_EVprotein_ED_associations.csv"))
write_csv(ev_cmd_results, file.path(out_dir, "Figure2_EVprotein_CMD_associations.csv"))
write_csv(plasma_ed_results, file.path(out_dir, "Figure2_plasma_protein_ED_associations.csv"))
write_csv(plasma_cmd_results, file.path(out_dir, "Figure2_plasma_protein_CMD_associations.csv"))

ed_merge <- merge_discovery_replication(ev_ed_results, c("Protein_type", "Trait", "protein"))
cmd_merge <- merge_discovery_replication(ev_cmd_results, c("Protein_type", "Trait", "protein"))

write_csv(ed_merge, file.path(out_dir, "Figure2_EVprotein_ED_discovery_replication_merged.csv"))
write_csv(cmd_merge, file.path(out_dir, "Figure2_EVprotein_CMD_discovery_replication_merged.csv"))
write_csv(ed_merge %>% filter(FDR_by_trait_discovery < 0.05), file.path(out_dir, "Figure2_EVprotein_ED_discovery_FDR005_merged.csv"))
write_csv(cmd_merge %>% filter(FDR_by_trait_discovery < 0.05), file.path(out_dir, "Figure2_EVprotein_CMD_discovery_FDR005_merged.csv"))

shared_ev_proteins <- ed_merge %>%
  filter(FDR_by_trait_discovery < 0.05) %>%
  transmute(ED = Trait, protein) %>%
  inner_join(cmd_merge %>% filter(FDR_by_trait_discovery < 0.05) %>% transmute(CMD = Trait, protein), by = "protein") %>%
  distinct()
write_csv(shared_ev_proteins, file.path(out_dir, "Figure2_ED_CMD_shared_significant_EVproteins.csv"))

if (requireNamespace("clusterProfiler", quietly = TRUE) && requireNamespace("ReactomePA", quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  gene_symbol <- unique(shared_ev_proteins$protein)
  gene_entrez <- clusterProfiler::bitr(gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  reactome_res <- ReactomePA::enrichPathway(gene = unique(gene_entrez$ENTREZID), organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.20, readable = TRUE)
  write_csv(as.data.frame(reactome_res), file.path(out_dir, "Figure2_shared_EVproteins_Reactome_enrichment.csv"))
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

write_csv(ed_heterogeneity, file.path(out_dir, "Figure2_EV_vs_plasma_ED_effect_heterogeneity.csv"))
write_csv(cmd_heterogeneity, file.path(out_dir, "Figure2_EV_vs_plasma_CMD_effect_heterogeneity.csv"))
