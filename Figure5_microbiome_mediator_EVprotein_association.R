# Figure 5: Association between gut microbiome features and mediator EV proteins.
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

mediator_file <- file.path(project_dir, "results", "Figure4", "Figure4_significant_forward_mediator_EVproteins.csv")

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

write_csv(species_results, file.path(out_dir, "Figure5_species_mediator_EVprotein_associations.csv"))
write_csv(pathway_results, file.path(out_dir, "Figure5_pathway_mediator_EVprotein_associations.csv"))

# Merge discovery and replication microbe-protein association results.
merge_discovery_replication <- function(res) {
  disc <- res %>% filter(Cohort == "Discovery") %>% rename_with(~ paste0(.x, "_discovery"), .cols = -c(Microbe_type, microbe, protein))
  rep <- res %>% filter(Cohort == "Replication") %>% rename_with(~ paste0(.x, "_replication"), .cols = -c(Microbe_type, microbe, protein))
  inner_join(disc, rep, by = c("Microbe_type", "microbe", "protein")) %>%
    mutate(replicated_nominal = sign(R_discovery) == sign(R_replication) & P_replication < 0.05)
}

species_merged <- merge_discovery_replication(species_results)
pathway_merged <- merge_discovery_replication(pathway_results)

write_csv(species_merged, file.path(out_dir, "Figure5_species_mediator_EVprotein_discovery_replication_merged.csv"))
write_csv(pathway_merged, file.path(out_dir, "Figure5_pathway_mediator_EVprotein_discovery_replication_merged.csv"))
write_csv(species_merged %>% filter(FDR_by_microbe_discovery < 0.05, replicated_nominal), file.path(out_dir, "Figure5_species_mediator_EVprotein_replicated_hits.csv"))
write_csv(pathway_merged %>% filter(FDR_by_microbe_discovery < 0.05, replicated_nominal), file.path(out_dir, "Figure5_pathway_mediator_EVprotein_replicated_hits.csv"))
