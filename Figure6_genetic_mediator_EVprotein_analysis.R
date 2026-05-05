# Figure 6: Genetic regulation of mediator EV proteins.
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
  Discovery = file.path(project_dir, "results", "Figure5", "Figure5_Discovery_Species_adjusted_residuals.csv"),
  Validation = file.path(project_dir, "results", "Figure5", "Figure5_Validation_Species_adjusted_residuals.csv")
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
