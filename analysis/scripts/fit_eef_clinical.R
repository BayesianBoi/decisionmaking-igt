# Fit EEF model to clinical populations (Ahn 2014 + Fridberg 2010)
# Research Question: Do substance users show elevated forgetting rates (lambda)?
#
# Run: Rscript analysis/scripts/fit_eef_clinical.R

library(rjags)
library(coda)

# Source utility functions
source("analysis/utils/load_data.R")
source("analysis/utils/prepare_eef_data.R")

#==============================================================================
# CONFIGURATION
#==============================================================================

# MCMC Settings - CRITICAL DECISIONS
# -----------------------------------
# Based on hierarchical Bayesian literature (Gelman et al., 2013):
# - Small N (N=38-48 per group) requires MORE iterations, not fewer
# - 4 parameters per subject = moderate complexity
# - Group-level parameters require adequate exploration
#
# RATIONALE for these settings:
# 1. n_adapt = 5000: Let JAGS tune sampler for 5K iterations
#    - Standard is 1K, but we have group structure (need more)
# 2. n_burnin = 10000: Discard first 10K samples
#    - Ensures convergence from arbitrary starting points
#    - Hierarchical models need longer burn-in
# 3. n_iter = 20000: Collect 20K samples per chain
#    - Yields 80K total samples (4 chains)
#    - Effective sample size should be >10K for stable inference
# 4. n_chains = 4: Standard for MCMC diagnostics
#    - Allows Gelman-Rubin R-hat convergence assessment
# 5. thin = 2: Store every 2nd sample
#    - Reduces autocorrelation
#    - Final: 40K samples (80K / 2)
#
# TOTAL RUNTIME ESTIMATE: 4-8 hours on M1 Mac (173 subjects, 100 trials each)

config <- list(
  # MCMC settings
  n_adapt = 5000,      # Adaptation phase
  n_burnin = 10000,    # Burn-in samples (discarded)
  n_iter = 20000,      # Sampling iterations per chain
  n_chains = 4,        # Number of chains
  thin = 2,            # Thinning interval

  # Convergence criteria
  rhat_threshold = 1.1,   # Gelman-Rubin R-hat (< 1.1 = good convergence)
  n_eff_min = 1000,       # Minimum effective sample size per parameter

  # Monitoring
  parameters_to_monitor = c(
    # Group-level means (primary research focus)
    "mu_theta",           # Outcome sensitivity
    "mu_lambda_forget",   # FORGETTING RATE (key parameter!)
    "mu_phi",             # Exploration incentive
    "mu_cons",            # Choice consistency

    # Group-level variability
    "sigma_theta",
    "sigma_lambda_forget",
    "sigma_phi",
    "sigma_cons",

    # Subject-level parameters (for individual differences & diagnostics)
    "theta",
    "lambda_forget",      # FORGETTING RATE (individual-level)
    "phi",
    "cons"
  )
)

# Output settings
output_dir <- "results/eef_clinical"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#==============================================================================
# LOAD AND PREPARE DATA
#==============================================================================

message("=== LOADING DATA ===\n")

# Load all IGT data
dat_all <- load_all_igt_data()

# Filter to clinical populations only
# Ahn 2014: HC, Amph, Heroin
# Fridberg 2010: HC, Cannabis
clinical_studies <- c("Ahn2014_HC", "Ahn2014_Amph", "Ahn2014_Hero",
                      "Fridberg2010_HC", "Fridberg2010_Cbis")

dat_clinical <- dat_all[dat_all$study %in% clinical_studies, ]

# Create group variable (collapse HC groups, separate substance types)
dat_clinical$group <- dat_clinical$study
dat_clinical$group[dat_clinical$group %in% c("Ahn2014_HC", "Fridberg2010_HC")] <- "HC"
dat_clinical$group[dat_clinical$group == "Ahn2014_Amph"] <- "Amphetamine"
dat_clinical$group[dat_clinical$group == "Ahn2014_Hero"] <- "Heroin"
dat_clinical$group[dat_clinical$group == "Fridberg2010_Cbis"] <- "Cannabis"

message(sprintf("Total subjects: %d", length(unique(dat_clinical$subj_unique))))
message(sprintf("Total trials: %d\n", nrow(dat_clinical)))

# Check first-choice patterns by group
message("=== FIRST-CHOICE ANALYSIS ===\n")
first_choice_props <- compare_first_choices_by_group(dat_clinical, group_var = "group")

# Prepare JAGS data with group structure
message("\n=== PREPARING JAGS DATA ===\n")
jags_data <- prepare_eef_jags_data(
  dat = dat_clinical,
  group_var = "group",
  reference_group = "HC"
)

# Validate data
check_eef_jags_data(jags_data)

# Save prepared data
saveRDS(jags_data, file.path(output_dir, "jags_data.rds"))
message(sprintf("\nJAGS data saved to: %s\n", file.path(output_dir, "jags_data.rds")))

#==============================================================================
# FIT MODEL
#==============================================================================

message("=== FITTING EEF MODEL ===\n")
message("This will take several hours. Progress will be displayed.\n")

# Load JAGS model
model_file <- "analysis/models/eef_clinical.jags"
if (!file.exists(model_file)) {
  stop(sprintf("Model file not found: %s", model_file))
}

# Initialize JAGS model
message(sprintf("Initializing model with %d chains...", config$n_chains))
start_time <- Sys.time()

jags_model <- jags.model(
  file = model_file,
  data = jags_data,
  n.chains = config$n_chains,
  n.adapt = config$n_adapt,
  quiet = FALSE  # Show adaptation progress
)

adapt_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
message(sprintf("Adaptation complete (%.1f minutes)", adapt_time))

# Burn-in
message(sprintf("\nBurn-in: %d iterations...", config$n_burnin))
update(jags_model, n.iter = config$n_burnin, progress.bar = "text")

burnin_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins")) - adapt_time
message(sprintf("Burn-in complete (%.1f minutes)", burnin_time))

# Sampling
message(sprintf("\nSampling: %d iterations x %d chains (thin=%d)...",
                config$n_iter, config$n_chains, config$thin))

samples <- coda.samples(
  model = jags_model,
  variable.names = config$parameters_to_monitor,
  n.iter = config$n_iter,
  thin = config$thin,
  progress.bar = "text"
)

total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
message(sprintf("\nSampling complete!"))
message(sprintf("Total runtime: %.1f minutes (%.1f hours)", total_time, total_time/60))

# Save raw samples
saveRDS(samples, file.path(output_dir, "mcmc_samples.rds"))
message(sprintf("MCMC samples saved to: %s", file.path(output_dir, "mcmc_samples.rds")))

#==============================================================================
# DIAGNOSTICS
#==============================================================================

message("\n=== CONVERGENCE DIAGNOSTICS ===\n")

# Gelman-Rubin R-hat
rhat <- gelman.diag(samples, multivariate = FALSE)
rhat_values <- rhat$psrf[, "Point est."]

message("R-hat statistics (want < 1.1):")
message(sprintf("  Range: [%.3f, %.3f]", min(rhat_values), max(rhat_values)))
message(sprintf("  Median: %.3f", median(rhat_values)))

n_converged <- sum(rhat_values < config$rhat_threshold)
n_total <- length(rhat_values)
message(sprintf("  Converged: %d/%d parameters (%.1f%%)",
                n_converged, n_total, 100*n_converged/n_total))

if (n_converged < n_total) {
  warning(sprintf("%d parameters did not converge (R-hat >= %.2f)",
                  n_total - n_converged, config$rhat_threshold))
  # Show worst offenders
  worst_idx <- order(rhat_values, decreasing = TRUE)[1:min(10, n_total)]
  message("\nWorst R-hat values:")
  print(sort(rhat_values[worst_idx], decreasing = TRUE))
}

# Effective sample size
eff_size <- effectiveSize(samples)
message(sprintf("\nEffective sample sizes:"))
message(sprintf("  Range: [%.0f, %.0f]", min(eff_size), max(eff_size)))
message(sprintf("  Median: %.0f", median(eff_size)))

n_adequate <- sum(eff_size >= config$n_eff_min)
message(sprintf("  Adequate (>=%d): %d/%d parameters (%.1f%%)",
                config$n_eff_min, n_adequate, n_total, 100*n_adequate/n_total))

# Save diagnostics
diagnostics <- list(
  rhat = rhat,
  eff_size = eff_size,
  runtime_minutes = total_time,
  config = config,
  timestamp = Sys.time()
)
saveRDS(diagnostics, file.path(output_dir, "diagnostics.rds"))

#==============================================================================
# PARAMETER SUMMARIES
#==============================================================================

message("\n=== PARAMETER ESTIMATES ===\n")

# Extract posterior means
posterior_summary <- summary(samples)
param_means <- posterior_summary$statistics[, "Mean"]
param_sds <- posterior_summary$statistics[, "SD"]

# Focus on group-level parameters
group_params <- c("mu_theta", "mu_lambda_forget", "mu_phi", "mu_cons")
message("Group-level means (posterior mean ± SD):")
for (p in group_params) {
  if (p %in% names(param_means)) {
    message(sprintf("  %s: %.3f ± %.3f", p, param_means[p], param_sds[p]))
  }
}

# KEY RESULT: Forgetting rate
lambda_mean <- param_means["mu_lambda_forget"]
lambda_sd <- param_sds["mu_lambda_forget"]
message(sprintf("\n*** KEY FINDING: Population forgetting rate = %.3f ± %.3f ***",
                lambda_mean, lambda_sd))

# Save summaries
saveRDS(posterior_summary, file.path(output_dir, "parameter_summary.rds"))

#==============================================================================
# VISUALIZATIONS
#==============================================================================

message("\n=== CREATING DIAGNOSTIC PLOTS ===\n")

# Trace plots for group-level parameters
pdf(file.path(output_dir, "trace_plots.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2))
for (p in group_params) {
  traceplot(samples[, p], main = p)
}
dev.off()

# Density plots
pdf(file.path(output_dir, "density_plots.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2))
for (p in group_params) {
  densplot(samples[, p], main = p)
}
dev.off()

message(sprintf("Diagnostic plots saved to: %s", output_dir))

#==============================================================================
# DONE
#==============================================================================

message("\n=== FITTING COMPLETE ===\n")
message(sprintf("Output directory: %s", output_dir))
message("Files created:")
message("  - jags_data.rds (prepared data)")
message("  - mcmc_samples.rds (posterior samples)")
message("  - parameter_summary.rds (posterior statistics)")
message("  - diagnostics.rds (convergence checks)")
message("  - trace_plots.pdf (MCMC trace plots)")
message("  - density_plots.pdf (posterior densities)")
message("\nNext steps:")
message("  1. Check diagnostics (R-hat, effective sample size)")
message("  2. Run group comparisons (analysis/scripts/compare_groups.R)")
message("  3. Run model comparison vs PVL-Delta and VSE")
