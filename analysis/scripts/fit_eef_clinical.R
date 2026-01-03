# ===========================================================================
# Fit EEF Model to Clinical Populations
# ===========================================================================
#
# Fits the Exploitation-Exploration with Forgetting (EEF) model to Iowa 
# Gambling Task data from substance use disorder populations.
#
# Reference: Yang, X., et al. (2025). Exploitation and Exploration with 
#            Forgetting. Frontiers in Psychology.
#
# Data: Ahn et al. (2014), Fridberg et al. (2010)
#
# Usage: Rscript analysis/scripts/fit_eef_clinical.R
#
# ===========================================================================

library(rjags)
library(coda)

source("analysis/utils/load_data.R")
source("analysis/utils/prepare_eef_data.R")

# ===========================================================================
# Configuration
# ===========================================================================

# MCMC settings based on Gelman et al. (2013) recommendations for 
# hierarchical Bayesian models with small group sizes (N=38-48 per group)

config <- list(
  n_adapt = 5000,       # Adaptation iterations for sampler tuning
  n_burnin = 10000,     # Burn-in iterations (discarded)
  n_iter = 20000,       # Sampling iterations per chain
  n_chains = 4,         # Number of parallel chains
  thin = 2,             # Thinning interval to reduce autocorrelation
  
  rhat_threshold = 1.1, # Gelman-Rubin convergence criterion
  n_eff_min = 1000,     # Minimum effective sample size
  
  parameters_to_monitor = c(
    # Group-level means
    "mu_theta",
    "mu_lambda_forget",
    "mu_phi",
    "mu_cons",
    # Group-level standard deviations
    "sigma_theta",
    "sigma_lambda_forget",
    "sigma_phi",
    "sigma_cons",
    # Subject-level parameters
    "theta",
    "lambda_forget",
    "phi",
    "cons"
  )
)

output_dir <- "results/eef_clinical"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# Load and Prepare Data
# ===========================================================================

message("Loading data...")

dat_all <- load_all_igt_data()

# Filter to clinical populations
clinical_studies <- c("Ahn2014_HC", "Ahn2014_Amph", "Ahn2014_Hero",
                      "Fridberg2010_HC", "Fridberg2010_Cbis")

dat_clinical <- dat_all[dat_all$study %in% clinical_studies, ]

# Create group variable
dat_clinical$group <- dat_clinical$study
dat_clinical$group[dat_clinical$group %in% c("Ahn2014_HC", "Fridberg2010_HC")] <- "HC"
dat_clinical$group[dat_clinical$group == "Ahn2014_Amph"] <- "Amphetamine"
dat_clinical$group[dat_clinical$group == "Ahn2014_Hero"] <- "Heroin"
dat_clinical$group[dat_clinical$group == "Fridberg2010_Cbis"] <- "Cannabis"

message(sprintf("Subjects: %d", length(unique(dat_clinical$subj_unique))))
message(sprintf("Trials: %d", nrow(dat_clinical)))

# Analyze first-choice patterns by group
message("\nFirst-choice analysis:")
first_choice_props <- compare_first_choices_by_group(dat_clinical, group_var = "group")

# Prepare JAGS data with group structure
message("\nPreparing JAGS data...")
jags_data <- prepare_eef_jags_data(
  dat = dat_clinical,
  group_var = "group",
  reference_group = "HC"
)

check_eef_jags_data(jags_data)

saveRDS(jags_data, file.path(output_dir, "jags_data.rds"))
message(sprintf("JAGS data saved: %s", file.path(output_dir, "jags_data.rds")))

# ===========================================================================
# Fit Model
# ===========================================================================

message("\nFitting EEF model...")
message("Estimated runtime: 4-8 hours")

model_file <- "analysis/models/eef_clinical.jags"
if (!file.exists(model_file)) {
  stop(sprintf("Model file not found: %s", model_file))
}

message(sprintf("Initializing %d chains...", config$n_chains))
start_time <- Sys.time()

jags_model <- jags.model(
  file = model_file,
  data = jags_data,
  n.chains = config$n_chains,
  n.adapt = config$n_adapt,
  quiet = FALSE
)

adapt_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
message(sprintf("Adaptation complete (%.1f min)", adapt_time))

message(sprintf("Burn-in: %d iterations...", config$n_burnin))
update(jags_model, n.iter = config$n_burnin, progress.bar = "text")

burnin_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins")) - adapt_time
message(sprintf("Burn-in complete (%.1f min)", burnin_time))

message(sprintf("Sampling: %d iterations x %d chains...", config$n_iter, config$n_chains))

samples <- coda.samples(
  model = jags_model,
  variable.names = config$parameters_to_monitor,
  n.iter = config$n_iter,
  thin = config$thin,
  progress.bar = "text"
)

total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
message(sprintf("Sampling complete (%.1f min total)", total_time))

saveRDS(samples, file.path(output_dir, "mcmc_samples.rds"))

# ===========================================================================
# Convergence Diagnostics
# ===========================================================================

message("\nConvergence diagnostics:")

rhat <- gelman.diag(samples, multivariate = FALSE)
rhat_values <- rhat$psrf[, "Point est."]

message(sprintf("R-hat range: [%.3f, %.3f]", min(rhat_values), max(rhat_values)))
message(sprintf("R-hat median: %.3f", median(rhat_values)))

n_converged <- sum(rhat_values < config$rhat_threshold)
n_total <- length(rhat_values)
message(sprintf("Converged: %d/%d (%.1f%%)", n_converged, n_total, 100*n_converged/n_total))

if (n_converged < n_total) {
  warning(sprintf("%d parameters did not converge", n_total - n_converged))
}

eff_size <- effectiveSize(samples)
message(sprintf("ESS range: [%.0f, %.0f]", min(eff_size), max(eff_size)))
message(sprintf("ESS median: %.0f", median(eff_size)))

diagnostics <- list(
  rhat = rhat,
  eff_size = eff_size,
  runtime_minutes = total_time,
  config = config,
  timestamp = Sys.time()
)
saveRDS(diagnostics, file.path(output_dir, "diagnostics.rds"))

# ===========================================================================
# Parameter Summaries
# ===========================================================================

message("\nParameter estimates:")

posterior_summary <- summary(samples)
param_means <- posterior_summary$statistics[, "Mean"]
param_sds <- posterior_summary$statistics[, "SD"]

group_params <- c("mu_theta", "mu_lambda_forget", "mu_phi", "mu_cons")
for (p in group_params) {
  if (p %in% names(param_means)) {
    message(sprintf("  %s: %.3f (SD=%.3f)", p, param_means[p], param_sds[p]))
  }
}

saveRDS(posterior_summary, file.path(output_dir, "parameter_summary.rds"))

# ===========================================================================
# Diagnostic Plots
# ===========================================================================

message("\nGenerating diagnostic plots...")

pdf(file.path(output_dir, "trace_plots.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2))
for (p in group_params) {
  traceplot(samples[, p], main = p)
}
dev.off()

pdf(file.path(output_dir, "density_plots.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2))
for (p in group_params) {
  densplot(samples[, p], main = p)
}
dev.off()

message(sprintf("Plots saved to: %s", output_dir))

# ===========================================================================
# Summary
# ===========================================================================

message("\nFitting complete.")
message(sprintf("Output: %s", output_dir))
message("Files: jags_data.rds, mcmc_samples.rds, parameter_summary.rds, diagnostics.rds")
