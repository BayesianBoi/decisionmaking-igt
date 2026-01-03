# ===========================================================================
# Fit VSE Model to Clinical Populations (Parallel Chains)
# ===========================================================================
#
# Reference: Worthy, D.A., Pang, B., & Byrne, K.A. (2013). Decomposing the
#            roles of perseveration and expected value. Frontiers in Psychology.
#
# Usage: Rscript analysis/scripts/fit_vse.R
#
# ===========================================================================

library(rjags)
library(coda)
library(parallel)

source("analysis/utils/load_data.R")
source("analysis/utils/prepare_jags_data.R")

# ===========================================================================
# Configuration
# ===========================================================================

config <- list(
  n_adapt = 5000,
  n_burnin = 10000,
  n_iter = 20000,
  n_chains = 4,
  thin = 2,
  
  rhat_threshold = 1.1,
  n_eff_min = 1000,
  
  parameters_to_monitor = c(
    "mu_A", "mu_alpha", "mu_cons", "mu_lambda",
    "mu_epP", "mu_epN", "mu_K", "mu_w",
    "sigma_A", "sigma_alpha", "sigma_cons", "sigma_lambda",
    "sigma_epP", "sigma_epN", "sigma_K", "sigma_w",
    "A", "alpha", "cons", "lambda",
    "epP", "epN", "K", "w"
  )
)

output_dir <- "results/vse"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# Load and Prepare Data
# ===========================================================================

message("Loading data...")

dat_all <- load_all_igt_data()

clinical_studies <- c("Ahn2014_HC", "Ahn2014_Amph", "Ahn2014_Hero",
                      "Fridberg2010_HC", "Fridberg2010_Cbis")

dat_clinical <- dat_all[dat_all$study %in% clinical_studies, ]

message(sprintf("Subjects: %d", length(unique(dat_clinical$subj_unique))))
message(sprintf("Trials: %d", nrow(dat_clinical)))

message("Preparing JAGS data...")
jags_data <- prepare_jags_data(dat_clinical)
check_jags_data(jags_data)

saveRDS(jags_data, file.path(output_dir, "jags_data.rds"))

# ===========================================================================
# Fit Model with Parallel Chains
# ===========================================================================

message("\nFitting VSE model with parallel chains...")

model_file <- "analysis/models/vse.jags"
if (!file.exists(model_file)) {
  stop(sprintf("Model file not found: %s", model_file))
}

start_time <- Sys.time()

# Function to run one chain
run_chain <- function(chain_id, model_file, jags_data, config) {
  library(rjags)
  library(coda)
  
  set.seed(chain_id * 12345)
  
  model <- jags.model(
    file = model_file,
    data = jags_data,
    n.chains = 1,
    n.adapt = config$n_adapt,
    quiet = TRUE
  )
  
  update(model, n.iter = config$n_burnin)
  
  samples <- coda.samples(
    model = model,
    variable.names = config$parameters_to_monitor,
    n.iter = config$n_iter,
    thin = config$thin
  )
  
  return(samples[[1]])
}

# Detect available cores
n_cores <- min(config$n_chains, detectCores() - 1)
message(sprintf("Running %d chains on %d cores...", config$n_chains, n_cores))

# Run chains in parallel
cl <- makeCluster(n_cores)
clusterExport(cl, c("model_file", "jags_data", "config"))

chain_results <- parLapply(cl, 1:config$n_chains, function(i) {
  run_chain(i, model_file, jags_data, config)
})

stopCluster(cl)

# Combine into mcmc.list
samples <- mcmc.list(chain_results)

total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
message(sprintf("Sampling complete (%.1f min)", total_time))

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

eff_size <- effectiveSize(samples)
message(sprintf("ESS range: [%.0f, %.0f]", min(eff_size), max(eff_size)))

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

group_params <- c("mu_A", "mu_alpha", "mu_cons", "mu_lambda",
                  "mu_epP", "mu_epN", "mu_K", "mu_w")
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

pdf(file.path(output_dir, "trace_plots.pdf"), width = 12, height = 10)
par(mfrow = c(3, 3))
for (p in group_params) {
  traceplot(samples[, p], main = p)
}
dev.off()

pdf(file.path(output_dir, "density_plots.pdf"), width = 12, height = 10)
par(mfrow = c(3, 3))
for (p in group_params) {
  densplot(samples[, p], main = p)
}
dev.off()

message(sprintf("Output saved to: %s", output_dir))
message("Fitting complete.")
