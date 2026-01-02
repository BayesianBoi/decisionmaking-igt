#!/usr/bin/env Rscript
# Fit a single IGT model
# Usage: Rscript analysis/fit_single_model.R pvl_delta
# Usage: Rscript analysis/fit_single_model.R vse
# Usage: Rscript analysis/fit_single_model.R orl

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  cat("Error: Please specify a model name\n")
  cat("Usage: Rscript analysis/fit_single_model.R <model_name>\n")
  cat("Available models: pvl_delta, vse, orl\n")
  quit(status = 1)
}

model_to_fit <- args[1]

# Validate model name
valid_models <- c("pvl_delta", "vse", "orl")
if (!model_to_fit %in% valid_models) {
  cat(sprintf("Error: Invalid model '%s'\n", model_to_fit))
  cat(sprintf("Available models: %s\n", paste(valid_models, collapse = ", ")))
  quit(status = 1)
}

# Load required libraries
required_packages <- c("rjags", "dplyr", "coda", "parallel")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Install with: install.packages('%s')", pkg, pkg))
  }
}

# Load utility functions
source("analysis/utils/load_data.R")
source("analysis/utils/prepare_jags_data.R")

# Set random seed
set.seed(12345)

# Configuration
config <- list(
  n_chains = 4,
  n_adapt = 2000,
  n_burnin = 2000,
  n_iter = 5000,
  thin = 1,
  parallel = TRUE,
  n_cores = NULL  # Auto-detect
)

# Setup parallel execution
if (config$parallel) {
  n_available <- parallel::detectCores()
  if (is.null(config$n_cores)) {
    config$n_cores <- min(n_available, config$n_chains)
  } else {
    config$n_cores <- min(config$n_cores, n_available, config$n_chains)
  }
  cat(sprintf("Parallel mode: using %d cores (out of %d available)\n\n",
              config$n_cores, n_available))
}

cat(sprintf("=================================================\n"))
cat(sprintf("Fitting single model: %s\n", toupper(model_to_fit)))
cat(sprintf("=================================================\n"))
cat(sprintf("Start time: %s\n\n", Sys.time()))

# Check if already fitted
output_file <- sprintf("analysis/outputs/%s_fit.rds", model_to_fit)
if (file.exists(output_file)) {
  cat(sprintf("⚠ Model already fitted!\n"))
  cat(sprintf("  File: %s\n", output_file))
  cat(sprintf("  To refit, delete this file first.\n"))
  quit(status = 0)
}

# Load data (with caching)
data_cache_file <- "analysis/outputs/cached_data.rds"

if (file.exists(data_cache_file)) {
  cat("Loading cached data...\n")
  cached <- readRDS(data_cache_file)
  all_data <- cached$all_data
  jags_data <- cached$jags_data
  cat(sprintf("  ✓ Loaded: %d subjects, %d trials\n\n",
              jags_data$N, sum(jags_data$Tsubj)))
} else {
  cat("Loading and validating data...\n")
  all_data <- load_all_igt_data()
  validate_igt_data(all_data)
  jags_data <- prepare_jags_data(all_data)
  check_jags_data(jags_data)

  cat("Caching data...\n")
  saveRDS(list(all_data = all_data, jags_data = jags_data), file = data_cache_file)
  cat(sprintf("  ✓ Cached to: %s\n\n", data_cache_file))
}

# Prepare model-specific data
cat("Step 1/4: Preparing model-specific data...\n")
model_file <- sprintf("analysis/models/%s_v2.jags", model_to_fit)

if (!file.exists(model_file)) {
  stop(sprintf("Model file not found: %s", model_file))
}

model_data <- prepare_jags_data_for_model(all_data, model_to_fit, study_filter = NULL)
cat(sprintf("  ✓ Data prepared\n\n"))

# Parameters to monitor
if (model_to_fit == "pvl_delta") {
  params <- c("mu_A", "mu_alpha", "mu_cons", "mu_lambda",
              "sigma_A", "sigma_alpha", "sigma_cons", "sigma_lambda",
              "A", "alpha", "cons", "lambda")
} else if (model_to_fit == "vse") {
  params <- c("mu_A", "mu_alpha", "mu_cons", "mu_lambda", "mu_epP", "mu_epN", "mu_K", "mu_w",
              "sigma_A", "sigma_alpha", "sigma_cons", "sigma_lambda", "sigma_epP", "sigma_epN", "sigma_K", "sigma_w",
              "A", "alpha", "cons", "lambda", "epP", "epN", "K", "w")
} else if (model_to_fit == "orl") {
  params <- c("mu_Arew", "mu_Apun", "mu_K", "mu_betaF", "mu_betaP",
              "sigma_Arew", "sigma_Apun", "sigma_K", "sigma_betaF", "sigma_betaP",
              "Arew", "Apun", "K", "betaF", "betaP")
}

# Run parallel MCMC
cat("Step 2/4: Initializing parallel execution...\n")
cat(sprintf("  • Chains: %d\n", config$n_chains))
cat(sprintf("  • Cores: %d\n", config$n_cores))
cat(sprintf("  • Adaptation: %d iterations\n", config$n_adapt))
cat(sprintf("  • Burn-in: %d iterations\n", config$n_burnin))
cat(sprintf("  • Sampling: %d iterations\n\n", config$n_iter))

fit_chain <- function(chain_id) {
  chain_model <- jags.model(
    file = model_file,
    data = model_data,
    n.chains = 1,
    n.adapt = config$n_adapt,
    quiet = TRUE
  )

  update(chain_model, n.iter = config$n_burnin)

  chain_samples <- coda.samples(
    model = chain_model,
    variable.names = params,
    n.iter = config$n_iter,
    thin = config$thin
  )

  return(chain_samples)
}

cat("Step 3/4: Running MCMC chains in parallel...\n")
chain_start <- Sys.time()

cl <- parallel::makeCluster(config$n_cores)
parallel::clusterExport(cl, c("model_file", "model_data", "params", "config"),
                       envir = environment())
parallel::clusterEvalQ(cl, {
  library(rjags)
  library(coda)
})

cat(sprintf("  • Sampling started at: %s\n", Sys.time()))
cat("  • This may take 15-30 minutes...\n")
flush.console()

chain_list <- parallel::parLapply(cl, 1:config$n_chains, fit_chain)

chain_end <- Sys.time()
chain_duration <- as.numeric(difftime(chain_end, chain_start, units = "mins"))

cat(sprintf("  ✓ Sampling complete! Duration: %.1f minutes\n\n", chain_duration))

samples <- as.mcmc.list(lapply(chain_list, function(x) x[[1]]))
parallel::stopCluster(cl)

# Save results
cat("Step 4/4: Saving results...\n")
fit_result <- list(
  samples = samples,
  model_file = model_file,
  data = model_data,
  config = config,
  timestamp = Sys.time()
)

saveRDS(fit_result, file = output_file)
file_size_mb <- file.size(output_file) / 1024^2
cat(sprintf("  ✓ Saved to: %s (%.1f MB)\n\n", output_file, file_size_mb))

# Quick convergence check
cat("Quick convergence check (group-level parameters):\n")
gelman_result <- try(gelman.diag(samples), silent = TRUE)
if (!inherits(gelman_result, "try-error")) {
  mu_params <- grep("^mu_", rownames(gelman_result$psrf), value = TRUE)
  if (length(mu_params) > 0) {
    rhat_vals <- gelman_result$psrf[mu_params, "Point est."]
    max_rhat <- max(rhat_vals, na.rm = TRUE)
    if (max_rhat < 1.1) {
      cat(sprintf("  ✓ Convergence looks good! Max R-hat = %.3f\n", max_rhat))
    } else {
      cat(sprintf("  ⚠ Warning: Max R-hat = %.3f (target < 1.1)\n", max_rhat))
    }
  }
}

cat(sprintf("\n=================================================\n"))
cat(sprintf("Model complete at: %s\n", Sys.time()))
cat(sprintf("Total duration: %.1f minutes\n", chain_duration))
cat(sprintf("=================================================\n"))
