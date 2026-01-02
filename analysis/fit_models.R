#!/usr/bin/env Rscript
# Main fitting pipeline for IGT models
# Run: Rscript analysis/fit_models.R

# Load required libraries
required_packages <- c("rjags", "dplyr", "coda")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Install with: install.packages('%s')", pkg, pkg))
  }
}

# Load utility functions
source("analysis/utils/load_data.R")
source("analysis/utils/prepare_jags_data.R")

# Set random seed for reproducibility
set.seed(12345)

# Configuration
config <- list(
  n_chains = 4,
  n_adapt = 2000,
  n_burnin = 2000,
  n_iter = 5000,
  thin = 1,
  models = c("pvl_delta", "vse", "orl"),
  fit_all_studies = TRUE,  # If FALSE, fits only first study for testing
  parallel = FALSE         # Set to TRUE for parallel chains (requires parallel package)
)

cat("=== IGT Model Fitting Pipeline ===\n\n")

# Step 1: Load and validate data
cat("Step 1: Loading and validating data...\n")
all_data <- load_all_igt_data()
validate_igt_data(all_data)

# Prepare JAGS data
if (config$fit_all_studies) {
  cat("\nPreparing combined dataset for model fitting...\n")
  jags_data <- prepare_jags_data(all_data)
} else {
  cat("\nPreparing single study (Ahn2014_HC) for testing...\n")
  jags_data <- prepare_jags_data(all_data, study_filter = "Ahn2014_HC")
}

check_jags_data(jags_data)

# Step 2: Fit models
cat("\n=== Step 2: Fitting Models ===\n")

fit_results <- list()

for (model_name in config$models) {
  cat(sprintf("\n--- Fitting %s model ---\n", toupper(model_name)))

  # Model file path (use v2 versions for numerical stability)
  model_file <- sprintf("analysis/models/%s_v2.jags", model_name)

  if (!file.exists(model_file)) {
    warning(sprintf("Model file not found: %s. Skipping.", model_file))
    next
  }

  # Prepare model-specific data
  model_data <- prepare_jags_data_for_model(all_data, model_name,
                                            study_filter = if(!config$fit_all_studies) "Ahn2014_HC" else NULL)

  # Initialize JAGS model
  cat("Initializing JAGS model...\n")
  jags_model <- jags.model(
    file = model_file,
    data = model_data,
    n.chains = config$n_chains,
    n.adapt = config$n_adapt,
    quiet = FALSE
  )

  # Burn-in
  cat(sprintf("Running burn-in (%d iterations)...\n", config$n_burnin))
  update(jags_model, n.iter = config$n_burnin)

  # Parameters to monitor
  if (model_name == "pvl_delta") {
    params <- c("mu_A", "mu_alpha", "mu_cons", "mu_lambda",
                "sigma_A", "sigma_alpha", "sigma_cons", "sigma_lambda",
                "A", "alpha", "cons", "lambda")
  } else if (model_name == "vse") {
    params <- c("mu_A", "mu_alpha", "mu_cons", "mu_lambda", "mu_epP", "mu_epN", "mu_K", "mu_w",
                "sigma_A", "sigma_alpha", "sigma_cons", "sigma_lambda", "sigma_epP", "sigma_epN", "sigma_K", "sigma_w",
                "A", "alpha", "cons", "lambda", "epP", "epN", "K", "w")
  } else if (model_name == "orl") {
    params <- c("mu_Arew", "mu_Apun", "mu_K", "mu_betaF", "mu_betaP",
                "sigma_Arew", "sigma_Apun", "sigma_K", "sigma_betaF", "sigma_betaP",
                "Arew", "Apun", "K", "betaF", "betaP")
  }

  # Sample from posterior
  cat(sprintf("Sampling from posterior (%d iterations)...\n", config$n_iter))
  samples <- coda.samples(
    model = jags_model,
    variable.names = params,
    n.iter = config$n_iter,
    thin = config$thin
  )

  # Store results
  fit_results[[model_name]] <- list(
    samples = samples,
    model_file = model_file,
    data = model_data,
    config = config
  )

  # Save individual model results
  output_file <- sprintf("analysis/outputs/%s_fit.rds", model_name)
  saveRDS(fit_results[[model_name]], file = output_file)
  cat(sprintf("Saved results to: %s\n", output_file))
}

# Step 3: Quick summary
cat("\n=== Step 3: Summary ===\n")
for (model_name in names(fit_results)) {
  cat(sprintf("\n%s model:\n", toupper(model_name)))

  # Get summary statistics
  smry <- summary(fit_results[[model_name]]$samples)

  # Print group-level parameters
  cat("Group-level parameters (mu):\n")
  mu_params <- grep("^mu_", rownames(smry$statistics), value = TRUE)
  if (length(mu_params) > 0) {
    print(round(smry$statistics[mu_params, c("Mean", "SD")], 3))
  }
}

cat("\n=== Fitting Complete ===\n")
cat("Results saved to: analysis/outputs/\n")
cat("Next steps:\n")
cat("  1. Run diagnostics: source('analysis/utils/diagnostics.R')\n")
cat("  2. Run posterior predictive checks: source('analysis/utils/ppc.R')\n")
cat("  3. Compare models: source('analysis/utils/model_comparison.R')\n")
