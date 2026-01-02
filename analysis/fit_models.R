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
  parallel = TRUE,         # Run chains in parallel for faster execution
  n_cores = NULL           # NULL = auto-detect and use all available cores
)

# Setup parallel execution if requested
if (config$parallel) {
  if (!require("parallel", quietly = TRUE)) {
    warning("parallel package not available. Falling back to sequential execution.")
    config$parallel <- FALSE
  } else {
    # Detect available cores
    n_available <- parallel::detectCores()

    # Set n_cores: use all available if NULL, otherwise use specified value
    if (is.null(config$n_cores)) {
      config$n_cores <- min(n_available, config$n_chains)
    } else {
      config$n_cores <- min(config$n_cores, n_available, config$n_chains)
    }

    cat(sprintf("Parallel mode enabled: using %d cores (out of %d available)\n",
                config$n_cores, n_available))
  }
}

cat("=== IGT Model Fitting Pipeline ===\n")
cat(sprintf("Start time: %s\n\n", Sys.time()))

# Step 1: Load and validate data (with caching)
data_cache_file <- "analysis/outputs/cached_data.rds"

if (file.exists(data_cache_file)) {
  cat("Step 1: Loading cached data...\n")
  cached <- readRDS(data_cache_file)
  all_data <- cached$all_data
  jags_data <- cached$jags_data
  cat(sprintf("  ✓ Loaded from cache: %d subjects, %d trials\n",
              jags_data$N, sum(jags_data$Tsubj)))
} else {
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

  # Cache the data for future runs
  cat("\nCaching data for future runs...\n")
  saveRDS(list(all_data = all_data, jags_data = jags_data), file = data_cache_file)
  cat(sprintf("  ✓ Saved to: %s\n", data_cache_file))
}

# Step 2: Fit models
cat("\n=== Step 2: Fitting Models ===\n")

fit_results <- list()

for (model_name in config$models) {
  cat(sprintf("\n=================================================\n"))
  cat(sprintf("MODEL: %s\n", toupper(model_name)))
  cat(sprintf("=================================================\n"))
  cat(sprintf("Start time: %s\n\n", Sys.time()))

  # Check if already fitted
  output_file <- sprintf("analysis/outputs/%s_fit.rds", model_name)
  if (file.exists(output_file)) {
    cat(sprintf("  ⚠ Model already fitted. Loading existing results from:\n"))
    cat(sprintf("    %s\n", output_file))
    cat(sprintf("  To refit, delete this file and re-run the pipeline.\n\n"))
    fit_results[[model_name]] <- readRDS(output_file)
    next
  }

  # Model file path (use v2 versions for numerical stability)
  model_file <- sprintf("analysis/models/%s_v2.jags", model_name)

  if (!file.exists(model_file)) {
    warning(sprintf("Model file not found: %s. Skipping.", model_file))
    next
  }

  cat(sprintf("Step 1/4: Preparing model-specific data...\n"))
  # Prepare model-specific data
  model_data <- prepare_jags_data_for_model(all_data, model_name,
                                            study_filter = if(!config$fit_all_studies) "Ahn2014_HC" else NULL)
  cat(sprintf("  ✓ Data prepared for %s model\n\n", model_name))

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

  if (config$parallel) {
    # Parallel execution: run each chain as a separate process
    cat(sprintf("Step 2/4: Initializing parallel execution...\n"))
    cat(sprintf("  • Chains: %d\n", config$n_chains))
    cat(sprintf("  • Cores: %d\n", config$n_cores))
    cat(sprintf("  • Adaptation: %d iterations\n", config$n_adapt))
    cat(sprintf("  • Burn-in: %d iterations\n", config$n_burnin))
    cat(sprintf("  • Sampling: %d iterations\n\n", config$n_iter))

    # Function to fit a single chain with progress logging
    fit_chain <- function(chain_id) {
      # Initialize model for this chain
      chain_model <- jags.model(
        file = model_file,
        data = model_data,
        n.chains = 1,  # One chain per process
        n.adapt = config$n_adapt,
        quiet = TRUE
      )

      # Burn-in
      update(chain_model, n.iter = config$n_burnin)

      # Sample
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

    # Run chains in parallel
    cl <- parallel::makeCluster(config$n_cores)
    on.exit(parallel::stopCluster(cl))

    # Export necessary objects to cluster
    parallel::clusterExport(cl, c("model_file", "model_data", "params", "config"),
                           envir = environment())

    # Load rjags on each worker
    cat("  • Loading JAGS on worker processes...\n")
    parallel::clusterEvalQ(cl, {
      library(rjags)
      library(coda)
    })

    # Run chains
    cat(sprintf("  • Sampling started at: %s\n", Sys.time()))
    cat("  • This may take 15-30 minutes depending on model complexity...\n")
    flush.console()  # Force output to appear immediately

    chain_list <- parallel::parLapply(cl, 1:config$n_chains, fit_chain)

    chain_end <- Sys.time()
    chain_duration <- as.numeric(difftime(chain_end, chain_start, units = "mins"))

    cat(sprintf("  ✓ Sampling complete! Duration: %.1f minutes\n\n", chain_duration))

    # Combine chains into mcmc.list
    cat("Step 4/4: Combining chains and saving results...\n")
    samples <- as.mcmc.list(lapply(chain_list, function(x) x[[1]]))

    parallel::stopCluster(cl)

  } else {
    # Sequential execution: original approach
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

    # Sample from posterior
    cat(sprintf("Sampling from posterior (%d iterations)...\n", config$n_iter))
    samples <- coda.samples(
      model = jags_model,
      variable.names = params,
      n.iter = config$n_iter,
      thin = config$thin
    )
  }

  # Store results
  fit_results[[model_name]] <- list(
    samples = samples,
    model_file = model_file,
    data = model_data,
    config = config,
    timestamp = Sys.time()
  )

  # Save individual model results
  cat(sprintf("  • Saving results to: %s\n", output_file))
  saveRDS(fit_results[[model_name]], file = output_file)
  file_size_mb <- file.size(output_file) / 1024^2
  cat(sprintf("  ✓ Saved successfully (%.1f MB)\n", file_size_mb))

  # Quick convergence check
  cat("\n  Quick convergence check (group-level parameters):\n")
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

  model_end_time <- Sys.time()
  cat(sprintf("\n  Model complete at: %s\n", model_end_time))
  cat(sprintf("=================================================\n"))
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
