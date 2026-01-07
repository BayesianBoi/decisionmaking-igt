# MCMC Diagnostics for JAGS model fits
# Run: source("utils/diagnostics.R")

library(coda)

#' Compute R-hat (Gelman-Rubin diagnostic) for all parameters
#' @param mcmc_samples MCMC list object from coda.samples
#' @return Data frame with R-hat values
compute_rhat <- function(mcmc_samples) {
  gr_diag <- gelman.diag(mcmc_samples, autoburnin = FALSE, multivariate = FALSE)
  rhat_df <- data.frame(
    parameter = rownames(gr_diag$psrf),
    rhat = gr_diag$psrf[, "Point est."],
    rhat_upper = gr_diag$psrf[, "Upper C.I."]
  )
  return(rhat_df)
}

#' Compute effective sample size for all parameters
#' @param mcmc_samples MCMC list object from coda.samples
#' @return Data frame with ESS values
compute_ess <- function(mcmc_samples) {
  ess_vals <- effectiveSize(mcmc_samples)
  ess_df <- data.frame(
    parameter = names(ess_vals),
    ess = as.numeric(ess_vals),
    stringsAsFactors = FALSE
  )
  return(ess_df)
}

#' Check convergence based on R-hat threshold
#' @param rhat_df Data frame from compute_rhat
#' @param threshold R-hat threshold (default 1.1)
#' @return List with convergence status and problematic parameters
check_convergence <- function(rhat_df, threshold = 1.1) {
  problematic <- rhat_df[rhat_df$rhat > threshold, ]

  converged <- nrow(problematic) == 0

  result <- list(
    converged = converged,
    threshold = threshold,
    n_total = nrow(rhat_df),
    n_problematic = nrow(problematic),
    problematic_params = problematic
  )

  return(result)
}

#' Generate diagnostic summary for a model fit
#' @param fit_result Model fit object with $samples component
#' @param model_name Name of the model
#' @return List with diagnostic results
generate_diagnostics <- function(fit_result, model_name = "model") {
  cat(sprintf("\n=== Diagnostics for %s ===\n", model_name))

  # Extract samples
  samples <- fit_result$samples

  # Compute R-hat
  cat("Computing R-hat...\n")
  rhat_df <- compute_rhat(samples)

  # Compute ESS
  cat("Computing effective sample size...\n")
  ess_df <- compute_ess(samples)

  # Merge diagnostics
  diag_df <- merge(rhat_df, ess_df, by = "parameter")

  # Check convergence
  conv_check <- check_convergence(rhat_df)

  # Print summary
  cat(sprintf("\nConvergence check (R-hat < %.2f):\n", conv_check$threshold))
  if (conv_check$converged) {
    cat("  Status: PASSED - All parameters converged\n")
  } else {
    cat(sprintf("  Status: WARNING - %d/%d parameters did not converge\n",
                conv_check$n_problematic, conv_check$n_total))
    cat("\nProblematic parameters:\n")
    print(conv_check$problematic_params)
  }

  # ESS summary
  cat(sprintf("\nEffective sample size:\n"))
  cat(sprintf("  Mean: %.0f\n", mean(diag_df$ess)))
  cat(sprintf("  Median: %.0f\n", median(diag_df$ess)))
  cat(sprintf("  Min: %.0f (%s)\n",
              min(diag_df$ess),
              diag_df$parameter[which.min(diag_df$ess)]))

  # Return results
  result <- list(
    diagnostics = diag_df,
    convergence = conv_check,
    summary_stats = summary(samples)
  )

  return(result)
}

#' Save trace plots for key parameters
#' @param mcmc_samples MCMC list object
#' @param output_file Path to save plot
#' @param params Parameters to plot (if NULL, plots group-level params)
save_trace_plots <- function(mcmc_samples, output_file = "analysis/outputs/trace_plots.pdf",
                             params = NULL) {

  if (is.null(params)) {
    # Default to group-level parameters (mu and sigma)
    all_params <- colnames(mcmc_samples[[1]])
    params <- grep("^(mu_|sigma_)", all_params, value = TRUE)
  }

  # Filter samples
  if (length(params) > 0) {
    samples_subset <- mcmc_samples[, params]

    pdf(output_file, width = 10, height = 8)
    plot(samples_subset, ask = FALSE)
    dev.off()

    cat(sprintf("Trace plots saved to: %s\n", output_file))
  } else {
    warning("No parameters to plot")
  }
}

#' Run full diagnostic pipeline on all fitted models
#' @param models_dir Directory containing fitted model .rds files
#' @return List of diagnostic results for each model
run_full_diagnostics <- function(models_dir = "analysis/outputs") {
  # Find all fitted models
  fit_files <- list.files(models_dir, pattern = "_fit\\.rds$", full.names = TRUE)

  if (length(fit_files) == 0) {
    stop("No fitted models found in: ", models_dir)
  }

  cat(sprintf("Found %d fitted models\n", length(fit_files)))

  all_diagnostics <- list()

  for (fit_file in fit_files) {
    # Extract model name
    model_name <- gsub("_fit\\.rds$", "", basename(fit_file))

    # Load fit
    cat(sprintf("\nLoading %s...\n", model_name))
    fit <- readRDS(fit_file)

    # Generate diagnostics
    diag_result <- generate_diagnostics(fit, model_name)

    # Save trace plots
    trace_file <- file.path(models_dir, sprintf("%s_trace_plots.pdf", model_name))
    save_trace_plots(fit$samples, trace_file)

    # Store results
    all_diagnostics[[model_name]] <- diag_result
  }

  # Save all diagnostics
  output_file <- file.path(models_dir, "all_diagnostics.rds")
  saveRDS(all_diagnostics, output_file)
  cat(sprintf("\nAll diagnostics saved to: %s\n", output_file))

  return(all_diagnostics)
}

# If running as script, run diagnostics on all models
if (!interactive()) {
  diagnostics <- run_full_diagnostics()
}
