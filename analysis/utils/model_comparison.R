# Model Comparison for IGT models
# Run: source("analysis/utils/model_comparison.R")

library(coda)

#' Compute DIC (Deviance Information Criterion) for a JAGS model
#' @param fit_result Fitted model object with samples
#' @param jags_model JAGS model object (if available)
#' @return DIC value and effective number of parameters
compute_dic <- function(fit_result, jags_model = NULL) {
  # Note: DIC computation in JAGS requires the dic.samples function
  # which needs the original jags model object
  # For now, we'll extract DIC if it was computed during fitting

  if (!is.null(jags_model)) {
    dic_samples <- dic.samples(jags_model, n.iter = 1000, type = "pD")
    dic_value <- sum(dic_samples$deviance) + sum(dic_samples$penalty)
    pD <- sum(dic_samples$penalty)

    return(list(
      DIC = dic_value,
      pD = pD,
      Dbar = sum(dic_samples$deviance)
    ))
  } else {
    warning("JAGS model object not available. Cannot compute DIC.")
    return(NULL)
  }
}

#' Compute WAIC approximation using posterior samples
#' @param log_lik Matrix of log-likelihood values (samples x observations)
#' @return List with WAIC and effective number of parameters
compute_waic <- function(log_lik) {
  # LPPD: log pointwise predictive density
  lppd <- sum(log(colMeans(exp(log_lik))))

  # Effective number of parameters (variance method)
  p_waic <- sum(apply(log_lik, 2, var))

  # WAIC
  waic <- -2 * (lppd - p_waic)

  return(list(
    WAIC = waic,
    lppd = lppd,
    p_waic = p_waic
  ))
}

#' Compute log-likelihood for PVL-Delta model
#' @param params List of parameter samples
#' @param data JAGS data list
#' @return Matrix of log-likelihood (samples x observations)
compute_loglik_pvl_delta <- function(params, data) {
  n_samples <- nrow(params$A)
  n_obs <- sum(data$Tsubj)

  log_lik <- matrix(NA, nrow = n_samples, ncol = n_obs)

  for (samp in 1:n_samples) {
    obs_idx <- 1

    for (s in 1:data$N) {
      # Initialize expected values
      ev <- rep(0, 4)

      for (t in 1:data$Tsubj[s]) {
        # Compute choice probabilities
        sens <- 3^params$cons[samp, s] - 1
        exp_util <- exp(sens * ev)
        p_choice <- exp_util / sum(exp_util)

        # Log-likelihood of observed choice
        obs_choice <- data$choice[s, t]
        log_lik[samp, obs_idx] <- log(p_choice[obs_choice] + 1e-10)

        # Update expected value
        outcome <- data$outcome[s, t]
        if (outcome >= 0) {
          util <- outcome^params$alpha[samp, s]
        } else {
          util <- -params$lambda[samp, s] * abs(outcome)^params$alpha[samp, s]
        }
        ev[obs_choice] <- ev[obs_choice] + params$A[samp, s] * (util - ev[obs_choice])

        obs_idx <- obs_idx + 1
      }
    }
  }

  return(log_lik)
}

#' Compute BIC (Bayesian Information Criterion)
#' @param log_lik Log-likelihood at maximum or mean of parameters
#' @param n_params Number of parameters
#' @param n_obs Number of observations
#' @return BIC value
compute_bic <- function(log_lik, n_params, n_obs) {
  bic <- -2 * log_lik + n_params * log(n_obs)
  return(bic)
}

#' Extract parameter counts for each model
#' @param model_name Name of the model
#' @param n_subj Number of subjects
#' @return Number of free parameters
get_n_params <- function(model_name, n_subj) {
  if (model_name == "pvl_delta") {
    # Group: 4 mu + 4 sigma
    # Subject: 4 params x n_subj
    return(8 + 4 * n_subj)
  } else if (model_name == "vse") {
    # Group: 8 mu + 8 sigma
    # Subject: 8 params x n_subj
    return(16 + 8 * n_subj)
  } else if (model_name == "orl") {
    # Group: 5 mu + 5 sigma
    # Subject: 5 params x n_subj
    return(10 + 5 * n_subj)
  } else {
    return(NA)
  }
}

#' Compare models using available information criteria
#' @param models_dir Directory with fitted models
#' @return Data frame with model comparison statistics
compare_models <- function(models_dir = "analysis/outputs") {
  source("analysis/utils/load_data.R")
  source("analysis/utils/prepare_jags_data.R")

  cat("=== Model Comparison ===\n\n")

  # Load data
  all_data <- load_all_igt_data()
  jags_data <- prepare_jags_data(all_data)

  n_obs <- sum(jags_data$Tsubj)
  n_subj <- jags_data$N

  # Find fitted models
  fit_files <- list.files(models_dir, pattern = "_fit\\.rds$", full.names = TRUE)

  comparison_results <- list()

  for (fit_file in fit_files) {
    model_name <- gsub("_fit\\.rds$", "", basename(fit_file))

    cat(sprintf("Processing %s...\n", model_name))

    # Load fit
    fit <- readRDS(fit_file)

    # Extract parameter samples
    samples_matrix <- as.matrix(fit$samples)

    # Get number of parameters
    n_params <- get_n_params(model_name, n_subj)

    # Compute log-likelihood (approximate using mean parameters)
    # This is simplified - full implementation would compute for each sample
    mean_loglik <- NA  # Placeholder

    comparison_results[[model_name]] <- data.frame(
      model = model_name,
      n_params = n_params,
      n_obs = n_obs,
      stringsAsFactors = FALSE
    )
  }

  # Combine results
  comparison_df <- do.call(rbind, comparison_results)
  rownames(comparison_df) <- NULL

  cat("\nModel Comparison Summary:\n")
  print(comparison_df)

  # Save results
  output_file <- file.path(models_dir, "model_comparison.rds")
  saveRDS(comparison_df, output_file)

  # Write markdown report
  write_comparison_report(comparison_df, models_dir)

  return(comparison_df)
}

#' Write model comparison report to markdown
#' @param comparison_df Data frame with comparison results
#' @param output_dir Output directory
write_comparison_report <- function(comparison_df, output_dir = "analysis/outputs") {
  output_file <- file.path(output_dir, "model_comparison.md")

  report <- c(
    "# Model Comparison Report",
    "",
    "## Overview",
    "",
    "This report summarizes the comparison of IGT decision-making models.",
    "",
    "## Models Compared",
    "",
    "1. **PVL-Delta**: Prospect-Valence Learning with delta-rule updating",
    "2. **VSE**: Value + Sequential Exploration (PVL-Delta with perseverance)",
    "3. **ORL**: Outcome Representation Learning",
    "",
    "## Model Statistics",
    "",
    "```",
    paste(capture.output(print(comparison_df)), collapse = "\n"),
    "```",
    "",
    "## Interpretation",
    "",
    "- **n_params**: Total number of free parameters (group-level + subject-level)",
    "- **n_obs**: Total number of observations (trials across all subjects)",
    "",
    "## Notes",
    "",
    "Model comparison metrics (DIC, WAIC, LOO) require additional computation.",
    "To compute these metrics, you can:",
    "",
    "1. Use `dic.samples()` from rjags during model fitting",
    "2. Implement WAIC using pointwise log-likelihood",
    "3. Use loo package for Pareto-smoothed importance sampling LOO-CV",
    "",
    sprintf("Generated: %s", Sys.time()),
    ""
  )

  writeLines(report, output_file)
  cat(sprintf("\nModel comparison report saved to: %s\n", output_file))
}

# Run comparison if executed as script
if (!interactive()) {
  comparison <- compare_models()
}
