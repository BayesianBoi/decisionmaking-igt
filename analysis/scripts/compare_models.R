# Compare PVL-Delta, ORL, and EEF models
# Determines which model best accounts for the clinical IGT data
#
# Run: Rscript analysis/scripts/compare_models.R

library(coda)
library(loo)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Model directories
model_dirs <- c(
  pvl_delta = "results/pvl_delta",
  orl = "results/orl",
  eef = "results/eef"
)

# Output directory
output_dir <- "results/model_comparison"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Compute log-likelihood for each observation
#' @param samples MCMC samples object
#' @param jags_data Data used for fitting
#' @return Matrix of log-likelihoods (iterations x observations)
compute_log_lik <- function(samples, jags_data) {
  # Extract choice probabilities from samples
  # For hierarchical models, we need to reconstruct p[s,t,choice[s,t]]

  # Get samples as matrix
  samples_mat <- as.matrix(samples)
  n_iter <- nrow(samples_mat)

  # Count total observations
  n_obs <- sum(jags_data$Tsubj)

  # Initialize log-likelihood matrix
  log_lik <- matrix(NA, nrow = n_iter, ncol = n_obs)

  message(sprintf(
    "Computing log-likelihood for %d observations across %d iterations",
    n_obs, n_iter
  ))
  message("This may take several minutes...")

  # For each iteration, compute log P(choice | parameters)
  # This requires re-running the model forward pass for each posterior sample
  # Since JAGS doesn't save trial-by-trial likelihoods, we need to compute them

  # NOTE: This is computationally expensive. For a simpler approach,
  # we can use DIC from JAGS directly, or implement WAIC/LOO using
  # approximations based on deviance

  warning("Full log-likelihood computation not yet implemented.")
  warning("Using DIC from JAGS diagnostics instead.")

  return(NULL)
}

#' Compute WAIC and other metrics
#' @param model_dirs Named vector of model directories
#' @return Data frame with comparison metrics
compare_models_waic <- function(model_dirs) {
  results <- data.frame(
    model = names(model_dirs),
    n_params = NA,
    mean_rhat = NA,
    median_ess = NA,
    waic = NA,
    waic_se = NA,
    loo = NA,
    loo_se = NA,
    p_loo = NA,
    pareto_k_good = NA, # Proportion of k < 0.7
    dic = NA,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(model_dirs)) {
    model_name <- names(model_dirs)[i]
    model_dir <- model_dirs[i]

    message(sprintf("\nProcessing %s...", model_name))

    # Load diagnostics
    diag_file <- file.path(model_dir, "diagnostics.rds")
    if (file.exists(diag_file)) {
      diag <- readRDS(diag_file)
      rhat_vals <- diag$rhat$psrf[, "Point est."]
      ess_vals <- diag$eff_size

      results$mean_rhat[i] <- mean(rhat_vals, na.rm = TRUE)
      results$median_ess[i] <- median(ess_vals, na.rm = TRUE)
    }

    # Load samples for WAIC and LOO
    samples_file <- file.path(model_dir, "mcmc_samples.rds")
    if (file.exists(samples_file)) {
      samples <- readRDS(samples_file)
      results$n_params[i] <- ncol(as.matrix(samples[, !grepl("^log_lik", colnames(samples[[1]]))]))

      # Extract log-likelihoods
      message("  Extracting log-likelihoods...")
      mat <- as.matrix(samples)
      ll_cols <- grep("^log_lik", colnames(mat))

      if (length(ll_cols) > 0) {
        log_lik_mat <- mat[, ll_cols]
        message(sprintf("  Computing WAIC and LOO for %d obs x %d draws...", ncol(log_lik_mat), nrow(log_lik_mat)))

        # =====================================================================
        # WAIC (Widely Applicable Information Criterion)
        # Lower is better. SE allows pairwise comparison confidence.
        # =====================================================================
        waic_res <- tryCatch(
          {
            loo::waic(log_lik_mat)
          },
          error = function(e) {
            warning(sprintf("WAIC computation failed for %s: %s", model_name, e$message))
            return(NULL)
          }
        )

        if (!is.null(waic_res)) {
          results$waic[i] <- waic_res$estimates["waic", "Estimate"]
          results$waic_se[i] <- waic_res$estimates["waic", "SE"]
          message(sprintf("  WAIC: %.2f (SE: %.2f)", results$waic[i], results$waic_se[i]))
        }

        # =====================================================================
        # LOO-CV (Leave-One-Out Cross-Validation via PSIS)
        # More robust than WAIC for small effective sample sizes.
        # Pareto-k diagnostics indicate reliability:
        #   k < 0.5: Excellent
        #   0.5 < k < 0.7: Good
        #   k > 0.7: Problematic (observation is influential)
        # =====================================================================
        loo_res <- tryCatch(
          {
            loo::loo(log_lik_mat)
          },
          error = function(e) {
            warning(sprintf("LOO computation failed for %s: %s", model_name, e$message))
            return(NULL)
          }
        )

        if (!is.null(loo_res)) {
          results$loo[i] <- loo_res$estimates["looic", "Estimate"]
          results$loo_se[i] <- loo_res$estimates["looic", "SE"]
          results$p_loo[i] <- loo_res$estimates["p_loo", "Estimate"]

          # Pareto-k diagnostics
          pareto_k <- loo_res$diagnostics$pareto_k
          prop_good <- mean(pareto_k < 0.7, na.rm = TRUE)
          results$pareto_k_good[i] <- prop_good

          message(sprintf("  LOO: %.2f (SE: %.2f)", results$loo[i], results$loo_se[i]))
          message(sprintf("  Pareto-k < 0.7: %.1f%% of observations", prop_good * 100))

          # Warn if many problematic observations
          if (prop_good < 0.95) {
            warning(sprintf(
              "  %s: %.1f%% of Pareto-k values > 0.7 (may indicate model misfit)",
              model_name, (1 - prop_good) * 100
            ))
          }
        }
      } else {
        warning("No 'log_lik' parameters found in samples.")
      }
    }

    # DIC (Legacy check)
    results$dic[i] <- NA
  }

  return(results)
}

#' Extract parameter estimates for comparison
#' @param model_dirs Named vector of model directories
#' @return Data frame with parameter estimates
extract_parameter_estimates <- function(model_dirs) {
  param_list <- list()

  for (model_name in names(model_dirs)) {
    model_dir <- model_dirs[model_name]

    summary_file <- file.path(model_dir, "parameter_summary.rds")
    if (!file.exists(summary_file)) {
      warning(sprintf("Parameter summary not found for %s", model_name))
      next
    }

    summary_obj <- readRDS(summary_file)

    # Extract group-level parameters only
    stats <- summary_obj$statistics
    group_params <- grep("^mu_", rownames(stats), value = TRUE)

    param_df <- data.frame(
      model = model_name,
      parameter = group_params,
      mean = stats[group_params, "Mean"],
      sd = stats[group_params, "SD"],
      stringsAsFactors = FALSE
    )

    param_list[[model_name]] <- param_df
  }

  all_params <- do.call(rbind, param_list)
  rownames(all_params) <- NULL

  return(all_params)
}

# ==============================================================================
# RUN COMPARISON
# ==============================================================================

message("=== MODEL COMPARISON ===\n")

# Check that all model fits exist
for (model_name in names(model_dirs)) {
  model_dir <- model_dirs[model_name]
  samples_file <- file.path(model_dir, "mcmc_samples.rds")

  if (!file.exists(samples_file)) {
    stop(sprintf(
      "Model fit not found for %s. Run fit_%s.R first.",
      model_name, model_name
    ))
  }
}

message("All model fits found.\n")

# ==============================================================================
# CONVERGENCE COMPARISON
# ==============================================================================

message("=== CONVERGENCE STATISTICS ===\n")

comparison_df <- compare_models_waic(model_dirs)

print(comparison_df)

# Save comparison table
write.csv(comparison_df, file.path(output_dir, "model_comparison_table.csv"),
  row.names = FALSE
)
message(sprintf(
  "\nComparison table saved to: %s",
  file.path(output_dir, "model_comparison_table.csv")
))

# ==============================================================================
# PARAMETER ESTIMATES
# ==============================================================================

message("\n=== PARAMETER ESTIMATES ===\n")

param_estimates <- extract_parameter_estimates(model_dirs)

print(param_estimates)

# Save parameter estimates
write.csv(param_estimates, file.path(output_dir, "parameter_estimates.csv"),
  row.names = FALSE
)
message(sprintf(
  "\nParameter estimates saved to: %s",
  file.path(output_dir, "parameter_estimates.csv")
))

# ==============================================================================
# VISUALIZATIONS
# ==============================================================================

message("\n=== CREATING COMPARISON PLOTS ===\n")

# Convergence comparison plot
pdf(file.path(output_dir, "convergence_comparison.pdf"), width = 10, height = 6)
par(mfrow = c(1, 2))

# R-hat comparison
barplot(comparison_df$mean_rhat,
  names.arg = comparison_df$model,
  main = "Mean R-hat by Model",
  ylab = "Mean R-hat",
  ylim = c(0.9, 1.2),
  col = "steelblue"
)
abline(h = 1.1, col = "red", lty = 2, lwd = 2)
text(x = 1.5, y = 1.15, labels = "Threshold = 1.1", col = "red")

# ESS comparison
barplot(comparison_df$median_ess,
  names.arg = comparison_df$model,
  main = "Median Effective Sample Size",
  ylab = "Median ESS",
  col = "steelblue"
)
abline(h = 1000, col = "red", lty = 2, lwd = 2)
text(x = 1.5, y = 1500, labels = "Target = 1000", col = "red")

dev.off()

message(sprintf(
  "Convergence plots saved to: %s",
  file.path(output_dir, "convergence_comparison.pdf")
))

# Parameter comparison plot (for common parameters)
common_params <- c("mu_cons") # cons is in all three models

pdf(file.path(output_dir, "parameter_comparison.pdf"), width = 10, height = 6)

for (param in common_params) {
  param_subset <- param_estimates[param_estimates$parameter == param, ]

  if (nrow(param_subset) > 0) {
    # Create barplot with error bars
    bp <- barplot(param_subset$mean,
      names.arg = param_subset$model,
      main = sprintf("Comparison: %s", param),
      ylab = "Posterior Mean",
      ylim = range(c(
        param_subset$mean - param_subset$sd,
        param_subset$mean + param_subset$sd
      )),
      col = "steelblue"
    )

    # Add error bars
    arrows(
      x0 = bp,
      y0 = param_subset$mean - param_subset$sd,
      x1 = bp,
      y1 = param_subset$mean + param_subset$sd,
      angle = 90, code = 3, length = 0.1, lwd = 2
    )
  }
}

dev.off()

message(sprintf(
  "Parameter plots saved to: %s",
  file.path(output_dir, "parameter_comparison.pdf")
))

# ==============================================================================
# SUMMARY REPORT
# ==============================================================================

message("\n=== SUMMARY ===\n")
message("Model comparison complete.\n")
message("Files created:")
message("  - model_comparison_table.csv")
message("  - parameter_estimates.csv")
message("  - convergence_comparison.pdf")
message("  - parameter_comparison.pdf")
message("\nNOTE: For full model comparison with WAIC/LOO-CV:")
message("  1. Modify JAGS models to save log-likelihoods")
message("  2. Re-run fits with dic.samples() for DIC")
message("  3. Use loo package for WAIC/LOO computation")
message("\nCurrent comparison is based on:")
message("  - Convergence diagnostics (R-hat, ESS)")
message("  - Parameter count (complexity)")
message("  - Parameter estimates")
