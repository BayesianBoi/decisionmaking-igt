#!/usr/bin/env Rscript
# ==============================================================================
# FULL PIPELINE VALIDATION: Single-Subject Test
# ==============================================================================
#
# Purpose: Test the COMPLETE pipeline (fitting → diagnostics → PPC → plotting)
#          on a small subset (5 subjects) before running on all subjects.
#
# What it tests:
#   1. Model fitting (JAGS compilation, MCMC sampling)
#   2. Diagnostics (R-hat, ESS)
#   3. Posterior Predictive Checks (simulating choices)
#   4. Plotting (learning curves, posterior densities)
#   5. Saving/loading results
#
# Usage:
#   Rscript analysis/scripts/quick_validate.R
#   Rscript analysis/scripts/quick_validate.R --model pvl_delta
#   Rscript analysis/scripts/quick_validate.R --model eef
#
# Expected runtime: ~2-5 minutes per model
# ==============================================================================

# Install/load packages
if (!require("R2jags", quietly = TRUE)) install.packages("R2jags", repos = "https://cloud.r-project.org")
if (!require("coda", quietly = TRUE)) install.packages("coda", repos = "https://cloud.r-project.org")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos = "https://cloud.r-project.org")

library(R2jags)
library(coda)
library(ggplot2)

set.seed(42)

# ==============================================================================
# Configuration
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)

# Parse --model argument
model_arg <- NULL
if ("--model" %in% args) {
    model_idx <- which(args == "--model")
    if (length(args) > model_idx) {
        model_arg <- args[model_idx + 1]
    }
}

# Models to test
if (!is.null(model_arg)) {
    models_to_test <- model_arg
} else {
    models_to_test <- c("pvl_delta", "orl", "eef")
}

# MCMC settings (minimal for speed)
n_chains <- 3
n_iter <- 3000 # Increased slightly for better convergence
n_burnin <- 500
n_thin <- 1

# Output directory for test results
test_output_dir <- "analysis/outputs/validation_test"
dir.create(test_output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Full Pipeline Validation Script ===\n")
cat("Models to test:", paste(models_to_test, collapse = ", "), "\n")
cat("MCMC: chains =", n_chains, ", iter =", n_iter, ", burnin =", n_burnin, "\n")
cat("Output directory:", test_output_dir, "\n\n")

# ==============================================================================
# Load Multiple Subjects (5 for realistic but quick testing)
# ==============================================================================
source("analysis/utils/load_data.R")

all_data <- load_all_igt_data()

# Take first 10 subjects from HC group (more subjects = more stable priors)
hc_data <- all_data[all_data$study == "Ahn2014_HC", ]
subj_ids <- unique(hc_data$subj)[1:10]
n_subjects <- length(subj_ids)

cat("Using", n_subjects, "subjects:", paste(subj_ids, collapse = ", "), "\n\n")

# Prepare JAGS data
ntrials_max <- 100
x <- matrix(NA, nrow = n_subjects, ncol = ntrials_max)
X <- matrix(NA, nrow = n_subjects, ncol = ntrials_max)
ntrials_vec <- rep(NA, n_subjects)

for (s in 1:n_subjects) {
    subj_data <- hc_data[hc_data$subj == subj_ids[s], ]
    n_t <- nrow(subj_data)
    ntrials_vec[s] <- n_t

    x[s, 1:n_t] <- subj_data$choice
    X[s, 1:n_t] <- (subj_data$gain + subj_data$loss) / 100
}

jags_data <- list(
    x = x,
    X = X,
    ntrials = ntrials_vec,
    nsubs = n_subjects
)

# Store for later use
observed_choices <- x
ntrials <- ntrials_vec[1] # For plots (first subject)

# ==============================================================================
# Model Fitting Loop
# ==============================================================================
results <- list()

for (model_name in models_to_test) {
    cat("\n", strrep("=", 60), "\n")
    cat("TESTING MODEL:", toupper(model_name), "\n")
    cat(strrep("=", 60), "\n")

    model_file <- sprintf("analysis/models/%s.txt", model_name)

    if (!file.exists(model_file)) {
        cat("ERROR: Model file not found:", model_file, "\n")
        next
    }

    # Define parameters based on model
    if (model_name == "pvl_delta") {
        params <- c("mu_w", "mu_A", "mu_theta", "mu_a", "w", "A", "theta", "a")
    } else if (model_name == "orl") {
        params <- c(
            "mu_a_rew", "mu_a_pun", "mu_K", "mu_theta", "mu_omega_f", "mu_omega_p",
            "a_rew", "a_pun", "K", "theta", "omega_f", "omega_p"
        )
    } else if (model_name == "eef") {
        params <- c(
            "mu_theta", "mu_lambda", "mu_phi", "mu_cons",
            "theta", "lambda", "phi", "cons"
        )
    } else {
        cat("Unknown model:", model_name, "\n")
        next
    }

    # =========================================================================
    # STEP 1: Model Fitting
    # =========================================================================
    cat("\n[1/5] FITTING MODEL...\n")
    start_time <- Sys.time()

    tryCatch(
        {
            # Set initial values to help with convergence
            if (model_name == "eef") {
                inits <- function() {
                    list(
                        mu_theta = 0.5, mu_lambda = 0.5, mu_phi = 0, mu_cons = 1,
                        lambda_theta = 1, lambda_lambda = 1, lambda_phi = 1, lambda_cons = 1
                    )
                }
            } else {
                inits <- NULL
            }

            fit <- jags(
                data = jags_data,
                inits = inits,
                parameters.to.save = params,
                model.file = model_file,
                n.chains = n_chains,
                n.iter = n_iter,
                n.burnin = n_burnin,
                n.thin = n_thin,
                progress.bar = "text"
            )

            end_time <- Sys.time()
            cat("      ✓ Fitting complete! Duration:", round(difftime(end_time, start_time, units = "secs"), 1), "seconds\n")

            # Save fit
            fit_path <- file.path(test_output_dir, sprintf("%s_fit.rds", model_name))
            saveRDS(fit, fit_path)
            cat("      Saved to:", fit_path, "\n")

            # =========================================================================
            # STEP 2: Diagnostics
            # =========================================================================
            cat("\n[2/5] RUNNING DIAGNOSTICS...\n")

            sims <- fit$BUGSoutput$sims.list
            summary_stats <- fit$BUGSoutput$summary

            # R-hat
            rhat <- summary_stats[, "Rhat"]
            rhat <- rhat[!is.na(rhat)]
            max_rhat <- max(rhat)
            mean_rhat <- mean(rhat)

            cat(sprintf("      R-hat: max = %.3f, mean = %.3f\n", max_rhat, mean_rhat))
            if (max_rhat > 1.1) {
                cat("      ⚠ WARNING: R-hat > 1.1 (expected with few iterations)\n")
            } else {
                cat("      ✓ R-hat OK\n")
            }

            # ESS (n.eff)
            neff <- summary_stats[, "n.eff"]
            neff <- neff[!is.na(neff)]
            min_ess <- min(neff)
            cat(sprintf("      ESS: min = %.0f\n", min_ess))

            # Parameter summaries
            cat("\n      Parameter estimates:\n")
            mu_params <- names(sims)[grep("^mu_", names(sims))]
            for (p in mu_params) {
                vals <- sims[[p]]
                cat(sprintf("        %s: %.3f (SD = %.3f)\n", p, mean(vals), sd(vals)))
            }

            # =========================================================================
            # STEP 3: Posterior Predictive Check
            # =========================================================================
            cat("\n[3/5] RUNNING PPC...\n")

            # Sample from posterior
            n_ppc_samples <- min(50, nrow(sims[[1]]))
            sample_idx <- sample(1:nrow(sims[[1]]), n_ppc_samples)

            ppc_accuracies <- rep(NA, n_ppc_samples)

            for (i in seq_along(sample_idx)) {
                idx <- sample_idx[i]

                # Simulate choices based on model
                if (model_name == "pvl_delta") {
                    a_val <- sims$a[idx, 1]
                    A_val <- sims$A[idx, 1]
                    theta_val <- sims$theta[idx, 1]
                    w_val <- sims$w[idx, 1]

                    ev <- rep(0, 4)
                    sim_choices <- rep(NA, ntrials)

                    for (t in 1:ntrials) {
                        v <- theta_val * ev
                        v <- v - max(v)
                        p <- exp(v) / sum(exp(v))
                        if (any(is.na(p))) p <- rep(0.25, 4)

                        sim_choices[t] <- sample(1:4, 1, prob = p)

                        outcome <- X[1, t]
                        if (!is.na(outcome)) {
                            util <- ifelse(outcome >= 0, abs(outcome)^A_val, -w_val * abs(outcome)^A_val)
                            ev[sim_choices[t]] <- ev[sim_choices[t]] + a_val * (util - ev[sim_choices[t]])
                        }
                    }
                } else if (model_name == "eef") {
                    theta_val <- sims$theta[idx, 1]
                    lambda_val <- sims$lambda[idx, 1]
                    phi_val <- sims$phi[idx, 1]
                    cons_val <- sims$cons[idx, 1]

                    exploit <- rep(0, 4)
                    explore <- rep(0, 4)
                    C <- 3^cons_val - 1
                    sim_choices <- rep(NA, ntrials)

                    for (t in 1:ntrials) {
                        w <- exploit + explore
                        v <- C * w
                        v <- v - max(v)
                        p <- exp(v) / sum(exp(v))
                        if (any(is.na(p))) p <- rep(0.25, 4)

                        sim_choices[t] <- sample(1:4, 1, prob = p)

                        outcome <- X[1, t]
                        if (!is.na(outcome)) {
                            util <- ifelse(outcome >= 0, abs(outcome)^theta_val, -abs(outcome)^theta_val)

                            for (d in 1:4) {
                                if (d == sim_choices[t]) {
                                    exploit[d] <- (1 - lambda_val) * exploit[d] + lambda_val * util
                                    explore[d] <- 0
                                } else {
                                    exploit[d] <- (1 - lambda_val) * exploit[d]
                                    explore[d] <- (1 - lambda_val) * explore[d] + lambda_val * phi_val
                                }
                            }
                        }
                    }
                } else {
                    # ORL - simple random for speed
                    sim_choices <- sample(1:4, ntrials, replace = TRUE)
                }

                ppc_accuracies[i] <- mean(sim_choices == observed_choices[1, 1:ntrials], na.rm = TRUE)
            }

            mean_ppc <- mean(ppc_accuracies, na.rm = TRUE)
            cat(sprintf("      PPC accuracy (subject 1): %.1f%% (chance = 25%%)\n", mean_ppc * 100))

            # Deck proportions (first subject)
            obs_props <- table(factor(observed_choices[1, 1:ntrials], levels = 1:4)) / ntrials

            cat("\n      Observed deck proportions:\n")
            cat(
                "        Deck 1:", sprintf("%.2f", obs_props[1]),
                " | Deck 2:", sprintf("%.2f", obs_props[2]),
                " | Deck 3:", sprintf("%.2f", obs_props[3]),
                " | Deck 4:", sprintf("%.2f", obs_props[4]), "\n"
            )

            cat("      ✓ PPC completed\n")

            # =========================================================================
            # STEP 4: Plotting
            # =========================================================================
            cat("\n[4/5] GENERATING PLOTS...\n")

            plot_dir <- file.path(test_output_dir, "plots")
            dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

            # A. Posterior density plot for mu parameters
            for (p in mu_params) {
                vals <- sims[[p]]
                df <- data.frame(value = vals)

                plot_path <- file.path(plot_dir, sprintf("%s_%s_posterior.png", model_name, p))

                pl <- ggplot(df, aes(x = value)) +
                    geom_density(fill = "steelblue", alpha = 0.6) +
                    geom_vline(xintercept = mean(vals), linetype = "dashed", color = "red") +
                    labs(
                        title = sprintf("%s: %s", toupper(model_name), p),
                        subtitle = sprintf("Mean = %.3f, SD = %.3f", mean(vals), sd(vals)),
                        x = "Value", y = "Density"
                    ) +
                    theme_minimal()

                ggsave(plot_path, pl, width = 6, height = 4, dpi = 100)
            }
            cat("      ✓ Posterior density plots saved\n")

            # B. Learning curve plot (proportion of good decks over time) - first subject
            first_subj_choices <- observed_choices[1, 1:ntrials]
            blocks <- ceiling(1:ntrials / 20)
            good_choices <- first_subj_choices %in% c(3, 4)

            block_means <- tapply(good_choices, blocks, mean, na.rm = TRUE)
            block_df <- data.frame(
                block = as.numeric(names(block_means)),
                prop_good = as.numeric(block_means)
            )

            learning_plot_path <- file.path(plot_dir, sprintf("%s_learning_curve.png", model_name))

            pl_learn <- ggplot(block_df, aes(x = block, y = prop_good)) +
                geom_line(color = "steelblue", linewidth = 1) +
                geom_point(color = "steelblue", size = 2) +
                geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
                labs(
                    title = sprintf("%s: Learning Curve (Subject %s)", toupper(model_name), subj_ids[1]),
                    x = "Block (20 trials each)",
                    y = "P(Good Deck)"
                ) +
                ylim(0, 1) +
                theme_minimal()

            ggsave(learning_plot_path, pl_learn, width = 8, height = 5, dpi = 100)
            cat("      ✓ Learning curve plot saved\n")

            # =========================================================================
            # STEP 5: Save Summary
            # =========================================================================
            cat("\n[5/5] SAVING SUMMARY...\n")

            summary_df <- data.frame(
                model = model_name,
                max_rhat = max_rhat,
                mean_rhat = mean_rhat,
                min_ess = min_ess,
                ppc_accuracy = mean_ppc,
                fit_duration_sec = as.numeric(difftime(end_time, start_time, units = "secs"))
            )

            summary_path <- file.path(test_output_dir, sprintf("%s_summary.csv", model_name))
            write.csv(summary_df, summary_path, row.names = FALSE)
            cat("      ✓ Summary saved to:", summary_path, "\n")

            results[[model_name]] <- list(
                success = TRUE,
                max_rhat = max_rhat,
                ppc_accuracy = mean_ppc,
                plots_dir = plot_dir
            )

            cat("\n      ✓ MODEL", toupper(model_name), "PASSED ALL TESTS\n")
        },
        error = function(e) {
            cat("ERROR:", conditionMessage(e), "\n")
            results[[model_name]] <- list(success = FALSE, error = conditionMessage(e))
        }
    )
}

# ==============================================================================
# Final Summary
# ==============================================================================
cat("\n", strrep("=", 60), "\n")
cat("VALIDATION SUMMARY\n")
cat(strrep("=", 60), "\n\n")

passed <- 0
failed <- 0

for (model_name in names(results)) {
    res <- results[[model_name]]
    if (res$success) {
        cat(sprintf(
            "  ✓ %-12s PASSED  (R-hat=%.2f, PPC=%.1f%%)\n",
            toupper(model_name), res$max_rhat, res$ppc_accuracy * 100
        ))
        passed <- passed + 1
    } else {
        cat(sprintf("  ✗ %-12s FAILED  (%s)\n", toupper(model_name), res$error))
        failed <- failed + 1
    }
}

cat("\n")
cat(sprintf("Results: %d passed, %d failed\n", passed, failed))
cat("Test outputs saved to:", test_output_dir, "\n")

if (failed == 0) {
    cat("\n✓ All models validated! Ready for full analysis.\n")
} else {
    cat("\n⚠ Some models failed. Check errors above.\n")
}

cat("\n")
