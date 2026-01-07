# LOO-CV and Model Comparison for IGT Models
# Run: source("utils/loo_cv.R")
#
# Implements Leave-One-Out Cross-Validation using the 'loo' package
# for robust model comparison in hierarchical Bayesian models.

if (!require("pacman")) install.packages("pacman")
pacman::p_load(loo, coda)

# ==============================================================================
# POINTWISE LOG-LIKELIHOOD FUNCTIONS
# ==============================================================================

#' Compute pointwise log-likelihood for PVL-Delta model
#' @param fit JAGS fit object
#' @param data JAGS data list with x, X, ntrials, nsubs
#' @return Matrix of log-likelihood (n_samples x n_obs)
compute_loglik_pvl_delta <- function(fit, data) {
    samples <- fit$BUGSoutput$sims.list
    n_samples <- length(samples$mu_w)

    # Get subject-level parameters
    w <- samples$w # [n_samples, nsubs]
    A <- samples$A # [n_samples, nsubs]
    a <- samples$a # [n_samples, nsubs]
    theta <- samples$theta # [n_samples, nsubs]

    nsubs <- data$nsubs

    # Count total observations
    n_obs <- sum(data$ntrials)

    log_lik <- matrix(NA, nrow = n_samples, ncol = n_obs)

    for (samp in 1:n_samples) {
        obs_idx <- 1

        for (s in 1:nsubs) {
            # Initialize EV
            ev <- rep(0, 4)

            for (t in 1:data$ntrials[s]) {
                # Choice probabilities
                v <- theta[samp, s] * ev
                v <- v - max(v) # Numerical stability
                p_choice <- exp(v) / sum(exp(v))

                # Log-likelihood of observed choice
                obs_choice <- data$x[s, t]
                if (!is.na(obs_choice)) {
                    log_lik[samp, obs_idx] <- log(max(p_choice[obs_choice], 1e-10))
                    obs_idx <- obs_idx + 1
                }

                # Update EV (for next trial)
                if (t < data$ntrials[s]) {
                    outcome <- data$X[s, t]
                    if (!is.na(outcome)) {
                        # Utility
                        if (outcome >= 0) {
                            util <- abs(outcome)^A[samp, s]
                        } else {
                            util <- -w[samp, s] * abs(outcome)^A[samp, s]
                        }
                        # Delta learning rule
                        ev[obs_choice] <- ev[obs_choice] + a[samp, s] * (util - ev[obs_choice])
                    }
                }
            }
        }
    }

    return(log_lik[, 1:(obs_idx - 1)])
}

#' Compute pointwise log-likelihood for ORL model
#' @param fit JAGS fit object
#' @param data JAGS data list
#' @return Matrix of log-likelihood (n_samples x n_obs)
compute_loglik_orl <- function(fit, data) {
    samples <- fit$BUGSoutput$sims.list
    n_samples <- nrow(samples$a_rew)

    # Get subject-level parameters
    a_rew <- samples$a_rew
    a_pun <- samples$a_pun
    K <- samples$K
    theta <- samples$theta
    omega_f <- samples$omega_f
    omega_p <- samples$omega_p

    nsubs <- data$nsubs
    n_obs <- sum(data$ntrials)

    log_lik <- matrix(NA, nrow = n_samples, ncol = n_obs)

    for (samp in 1:n_samples) {
        obs_idx <- 1

        for (s in 1:nsubs) {
            # Initialize
            ev <- rep(0, 4)
            ef <- rep(0, 4)
            ps <- rep(0, 4)

            for (t in 1:data$ntrials[s]) {
                # Combined utility
                V <- ev + ef * omega_f[samp, s] + ps * omega_p[samp, s]

                # Choice probabilities
                v <- theta[samp, s] * V
                v <- v - max(v)
                p_choice <- exp(v) / sum(exp(v))

                # Log-likelihood
                obs_choice <- data$x[s, t]
                if (!is.na(obs_choice)) {
                    log_lik[samp, obs_idx] <- log(max(p_choice[obs_choice], 1e-10))
                    obs_idx <- obs_idx + 1
                }

                # Updates for next trial
                if (t < data$ntrials[s]) {
                    outcome <- data$X[s, t]
                    if (!is.na(outcome)) {
                        sign_out <- ifelse(outcome >= 0, 1, -1)
                        lr <- ifelse(outcome >= 0, a_rew[samp, s], a_pun[samp, s])
                        lr_rev <- ifelse(outcome >= 0, a_pun[samp, s], a_rew[samp, s])

                        # Update chosen
                        ev[obs_choice] <- ev[obs_choice] + lr * (outcome - ev[obs_choice])
                        ef[obs_choice] <- ef[obs_choice] + lr * (sign_out - ef[obs_choice])

                        # Update unchosen (fictive)
                        for (d in setdiff(1:4, obs_choice)) {
                            ef[d] <- ef[d] + lr_rev * (-sign_out / 3 - ef[d])
                        }

                        # Perseverance
                        ps[obs_choice] <- 1 / (1 + K[samp, s])
                        for (d in setdiff(1:4, obs_choice)) {
                            ps[d] <- ps[d] / (1 + K[samp, s])
                        }
                    }
                }
            }
        }
    }

    return(log_lik[, 1:(obs_idx - 1)])
}

#' Compute pointwise log-likelihood for EEF model
#' @param fit JAGS fit object
#' @param data JAGS data list
#' @return Matrix of log-likelihood (n_samples x n_obs)
compute_loglik_eef <- function(fit, data) {
    samples <- fit$BUGSoutput$sims.list
    n_samples <- nrow(samples$theta)

    # Get subject-level parameters
    p_theta <- samples$theta
    p_lambda <- samples$lambda
    p_phi <- samples$phi
    p_cons <- samples$cons

    nsubs <- data$nsubs
    n_obs <- sum(data$ntrials)

    log_lik <- matrix(NA, nrow = n_samples, ncol = n_obs)

    for (samp in 1:n_samples) {
        obs_idx <- 1

        for (s in 1:nsubs) {
            # Initialize
            exploit <- rep(0, 4)
            explore <- rep(0, 4)

            # Consistency transformation
            C <- 3^p_cons[samp, s] - 1

            for (t in 1:data$ntrials[s]) {
                # Value
                w <- exploit + explore
                v <- C * w
                v <- v - max(v)
                p_choice <- exp(v) / sum(exp(v))

                # Log-likelihood
                obs_choice <- data$x[s, t]
                if (!is.na(obs_choice)) {
                    log_lik[samp, obs_idx] <- log(max(p_choice[obs_choice], 1e-10))
                    obs_idx <- obs_idx + 1
                }

                # Updates for next trial
                if (t < data$ntrials[s]) {
                    outcome <- data$X[s, t]
                    if (!is.na(outcome)) {
                        # Utility (power function)
                        util <- ifelse(outcome >= 0,
                            abs(outcome)^p_theta[samp, s],
                            -abs(outcome)^p_theta[samp, s]
                        )

                        # Update exploit/explore
                        for (d in 1:4) {
                            if (d == obs_choice) {
                                exploit[d] <- (1 - p_lambda[samp, s]) * exploit[d] + p_lambda[samp, s] * util
                                explore[d] <- 0 # Reset
                            } else {
                                exploit[d] <- (1 - p_lambda[samp, s]) * exploit[d]
                                explore[d] <- (1 - p_lambda[samp, s]) * explore[d] + p_lambda[samp, s] * p_phi[samp, s]
                            }
                        }
                    }
                }
            }
        }
    }

    return(log_lik[, 1:(obs_idx - 1)])
}

# ==============================================================================
# LOO-CV COMPUTATION
# ==============================================================================

#' Compute LOO-CV for a fitted model
#' @param fit JAGS fit object
#' @param data JAGS data list
#' @param model_name One of "pvl_delta", "orl", "eef"
#' @return loo object from loo package
compute_loo <- function(fit, data, model_name) {
    cat(sprintf("\nComputing LOO-CV for %s...\n", model_name))

    # Compute log-likelihood
    if (model_name == "pvl_delta") {
        log_lik <- compute_loglik_pvl_delta(fit, data)
    } else if (model_name == "orl") {
        log_lik <- compute_loglik_orl(fit, data)
    } else if (model_name == "eef") {
        log_lik <- compute_loglik_eef(fit, data)
    } else {
        stop("Unknown model: ", model_name)
    }

    # Compute relative effective sample size
    r_eff <- loo::relative_eff(exp(log_lik))

    # Compute LOO
    loo_result <- loo::loo(log_lik, r_eff = r_eff)

    cat("\nLOO-CV Summary:\n")
    print(loo_result)

    return(loo_result)
}

#' Compare multiple models using LOO-CV
#' @param fits Named list of JAGS fit objects
#' @param data JAGS data list
#' @param model_names Vector of model names matching fits
#' @return loo_compare object
compare_models_loo <- function(fits, data, model_names = names(fits)) {
    cat("=== LOO-CV Model Comparison ===\n")

    loo_results <- list()

    for (i in seq_along(fits)) {
        model_name <- model_names[i]
        loo_results[[model_name]] <- compute_loo(fits[[i]], data, model_name)
    }

    # Compare
    cat("\n=== Model Comparison ===\n")
    comparison <- loo::loo_compare(loo_results)
    print(comparison)

    # Create summary table
    summary_df <- data.frame(
        model = rownames(comparison),
        elpd_diff = comparison[, "elpd_diff"],
        se_diff = comparison[, "se_diff"],
        elpd_loo = comparison[, "elpd_loo"],
        se_elpd_loo = comparison[, "se_elpd_loo"],
        p_loo = comparison[, "p_loo"],
        looic = -2 * comparison[, "elpd_loo"] # LOOIC for AIC-like interpretation
    )
    rownames(summary_df) <- NULL

    return(list(
        loo_results = loo_results,
        comparison = comparison,
        summary = summary_df
    ))
}

#' Run LOO-CV for all fitted models in a directory
#' @param models_dir Directory with fitted models (*_fit.rds)
#' @param output_file Output file for results
#' @return Model comparison results
run_all_loo <- function(models_dir = "analysis/outputs",
                        output_file = "analysis/outputs/loo_comparison.rds") {
    source("utils/load_data.R")

    # Load data
    all_data <- load_all_igt_data()

    # Find model fits
    fit_files <- list.files(models_dir, pattern = "_fit_.*\\.rds$", full.names = TRUE)

    if (length(fit_files) == 0) {
        cat("No fitted models found in", models_dir, "\n")
        cat("Looking for files matching: *_fit_*.rds\n")
        return(NULL)
    }

    # Group by model type
    model_groups <- list()

    for (fit_file in fit_files) {
        # Extract model name and group
        basename <- gsub("\\.rds$", "", basename(fit_file))
        parts <- strsplit(basename, "_fit_")[[1]]
        model_name <- parts[1]
        group_name <- parts[2]

        cat(sprintf("Found: %s (group: %s)\n", model_name, group_name))

        # Load and compute for each group separately
        fit <- readRDS(fit_file)

        # Prepare data for this group
        group_data <- all_data[all_data$study == paste0("Ahn2014_", group_name), ]
        jags_data <- prepare_group_jags_data(group_data)

        key <- paste0(model_name, "_", group_name)
        model_groups[[key]] <- list(
            fit = fit,
            data = jags_data,
            model = model_name,
            group = group_name
        )
    }

    # Save results
    saveRDS(model_groups, output_file)
    cat("\nResults saved to:", output_file, "\n")

    return(model_groups)
}

#' Prepare JAGS data for a group (helper)
#' @param group_data Data frame for one group
#' @return JAGS-formatted data list
prepare_group_jags_data <- function(group_data) {
    subIDs <- unique(group_data$subj)
    nsubs <- length(subIDs)
    ntrials_max <- 100

    ntrials_all <- array(0, c(nsubs))
    x_all <- array(0, c(nsubs, ntrials_max))
    X_all <- array(0, c(nsubs, ntrials_max))

    for (s in 1:nsubs) {
        subj_df <- group_data[group_data$subj == subIDs[s], ]
        ntrials_all[s] <- nrow(subj_df)

        x_sub <- subj_df$choice
        length(x_sub) <- ntrials_max
        x_all[s, ] <- x_sub

        X_raw <- subj_df$gain + subj_df$loss
        X_sub <- X_raw / 100
        length(X_sub) <- ntrials_max
        X_all[s, ] <- X_sub
    }

    return(list(
        x = x_all,
        X = X_all,
        ntrials = ntrials_all,
        nsubs = nsubs
    ))
}

# ==============================================================================
# USAGE EXAMPLE
# ==============================================================================
#
# # After fitting models:
# source("utils/loo_cv.R")
#
# # For single model
# fit <- readRDS("analysis/outputs/pvl_delta/pvl_delta_fit_HC.rds")
# data <- prepare_group_jags_data(hc_data)
# loo_result <- compute_loo(fit, data, "pvl_delta")
#
# # For model comparison (same group, different models)
# fits <- list(
#   pvl_delta = readRDS("analysis/outputs/pvl_delta/pvl_delta_fit_HC.rds"),
#   orl = readRDS("analysis/outputs/orl/orl_fit_HC.rds"),
#   eef = readRDS("analysis/outputs/eef/eef_fit_HC.rds")
# )
# comparison <- compare_models_loo(fits, data, c("pvl_delta", "orl", "eef"))
