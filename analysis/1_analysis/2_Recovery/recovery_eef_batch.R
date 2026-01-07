# ==============================================================================
# EEF Parameter Recovery - Batch Worker Script
# ==============================================================================
# This script runs a subset of recovery iterations.
# It is designed to be run in parallel via tmux or another scheduler.
#
# Usage:
#   Rscript recovery_eef_batch.R --seed 123 --iter 10 --output path/to/save.rds
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Parse Command Line Arguments
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# Default values
seed <- NULL
n_iter <- 5
output_file <- "recovery_batch.rds"

# Simple argument parser
for (i in seq_along(args)) {
    if (args[i] == "--seed" && i < length(args)) seed <- as.integer(args[i + 1])
    if (args[i] == "--iter" && i < length(args)) n_iter <- as.integer(args[i + 1])
    if (args[i] == "--output" && i < length(args)) output_file <- args[i + 1]
}

if (is.null(seed)) stop("Error: --seed argument is required.")

cat("----------------------------------------------------------------\n")
cat("EEF Recovery Batch Worker\n")
cat("Seed:", seed, "\n")
cat("Iterations:", n_iter, "\n")
cat("Output:", output_file, "\n")
cat("----------------------------------------------------------------\n")

# ------------------------------------------------------------------------------
# 2. Dependencies & Setup
# ------------------------------------------------------------------------------
# Load required packages (Base R style)
required_packages <- c("R2jags", "parallel", "extraDistr", "truncnorm")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")
invisible(lapply(required_packages, library, character.only = TRUE))

# Source model and helper functions
# (Assuming script is run from project root)
source("analysis/1_analysis/2_Recovery/simulation_eef.R")
source("analysis/utils/payoff_scheme.R")

# helper helper
MPD <- function(x) {
    if (all(is.na(x))) {
        return(NA)
    }
    density(x)$x[which.max(density(x)$y)]
}

# ------------------------------------------------------------------------------
# 3. Configuration
# ------------------------------------------------------------------------------
set.seed(seed)
nsubs <- 24 # Subjects per group (same as main recovery)
ntrials_all <- rep(100, nsubs)
payoff_struct <- generate_modified_igt_payoff(ntrials = 100, scale = TRUE)

# ------------------------------------------------------------------------------
# 4. Run Iterations
# ------------------------------------------------------------------------------
results_list <- list()

for (i in 1:n_iter) {
    cat(sprintf("[%d/%d] Generating and fitting...\n", i, n_iter))

    # A. Generate True Parameters
    mu_theta <- runif(1, 0.2, 0.8)
    sigma_theta <- runif(1, 0.05, 0.15)

    mu_lambda <- runif(1, 0.2, 0.8)
    sigma_lambda <- runif(1, 0.05, 0.15)

    mu_phi <- runif(1, -2, 2)
    sigma_phi <- runif(1, 0.1, 0.5)

    mu_cons <- runif(1, 1, 3)
    sigma_cons <- runif(1, 0.1, 0.5)

    # B. Simulate Data
    sim_data <- tryCatch(
        {
            simulation_eef(
                payoff_struct = payoff_struct,
                nsubs = nsubs,
                ntrials = ntrials_all,
                mu_theta = mu_theta, mu_lambda = mu_lambda,
                mu_phi = mu_phi, mu_cons = mu_cons,
                sigma_theta = sigma_theta, sigma_lambda = sigma_lambda,
                sigma_phi = sigma_phi, sigma_cons = sigma_cons
            )
        },
        error = function(e) {
            cat("  Simulation failed:", conditionMessage(e), "\n")
            return(NULL)
        }
    )

    if (is.null(sim_data)) next

    # C. Fit JAGS Model
    jags_data <- list(
        "x" = sim_data$x,
        "X" = sim_data$X,
        "ntrials" = ntrials_all,
        "nsubs" = nsubs
    )

    params <- c(
        "mu_theta", "mu_lambda", "mu_phi", "mu_cons",
        "lambda_theta", "lambda_lambda", "lambda_phi", "lambda_cons"
    )

    # Use tryCatch for JAGS fitting
    fit_result <- tryCatch(
        {
            jags(
                data = jags_data,
                inits = NULL,
                parameters.to.save = params,
                model.file = "analysis/models/eef.txt",
                n.chains = 3,
                n.iter = 3000,
                n.burnin = 1000,
                n.thin = 1,
                progress.bar = "none" # keep log clean
            )
        },
        error = function(e) {
            cat("  JAGS fit failed:", conditionMessage(e), "\n")
            return(NULL)
        }
    )

    if (!is.null(fit_result)) {
        # D. Extract Estimates
        Y <- fit_result$BUGSoutput$sims.list

        # Helper to extract quantiles
        get_q <- function(param_name, q_col) {
            if (param_name %in% rownames(fit_result$BUGSoutput$summary)) {
                return(fit_result$BUGSoutput$summary[param_name, q_col])
            } else {
                return(NA)
            }
        }

        res <- list(
            # True Means
            true_mu_theta = mu_theta, true_mu_lambda = mu_lambda,
            true_mu_phi = mu_phi, true_mu_cons = mu_cons,

            # True Sigmas
            true_sigma_theta = sigma_theta, true_sigma_lambda = sigma_lambda,
            true_sigma_phi = sigma_phi, true_sigma_cons = sigma_cons,

            # Inferred Means
            infer_mu_theta = MPD(Y$mu_theta), infer_mu_lambda = MPD(Y$mu_lambda),
            infer_mu_phi = MPD(Y$mu_phi), infer_mu_cons = MPD(Y$mu_cons),

            # Inferred Sigmas (converted from precision lambda)
            infer_sigma_theta = MPD(1 / sqrt(Y$lambda_theta)),
            infer_sigma_lambda = MPD(1 / sqrt(Y$lambda_lambda)),
            infer_sigma_phi = MPD(1 / sqrt(Y$lambda_phi)),
            infer_sigma_cons = MPD(1 / sqrt(Y$lambda_cons)),

            # 95% CI Lower (2.5%)
            lower_mu_theta = get_q("mu_theta", "2.5%"),
            lower_mu_lambda = get_q("mu_lambda", "2.5%"),
            lower_mu_phi = get_q("mu_phi", "2.5%"),
            lower_mu_cons = get_q("mu_cons", "2.5%"),

            # 95% CI Upper (97.5%)
            upper_mu_theta = get_q("mu_theta", "97.5%"),
            upper_mu_lambda = get_q("mu_lambda", "97.5%"),
            upper_mu_phi = get_q("mu_phi", "97.5%"),
            upper_mu_cons = get_q("mu_cons", "97.5%")
        )

        results_list[[length(results_list) + 1]] <- res
        cat("  Success.\n")
    }
}

# ------------------------------------------------------------------------------
# 5. Save Results
# ------------------------------------------------------------------------------
# Ensure directory exists
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

if (length(results_list) > 0) {
    saveRDS(results_list, file = output_file)
    cat("Saved", length(results_list), "results to", output_file, "\n")
} else {
    cat("No successful iterations to save.\n")
}
