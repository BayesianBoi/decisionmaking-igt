# PVL-Delta Parameter Recovery - Batch Worker Script
# This script runs a subset of recovery iterations for the PVL-Delta model.
# It is designed to be run in parallel via tmux or another scheduler.
#
# Usage:
#   Rscript recovery_pvl_delta_batch.R --seed 123 --iter 10 --output path/to/save.rds

# parse args
args <- commandArgs(trailingOnly = TRUE)

# Default values
seed <- NULL
n_iter <- 5
output_file <- "recovery_pvl_delta_batch.rds"

# Simple argument parser
for (i in seq_along(args)) {
    if (args[i] == "--seed" && i < length(args)) seed <- as.integer(args[i + 1])
    if (args[i] == "--iter" && i < length(args)) n_iter <- as.integer(args[i + 1])
    if (args[i] == "--output" && i < length(args)) output_file <- args[i + 1]
}

if (is.null(seed)) stop("Error: --seed argument is required.")

cat("----------------------------------------------------------------\n")
cat("PVL-Delta Recovery Batch Worker\n")
cat("Seed:", seed, "\n")
cat("Iterations:", n_iter, "\n")
cat("Output:", output_file, "\n")
cat("----------------------------------------------------------------\n")

# dependencies
# Load required packages (Base R style)
required_packages <- c("R2jags", "parallel", "extraDistr", "truncnorm")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")
invisible(lapply(required_packages, library, character.only = TRUE))

# Source model and helper functions
# (Assuming script is run from project root)
source("scripts/recovery/simulation_pvl_delta.R")
source("utils/payoff_scheme.R")

# MPD helper (Maximum Posterior Density)
MPD <- function(x) {
    if (all(is.na(x))) {
        return(NA)
    }
    density(x)$x[which.max(density(x)$y)]
}

# config
set.seed(seed)
nsubs <- 24 # Subjects per group (same as main recovery)
ntrials_all <- rep(100, nsubs)
payoff_struct <- generate_modified_igt_payoff(ntrials = 100, scale = TRUE)

# main loop
results_list <- list()

for (i in 1:n_iter) {
    cat(sprintf("[%d/%d] Generating and fitting...\n", i, n_iter))

    # A. Generate True Parameters
    mu_w <- runif(1, 1, 3)
    sigma_w <- runif(1, 0.2, 0.5)

    mu_A <- runif(1, 0.3, 0.8)
    sigma_A <- runif(1, 0.1, 0.3)

    mu_a <- runif(1, 0.1, 0.5)
    sigma_a <- runif(1, 0.05, 0.15)

    mu_theta <- runif(1, 0.5, 2.0)
    sigma_theta <- runif(1, 0.1, 0.3)

    # B. Simulate Data
    sim_data <- tryCatch(
        {
            hier_PVL_sim(
                payoff_struct = payoff_struct,
                nsubs = nsubs,
                ntrials = ntrials_all,
                mu_w = mu_w, mu_A = mu_A, mu_a = mu_a, mu_theta = mu_theta,
                sigma_w = sigma_w, sigma_A = sigma_A, sigma_a = sigma_a, sigma_theta = sigma_theta
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
        "mu_w", "mu_A", "mu_a", "mu_theta",
        "lambda_w", "lambda_A", "lambda_a", "lambda_theta"
    )

    # Use tryCatch for JAGS fitting
    fit_result <- tryCatch(
        {
            jags(
                data = jags_data,
                inits = NULL,
                parameters.to.save = params,
                model.file = "models/pvl_delta.txt",
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
            true_mu_w = mu_w, true_mu_A = mu_A,
            true_mu_a = mu_a, true_mu_theta = mu_theta,

            # True Sigmas
            true_sigma_w = sigma_w, true_sigma_A = sigma_A,
            true_sigma_a = sigma_a, true_sigma_theta = sigma_theta,

            # Inferred Means (MPD)
            infer_mu_w = MPD(Y$mu_w), infer_mu_A = MPD(Y$mu_A),
            infer_mu_a = MPD(Y$mu_a), infer_mu_theta = MPD(Y$mu_theta),

            # Inferred Sigmas (converted from precision lambda)
            infer_sigma_w = MPD(1 / sqrt(Y$lambda_w)),
            infer_sigma_A = MPD(1 / sqrt(Y$lambda_A)),
            infer_sigma_a = MPD(1 / sqrt(Y$lambda_a)),
            infer_sigma_theta = MPD(1 / sqrt(Y$lambda_theta)),

            # 95% CI Lower (2.5%)
            lower_mu_w = get_q("mu_w", "2.5%"),
            lower_mu_A = get_q("mu_A", "2.5%"),
            lower_mu_a = get_q("mu_a", "2.5%"),
            lower_mu_theta = get_q("mu_theta", "2.5%"),

            # 95% CI Upper (97.5%)
            upper_mu_w = get_q("mu_w", "97.5%"),
            upper_mu_A = get_q("mu_A", "97.5%"),
            upper_mu_a = get_q("mu_a", "97.5%"),
            upper_mu_theta = get_q("mu_theta", "97.5%")
        )

        results_list[[length(results_list) + 1]] <- res
        cat("  Success.\n")
    }
}

# save results
# Ensure directory exists
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

if (length(results_list) > 0) {
    saveRDS(results_list, file = output_file)
    cat("Saved", length(results_list), "results to", output_file, "\n")
} else {
    cat("No successful iterations to save.\n")
}
