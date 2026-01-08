# ==============================================================================
# ORL Parameter Recovery - Batch Worker Script
# ==============================================================================
# This script runs a subset of recovery iterations for the ORL model.
# It is designed to be run in parallel via tmux or another scheduler.
#
# Usage:
#   Rscript recovery_orl_batch.R --seed 123 --iter 10 --output path/to/save.rds
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Parse Command Line Arguments
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# Default values
seed <- NULL
n_iter <- 5
output_file <- "recovery_orl_batch.rds"

# Simple argument parser
for (i in seq_along(args)) {
    if (args[i] == "--seed" && i < length(args)) seed <- as.integer(args[i + 1])
    if (args[i] == "--iter" && i < length(args)) n_iter <- as.integer(args[i + 1])
    if (args[i] == "--output" && i < length(args)) output_file <- args[i + 1]
}

if (is.null(seed)) stop("Error: --seed argument is required.")

cat("----------------------------------------------------------------\n")
cat("ORL Recovery Batch Worker\n")
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
source("scripts/recovery/simulation_orl.R")
source("utils/payoff_scheme.R")

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
    mu_a_rew <- runif(1, 0.1, 0.5)
    mu_a_pun <- runif(1, 0.1, 0.5)
    sigma_a_rew <- runif(1, 0.05, 0.15)
    sigma_a_pun <- runif(1, 0.05, 0.15)

    mu_K <- runif(1, 0.5, 2.0)
    sigma_K <- runif(1, 0.1, 0.3)

    # mu_theta removed (fixed to 1)

    mu_omega_f <- runif(1, -1, 1)
    mu_omega_p <- runif(1, -1, 1)
    sigma_omega_f <- runif(1, 0.1, 0.5)
    sigma_omega_p <- runif(1, 0.1, 0.5)

    # B. Simulate Data
    sim_data <- tryCatch(
        {
            hier_ORL_sim(
                payoff_struct = payoff_struct,
                nsubs = nsubs,
                ntrials = ntrials_all,
                mu_a_rew = mu_a_rew, mu_a_pun = mu_a_pun,
                mu_K = mu_K,
                mu_omega_f = mu_omega_f, mu_omega_p = mu_omega_p,
                sigma_a_rew = sigma_a_rew, sigma_a_pun = sigma_a_pun,
                sigma_K = sigma_K,
                sigma_omega_f = sigma_omega_f, sigma_omega_p = sigma_omega_p
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
        "mu_a_rew", "mu_a_pun", "mu_K", "mu_omega_f", "mu_omega_p",
        "lambda_a_rew", "lambda_a_pun", "lambda_K",
        "lambda_omega_f", "lambda_omega_p"
    )

    # Use tryCatch for JAGS fitting
    fit_result <- tryCatch(
        {
            jags(
                data = jags_data,
                inits = NULL,
                parameters.to.save = params,
                model.file = "models/orl.txt",
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
        # Helper to extract quantiles
        get_q <- function(param_name, q_col) {
            if (param_name %in% rownames(fit_result$BUGSoutput$summary)) {
                return(fit_result$BUGSoutput$summary[param_name, q_col])
            } else {
                return(NA)
            }
        }

        # Extract posterior samples
        Y <- fit_result$BUGSoutput$sims.list

        res <- list(
            # True Means
            true_mu_a_rew = mu_a_rew, true_mu_a_pun = mu_a_pun,
            true_mu_K = mu_K,
            true_mu_omega_f = mu_omega_f, true_mu_omega_p = mu_omega_p,

            # True Sigmas
            true_sigma_a_rew = sigma_a_rew, true_sigma_a_pun = sigma_a_pun,
            true_sigma_K = sigma_K,
            true_sigma_omega_f = sigma_omega_f, true_sigma_omega_p = sigma_omega_p,

            # Inferred Means
            infer_mu_a_rew = MPD(Y$mu_a_rew), infer_mu_a_pun = MPD(Y$mu_a_pun),
            infer_mu_K = MPD(Y$mu_K),
            infer_mu_omega_f = MPD(Y$mu_omega_f), infer_mu_omega_p = MPD(Y$mu_omega_p),

            # Inferred Sigmas
            infer_sigma_a_rew = MPD(1 / sqrt(Y$lambda_a_rew)),
            infer_sigma_a_pun = MPD(1 / sqrt(Y$lambda_a_pun)),
            infer_sigma_K = MPD(1 / sqrt(Y$lambda_K)),
            infer_sigma_omega_f = MPD(1 / sqrt(Y$lambda_omega_f)),
            infer_sigma_omega_p = MPD(1 / sqrt(Y$lambda_omega_p)),

            # 95% CI Lower (2.5%)
            lower_mu_a_rew = get_q("mu_a_rew", "2.5%"),
            lower_mu_a_pun = get_q("mu_a_pun", "2.5%"),
            lower_mu_K = get_q("mu_K", "2.5%"),
            lower_mu_omega_f = get_q("mu_omega_f", "2.5%"),
            lower_mu_omega_p = get_q("mu_omega_p", "2.5%"),

            # 95% CI Upper (97.5%)
            upper_mu_a_rew = get_q("mu_a_rew", "97.5%"),
            upper_mu_a_pun = get_q("mu_a_pun", "97.5%"),
            upper_mu_K = get_q("mu_K", "97.5%"),
            upper_mu_omega_f = get_q("mu_omega_f", "97.5%"),
            upper_mu_omega_p = get_q("mu_omega_p", "97.5%")
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
