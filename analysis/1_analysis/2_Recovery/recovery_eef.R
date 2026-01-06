# ==============================================================================
# Dependencies
if (!require("pacman")) install.packages("pacman")
pacman::p_load(R2jags, parallel, ggpubr, extraDistr, truncnorm)
# ==============================================================================
# Parameter Recovery: EEF (Explore-Exploit with Forgetting) Model
# ==============================================================================
#
# Purpose: Validate that our JAGS implementation can accurately recover
# known parameter values from simulated data.
#
# The EEF model has 4 parameters:
#   - theta: Outcome sensitivity (0-1)
#   - lambda: Learning/forgetting rate (0-1)
#   - phi: Exploration bonus (-5 to 5)
#   - cons: Consistency, transformed to C = 3^cons - 1 (0-5)
#
# Success Criteria:
#   - Strong correlations between true and recovered parameters (r > 0.7)
# ==============================================================================

set.seed(69420)

# ==============================================================================
# Helper Functions
# ==============================================================================
MPD <- function(x) {
    density(x)$x[which(density(x)$y == max(density(x)$y))]
}

precision_to_sd <- function(lambda) {
    1 / sqrt(lambda)
}

source("analysis/1_analysis/2_Recovery/simulation_eef.R")
source("analysis/utils/payoff_scheme.R")
source("analysis/utils/plotting_utils.R")

# ==============================================================================
# Task Setup
# ==============================================================================
ntrials <- 100
payoff_struct <- generate_modified_igt_payoff(ntrials, scale = TRUE)

cat("Modified IGT Payoff (Scaled /100) loaded.\n")
cat("Payoff structure contains separate gain and loss matrices.\n\n")

# ==============================================================================
# Recovery Configuration
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)
if ("--test" %in% args) {
    niterations <- 2
    cat(">>> RUNNING IN TEST MODE (2 iterations) <<<\n")
} else {
    niterations <- 100 # Standard recovery count
}

nsubs <- 24 # Reduced for computational feasibility; still valid for recovery
ntrials_all <- rep(100, nsubs)

# ==============================================================================
# Main Recovery Loop (Parallelized)
# ==============================================================================
start_time <- Sys.time()
n_cores <- min(40, parallel::detectCores() - 1)
if (n_cores < 1) n_cores <- 1

cat("Starting EEF Parameter Recovery...\n")
cat("Configuration:", niterations, "iterations x", nsubs, "subjects\n")
cat("Running on", n_cores, "cores\n\n")

run_iteration <- function(i) {
    # -------------------------------------------------------------------------
    # Step 1: Generate True Parameters
    # -------------------------------------------------------------------------
    mu_theta <- runif(1, 0.2, 0.8)
    sigma_theta <- runif(1, 0.05, 0.15)

    mu_lambda <- runif(1, 0.2, 0.8)
    sigma_lambda <- runif(1, 0.05, 0.15)

    mu_phi <- runif(1, -2, 2)
    sigma_phi <- runif(1, 0.1, 0.5)

    mu_cons <- runif(1, 1, 3)
    sigma_cons <- runif(1, 0.1, 0.5)

    # -------------------------------------------------------------------------
    # Step 2: Simulate Data
    # -------------------------------------------------------------------------
    sim_data <- simulation_eef(
        payoff_struct = payoff_struct,
        nsubs = nsubs,
        ntrials = ntrials_all,
        mu_theta = mu_theta, mu_lambda = mu_lambda,
        mu_phi = mu_phi, mu_cons = mu_cons,
        sigma_theta = sigma_theta, sigma_lambda = sigma_lambda,
        sigma_phi = sigma_phi, sigma_cons = sigma_cons
    )

    # -------------------------------------------------------------------------
    # Step 3: Fit JAGS Model
    # -------------------------------------------------------------------------
    jags_data <- list("x" = sim_data$x, "X" = sim_data$X, "ntrials" = ntrials_all, "nsubs" = nsubs)

    params <- c(
        "mu_theta", "mu_lambda", "mu_phi", "mu_cons",
        "lambda_theta", "lambda_lambda", "lambda_phi", "lambda_cons"
    )

    samples <- jags(
        data = jags_data,
        inits = NULL,
        parameters.to.save = params,
        model.file = "analysis/models/eef.txt",
        n.chains = 3,
        n.iter = 3000,
        n.burnin = 1000,
        n.thin = 1,
        progress.bar = "none"
    )

    # -------------------------------------------------------------------------
    # Step 4: Extract Point Estimates
    # -------------------------------------------------------------------------
    Y <- samples$BUGSoutput$sims.list

    list(
        # True
        true_mu_theta = mu_theta, true_mu_lambda = mu_lambda,
        true_mu_phi = mu_phi, true_mu_cons = mu_cons,

        # Inferred ( Means)
        infer_mu_theta = MPD(Y$mu_theta), infer_mu_lambda = MPD(Y$mu_lambda),
        infer_mu_phi = MPD(Y$mu_phi), infer_mu_cons = MPD(Y$mu_cons)
    )
}

results_list <- mclapply(1:niterations, run_iteration, mc.cores = n_cores, mc.set.seed = TRUE)

cat("Simulation complete. Processing results...\n")

# Unpack results
true_mu_theta <- sapply(results_list, function(x) x$true_mu_theta)
infer_mu_theta <- sapply(results_list, function(x) x$infer_mu_theta)

true_mu_lambda <- sapply(results_list, function(x) x$true_mu_lambda)
infer_mu_lambda <- sapply(results_list, function(x) x$infer_mu_lambda)

true_mu_phi <- sapply(results_list, function(x) x$true_mu_phi)
infer_mu_phi <- sapply(results_list, function(x) x$infer_mu_phi)

true_mu_cons <- sapply(results_list, function(x) x$true_mu_cons)
infer_mu_cons <- sapply(results_list, function(x) x$infer_mu_cons)


end_time <- Sys.time()
cat("Total time:", round(difftime(end_time, start_time, units = "mins"), 1), "minutes\n")

# ==============================================================================
# Save Results (RDS & CSV)
# ==============================================================================
output_dir_data <- "analysis/outputs/recovery"
dir.create(output_dir_data, recursive = TRUE, showWarnings = FALSE)

# Save full list
saveRDS(results_list, file.path(output_dir_data, "recovery_eef.rds"))

# Create summary CSV for means
df_summary <- data.frame(
    true_mu_theta = true_mu_theta, infer_mu_theta = infer_mu_theta,
    true_mu_lambda = true_mu_lambda, infer_mu_lambda = infer_mu_lambda,
    true_mu_phi = true_mu_phi, infer_mu_phi = infer_mu_phi,
    true_mu_cons = true_mu_cons, infer_mu_cons = infer_mu_cons
)
write.csv(df_summary, file.path(output_dir_data, "recovery_eef.csv"), row.names = FALSE)

cat("Results saved to:", output_dir_data, "\n")

# ==============================================================================
# Plotting Results
# ==============================================================================
pl1 <- plot_recovery(true_mu_theta, infer_mu_theta, "mu_theta (Outcome Sensitivity)")
pl2 <- plot_recovery(true_mu_lambda, infer_mu_lambda, "mu_lambda (Learning Rate)")
pl3 <- plot_recovery(true_mu_phi, infer_mu_phi, "mu_phi (Exploration Bonus)")
pl4 <- plot_recovery(true_mu_cons, infer_mu_cons, "mu_cons (Consistency)")

final_plot <- ggarrange(pl1, pl2, pl3, pl4, nrow = 2, ncol = 2)

output_dir <- "analysis/plots/recovery"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output_dir, "recovery_eef.png"), final_plot,
    width = 12, height = 10, dpi = 150
)

cat("Recovery plot saved to:", file.path(output_dir, "recovery_eef.png"), "\n")

# Print correlation summary
cat("\n=== Recovery Correlations ===\n")
cat("mu_theta:  r =", round(cor(true_mu_theta, infer_mu_theta), 3), "\n")
cat("mu_lambda: r =", round(cor(true_mu_lambda, infer_mu_lambda), 3), "\n")
cat("mu_phi:    r =", round(cor(true_mu_phi, infer_mu_phi), 3), "\n")
cat("mu_cons:   r =", round(cor(true_mu_cons, infer_mu_cons), 3), "\n")
