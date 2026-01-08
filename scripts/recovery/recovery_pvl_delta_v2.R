# ==============================================================================
# Dependencies
if (!require("pacman")) install.packages("pacman")
pacman::p_load(R2jags, parallel, ggpubr, extraDistr, truncnorm)
# ==============================================================================
# Parameter Recovery: PVL-Delta Model
# ==============================================================================
#
# Purpose: Validate that our JAGS implementation can accurately recover
# known parameter values from simulated data.
#
# Method:
#   1. Generate "true" group-level parameters (random draws)
#   2. Simulate choices from the PVL-Delta model with these parameters
#   3. Fit the JAGS model to recover parameters
#   4. Compare true vs recovered (correlation should be high, r > 0.7)
#
# Success Criteria:
#   - Strong correlations between true and recovered parameters
#   - No systematic bias (regression slope near 1)
# ==============================================================================

set.seed(69420)

# ==============================================================================
# Helper Functions
# ==============================================================================

# Maximum Posterior Density
MPD <- function(x) {
    density(x)$x[which(density(x)$y == max(density(x)$y))]
}

# Convert JAGS precision (lambda) to standard deviation (sigma)
precision_to_sd <- function(lambda) {
    1 / sqrt(lambda)
}

# Source dependencies
source("scripts/recovery/simulation_pvl_delta_v2.R")
source("utils/payoff_scheme.R")
source("scripts/plotting/plotting_utils.R")

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
# Allow command line args for testing (e.g., --test)
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

cat("Starting PVL-Delta Parameter Recovery...\n")
cat("Configuration: ", niterations, " iterations x ", nsubs, " subjects\n")
cat("Running on", n_cores, "cores\n\n")

# Wrapper function for a single recovery iteration
run_iteration <- function(i) {
    # -------------------------------------------------------------------------
    # Step 1: Generate True Parameters
    # -------------------------------------------------------------------------
    mu_w <- runif(1, 0.5, 2.5) # Loss aversion
    sigma_w <- runif(1, 0.1, 0.3)

    mu_A <- runif(1, 0.2, 0.8) # Outcome sensitivity
    sigma_A <- runif(1, 0.05, 0.15)

    mu_theta <- runif(1, 0.5, 2.0) # Inverse temperature
    sigma_theta <- runif(1, 0.1, 0.3)

    mu_a <- runif(1, 0.1, 0.5) # Learning rate
    sigma_a <- runif(1, 0.05, 0.15)

    # -------------------------------------------------------------------------
    # Step 2: Simulate Data
    # -------------------------------------------------------------------------
    sim_data <- hier_PVL_sim(
        payoff_struct = payoff_struct,
        nsubs = nsubs,
        ntrials = ntrials_all,
        mu_w = mu_w, mu_A = mu_A, mu_a = mu_a, mu_theta = mu_theta,
        sigma_w = sigma_w, sigma_A = sigma_A, sigma_a = sigma_a, sigma_theta = sigma_theta
    )

    # -------------------------------------------------------------------------
    # Step 3: Fit JAGS Model
    # -------------------------------------------------------------------------
    jags_data <- list("x" = sim_data$x, "X" = sim_data$X, "ntrials" = ntrials_all, "nsubs" = nsubs)
    params <- c(
        "mu_w", "mu_A", "mu_theta", "mu_a",
        "lambda_w", "lambda_A", "lambda_theta", "lambda_a"
    )

    # Use standard jags() inside parallel worker to avoid nested parallel conflicts
    samples <- jags(
        data = jags_data,
        inits = NULL,
        parameters.to.save = params,
        model.file = "models/pvl_delta_v2.txt",
        n.chains = 3,
        n.iter = 3000,
        n.burnin = 1000,
        n.thin = 1,
        progress.bar = "none" # Suppress output
    )

    # -------------------------------------------------------------------------
    # Step 4: Extract Point Estimates
    # -------------------------------------------------------------------------
    Y <- samples$BUGSoutput$sims.list

    # Return list of results
    list(
        # True
        true_mu_w = mu_w, true_mu_A = mu_A, true_mu_theta = mu_theta, true_mu_a = mu_a,
        true_sigma_w = sigma_w, true_sigma_A = sigma_A, true_sigma_theta = sigma_theta, true_sigma_a = sigma_a,

        # Inferred (Means)
        infer_mu_w = MPD(Y$mu_w), infer_mu_A = MPD(Y$mu_A),
        infer_mu_theta = MPD(Y$mu_theta), infer_mu_a = MPD(Y$mu_a),

        # Inferred (SDs)
        infer_sigma_w = precision_to_sd(MPD(Y$lambda_w)),
        infer_sigma_A = precision_to_sd(MPD(Y$lambda_A)),
        infer_sigma_theta = precision_to_sd(MPD(Y$lambda_theta)),
        infer_sigma_a = precision_to_sd(MPD(Y$lambda_a))
    )
}

# Run Parallel Loop
results_list <- mclapply(1:niterations, run_iteration, mc.cores = n_cores, mc.set.seed = TRUE)

cat("Simulation complete. Processing results...\n")

# Unpack results
true_mu_w <- sapply(results_list, function(x) x$true_mu_w)
true_mu_A <- sapply(results_list, function(x) x$true_mu_A)
true_mu_theta <- sapply(results_list, function(x) x$true_mu_theta)
true_mu_a <- sapply(results_list, function(x) x$true_mu_a)

infer_mu_w <- sapply(results_list, function(x) x$infer_mu_w)
infer_mu_A <- sapply(results_list, function(x) x$infer_mu_A)
infer_mu_theta <- sapply(results_list, function(x) x$infer_mu_theta)
infer_mu_a <- sapply(results_list, function(x) x$infer_mu_a)

end_time <- Sys.time()
cat("Total time:", round(difftime(end_time, start_time, units = "mins"), 1), "minutes\n")

# ==============================================================================
# Save Results (RDS & CSV)
# ==============================================================================
output_dir_data <- "outputs/recovery"
dir.create(output_dir_data, recursive = TRUE, showWarnings = FALSE)

# Save full list (includes all parameters and true values)
saveRDS(results_list, file.path(output_dir_data, "recovery_pvl_delta.rds"))

# Create summary CSV for means
df_summary <- data.frame(
    true_mu_w = true_mu_w, infer_mu_w = infer_mu_w,
    true_mu_A = true_mu_A, infer_mu_A = infer_mu_A,
    true_mu_theta = true_mu_theta, infer_mu_theta = infer_mu_theta,
    true_mu_a = true_mu_a, infer_mu_a = infer_mu_a
)
write.csv(df_summary, file.path(output_dir_data, "recovery_pvl_delta.csv"), row.names = FALSE)

cat("Results saved to:", output_dir_data, "\n")

# ==============================================================================
# Plotting Results
# ==============================================================================
pl1 <- plot_recovery(true_mu_w, infer_mu_w, "mu_w (Loss Aversion)")
pl2 <- plot_recovery(true_mu_A, infer_mu_A, "mu_A (Outcome Sensitivity)")
pl3 <- plot_recovery(true_mu_theta, infer_mu_theta, "mu_theta (Inv. Temperature)")
pl4 <- plot_recovery(true_mu_a, infer_mu_a, "mu_a (Learning Rate)")

final_plot <- ggarrange(pl1, pl2, pl3, pl4, nrow = 2, ncol = 2)

output_dir <- "analysis/plots/recovery"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output_dir, "recovery_pvl_delta.png"), final_plot,
    width = 12, height = 10, dpi = 150
)

cat("Recovery plot saved to:", file.path(output_dir, "recovery_pvl_delta.png"), "\n")

# Print correlation summary
cat("\n=== Recovery Correlations ===\n")
cat("mu_w:     r =", round(cor(true_mu_w, infer_mu_w), 3), "\n")
cat("mu_A:     r =", round(cor(true_mu_A, infer_mu_A), 3), "\n")
cat("mu_theta: r =", round(cor(true_mu_theta, infer_mu_theta), 3), "\n")
cat("mu_a:     r =", round(cor(true_mu_a, infer_mu_a), 3), "\n")
