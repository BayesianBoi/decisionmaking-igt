# ==============================================================================
# Dependencies
if (!require("pacman")) install.packages("pacman")
pacman::p_load(R2jags, parallel, ggpubr, extraDistr, truncnorm)
# ==============================================================================
# ==============================================================================
# Parameter Recovery: ORL Model
# ==============================================================================
# Purpose: Validate parameter estimation accuracy using simulated data.
# Validation Criteria:
#   - Correlation between true and recovered parameters (r > 0.7)
#   - Absence of systematic bias
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

source("scripts/recovery/simulation_orl.R")
source("utils/payoff_scheme.R")
source("analysis/2_plotting/plotting_utils.R")

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

nsubs <- 48 # Number of simulated subjects (matches Ahn 2014 HC group size)
ntrials_all <- rep(100, nsubs)

# ==============================================================================
# Main Recovery Loop (Parallelized)
# ==============================================================================
start_time <- Sys.time()
n_cores <- parallel::detectCores() - 1
if (n_cores < 1) n_cores <- 1

cat("Starting ORL Parameter Recovery...\n")
cat("Configuration:", niterations, "iterations x", nsubs, "subjects\n")
cat("Running on", n_cores, "cores\n\n")

run_iteration <- function(i) {
    # -------------------------------------------------------------------------
    # Step 1: Generate True Parameters
    # -------------------------------------------------------------------------
    mu_a_rew <- runif(1, 0.1, 0.5)
    mu_a_pun <- runif(1, 0.1, 0.5)
    sigma_a_rew <- runif(1, 0.05, 0.15)
    sigma_a_pun <- runif(1, 0.05, 0.15)

    mu_K <- runif(1, 0.5, 2.0)
    sigma_K <- runif(1, 0.1, 0.3)

    mu_theta <- runif(1, 0.5, 2.0)
    sigma_theta <- runif(1, 0.1, 0.3)

    mu_omega_f <- runif(1, -1, 1)
    mu_omega_p <- runif(1, -1, 1)
    sigma_omega_f <- runif(1, 0.1, 0.5)
    sigma_omega_p <- runif(1, 0.1, 0.5)

    # -------------------------------------------------------------------------
    # Step 2: Simulate Data
    # -------------------------------------------------------------------------
    sim_data <- hier_ORL_sim(
        payoff_struct = payoff_struct,
        nsubs = nsubs,
        ntrials = ntrials_all,
        mu_a_rew = mu_a_rew, mu_a_pun = mu_a_pun,
        mu_K = mu_K, mu_theta = mu_theta,
        mu_omega_f = mu_omega_f, mu_omega_p = mu_omega_p,
        sigma_a_rew = sigma_a_rew, sigma_a_pun = sigma_a_pun,
        sigma_K = sigma_K, sigma_theta = sigma_theta,
        sigma_omega_f = sigma_omega_f, sigma_omega_p = sigma_omega_p
    )

    # -------------------------------------------------------------------------
    # Step 3: Fit JAGS Model
    # -------------------------------------------------------------------------
    jags_data <- list("x" = sim_data$x, "X" = sim_data$X, "ntrials" = ntrials_all, "nsubs" = nsubs)

    params <- c(
        "mu_a_rew", "mu_a_pun", "mu_K", "mu_theta", "mu_omega_f", "mu_omega_p",
        "lambda_a_rew", "lambda_a_pun", "lambda_K", "lambda_theta",
        "lambda_omega_f", "lambda_omega_p"
    )

    # Sequential chains inside parallel worker
    samples <- jags(
        data = jags_data,
        inits = NULL,
        parameters.to.save = params,
        model.file = "models/orl.txt",
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
        true_mu_a_rew = mu_a_rew, true_mu_a_pun = mu_a_pun,
        true_mu_K = mu_K, true_mu_theta = mu_theta,
        true_mu_omega_f = mu_omega_f, true_mu_omega_p = mu_omega_p,

        # Inferred ( Means)
        infer_mu_a_rew = MPD(Y$mu_a_rew), infer_mu_a_pun = MPD(Y$mu_a_pun),
        infer_mu_K = MPD(Y$mu_K), infer_mu_theta = MPD(Y$mu_theta),
        infer_mu_omega_f = MPD(Y$mu_omega_f), infer_mu_omega_p = MPD(Y$mu_omega_p)
    )
}

results_list <- mclapply(1:niterations, run_iteration, mc.cores = n_cores, mc.set.seed = TRUE)

cat("Simulation complete. Processing results...\n")

# Unpack results
true_mu_a_rew <- sapply(results_list, function(x) x$true_mu_a_rew)
infer_mu_a_rew <- sapply(results_list, function(x) x$infer_mu_a_rew)

true_mu_a_pun <- sapply(results_list, function(x) x$true_mu_a_pun)
infer_mu_a_pun <- sapply(results_list, function(x) x$infer_mu_a_pun)

true_mu_K <- sapply(results_list, function(x) x$true_mu_K)
infer_mu_K <- sapply(results_list, function(x) x$infer_mu_K)

true_mu_theta <- sapply(results_list, function(x) x$true_mu_theta)
infer_mu_theta <- sapply(results_list, function(x) x$infer_mu_theta)

true_mu_omega_f <- sapply(results_list, function(x) x$true_mu_omega_f)
infer_mu_omega_f <- sapply(results_list, function(x) x$infer_mu_omega_f)

true_mu_omega_p <- sapply(results_list, function(x) x$true_mu_omega_p)
infer_mu_omega_p <- sapply(results_list, function(x) x$infer_mu_omega_p)


end_time <- Sys.time()
cat("Total time:", round(difftime(end_time, start_time, units = "mins"), 1), "minutes\n")

# ==============================================================================
# Save Results (RDS & CSV)
# ==============================================================================
output_dir_data <- "outputs/recovery"
dir.create(output_dir_data, recursive = TRUE, showWarnings = FALSE)

# Save full list (includes all parameters and true values)
saveRDS(results_list, file.path(output_dir_data, "recovery_orl.rds"))

# Create summary CSV for means
df_summary <- data.frame(
    true_mu_a_rew = true_mu_a_rew, infer_mu_a_rew = infer_mu_a_rew,
    true_mu_a_pun = true_mu_a_pun, infer_mu_a_pun = infer_mu_a_pun,
    true_mu_K = true_mu_K, infer_mu_K = infer_mu_K,
    true_mu_theta = true_mu_theta, infer_mu_theta = infer_mu_theta,
    true_mu_omega_f = true_mu_omega_f, infer_mu_omega_f = infer_mu_omega_f,
    true_mu_omega_p = true_mu_omega_p, infer_mu_omega_p = infer_mu_omega_p
)
write.csv(df_summary, file.path(output_dir_data, "recovery_orl.csv"), row.names = FALSE)

cat("Results saved to:", output_dir_data, "\n")

# ==============================================================================
# Plotting Results
# ==============================================================================
pl1 <- plot_recovery(true_mu_a_rew, infer_mu_a_rew, "mu_a_rew (Reward Learning)")
pl2 <- plot_recovery(true_mu_a_pun, infer_mu_a_pun, "mu_a_pun (Punishment Learning)")
pl3 <- plot_recovery(true_mu_K, infer_mu_K, "mu_K (Perseverance Decay)")
pl4 <- plot_recovery(true_mu_theta, infer_mu_theta, "mu_theta (Inv. Temperature)")
pl5 <- plot_recovery(true_mu_omega_f, infer_mu_omega_f, "mu_omega_f (Frequency Weight)")
pl6 <- plot_recovery(true_mu_omega_p, infer_mu_omega_p, "mu_omega_p (Perseverance Weight)")

final_plot <- ggarrange(pl1, pl2, pl3, pl4, pl5, pl6, nrow = 2, ncol = 3)

output_dir <- "analysis/plots/recovery"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output_dir, "recovery_orl.png"), final_plot,
    width = 15, height = 10, dpi = 150
)

cat("Recovery plot saved to:", file.path(output_dir, "recovery_orl.png"), "\n")

# Print correlation summary
cat("\n=== Recovery Correlations ===\n")
cat("mu_a_rew:   r =", round(cor(true_mu_a_rew, infer_mu_a_rew), 3), "\n")
cat("mu_a_pun:   r =", round(cor(true_mu_a_pun, infer_mu_a_pun), 3), "\n")
cat("mu_K:       r =", round(cor(true_mu_K, infer_mu_K), 3), "\n")
cat("mu_theta:   r =", round(cor(true_mu_theta, infer_mu_theta), 3), "\n")
cat("mu_omega_f: r =", round(cor(true_mu_omega_f, infer_mu_omega_f), 3), "\n")
cat("mu_omega_p: r =", round(cor(true_mu_omega_p, infer_mu_omega_p), 3), "\n")
