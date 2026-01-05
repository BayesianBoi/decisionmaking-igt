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

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(R2jags, parallel, ggpubr, extraDistr, truncnorm)

set.seed(69420)

# ==============================================================================
# Helper Functions
# ==============================================================================

# Maximum Posterior Density (point estimate from posterior distribution)
# More robust than mean for skewed posteriors
MPD <- function(x) {
    density(x)$x[which(density(x)$y == max(density(x)$y))]
}

# Convert JAGS precision (lambda) to standard deviation (sigma)
# JAGS uses precision (1/variance), we want SD for comparison
precision_to_sd <- function(lambda) {
    1 / sqrt(lambda)
}

# Source dependencies
source("analysis/1_analysis/2_Recovery/simulation_pvl_delta.R")
source("analysis/utils/payoff_scheme.R")
source("analysis/utils/plotting_utils.R")

# ==============================================================================
# Task Setup - Using gain/loss structure for deck-based indexing
# ==============================================================================
ntrials <- 100
payoff_struct <- generate_modified_igt_payoff(ntrials, scale = TRUE)

cat("Modified IGT Payoff (Scaled /100) loaded.\n")
cat("Payoff structure contains separate gain and loss matrices.\n\n")

# ==============================================================================
# Recovery Configuration
# ==============================================================================
niterations <- 100 # Number of recovery iterations
nsubs <- 48 # Number of simulated subjects (matches Ahn 2014 HC group size)
ntrials_all <- rep(100, nsubs)

# ==============================================================================
# Storage Arrays for True and Recovered Parameters
# ==============================================================================
# Group-level means (mu)
true_mu_w <- array(NA, c(niterations))
true_mu_A <- array(NA, c(niterations))
true_mu_theta <- array(NA, c(niterations))
true_mu_a <- array(NA, c(niterations))

infer_mu_w <- array(NA, c(niterations))
infer_mu_A <- array(NA, c(niterations))
infer_mu_theta <- array(NA, c(niterations))
infer_mu_a <- array(NA, c(niterations))

# Group-level standard deviations (sigma)
# Note: JAGS estimates precision (lambda), we convert to sigma for comparison
true_sigma_w <- array(NA, c(niterations))
true_sigma_A <- array(NA, c(niterations))
true_sigma_theta <- array(NA, c(niterations))
true_sigma_a <- array(NA, c(niterations))

infer_sigma_w <- array(NA, c(niterations))
infer_sigma_A <- array(NA, c(niterations))
infer_sigma_theta <- array(NA, c(niterations))
infer_sigma_a <- array(NA, c(niterations))

# ==============================================================================
# Main Recovery Loop
# ==============================================================================
start_time <- Sys.time()
cat("Starting PVL-Delta Parameter Recovery...\n")
cat("Configuration: ", niterations, " iterations x ", nsubs, " subjects\n\n")

for (i in 1:niterations) {
    # -------------------------------------------------------------------------
    # Step 1: Generate True Parameters
    # -------------------------------------------------------------------------
    # Draw group-level means and SDs from plausible ranges
    # These ranges are informed by prior literature

    # Loss aversion (w): typically 1.5-2.5 in healthy adults
    mu_w <- runif(1, 0.5, 2.5)
    sigma_w <- runif(1, 0.1, 0.3)

    # Outcome sensitivity (A): typically 0.3-0.8
    mu_A <- runif(1, 0.2, 0.8)
    sigma_A <- runif(1, 0.05, 0.15)

    # Inverse temperature (theta): controls choice consistency
    mu_theta <- runif(1, 0.5, 2.0)
    sigma_theta <- runif(1, 0.1, 0.3)

    # Learning rate (a): bounded 0-1
    mu_a <- runif(1, 0.1, 0.5)
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

    x <- sim_data$x
    X <- sim_data$X

    # -------------------------------------------------------------------------
    # Step 3: Fit JAGS Model
    # -------------------------------------------------------------------------
    jags_data <- list("x" = x, "X" = X, "ntrials" = ntrials_all, "nsubs" = nsubs)
    params <- c(
        "mu_w", "mu_A", "mu_theta", "mu_a",
        "lambda_w", "lambda_A", "lambda_theta", "lambda_a"
    )

    samples <- jags.parallel(
        data = jags_data,
        inits = NULL,
        parameters.to.save = params,
        model.file = "analysis/models/pvl_delta.txt",
        n.chains = 3,
        n.iter = 3000,
        n.burnin = 1000,
        n.thin = 1
    )

    # -------------------------------------------------------------------------
    # Step 4: Extract Point Estimates
    # -------------------------------------------------------------------------
    Y <- samples$BUGSoutput$sims.list

    # Store true values
    true_mu_w[i] <- mu_w
    true_mu_A[i] <- mu_A
    true_mu_theta[i] <- mu_theta
    true_mu_a[i] <- mu_a

    true_sigma_w[i] <- sigma_w
    true_sigma_A[i] <- sigma_A
    true_sigma_theta[i] <- sigma_theta
    true_sigma_a[i] <- sigma_a

    # Extract MPD from posteriors (means)
    infer_mu_w[i] <- MPD(Y$mu_w)
    infer_mu_A[i] <- MPD(Y$mu_A)
    infer_mu_theta[i] <- MPD(Y$mu_theta)
    infer_mu_a[i] <- MPD(Y$mu_a)

    # Extract MPD and convert precision to SD
    infer_sigma_w[i] <- precision_to_sd(MPD(Y$lambda_w))
    infer_sigma_A[i] <- precision_to_sd(MPD(Y$lambda_A))
    infer_sigma_theta[i] <- precision_to_sd(MPD(Y$lambda_theta))
    infer_sigma_a[i] <- precision_to_sd(MPD(Y$lambda_a))

    cat("Iteration", i, "of", niterations, "complete\n")
}

end_time <- Sys.time()
cat("\nTotal time:", round(difftime(end_time, start_time, units = "mins"), 1), "minutes\n")

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

cat("\nRecovery plot saved to:", file.path(output_dir, "recovery_pvl_delta.png"), "\n")

# Print correlation summary
cat("\n=== Recovery Correlations ===\n")
cat("mu_w:     r =", round(cor(true_mu_w, infer_mu_w), 3), "\n")
cat("mu_A:     r =", round(cor(true_mu_A, infer_mu_A), 3), "\n")
cat("mu_theta: r =", round(cor(true_mu_theta, infer_mu_theta), 3), "\n")
cat("mu_a:     r =", round(cor(true_mu_a, infer_mu_a), 3), "\n")
