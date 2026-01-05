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
# Reference: Yang et al. (2025)
#
# Success Criteria:
#   - Strong correlations between true and recovered parameters (r > 0.7)
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(R2jags, parallel, ggpubr, extraDistr, truncnorm)

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
# Task Setup - Using gain/loss structure for deck-based indexing
# ==============================================================================
ntrials <- 100
payoff_struct <- generate_modified_igt_payoff(ntrials, scale = TRUE)

cat("Modified IGT Payoff (Scaled /100) loaded.\n")
cat("Payoff structure contains separate gain and loss matrices.\n\n")

# ==============================================================================
# Recovery Configuration
# ==============================================================================
niterations <- 100
nsubs <- 48 # Number of simulated subjects (matches Ahn 2014 HC group size)
ntrials_all <- rep(100, nsubs)

# ==============================================================================
# Storage Arrays
# ==============================================================================
true_mu_theta <- infer_mu_theta <- array(NA, c(niterations))
true_mu_lambda <- infer_mu_lambda <- array(NA, c(niterations))
true_mu_phi <- infer_mu_phi <- array(NA, c(niterations))
true_mu_cons <- infer_mu_cons <- array(NA, c(niterations))

true_sigma_theta <- infer_sigma_theta <- array(NA, c(niterations))
true_sigma_lambda <- infer_sigma_lambda <- array(NA, c(niterations))
true_sigma_phi <- infer_sigma_phi <- array(NA, c(niterations))
true_sigma_cons <- infer_sigma_cons <- array(NA, c(niterations))

# ==============================================================================
# Main Recovery Loop
# ==============================================================================
start_time <- Sys.time()
cat("Starting EEF Parameter Recovery...\n")
cat("Configuration:", niterations, "iterations x", nsubs, "subjects\n\n")

for (i in 1:niterations) {
    # -------------------------------------------------------------------------
    # Step 1: Generate True Parameters
    # -------------------------------------------------------------------------
    # Outcome sensitivity theta (0-1)
    mu_theta <- runif(1, 0.2, 0.8)
    sigma_theta <- runif(1, 0.05, 0.15)

    # Learning rate lambda (0-1)
    mu_lambda <- runif(1, 0.2, 0.8)
    sigma_lambda <- runif(1, 0.05, 0.15)

    # Exploration bonus phi (-5 to 5)
    # Can be negative (exploration aversion) or positive (curiosity)
    mu_phi <- runif(1, -2, 2)
    sigma_phi <- runif(1, 0.1, 0.5)

    # Consistency cons (0-5), transformed to C = 3^cons - 1
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

    x <- sim_data$x
    X <- sim_data$X

    # -------------------------------------------------------------------------
    # Step 3: Fit JAGS Model
    # -------------------------------------------------------------------------
    jags_data <- list("x" = x, "X" = X, "ntrials" = ntrials_all, "nsubs" = nsubs)

    params <- c(
        "mu_theta", "mu_lambda", "mu_phi", "mu_cons",
        "lambda_theta", "lambda_lambda", "lambda_phi", "lambda_cons"
    )

    samples <- jags.parallel(
        data = jags_data,
        inits = NULL,
        parameters.to.save = params,
        model.file = "analysis/models/eef.txt",
        n.chains = 3,
        n.iter = 3000,
        n.burnin = 1000,
        n.thin = 1
    )

    # -------------------------------------------------------------------------
    # Step 4: Extract Point Estimates
    # -------------------------------------------------------------------------
    Y <- samples$BUGSoutput$sims.list

    # True values
    true_mu_theta[i] <- mu_theta
    true_mu_lambda[i] <- mu_lambda
    true_mu_phi[i] <- mu_phi
    true_mu_cons[i] <- mu_cons

    true_sigma_theta[i] <- sigma_theta
    true_sigma_lambda[i] <- sigma_lambda
    true_sigma_phi[i] <- sigma_phi
    true_sigma_cons[i] <- sigma_cons

    # Inferred values (MPD)
    infer_mu_theta[i] <- MPD(Y$mu_theta)
    infer_mu_lambda[i] <- MPD(Y$mu_lambda)
    infer_mu_phi[i] <- MPD(Y$mu_phi)
    infer_mu_cons[i] <- MPD(Y$mu_cons)

    # Convert precision to SD
    infer_sigma_theta[i] <- precision_to_sd(MPD(Y$lambda_theta))
    infer_sigma_lambda[i] <- precision_to_sd(MPD(Y$lambda_lambda))
    infer_sigma_phi[i] <- precision_to_sd(MPD(Y$lambda_phi))
    infer_sigma_cons[i] <- precision_to_sd(MPD(Y$lambda_cons))

    cat("Iteration", i, "of", niterations, "complete\n")
}

end_time <- Sys.time()
cat("\nTotal time:", round(difftime(end_time, start_time, units = "mins"), 1), "minutes\n")

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

cat("\nRecovery plot saved to:", file.path(output_dir, "recovery_eef.png"), "\n")

# Print correlation summary
cat("\n=== Recovery Correlations ===\n")
cat("mu_theta:  r =", round(cor(true_mu_theta, infer_mu_theta), 3), "\n")
cat("mu_lambda: r =", round(cor(true_mu_lambda, infer_mu_lambda), 3), "\n")
cat("mu_phi:    r =", round(cor(true_mu_phi, infer_mu_phi), 3), "\n")
cat("mu_cons:   r =", round(cor(true_mu_cons, infer_mu_cons), 3), "\n")
