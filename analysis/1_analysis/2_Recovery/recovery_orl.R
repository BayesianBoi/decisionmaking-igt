# ==============================================================================
# Parameter Recovery: ORL (Outcome-Representation Learning) Model
# ==============================================================================
#
# Purpose: Validate that our JAGS implementation can accurately recover
# known parameter values from simulated data.
#
# The ORL model has 6 parameters:
#   - a_rew: Learning rate for rewards (0-1)
#   - a_pun: Learning rate for punishments (0-1)
#   - K: Perseverance decay parameter (>= 0)
#   - theta: Inverse temperature (>= 0)
#   - omega_f: Weight on expected frequency (unbounded)
#   - omega_p: Weight on perseverance (unbounded)
#
# Success Criteria:
#   - Strong correlations between true and recovered parameters (r > 0.7)
#   - No systematic bias
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

source("analysis/1_analysis/2_Recovery/simulation_orl.R")
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
# Group means
true_mu_a_rew <- infer_mu_a_rew <- array(NA, c(niterations))
true_mu_a_pun <- infer_mu_a_pun <- array(NA, c(niterations))
true_mu_K <- infer_mu_K <- array(NA, c(niterations))
true_mu_theta <- infer_mu_theta <- array(NA, c(niterations))
true_mu_omega_f <- infer_mu_omega_f <- array(NA, c(niterations))
true_mu_omega_p <- infer_mu_omega_p <- array(NA, c(niterations))

# Group SDs (converted from precision)
true_sigma_a_rew <- infer_sigma_a_rew <- array(NA, c(niterations))
true_sigma_a_pun <- infer_sigma_a_pun <- array(NA, c(niterations))
true_sigma_K <- infer_sigma_K <- array(NA, c(niterations))
true_sigma_theta <- infer_sigma_theta <- array(NA, c(niterations))
true_sigma_omega_f <- infer_sigma_omega_f <- array(NA, c(niterations))
true_sigma_omega_p <- infer_sigma_omega_p <- array(NA, c(niterations))

# ==============================================================================
# Main Recovery Loop
# ==============================================================================
start_time <- Sys.time()
cat("Starting ORL Parameter Recovery...\n")
cat("Configuration:", niterations, "iterations x", nsubs, "subjects\n\n")

for (i in 1:niterations) {
    # -------------------------------------------------------------------------
    # Step 1: Generate True Parameters
    # -------------------------------------------------------------------------
    # Learning rates (0-1)
    mu_a_rew <- runif(1, 0.1, 0.5)
    mu_a_pun <- runif(1, 0.1, 0.5)
    sigma_a_rew <- runif(1, 0.05, 0.15)
    sigma_a_pun <- runif(1, 0.05, 0.15)

    # Decay parameter K >= 0
    mu_K <- runif(1, 0.5, 2.0)
    sigma_K <- runif(1, 0.1, 0.3)

    # Inverse temperature theta >= 0
    mu_theta <- runif(1, 0.5, 2.0)
    sigma_theta <- runif(1, 0.1, 0.3)

    # Omega weights (can be negative!)
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

    x <- sim_data$x
    X <- sim_data$X

    # -------------------------------------------------------------------------
    # Step 3: Fit JAGS Model
    # -------------------------------------------------------------------------
    jags_data <- list("x" = x, "X" = X, "ntrials" = ntrials_all, "nsubs" = nsubs)

    params <- c(
        "mu_a_rew", "mu_a_pun", "mu_K", "mu_theta", "mu_omega_f", "mu_omega_p",
        "lambda_a_rew", "lambda_a_pun", "lambda_K", "lambda_theta",
        "lambda_omega_f", "lambda_omega_p"
    )

    samples <- jags.parallel(
        data = jags_data,
        inits = NULL,
        parameters.to.save = params,
        model.file = "analysis/models/orl.txt",
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
    true_mu_a_rew[i] <- mu_a_rew
    true_mu_a_pun[i] <- mu_a_pun
    true_mu_K[i] <- mu_K
    true_mu_theta[i] <- mu_theta
    true_mu_omega_f[i] <- mu_omega_f
    true_mu_omega_p[i] <- mu_omega_p

    true_sigma_a_rew[i] <- sigma_a_rew
    true_sigma_a_pun[i] <- sigma_a_pun
    true_sigma_K[i] <- sigma_K
    true_sigma_theta[i] <- sigma_theta
    true_sigma_omega_f[i] <- sigma_omega_f
    true_sigma_omega_p[i] <- sigma_omega_p

    # Inferred values (MPD)
    infer_mu_a_rew[i] <- MPD(Y$mu_a_rew)
    infer_mu_a_pun[i] <- MPD(Y$mu_a_pun)
    infer_mu_K[i] <- MPD(Y$mu_K)
    infer_mu_theta[i] <- MPD(Y$mu_theta)
    infer_mu_omega_f[i] <- MPD(Y$mu_omega_f)
    infer_mu_omega_p[i] <- MPD(Y$mu_omega_p)

    # Convert precision to SD
    infer_sigma_a_rew[i] <- precision_to_sd(MPD(Y$lambda_a_rew))
    infer_sigma_a_pun[i] <- precision_to_sd(MPD(Y$lambda_a_pun))
    infer_sigma_K[i] <- precision_to_sd(MPD(Y$lambda_K))
    infer_sigma_theta[i] <- precision_to_sd(MPD(Y$lambda_theta))
    infer_sigma_omega_f[i] <- precision_to_sd(MPD(Y$lambda_omega_f))
    infer_sigma_omega_p[i] <- precision_to_sd(MPD(Y$lambda_omega_p))

    cat("Iteration", i, "of", niterations, "complete\n")
}

end_time <- Sys.time()
cat("\nTotal time:", round(difftime(end_time, start_time, units = "mins"), 1), "minutes\n")

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

cat("\nRecovery plot saved to:", file.path(output_dir, "recovery_orl.png"), "\n")

# Print correlation summary
cat("\n=== Recovery Correlations ===\n")
cat("mu_a_rew:   r =", round(cor(true_mu_a_rew, infer_mu_a_rew), 3), "\n")
cat("mu_a_pun:   r =", round(cor(true_mu_a_pun, infer_mu_a_pun), 3), "\n")
cat("mu_K:       r =", round(cor(true_mu_K, infer_mu_K), 3), "\n")
cat("mu_theta:   r =", round(cor(true_mu_theta, infer_mu_theta), 3), "\n")
cat("mu_omega_f: r =", round(cor(true_mu_omega_f, infer_mu_omega_f), 3), "\n")
cat("mu_omega_p: r =", round(cor(true_mu_omega_p, infer_mu_omega_p), 3), "\n")
