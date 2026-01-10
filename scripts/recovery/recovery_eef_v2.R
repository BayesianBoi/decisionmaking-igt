# =============================================================================
# Parameter Recovery for EEF v2 (Yang et al. 2025)
# =============================================================================
# Tests whether the JAGS model can accurately recover known parameter values
# from simulated data. Runs 100 iterations with varying true parameters.
# =============================================================================

# Dependencies
required_packages <- c("R2jags", "parallel", "ggpubr", "extraDistr", "truncnorm")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")
lapply(required_packages, library, character.only = TRUE)

set.seed(69420)

# Helper functions
MPD <- function(x) {
    density(x)$x[which(density(x)$y == max(density(x)$y))]
}

precision_to_sd <- function(lambda) {
    1 / sqrt(lambda)
}

source("scripts/recovery/simulation_eef_v2.R")
source("utils/payoff_scheme.R")

# Plotting function
plot_recovery <- function(true_vals, infer_vals, param_name) {
    df <- data.frame(true = true_vals, inferred = infer_vals)
    r_val <- round(cor(true_vals, infer_vals, use = "complete.obs"), 3)
    rmse_val <- round(sqrt(mean((true_vals - infer_vals)^2, na.rm = TRUE)), 3)

    ggplot(df, aes(x = true, y = inferred)) +
        geom_point(alpha = 0.6, size = 2) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
        labs(
            title = param_name,
            subtitle = paste0("r = ", r_val, ", RMSE = ", rmse_val),
            x = "True Value",
            y = "Recovered Value"
        ) +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"))
}

# Task Setup
ntrials <- 100
payoff_struct <- generate_modified_igt_payoff(ntrials, scale = TRUE)

cat("Modified IGT Payoff (Scaled /100) loaded.\n")
cat("Payoff structure contains separate gain and loss matrices.\n\n")

# Configuration
args <- commandArgs(trailingOnly = TRUE)
if ("--test" %in% args) {
    niterations <- 2
    cat(">>> RUNNING IN TEST MODE (2 iterations) <<<\n")
} else {
    niterations <- 100
}

nsubs <- 24
ntrials_all <- rep(100, nsubs)

# Main Recovery Loop
start_time <- Sys.time()
n_cores <- min(40, parallel::detectCores() - 1)
if (n_cores < 1) n_cores <- 1

cat("Starting EEF v2 Parameter Recovery...\n")
cat("Configuration:", niterations, "iterations x", nsubs, "subjects\n")
cat("Running on", n_cores, "cores\n\n")

run_iteration <- function(i) {
    tryCatch(
        {
            # Sample true group-level parameters
            mu_theta <- runif(1, 0.2, 0.8)
            sigma_theta <- runif(1, 0.05, 0.15)

            mu_lambda <- runif(1, 0.2, 0.8)
            sigma_lambda <- runif(1, 0.05, 0.15)

            mu_phi <- runif(1, -2, 2)
            sigma_phi <- runif(1, 0.1, 0.5)

            mu_cons <- runif(1, 1, 3)
            sigma_cons <- runif(1, 0.1, 0.5)

            # Simulate data
            sim_data <- simulation_eef_v2(
                payoff_struct = payoff_struct,
                nsubs = nsubs,
                ntrials = ntrials_all,
                mu_theta = mu_theta, mu_lambda = mu_lambda,
                mu_phi = mu_phi, mu_cons = mu_cons,
                sigma_theta = sigma_theta, sigma_lambda = sigma_lambda,
                sigma_phi = sigma_phi, sigma_cons = sigma_cons
            )

            # Scale gain/loss for JAGS (already scaled in payoff_struct)
            jags_data <- list(
                "x" = sim_data$x,
                "Gain" = sim_data$Gain,
                "Loss" = sim_data$Loss,
                "ntrials" = ntrials_all,
                "nsubs" = nsubs
            )

            params <- c(
                "mu_theta", "mu_lambda", "mu_phi", "mu_cons",
                "lambda_theta", "lambda_lambda", "lambda_phi", "lambda_cons"
            )

            samples <- jags(
                data = jags_data,
                inits = NULL,
                parameters.to.save = params,
                model.file = "models/eef_v2.txt",
                n.chains = 3,
                n.iter = 3000,
                n.burnin = 1000,
                n.thin = 1,
                progress.bar = "none"
            )

            Y <- samples$BUGSoutput$sims.list

            list(
                true_mu_theta = mu_theta, true_mu_lambda = mu_lambda,
                true_mu_phi = mu_phi, true_mu_cons = mu_cons,
                infer_mu_theta = MPD(Y$mu_theta), infer_mu_lambda = MPD(Y$mu_lambda),
                infer_mu_phi = MPD(Y$mu_phi), infer_mu_cons = MPD(Y$mu_cons),
                true_sigma_theta = sigma_theta, true_sigma_lambda = sigma_lambda,
                true_sigma_phi = sigma_phi, true_sigma_cons = sigma_cons,
                infer_sigma_theta = MPD(precision_to_sd(Y$lambda_theta)),
                infer_sigma_lambda = MPD(precision_to_sd(Y$lambda_lambda)),
                infer_sigma_phi = MPD(precision_to_sd(Y$lambda_phi)),
                infer_sigma_cons = MPD(precision_to_sd(Y$lambda_cons))
            )
        },
        error = function(e) {
            cat("Error in iteration", i, ":", conditionMessage(e), "\n")
            return(NULL)
        }
    )
}

results_list <- mclapply(1:niterations, run_iteration, mc.cores = n_cores, mc.set.seed = TRUE)

# Filter failed iterations
results_list <- results_list[!sapply(results_list, is.null)]

if (length(results_list) == 0) {
    stop("All iterations failed!")
}

cat("Successful iterations:", length(results_list), "/", niterations, "\n")

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

# Save Results
output_dir_data <- "outputs/recovery"
dir.create(output_dir_data, recursive = TRUE, showWarnings = FALSE)

saveRDS(results_list, file.path(output_dir_data, "recovery_eef_v2.rds"))

df_summary <- data.frame(
    true_mu_theta = true_mu_theta, infer_mu_theta = infer_mu_theta,
    true_mu_lambda = true_mu_lambda, infer_mu_lambda = infer_mu_lambda,
    true_mu_phi = true_mu_phi, infer_mu_phi = infer_mu_phi,
    true_mu_cons = true_mu_cons, infer_mu_cons = infer_mu_cons
)
write.csv(df_summary, file.path(output_dir_data, "recovery_eef_v2.csv"), row.names = FALSE)

cat("Results saved to:", output_dir_data, "\n")

# Plotting
pl1 <- plot_recovery(true_mu_theta, infer_mu_theta, "mu_theta (Sensitivity)")
pl2 <- plot_recovery(true_mu_lambda, infer_mu_lambda, "mu_lambda (Forgetting)")
pl3 <- plot_recovery(true_mu_phi, infer_mu_phi, "mu_phi (Exploration)")
pl4 <- plot_recovery(true_mu_cons, infer_mu_cons, "mu_cons (Consistency)")

final_plot <- ggarrange(pl1, pl2, pl3, pl4, nrow = 2, ncol = 2)

output_dir <- "plots/recovery"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output_dir, "recovery_eef_v2.png"), final_plot,
    width = 12, height = 10, dpi = 150
)

cat("Recovery plot saved to:", file.path(output_dir, "recovery_eef_v2.png"), "\n")

# Print correlation summary
cat("\n=== Recovery Correlations ===\n")
cat("mu_theta:  r =", round(cor(true_mu_theta, infer_mu_theta), 3), "\n")
cat("mu_lambda: r =", round(cor(true_mu_lambda, infer_mu_lambda), 3), "\n")
cat("mu_phi:    r =", round(cor(true_mu_phi, infer_mu_phi), 3), "\n")
cat("mu_cons:   r =", round(cor(true_mu_cons, infer_mu_cons), 3), "\n")
