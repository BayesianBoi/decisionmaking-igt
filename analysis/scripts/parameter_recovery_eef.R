# Parameter recovery analysis for EEF model
# Simulates data with known parameters, fits model, checks recovery quality
#
# Run: Rscript analysis/scripts/parameter_recovery_eef.R

library(rjags)
library(coda)

#==============================================================================
# CONFIGURATION
#==============================================================================

# Simulation settings
n_subjects <- 50
n_trials <- 100
n_decks <- 4

# True parameter values (population means)
true_params <- list(
  theta = 0.35,
  lambda_forget = 0.50,
  phi = 0.60,
  cons = 1.50
)

# Individual variability (SD around population mean)
param_sd <- list(
  theta = 0.10,
  lambda_forget = 0.10,
  phi = 0.30,
  cons = 0.40
)

# MCMC settings (shorter for recovery test)
config <- list(
  n_adapt = 2000,
  n_burnin = 5000,
  n_iter = 10000,
  n_chains = 3,
  thin = 2,

  parameters_to_monitor = c(
    "mu_theta",
    "mu_lambda_forget",
    "mu_phi",
    "mu_cons",
    "theta",
    "lambda_forget",
    "phi",
    "cons"
  )
)

# Output directory
output_dir <- "results/parameter_recovery"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set seed for reproducibility
set.seed(42)

#==============================================================================
# DECK PAYOFFS (Standard IGT structure)
#==============================================================================

# Net outcomes per deck (simplified version)
# A: high immediate reward, higher frequent loss
# B: high immediate reward, low infrequent loss
# C: low immediate reward, higher frequent loss
# D: low immediate reward, low infrequent loss

# For simplicity, use deterministic outcomes
deck_outcomes <- list(
  A = rep(c(1.00, 1.00, 1.00, 1.00, 1.00, -2.50), length.out = n_trials),
  B = rep(c(1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, -12.50), length.out = n_trials),
  C = rep(c(0.50, 0.50, 0.50, 0.50, 0.50, -0.50), length.out = n_trials),
  D = rep(c(0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50), length.out = n_trials)
)

#==============================================================================
# SIMULATE DATA
#==============================================================================

message("=== SIMULATING DATA ===\n")

# Draw individual parameters from population distribution
subject_params <- data.frame(
  theta = pmax(0.01, pmin(0.99, rnorm(n_subjects, true_params$theta, param_sd$theta))),
  lambda_forget = pmax(0.01, pmin(0.99, rnorm(n_subjects, true_params$lambda_forget, param_sd$lambda_forget))),
  phi = pmax(-5, pmin(5, rnorm(n_subjects, true_params$phi, param_sd$phi))),
  cons = pmax(0.01, pmin(5, rnorm(n_subjects, true_params$cons, param_sd$cons)))
)

message(sprintf("Drew parameters for %d subjects", n_subjects))
message("Population means:")
message(sprintf("  theta: %.3f (true: %.3f)", mean(subject_params$theta), true_params$theta))
message(sprintf("  lambda_forget: %.3f (true: %.3f)", mean(subject_params$lambda_forget), true_params$lambda_forget))
message(sprintf("  phi: %.3f (true: %.3f)", mean(subject_params$phi), true_params$phi))
message(sprintf("  cons: %.3f (true: %.3f)", mean(subject_params$cons), true_params$cons))

# Initialize data structures
choice_mat <- matrix(NA, nrow = n_subjects, ncol = n_trials)
outcome_mat <- matrix(NA, nrow = n_subjects, ncol = n_trials)

# Simulate each subject
message("\nSimulating choices...")

for (s in 1:n_subjects) {

  # Extract subject parameters
  theta_s <- subject_params$theta[s]
  lambda_s <- subject_params$lambda_forget[s]
  phi_s <- subject_params$phi[s]
  cons_s <- subject_params$cons[s]

  # Initialize exploitation and exploration weights
  exploit <- rep(0, n_decks)
  explore <- rep(0, n_decks)

  for (t in 1:n_trials) {

    # Compute choice probabilities (softmax)
    weight <- exploit + explore
    v <- cons_s * weight
    exp_v <- exp(v - max(v))  # Subtract max for numerical stability
    p <- exp_v / sum(exp_v)

    # Sample choice
    choice <- sample(1:n_decks, size = 1, prob = p)
    choice_mat[s, t] <- choice

    # Get outcome for chosen deck
    outcome <- deck_outcomes[[choice]][t]
    outcome_mat[s, t] <- outcome

    # Compute utility
    gain <- max(outcome, 0)
    loss <- max(-outcome, 0)
    util <- (gain + 0.001)^theta_s - (loss + 0.001)^theta_s

    # Update exploitation weights
    for (d in 1:n_decks) {
      if (d == choice) {
        exploit[d] <- (1 - lambda_s) * exploit[d] + util
      } else {
        exploit[d] <- (1 - lambda_s) * exploit[d]
      }
    }

    # Update exploration weights
    for (d in 1:n_decks) {
      if (d == choice) {
        explore[d] <- 0
      } else {
        explore[d] <- lambda_s * explore[d] + (1 - lambda_s) * phi_s
      }
    }
  }

  if (s %% 10 == 0) {
    message(sprintf("  Simulated %d/%d subjects", s, n_subjects))
  }
}

message("Simulation complete.\n")

# Create JAGS data
w_ini <- rep(0, n_decks)  # Uniform initial weights for simulated data

jags_data <- list(
  N = n_subjects,
  Tsubj = rep(n_trials, n_subjects),
  choice = choice_mat,
  outcome = outcome_mat,
  w_ini = w_ini
)

# Save simulated data and true parameters
saveRDS(list(
  jags_data = jags_data,
  true_params = subject_params,
  true_population = true_params
), file.path(output_dir, "simulated_data.rds"))

#==============================================================================
# FIT MODEL TO SIMULATED DATA
#==============================================================================

message("=== FITTING EEF MODEL TO SIMULATED DATA ===\n")

# Load JAGS model
model_file <- "analysis/models/eef_clinical.jags"
if (!file.exists(model_file)) {
  stop(sprintf("Model file not found: %s", model_file))
}

# Initialize JAGS model
message(sprintf("Initializing model with %d chains...", config$n_chains))
start_time <- Sys.time()

jags_model <- jags.model(
  file = model_file,
  data = jags_data,
  n.chains = config$n_chains,
  n.adapt = config$n_adapt,
  quiet = FALSE
)

message("Adaptation complete.")

# Burn-in
message(sprintf("Burn-in: %d iterations...", config$n_burnin))
update(jags_model, n.iter = config$n_burnin, progress.bar = "text")

# Sampling
message(sprintf("Sampling: %d iterations x %d chains (thin=%d)...",
                config$n_iter, config$n_chains, config$thin))

samples <- coda.samples(
  model = jags_model,
  variable.names = config$parameters_to_monitor,
  n.iter = config$n_iter,
  thin = config$thin,
  progress.bar = "text"
)

total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
message(sprintf("Fitting complete (%.1f minutes)", total_time))

# Save samples
saveRDS(samples, file.path(output_dir, "recovery_samples.rds"))

#==============================================================================
# EVALUATE RECOVERY
#==============================================================================

message("\n=== PARAMETER RECOVERY EVALUATION ===\n")

# Extract posterior means
samples_mat <- as.matrix(samples)

# Subject-level recovery
param_names <- c("theta", "lambda_forget", "phi", "cons")
recovery_results <- list()

for (param in param_names) {

  # Extract recovered values (posterior means)
  param_cols <- grep(sprintf("^%s\\[", param), colnames(samples_mat), value = TRUE)
  recovered <- colMeans(samples_mat[, param_cols])

  # True values
  true <- subject_params[[param]]

  # Compute correlation
  corr <- cor(true, recovered)

  # Compute RMSE
  rmse <- sqrt(mean((true - recovered)^2))

  # Compute bias
  bias <- mean(recovered - true)

  recovery_results[[param]] <- list(
    true = true,
    recovered = recovered,
    correlation = corr,
    rmse = rmse,
    bias = bias
  )

  message(sprintf("%s recovery:", param))
  message(sprintf("  Correlation: %.3f", corr))
  message(sprintf("  RMSE: %.3f", rmse))
  message(sprintf("  Bias: %.3f", bias))
}

# Population-level recovery
message("\nPopulation parameter recovery:")
for (param in param_names) {
  mu_param <- sprintf("mu_%s", param)
  recovered_pop <- mean(samples_mat[, mu_param])
  true_pop <- true_params[[param]]

  message(sprintf("  %s: true=%.3f, recovered=%.3f, diff=%.3f",
                  param, true_pop, recovered_pop, recovered_pop - true_pop))
}

# Save recovery results
saveRDS(recovery_results, file.path(output_dir, "recovery_results.rds"))

#==============================================================================
# VISUALIZATIONS
#==============================================================================

message("\n=== CREATING RECOVERY PLOTS ===\n")

pdf(file.path(output_dir, "parameter_recovery.pdf"), width = 12, height = 10)
par(mfrow = c(2, 2))

for (param in param_names) {

  res <- recovery_results[[param]]

  # Scatter plot with identity line
  plot(res$true, res$recovered,
       xlab = sprintf("True %s", param),
       ylab = sprintf("Recovered %s", param),
       main = sprintf("%s Recovery (r=%.3f)", param, res$correlation),
       pch = 19, col = adjustcolor("steelblue", alpha = 0.6),
       xlim = range(c(res$true, res$recovered)),
       ylim = range(c(res$true, res$recovered)))

  # Add identity line
  abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)

  # Add best fit line
  fit <- lm(res$recovered ~ res$true)
  abline(fit, col = "blue", lwd = 2)

  # Add legend
  legend("topleft",
         legend = c(sprintf("r = %.3f", res$correlation),
                    sprintf("RMSE = %.3f", res$rmse),
                    sprintf("Bias = %.3f", res$bias)),
         bty = "n")
}

dev.off()

message(sprintf("Recovery plots saved to: %s",
                file.path(output_dir, "parameter_recovery.pdf")))

#==============================================================================
# SUMMARY
#==============================================================================

message("\n=== PARAMETER RECOVERY COMPLETE ===\n")
message(sprintf("Output directory: %s", output_dir))
message("Files created:")
message("  - simulated_data.rds (simulated data and true parameters)")
message("  - recovery_samples.rds (MCMC samples from recovery fit)")
message("  - recovery_results.rds (correlation, RMSE, bias)")
message("  - parameter_recovery.pdf (scatter plots)")
message("\nInterpretation:")
message("  - Correlation > 0.80: Excellent recovery")
message("  - Correlation 0.60-0.80: Good recovery")
message("  - Correlation < 0.60: Poor recovery (model may not be identifiable)")
