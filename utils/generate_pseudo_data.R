# ==============================================================================
# Generate Pseudo Recovery Data for Plot Testing
# ==============================================================================
# Run this script to create fake recovery results that can be used to
# test and refine plotting utilities without waiting for real recovery.
#
# Usage: Rscript utils/generate_pseudo_data.R
# ==============================================================================

set.seed(42)

# Define number of recovery iterations to simulate
n_iterations <- 100

# ==============================================================================
# PVL-Delta Pseudo Recovery
# ==============================================================================
cat("Generating PVL-Delta pseudo recovery data...\n")

# True values drawn from the simulation ranges
true_mu_w <- runif(n_iterations, 0.5, 2.5)
true_mu_A <- runif(n_iterations, 0.2, 0.8)
true_mu_theta <- runif(n_iterations, 0.5, 2.0)
true_mu_a <- runif(n_iterations, 0.1, 0.5)

# Simulate "recovered" values with noise + some correlation
# Realistic recovery: varies from r~0.3 (poor) to r~0.9 (good)
add_recovery_noise <- function(true_vals, r = 0.7, noise_sd = 0.15) {
    noise <- rnorm(length(true_vals), 0, noise_sd)
    recovered <- r * true_vals + (1 - r) * mean(true_vals) + noise
    return(recovered)
}

# PVL: Some parameters recover well, theta is often tricky
infer_mu_w <- add_recovery_noise(true_mu_w, r = 0.75, noise_sd = 0.25)
infer_mu_A <- add_recovery_noise(true_mu_A, r = 0.70, noise_sd = 0.08)
infer_mu_theta <- add_recovery_noise(true_mu_theta, r = 0.50, noise_sd = 0.20) # Poor
infer_mu_a <- add_recovery_noise(true_mu_a, r = 0.80, noise_sd = 0.05)

df_pvl <- data.frame(
    true_mu_w = true_mu_w, infer_mu_w = infer_mu_w,
    true_mu_A = true_mu_A, infer_mu_A = infer_mu_A,
    true_mu_theta = true_mu_theta, infer_mu_theta = infer_mu_theta,
    true_mu_a = true_mu_a, infer_mu_a = infer_mu_a
)

# ==============================================================================
# ORL Pseudo Recovery
# ==============================================================================
cat("Generating ORL pseudo recovery data...\n")

true_mu_a_rew <- runif(n_iterations, 0.1, 0.5)
true_mu_a_pun <- runif(n_iterations, 0.1, 0.5)
true_mu_K <- runif(n_iterations, 0.5, 2.0)
true_mu_theta_orl <- runif(n_iterations, 0.5, 2.0)
true_mu_omega_f <- runif(n_iterations, -1, 1)
true_mu_omega_p <- runif(n_iterations, -1, 1)

# ORL: More parameters, more variability (like in the Haines paper)
infer_mu_a_rew <- add_recovery_noise(true_mu_a_rew, r = 0.65, noise_sd = 0.06)
infer_mu_a_pun <- add_recovery_noise(true_mu_a_pun, r = 0.60, noise_sd = 0.06)
infer_mu_K <- add_recovery_noise(true_mu_K, r = 0.45, noise_sd = 0.20) # K is notoriously hard
infer_mu_theta_orl <- add_recovery_noise(true_mu_theta_orl, r = 0.40, noise_sd = 0.25) # Poor
infer_mu_omega_f <- add_recovery_noise(true_mu_omega_f, r = 0.75, noise_sd = 0.12)
infer_mu_omega_p <- add_recovery_noise(true_mu_omega_p, r = 0.70, noise_sd = 0.15)

df_orl <- data.frame(
    true_mu_a_rew = true_mu_a_rew, infer_mu_a_rew = infer_mu_a_rew,
    true_mu_a_pun = true_mu_a_pun, infer_mu_a_pun = infer_mu_a_pun,
    true_mu_K = true_mu_K, infer_mu_K = infer_mu_K,
    true_mu_theta = true_mu_theta_orl, infer_mu_theta = infer_mu_theta_orl,
    true_mu_omega_f = true_mu_omega_f, infer_mu_omega_f = infer_mu_omega_f,
    true_mu_omega_p = true_mu_omega_p, infer_mu_omega_p = infer_mu_omega_p
)

# ==============================================================================
# EEF Pseudo Recovery
# ==============================================================================
cat("Generating EEF pseudo recovery data...\n")

true_mu_theta_eef <- runif(n_iterations, 0.2, 0.8)
true_mu_lambda <- runif(n_iterations, 0.2, 0.8)
true_mu_phi <- runif(n_iterations, -2, 2)
true_mu_cons <- runif(n_iterations, 1, 3)

# EEF: Newer model, mixed recovery
infer_mu_theta_eef <- add_recovery_noise(true_mu_theta_eef, r = 0.70, noise_sd = 0.08)
infer_mu_lambda <- add_recovery_noise(true_mu_lambda, r = 0.55, noise_sd = 0.10) # Moderate
infer_mu_phi <- add_recovery_noise(true_mu_phi, r = 0.80, noise_sd = 0.20)
infer_mu_cons <- add_recovery_noise(true_mu_cons, r = 0.50, noise_sd = 0.25) # Poor

df_eef <- data.frame(
    true_mu_theta = true_mu_theta_eef, infer_mu_theta = infer_mu_theta_eef,
    true_mu_lambda = true_mu_lambda, infer_mu_lambda = infer_mu_lambda,
    true_mu_phi = true_mu_phi, infer_mu_phi = infer_mu_phi,
    true_mu_cons = true_mu_cons, infer_mu_cons = infer_mu_cons
)

# ==============================================================================
# Save Pseudo Data
# ==============================================================================
output_dir <- "analysis/outputs/recovery"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(df_pvl, file.path(output_dir, "pseudo_recovery_pvl_delta.csv"), row.names = FALSE)
write.csv(df_orl, file.path(output_dir, "pseudo_recovery_orl.csv"), row.names = FALSE)
write.csv(df_eef, file.path(output_dir, "pseudo_recovery_eef.csv"), row.names = FALSE)

cat("\n=== Pseudo Recovery Data Generated ===\n")
cat("Saved to:\n")
cat("  -", file.path(output_dir, "pseudo_recovery_pvl_delta.csv"), "\n")
cat("  -", file.path(output_dir, "pseudo_recovery_orl.csv"), "\n")
cat("  -", file.path(output_dir, "pseudo_recovery_eef.csv"), "\n")

cat("\nCorrelations (simulated):\n")
cat("PVL mu_w:     r =", round(cor(true_mu_w, infer_mu_w), 3), "\n")
cat("PVL mu_A:     r =", round(cor(true_mu_A, infer_mu_A), 3), "\n")
cat("PVL mu_theta: r =", round(cor(true_mu_theta, infer_mu_theta), 3), "\n")
cat("PVL mu_a:     r =", round(cor(true_mu_a, infer_mu_a), 3), "\n")
