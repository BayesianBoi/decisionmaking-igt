# EEF Parameter Recovery - Batch Results Combiner
# This script reads all partial .rds files from "outputs/recovery/parts",
# combines them, and generates the final recovery plots and summary CSV.
#
# Usage:
#   a1

# dependencies
# Base R package loading
required_packages <- c("ggplot2", "ggpubr", "dplyr", "tidyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")
invisible(lapply(required_packages, library, character.only = TRUE))

source("scripts/plotting/plotting_utils.R")

# load and combine results
parts_dir <- "data/processed/recovery/parts_eef"
files <- list.files(parts_dir, pattern = "\\.rds$", full.names = TRUE)

if (length(files) == 0) {
    stop(paste("No .rds files found in", parts_dir))
}

cat("Found", length(files), "batch files. Combining...\n")

combined_results <- list()

for (f in files) {
    res <- readRDS(f)
    if (is.list(res) && length(res) > 0) {
        combined_results <- c(combined_results, res)
    } else {
        warning(paste("File", f, "was empty or invalid."))
    }
}

N <- length(combined_results)
cat("Total valid iterations recovered:", N, "\n")

if (N == 0) stop("No valid data to process.")

# extract params
extract_param <- function(res_list, param_name) {
    # Robust extraction: returns NA if param missing (e.g. old batch files)
    sapply(res_list, function(x) if (!is.null(x[[param_name]])) x[[param_name]] else NA)
}

# Means
true_mu_theta <- extract_param(combined_results, "true_mu_theta")
infer_mu_theta <- extract_param(combined_results, "infer_mu_theta")

true_mu_lambda <- extract_param(combined_results, "true_mu_lambda")
infer_mu_lambda <- extract_param(combined_results, "infer_mu_lambda")

true_mu_phi <- extract_param(combined_results, "true_mu_phi")
infer_mu_phi <- extract_param(combined_results, "infer_mu_phi")

true_mu_cons <- extract_param(combined_results, "true_mu_cons")
infer_mu_cons <- extract_param(combined_results, "infer_mu_cons")

# Sigmas
true_sigma_theta <- extract_param(combined_results, "true_sigma_theta")
infer_sigma_theta <- extract_param(combined_results, "infer_sigma_theta")

true_sigma_lambda <- extract_param(combined_results, "true_sigma_lambda")
infer_sigma_lambda <- extract_param(combined_results, "infer_sigma_lambda")

true_sigma_phi <- extract_param(combined_results, "true_sigma_phi")
infer_sigma_phi <- extract_param(combined_results, "infer_sigma_phi")

true_sigma_cons <- extract_param(combined_results, "true_sigma_cons")
infer_sigma_cons <- extract_param(combined_results, "infer_sigma_cons")

# save csv
output_dir <- "data/processed/recovery/eef"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df_summary <- data.frame(
    # Means
    true_mu_theta = true_mu_theta, infer_mu_theta = infer_mu_theta,
    true_mu_lambda = true_mu_lambda, infer_mu_lambda = infer_mu_lambda,
    true_mu_phi = true_mu_phi, infer_mu_phi = infer_mu_phi,
    true_mu_cons = true_mu_cons, infer_mu_cons = infer_mu_cons,

    # Sigmas
    true_sigma_theta = true_sigma_theta, infer_sigma_theta = infer_sigma_theta,
    true_sigma_lambda = true_sigma_lambda, infer_sigma_lambda = infer_sigma_lambda,
    true_sigma_phi = true_sigma_phi, infer_sigma_phi = infer_sigma_phi,
    true_sigma_cons = true_sigma_cons, infer_sigma_cons = infer_sigma_cons
)

csv_path <- file.path(output_dir, "recovery_eef.csv")
write.csv(df_summary, csv_path, row.names = FALSE)
cat("Summary CSV saved to:", csv_path, "\n")

# plots
plot_dir <- "figures/recovery"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Means
p_theta <- plot_recovery(true_mu_theta, infer_mu_theta, "Value Sensitivity (theta)")
p_lambda <- plot_recovery(true_mu_lambda, infer_mu_lambda, "Forgetting (lambda)")
p_phi <- plot_recovery(true_mu_phi, infer_mu_phi, "Exploration (phi)")
p_cons <- plot_recovery(true_mu_cons, infer_mu_cons, "Consistency (c)")

# Sigmas
p_sigma_theta <- plot_recovery(true_sigma_theta, infer_sigma_theta, "Sigma Theta")
p_sigma_lambda <- plot_recovery(true_sigma_lambda, infer_sigma_lambda, "Sigma Lambda")
p_sigma_phi <- plot_recovery(true_sigma_phi, infer_sigma_phi, "Sigma Phi")
p_sigma_cons <- plot_recovery(true_sigma_cons, infer_sigma_cons, "Sigma Consistency")

# Combine all into one 2x4 grid
combined_plot <- ggarrange(
    p_theta, p_lambda, p_phi, p_cons,
    p_sigma_theta, p_sigma_lambda, p_sigma_phi, p_sigma_cons,
    ncol = 4, nrow = 2
)

ggsave(file.path(plot_dir, "recovery_eef_combined.png"), combined_plot, width = 16, height = 8)
print(paste("Combined plot saved to:", file.path(plot_dir, "recovery_eef_combined.png")))

# Mu-only plot (2x2 grid)
mu_plot <- ggarrange(
    p_theta, p_lambda, p_phi, p_cons,
    ncol = 2, nrow = 2
)
ggsave(file.path(plot_dir, "recovery_eef_mu.png"), mu_plot, width = 10, height = 10)
print(paste("Mu plot saved to:", file.path(plot_dir, "recovery_eef_mu.png")))

# Sigma-only plot (2x2 grid)
sigma_plot <- ggarrange(
    p_sigma_theta, p_sigma_lambda, p_sigma_phi, p_sigma_cons,
    ncol = 2, nrow = 2
)
ggsave(file.path(plot_dir, "recovery_eef_sigma.png"), sigma_plot, width = 10, height = 10)
print(paste("Sigma plot saved to:", file.path(plot_dir, "recovery_eef_sigma.png")))

# Print correlations

valid_idx <- !is.na(true_sigma_theta) # Re-evaluate valid_idx for correlations
if (sum(valid_idx) > 0) {
    # Print correlations
    cat("\n=== Recovery Correlations (Sigmas) ===\n")
    cat("sigma_theta:  r =", round(cor(true_sigma_theta[valid_idx], infer_sigma_theta[valid_idx]), 2), "\n")
    cat("sigma_lambda: r =", round(cor(true_sigma_lambda[valid_idx], infer_sigma_lambda[valid_idx]), 2), "\n")
    cat("sigma_phi:    r =", round(cor(true_sigma_phi[valid_idx], infer_sigma_phi[valid_idx]), 2), "\n")
    cat("sigma_cons:   r =", round(cor(true_sigma_cons[valid_idx], infer_sigma_cons[valid_idx]), 2), "\n")
} else {
    cat("No valid sigma data found (please re-run batch recovery).\n")
}

cat("\n=== Recovery Correlations (Means) ===\n")
cat("mu_theta:  r =", round(cor(true_mu_theta, infer_mu_theta, use = "pairwise.complete.obs"), 2), "\n")
cat("mu_lambda: r =", round(cor(true_mu_lambda, infer_mu_lambda, use = "pairwise.complete.obs"), 2), "\n")
cat("mu_phi:    r =", round(cor(true_mu_phi, infer_mu_phi, use = "pairwise.complete.obs"), 2), "\n")
cat("mu_cons:   r =", round(cor(true_mu_cons, infer_mu_cons, use = "pairwise.complete.obs"), 2), "\n")

# Extract CIs
lower_mu_theta <- extract_param(combined_results, "lower_mu_theta")
upper_mu_theta <- extract_param(combined_results, "upper_mu_theta")
lower_mu_lambda <- extract_param(combined_results, "lower_mu_lambda")
upper_mu_lambda <- extract_param(combined_results, "upper_mu_lambda")
lower_mu_phi <- extract_param(combined_results, "lower_mu_phi")
upper_mu_phi <- extract_param(combined_results, "upper_mu_phi")
lower_mu_cons <- extract_param(combined_results, "lower_mu_cons")
upper_mu_cons <- extract_param(combined_results, "upper_mu_cons")

# Helper to calc metrics
calc_metrics <- function(true, infer, lower, upper) {
    valid <- !is.na(true) & !is.na(lower) & !is.na(upper)
    true <- true[valid]
    infer <- infer[valid]
    lower <- lower[valid]
    upper <- upper[valid]

    rmse <- sqrt(mean((true - infer)^2))
    coverage <- mean(true >= lower & true <= upper) * 100
    r <- cor(true, infer)
    return(c(RMSE = rmse, Coverage = coverage, R = r))
}

# Print Coverage if available
if (any(!is.na(lower_mu_theta))) {
    cat("\n=== Recovery Performance Metrics (EEF) ===\n")
    cat("parameter     | RMSE   | Coverage | Correlation\n")
    cat("---------------------------------------------\n")
    m1 <- calc_metrics(true_mu_theta, infer_mu_theta, lower_mu_theta, upper_mu_theta)
    cat(sprintf("mu_theta      | %.3f  | %.1f%%   | %.3f\n", m1[1], m1[2], m1[3]))

    m2 <- calc_metrics(true_mu_lambda, infer_mu_lambda, lower_mu_lambda, upper_mu_lambda)
    cat(sprintf("mu_lambda     | %.3f  | %.1f%%   | %.3f\n", m2[1], m2[2], m2[3]))

    m3 <- calc_metrics(true_mu_phi, infer_mu_phi, lower_mu_phi, upper_mu_phi)
    cat(sprintf("mu_phi        | %.3f  | %.1f%%   | %.3f\n", m3[1], m3[2], m3[3]))

    m4 <- calc_metrics(true_mu_cons, infer_mu_cons, lower_mu_cons, upper_mu_cons)
    cat(sprintf("mu_cons       | %.3f  | %.1f%%   | %.3f\n", m4[1], m4[2], m4[3]))
    cat("---------------------------------------------\n")
}
