# ==============================================================================
# PVL-Delta Parameter Recovery - Batch Results Combiner
# ==============================================================================
# This script reads all partial .rds files from "analysis/outputs/recovery/parts_pvl_delta",
# combines them, and generates the final recovery plots and summary CSV.
#
# Usage:
#   Rscript combine_pvl_delta_batches.R
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Dependencies & Utils
# ------------------------------------------------------------------------------
required_packages <- c("ggplot2", "ggpubr", "dplyr", "tidyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")
invisible(lapply(required_packages, library, character.only = TRUE))

source("analysis/utils/plotting_utils.R")

# ------------------------------------------------------------------------------
# 2. Load and Combine Results
# ------------------------------------------------------------------------------
parts_dir <- "analysis/outputs/recovery/parts_pvl_delta"
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

# ------------------------------------------------------------------------------
# 3. Extract Vectors for Plotting (Means + Sigmas)
# ------------------------------------------------------------------------------
extract_param <- function(res_list, param_name) {
    sapply(res_list, function(x) if (!is.null(x[[param_name]])) x[[param_name]] else NA)
}

# Means
true_mu_w <- extract_param(combined_results, "true_mu_w")
infer_mu_w <- extract_param(combined_results, "infer_mu_w")
true_mu_A <- extract_param(combined_results, "true_mu_A")
infer_mu_A <- extract_param(combined_results, "infer_mu_A")
true_mu_a <- extract_param(combined_results, "true_mu_a")
infer_mu_a <- extract_param(combined_results, "infer_mu_a")
true_mu_theta <- extract_param(combined_results, "true_mu_theta")
infer_mu_theta <- extract_param(combined_results, "infer_mu_theta")

# Sigmas
true_sigma_w <- extract_param(combined_results, "true_sigma_w")
infer_sigma_w <- extract_param(combined_results, "infer_sigma_w")
true_sigma_A <- extract_param(combined_results, "true_sigma_A")
infer_sigma_A <- extract_param(combined_results, "infer_sigma_A")
true_sigma_a <- extract_param(combined_results, "true_sigma_a")
infer_sigma_a <- extract_param(combined_results, "infer_sigma_a")
true_sigma_theta <- extract_param(combined_results, "true_sigma_theta")
infer_sigma_theta <- extract_param(combined_results, "infer_sigma_theta")

# CIs
lower_mu_w <- extract_param(combined_results, "lower_mu_w")
upper_mu_w <- extract_param(combined_results, "upper_mu_w")
lower_mu_A <- extract_param(combined_results, "lower_mu_A")
upper_mu_A <- extract_param(combined_results, "upper_mu_A")
lower_mu_a <- extract_param(combined_results, "lower_mu_a")
upper_mu_a <- extract_param(combined_results, "upper_mu_a")
lower_mu_theta <- extract_param(combined_results, "lower_mu_theta")
upper_mu_theta <- extract_param(combined_results, "upper_mu_theta")

# ------------------------------------------------------------------------------
# 4. Save Summary CSV
# ------------------------------------------------------------------------------
output_dir <- "analysis/outputs/recovery"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df_summary <- data.frame(
    # Means
    true_mu_w = true_mu_w, infer_mu_w = infer_mu_w,
    true_mu_A = true_mu_A, infer_mu_A = infer_mu_A,
    true_mu_a = true_mu_a, infer_mu_a = infer_mu_a,
    true_mu_theta = true_mu_theta, infer_mu_theta = infer_mu_theta,

    # Sigmas
    true_sigma_w = true_sigma_w, infer_sigma_w = infer_sigma_w,
    true_sigma_A = true_sigma_A, infer_sigma_A = infer_sigma_A,
    true_sigma_a = true_sigma_a, infer_sigma_a = infer_sigma_a,
    true_sigma_theta = true_sigma_theta, infer_sigma_theta = infer_sigma_theta
)

csv_path <- file.path(output_dir, "recovery_pvl_delta.csv")
write.csv(df_summary, csv_path, row.names = FALSE)
cat("Summary CSV saved to:", csv_path, "\n")


# ------------------------------------------------------------------------------
# 5. Generate Plots
# ------------------------------------------------------------------------------
plot_dir <- "analysis/plots/recovery"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Means
p1 <- plot_recovery(true_mu_w, infer_mu_w, "Loss Aversion (w)")
p2 <- plot_recovery(true_mu_A, infer_mu_A, "Outcome Sensitivity (A)")
p3 <- plot_recovery(true_mu_a, infer_mu_a, "Learning Rate (a)")
p4 <- plot_recovery(true_mu_theta, infer_mu_theta, "Response Consistency (theta)")

# Sigmas
s1 <- plot_recovery(true_sigma_w, infer_sigma_w, "Sigma w")
s2 <- plot_recovery(true_sigma_A, infer_sigma_A, "Sigma A")
s3 <- plot_recovery(true_sigma_a, infer_sigma_a, "Sigma a")
s4 <- plot_recovery(true_sigma_theta, infer_sigma_theta, "Sigma theta")

# Combine all into one 2x4 grid
combined_plot <- ggarrange(
    p1, p2, p3, p4,
    s1, s2, s3, s4,
    ncol = 4, nrow = 2
)

ggsave(file.path(plot_dir, "recovery_pvl_delta_combined.png"), combined_plot, width = 16, height = 8)
print(paste("Combined plot saved to:", file.path(plot_dir, "recovery_pvl_delta_combined.png")))


# ------------------------------------------------------------------------------
# 6. Calculate and Print Metrics
# ------------------------------------------------------------------------------
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

cat("\n=== Recovery Correlations (Means) ===\n")
cat("mu_w:     r =", round(cor(true_mu_w, infer_mu_w, use = "pairwise.complete.obs"), 2), "\n")
cat("mu_A:     r =", round(cor(true_mu_A, infer_mu_A, use = "pairwise.complete.obs"), 2), "\n")
cat("mu_a:     r =", round(cor(true_mu_a, infer_mu_a, use = "pairwise.complete.obs"), 2), "\n")
cat("mu_theta: r =", round(cor(true_mu_theta, infer_mu_theta, use = "pairwise.complete.obs"), 2), "\n")

# Print Coverage if available
if (any(!is.na(lower_mu_w))) {
    cat("\n=== Recovery Performance Metrics ===\n")
    cat("parameter     | RMSE   | Coverage | Correlation\n")
    cat("---------------------------------------------\n")
    m1 <- calc_metrics(true_mu_w, infer_mu_w, lower_mu_w, upper_mu_w)
    cat(sprintf("mu_w          | %.3f  | %.1f%%   | %.3f\n", m1[1], m1[2], m1[3]))

    m2 <- calc_metrics(true_mu_A, infer_mu_A, lower_mu_A, upper_mu_A)
    cat(sprintf("mu_A          | %.3f  | %.1f%%   | %.3f\n", m2[1], m2[2], m2[3]))

    m3 <- calc_metrics(true_mu_a, infer_mu_a, lower_mu_a, upper_mu_a)
    cat(sprintf("mu_a          | %.3f  | %.1f%%   | %.3f\n", m3[1], m3[2], m3[3]))

    m4 <- calc_metrics(true_mu_theta, infer_mu_theta, lower_mu_theta, upper_mu_theta)
    cat(sprintf("mu_theta      | %.3f  | %.1f%%   | %.3f\n", m4[1], m4[2], m4[3]))
    cat("---------------------------------------------\n")
}
