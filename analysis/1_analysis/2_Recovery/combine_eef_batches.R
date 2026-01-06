# ==============================================================================
# EEF Parameter Recovery - Batch Results Combiner
# ==============================================================================
# This script reads all partial .rds files from "analysis/outputs/recovery/parts",
# combines them, and generates the final recovery plots and summary CSV.
#
# Usage:
#   Rscript combine_eef_batches.R
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Dependencies & Utils
# ------------------------------------------------------------------------------
# Base R package loading
required_packages <- c("ggplot2", "ggpubr", "dplyr", "tidyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")
invisible(lapply(required_packages, library, character.only = TRUE))

source("analysis/utils/plotting_utils.R")

# ------------------------------------------------------------------------------
# 2. Load and Combine Results
# ------------------------------------------------------------------------------
parts_dir <- "analysis/outputs/recovery/parts"
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
# 3. Extract Vectors for Plotting
# ------------------------------------------------------------------------------
true_mu_theta <- sapply(combined_results, function(x) x$true_mu_theta)
infer_mu_theta <- sapply(combined_results, function(x) x$infer_mu_theta)

true_mu_lambda <- sapply(combined_results, function(x) x$true_mu_lambda)
infer_mu_lambda <- sapply(combined_results, function(x) x$infer_mu_lambda)

true_mu_phi <- sapply(combined_results, function(x) x$true_mu_phi)
infer_mu_phi <- sapply(combined_results, function(x) x$infer_mu_phi)

true_mu_cons <- sapply(combined_results, function(x) x$true_mu_cons)
infer_mu_cons <- sapply(combined_results, function(x) x$infer_mu_cons)

# ------------------------------------------------------------------------------
# 4. Save Summary CSV
# ------------------------------------------------------------------------------
output_dir <- "analysis/outputs/recovery"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df_summary <- data.frame(
    true_mu_theta = true_mu_theta, infer_mu_theta = infer_mu_theta,
    true_mu_lambda = true_mu_lambda, infer_mu_lambda = infer_mu_lambda,
    true_mu_phi = true_mu_phi, infer_mu_phi = infer_mu_phi,
    true_mu_cons = true_mu_cons, infer_mu_cons = infer_mu_cons
)

csv_path <- file.path(output_dir, "recovery_eef.csv")
write.csv(df_summary, csv_path, row.names = FALSE)
cat("Summary CSV saved to:", csv_path, "\n")

# ------------------------------------------------------------------------------
# 5. Generate Plots
# ------------------------------------------------------------------------------
cat("Generating plots...\n")

p1 <- plot_recovery(true_mu_theta, infer_mu_theta, "Sensitivity (theta)")
p2 <- plot_recovery(true_mu_lambda, infer_mu_lambda, "Decay Rate (lambda)")
p3 <- plot_recovery(true_mu_phi, infer_mu_phi, "Exploration (phi)")
p4 <- plot_recovery(true_mu_cons, infer_mu_cons, "Consistency (cons)")

final_plot <- combine_plots(list(p1, p2, p3, p4), "EEF Parameter Recovery (Distributed)")

plot_dir <- "analysis/plots/recovery"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(plot_dir, "recovery_eef.png"), final_plot, width = 10, height = 8)

cat("Recovery plot saved to:", file.path(plot_dir, "recovery_eef.png"), "\n")

# Print correlations
cat("\n=== Recovery Correlations ===\n")
cat("mu_theta:  r =", round(cor(true_mu_theta, infer_mu_theta), 2), "\n")
cat("mu_lambda: r =", round(cor(true_mu_lambda, infer_mu_lambda), 2), "\n")
cat("mu_phi:    r =", round(cor(true_mu_phi, infer_mu_phi), 2), "\n")
cat("mu_cons:   r =", round(cor(true_mu_cons, infer_mu_cons), 2), "\n")
