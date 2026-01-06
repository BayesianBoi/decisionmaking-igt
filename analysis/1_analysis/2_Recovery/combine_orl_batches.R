# ==============================================================================
# ORL Parameter Recovery - Batch Results Combiner
# ==============================================================================
# This script reads all partial .rds files from "analysis/outputs/recovery/parts_orl",
# combines them, and generates the final recovery plots and summary CSV.
#
# Usage:
#   Rscript combine_orl_batches.R
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
parts_dir <- "analysis/outputs/recovery/parts_orl"
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
true_mu_a_rew <- sapply(combined_results, function(x) x$true_mu_a_rew)
infer_mu_a_rew <- sapply(combined_results, function(x) x$infer_mu_a_rew)

true_mu_a_pun <- sapply(combined_results, function(x) x$true_mu_a_pun)
infer_mu_a_pun <- sapply(combined_results, function(x) x$infer_mu_a_pun)

true_mu_K <- sapply(combined_results, function(x) x$true_mu_K)
infer_mu_K <- sapply(combined_results, function(x) x$infer_mu_K)

true_mu_theta <- sapply(combined_results, function(x) x$true_mu_theta)
infer_mu_theta <- sapply(combined_results, function(x) x$infer_mu_theta)

true_mu_omega_f <- sapply(combined_results, function(x) x$true_mu_omega_f)
infer_mu_omega_f <- sapply(combined_results, function(x) x$infer_mu_omega_f)

true_mu_omega_p <- sapply(combined_results, function(x) x$true_mu_omega_p)
infer_mu_omega_p <- sapply(combined_results, function(x) x$infer_mu_omega_p)

# ------------------------------------------------------------------------------
# 4. Save Summary CSV
# ------------------------------------------------------------------------------
output_dir <- "analysis/outputs/recovery"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df_summary <- data.frame(
    true_mu_a_rew = true_mu_a_rew, infer_mu_a_rew = infer_mu_a_rew,
    true_mu_a_pun = true_mu_a_pun, infer_mu_a_pun = infer_mu_a_pun,
    true_mu_K = true_mu_K, infer_mu_K = infer_mu_K,
    true_mu_theta = true_mu_theta, infer_mu_theta = infer_mu_theta,
    true_mu_omega_f = true_mu_omega_f, infer_mu_omega_f = infer_mu_omega_f,
    true_mu_omega_p = true_mu_omega_p, infer_mu_omega_p = infer_mu_omega_p
)

csv_path <- file.path(output_dir, "recovery_orl.csv")
write.csv(df_summary, csv_path, row.names = FALSE)
cat("Summary CSV saved to:", csv_path, "\n")

# ------------------------------------------------------------------------------
# 5. Generate Plots
# ------------------------------------------------------------------------------
cat("Generating plots...\n")

pl1 <- plot_recovery(true_mu_a_rew, infer_mu_a_rew, "mu_a_rew (Reward Learning)")
pl2 <- plot_recovery(true_mu_a_pun, infer_mu_a_pun, "mu_a_pun (Punishment Learning)")
pl3 <- plot_recovery(true_mu_K, infer_mu_K, "mu_K (Perseverance Decay)")
pl4 <- plot_recovery(true_mu_theta, infer_mu_theta, "mu_theta (Inv. Temperature)")
pl5 <- plot_recovery(true_mu_omega_f, infer_mu_omega_f, "mu_omega_f (Frequency Weight)")
pl6 <- plot_recovery(true_mu_omega_p, infer_mu_omega_p, "mu_omega_p (Perseverance Weight)")

# Use ggarrange directly for 6 plots (2x3 grid)
final_plot <- ggarrange(pl1, pl2, pl3, pl4, pl5, pl6, nrow = 2, ncol = 3)
final_plot <- annotate_figure(final_plot, top = text_grob("ORL Parameter Recovery (Distributed)", face = "bold", size = 14))

plot_dir <- "analysis/plots/recovery"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(plot_dir, "recovery_orl.png"), final_plot, width = 15, height = 10)

cat("Recovery plot saved to:", file.path(plot_dir, "recovery_orl.png"), "\n")

# Print correlations
cat("\n=== Recovery Correlations ===\n")
cat("mu_a_rew:   r =", round(cor(true_mu_a_rew, infer_mu_a_rew), 2), "\n")
cat("mu_a_pun:   r =", round(cor(true_mu_a_pun, infer_mu_a_pun), 2), "\n")
cat("mu_K:       r =", round(cor(true_mu_K, infer_mu_K), 2), "\n")
cat("mu_theta:   r =", round(cor(true_mu_theta, infer_mu_theta), 2), "\n")
cat("mu_omega_f: r =", round(cor(true_mu_omega_f, infer_mu_omega_f), 2), "\n")
cat("mu_omega_p: r =", round(cor(true_mu_omega_p, infer_mu_omega_p), 2), "\n")
