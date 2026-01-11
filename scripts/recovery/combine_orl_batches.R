# ORL Parameter Recovery - Batch Results Combiner
# This script reads all partial .rds files from "outputs/recovery/parts_orl",
# combines them, and generates the final recovery plots and summary CSV.
#
# Usage:
#   Rscript combine_orl_batches.R

# dependencies
required_packages <- c("ggplot2", "ggpubr", "dplyr", "tidyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")
invisible(lapply(required_packages, library, character.only = TRUE))

source("scripts/plotting/plotting_utils.R")

# load and combine results
parts_dir <- "data/processed/recovery/parts_orl"
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

# extract vectors for plotting
extract_param <- function(res_list, param_name) {
    sapply(res_list, function(x) if (!is.null(x[[param_name]])) x[[param_name]] else NA)
}

# Means
true_mu_a_rew <- extract_param(combined_results, "true_mu_a_rew")
infer_mu_a_rew <- extract_param(combined_results, "infer_mu_a_rew")
true_mu_a_pun <- extract_param(combined_results, "true_mu_a_pun")
infer_mu_a_pun <- extract_param(combined_results, "infer_mu_a_pun")
true_mu_K <- extract_param(combined_results, "true_mu_K")
infer_mu_K <- extract_param(combined_results, "infer_mu_K")
true_mu_omega_f <- extract_param(combined_results, "true_mu_omega_f")
infer_mu_omega_f <- extract_param(combined_results, "infer_mu_omega_f")
true_mu_omega_p <- extract_param(combined_results, "true_mu_omega_p")
infer_mu_omega_p <- extract_param(combined_results, "infer_mu_omega_p")

# Sigmas
true_sigma_a_rew <- extract_param(combined_results, "true_sigma_a_rew")
infer_sigma_a_rew <- extract_param(combined_results, "infer_sigma_a_rew")
true_sigma_a_pun <- extract_param(combined_results, "true_sigma_a_pun")
infer_sigma_a_pun <- extract_param(combined_results, "infer_sigma_a_pun")
true_sigma_K <- extract_param(combined_results, "true_sigma_K")
infer_sigma_K <- extract_param(combined_results, "infer_sigma_K")
true_sigma_omega_f <- extract_param(combined_results, "true_sigma_omega_f")
infer_sigma_omega_f <- extract_param(combined_results, "infer_sigma_omega_f")
true_sigma_omega_p <- extract_param(combined_results, "true_sigma_omega_p")
infer_sigma_omega_p <- extract_param(combined_results, "infer_sigma_omega_p")

# CIs
lower_mu_a_rew <- extract_param(combined_results, "lower_mu_a_rew")
upper_mu_a_rew <- extract_param(combined_results, "upper_mu_a_rew")
lower_mu_a_pun <- extract_param(combined_results, "lower_mu_a_pun")
upper_mu_a_pun <- extract_param(combined_results, "upper_mu_a_pun")
lower_mu_K <- extract_param(combined_results, "lower_mu_K")
upper_mu_K <- extract_param(combined_results, "upper_mu_K")
lower_mu_omega_f <- extract_param(combined_results, "lower_mu_omega_f")
upper_mu_omega_f <- extract_param(combined_results, "upper_mu_omega_f")
lower_mu_omega_p <- extract_param(combined_results, "lower_mu_omega_p")
upper_mu_omega_p <- extract_param(combined_results, "upper_mu_omega_p")


# save csv
output_dir <- "data/processed/recovery/orl"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

df_summary <- data.frame(
    # Means
    true_mu_a_rew = true_mu_a_rew, infer_mu_a_rew = infer_mu_a_rew,
    true_mu_a_pun = true_mu_a_pun, infer_mu_a_pun = infer_mu_a_pun,
    true_mu_K = true_mu_K, infer_mu_K = infer_mu_K,
    true_mu_omega_f = true_mu_omega_f, infer_mu_omega_f = infer_mu_omega_f,
    true_mu_omega_p = true_mu_omega_p, infer_mu_omega_p = infer_mu_omega_p,

    # Sigmas
    true_sigma_a_rew = true_sigma_a_rew, infer_sigma_a_rew = infer_sigma_a_rew,
    true_sigma_a_pun = true_sigma_a_pun, infer_sigma_a_pun = infer_sigma_a_pun,
    true_sigma_K = true_sigma_K, infer_sigma_K = infer_sigma_K,
    true_sigma_omega_f = true_sigma_omega_f, infer_sigma_omega_f = infer_sigma_omega_f,
    true_sigma_omega_p = true_sigma_omega_p, infer_sigma_omega_p = infer_sigma_omega_p
)

csv_path <- file.path(output_dir, "recovery_orl.csv")
write.csv(df_summary, csv_path, row.names = FALSE)
cat("Summary CSV saved to:", csv_path, "\n")


# generate plots
plot_dir <- "figures/recovery"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

cat("Generating combined plots...\n")
# Means with meaningful titles
p1 <- plot_recovery(true_mu_a_rew, infer_mu_a_rew, "Reward Learning (A_rew)")
p2 <- plot_recovery(true_mu_a_pun, infer_mu_a_pun, "Punishment Learning (A_pun)")
p3 <- plot_recovery(true_mu_K, infer_mu_K, "Decay (K)")
p4 <- plot_recovery(true_mu_omega_f, infer_mu_omega_f, "Frequency Weight (omega_f)")
p5 <- plot_recovery(true_mu_omega_p, infer_mu_omega_p, "Perseverance (omega_p)")

# Sigmas with meaningful titles
s1 <- plot_recovery(true_sigma_a_rew, infer_sigma_a_rew, "Sigma A_rew")
s2 <- plot_recovery(true_sigma_a_pun, infer_sigma_a_pun, "Sigma A_pun")
s3 <- plot_recovery(true_sigma_K, infer_sigma_K, "Sigma K")
s4 <- plot_recovery(true_sigma_omega_f, infer_sigma_omega_f, "Sigma Omega_F")
s5 <- plot_recovery(true_sigma_omega_p, infer_sigma_omega_p, "Sigma Omega_P")

# Combined plot: 5x2 grid (means left, sigmas right)
combined_plot <- ggarrange(
    p1, s1, # Row 1: A_rew
    p2, s2, # Row 2: A_pun
    p3, s3, # Row 3: K
    p4, s4, # Row 4: omega_f
    p5, s5, # Row 5: omega_p
    nrow = 5, ncol = 2
)

ggsave(file.path(plot_dir, "recovery_orl_combined.png"), combined_plot, width = 8, height = 20)
print(paste("Combined plot saved to:", file.path(plot_dir, "recovery_orl_combined.png")))

# Mu-only plot (3x2 grid)
mu_plot <- ggarrange(
    p1, p2, p3, p4, p5, NULL,
    ncol = 3, nrow = 2
)
ggsave(file.path(plot_dir, "recovery_orl_mu.png"), mu_plot, width = 12, height = 8)
print(paste("Mu plot saved to:", file.path(plot_dir, "recovery_orl_mu.png")))

# Sigma-only plot (3x2 grid)
sigma_plot <- ggarrange(
    s1, s2, s3, s4, s5, NULL,
    ncol = 3, nrow = 2
)
ggsave(file.path(plot_dir, "recovery_orl_sigma.png"), sigma_plot, width = 12, height = 8)
print(paste("Sigma plot saved to:", file.path(plot_dir, "recovery_orl_sigma.png")))


# Print major correlations (Means)
cat("\n=== Recovery Correlations ===\n")
cat("mu_a_rew:   r =", round(cor(true_mu_a_rew, infer_mu_a_rew, use = "pairwise.complete.obs"), 2), "\n")
cat("mu_a_pun:   r =", round(cor(true_mu_a_pun, infer_mu_a_pun, use = "pairwise.complete.obs"), 2), "\n")
cat("mu_K:       r =", round(cor(true_mu_K, infer_mu_K, use = "pairwise.complete.obs"), 2), "\n")
cat("mu_omega_f: r =", round(cor(true_mu_omega_f, infer_mu_omega_f, use = "pairwise.complete.obs"), 2), "\n")
cat("mu_omega_p: r =", round(cor(true_mu_omega_p, infer_mu_omega_p, use = "pairwise.complete.obs"), 2), "\n")

# Helper to calc metrics
calc_metrics <- function(true_val, infer_val, lower, upper) {
    if (length(true_val) == 0) {
        return(c(NA, NA, NA))
    }
    ci_width <- mean(upper - lower, na.rm = TRUE)
    coverage <- mean(true_val >= lower & true_val <= upper, na.rm = TRUE) * 100
    rmse <- sqrt(mean((true_val - infer_val)^2, na.rm = TRUE))
    return(c(rmse, coverage, ci_width))
}

# Print Coverage if available
if (any(!is.na(lower_mu_a_rew))) {
    cat("\n=== Parameter Metrics (RMSE | Coverage | Width) ===\n")
    cat("parameter   | RMSE   | Coverage | Width\n")
    cat("---------------------------------------------\n")
    m1 <- calc_metrics(true_mu_a_rew, infer_mu_a_rew, lower_mu_a_rew, upper_mu_a_rew)
    cat(sprintf("mu_a_rew    | %.3f  | %.1f%%   | %.3f\n", m1[1], m1[2], m1[3]))

    m2 <- calc_metrics(true_mu_a_pun, infer_mu_a_pun, lower_mu_a_pun, upper_mu_a_pun)
    cat(sprintf("mu_a_pun    | %.3f  | %.1f%%   | %.3f\n", m2[1], m2[2], m2[3]))

    m3 <- calc_metrics(true_mu_K, infer_mu_K, lower_mu_K, upper_mu_K)
    cat(sprintf("mu_K        | %.3f  | %.1f%%   | %.3f\n", m3[1], m3[2], m3[3]))

    m5 <- calc_metrics(true_mu_omega_f, infer_mu_omega_f, lower_mu_omega_f, upper_mu_omega_f)
    cat(sprintf("mu_omega_f  | %.3f  | %.1f%%   | %.3f\n", m5[1], m5[2], m5[3]))

    m6 <- calc_metrics(true_mu_omega_p, infer_mu_omega_p, lower_mu_omega_p, upper_mu_omega_p)
    cat("---------------------------------------------\n")
}
