#!/usr/bin/env Rscript
# Generates recovery scatter plots from saved RDS files.
# Run with: Rscript scripts/plotting/plot_recovery_from_rds.R <model>
# Example: Rscript scripts/plotting/plot_recovery_from_rds.R pvl_delta

library(ggplot2)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript plot_recovery_from_rds.R <model>\nModel: eef, pvl_delta, orl")
}

model <- tolower(args[1])

rds_file <- paste0("data/processed/recovery/", model, "/recovery_", model, ".rds")
if (!file.exists(rds_file)) {
    stop("Recovery file not found: ", rds_file)
}

cat("Loading:", rds_file, "\n")
rec <- readRDS(rds_file)

output_dir <- "figures/paper"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Simple scatter plot with correlation
plot_scatter <- function(true_vals, infer_vals, param_name) {
    df <- data.frame(True = true_vals, Inferred = infer_vals)
    r <- cor(true_vals, infer_vals, use = "complete.obs")
    rmse <- sqrt(mean((true_vals - infer_vals)^2, na.rm = TRUE))

    ggplot(df, aes(x = True, y = Inferred)) +
        geom_point(alpha = 0.5, size = 2, color = "#2C3E50") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
        geom_smooth(method = "lm", se = TRUE, color = "#3498DB", fill = "lightblue", alpha = 0.3) +
        labs(
            title = param_name,
            subtitle = sprintf("r = %.2f, RMSE = %.2f", r, rmse),
            x = "True", y = "Recovered"
        ) +
        theme_minimal(base_size = 11) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
            plot.subtitle = element_text(hjust = 0.5, size = 12)
        )
}

# parameter names and display labels with Greek symbols
if (model == "eef") {
    mu_params <- c("mu_theta", "mu_lambda", "mu_phi", "mu_cons")
    mu_labels <- c(
        expression(paste("Value Sensitivity (", theta, ")")),
        expression(paste("Forgetting Rate (", lambda, ")")),
        expression(paste("Exploration Bonus (", phi, ")")),
        expression(paste("Consistency (", beta, ")"))
    )
    sigma_params <- c("sigma_theta", "sigma_lambda", "sigma_phi", "sigma_cons")
    sigma_labels <- c(
        expression(paste(sigma[theta])),
        expression(paste(sigma[lambda])),
        expression(paste(sigma[phi])),
        expression(paste(sigma[beta]))
    )
} else if (model == "pvl_delta") {
    mu_params <- c("mu_A", "mu_a", "mu_w", "mu_theta")
    mu_labels <- c(
        expression(paste("Outcome Sensitivity (", A, ")")),
        expression(paste("Learning Rate (", a, ")")),
        expression(paste("Loss Aversion (", w, ")")),
        expression(paste("Inverse Temperature (", theta, ")"))
    )
    sigma_params <- c("sigma_A", "sigma_a", "sigma_w", "sigma_theta")
    sigma_labels <- c(
        expression(paste(sigma[A])),
        expression(paste(sigma[a])),
        expression(paste(sigma[w])),
        expression(paste(sigma[theta]))
    )
} else if (model == "orl") {
    mu_params <- c("mu_a_rew", "mu_a_pun", "mu_K", "mu_omega_f", "mu_omega_p")
    mu_labels <- c(
        expression(paste("Reward LR (", A[rew], ")")),
        expression(paste("Punishment LR (", A[pun], ")")),
        expression(paste("Decay (", K, ")")),
        expression(paste("Frequency Weight (", omega[F], ")")),
        expression(paste("Perseverance Weight (", omega[P], ")"))
    )
    sigma_params <- c("sigma_a_rew", "sigma_a_pun", "sigma_K", "sigma_omega_f", "sigma_omega_p")
    sigma_labels <- c(
        expression(paste(sigma[A[rew]])),
        expression(paste(sigma[A[pun]])),
        expression(paste(sigma[K])),
        expression(paste(sigma[omega[F]])),
        expression(paste(sigma[omega[P]]))
    )
} else {
    stop("Unknown model: ", model)
}

# Mu recovery plots
mu_plots <- list()
for (i in seq_along(mu_params)) {
    p <- mu_params[i]
    true_col <- paste0("true_", p)
    infer_col <- paste0("infer_", p)

    true_vals <- sapply(rec, function(x) x[[true_col]])
    infer_vals <- sapply(rec, function(x) x[[infer_col]])

    if (length(true_vals) > 0 && !all(is.na(true_vals))) {
        mu_plots[[i]] <- plot_scatter(true_vals, infer_vals, mu_labels[i])
    }
}

if (length(mu_plots) > 0) {
    ncol_mu <- min(length(mu_plots), 2)
    nrow_mu <- ceiling(length(mu_plots) / ncol_mu)
    combined_mu <- ggarrange(plotlist = mu_plots, ncol = ncol_mu, nrow = nrow_mu)

    mu_file <- file.path(output_dir, paste0("recovery_", model, "_mu.png"))
    ggsave(mu_file, combined_mu, width = 10, height = 8, dpi = 600)
    cat("Saved:", mu_file, "\n")
}

# Sigma recovery plots
sigma_plots <- list()
for (i in seq_along(sigma_params)) {
    p <- sigma_params[i]
    true_col <- paste0("true_", p)
    infer_col <- paste0("infer_", p)

    true_vals <- sapply(rec, function(x) x[[true_col]])
    infer_vals <- sapply(rec, function(x) x[[infer_col]])

    if (length(true_vals) > 0 && !all(is.na(true_vals))) {
        sigma_plots[[length(sigma_plots) + 1]] <- plot_scatter(true_vals, infer_vals, sigma_labels[i])
    }
}

if (length(sigma_plots) > 0) {
    ncol_s <- min(length(sigma_plots), 2)
    nrow_s <- ceiling(length(sigma_plots) / ncol_s)
    combined_sigma <- ggarrange(plotlist = sigma_plots, ncol = ncol_s, nrow = nrow_s)

    sigma_file <- file.path(output_dir, paste0("recovery_", model, "_sigma.png"))
    ggsave(sigma_file, combined_sigma, width = 10, height = 8, dpi = 600)
    cat("Saved:", sigma_file, "\n")
}

cat("Done!\n")
