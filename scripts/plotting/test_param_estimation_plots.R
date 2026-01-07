#!/usr/bin/env Rscript
# ==============================================================================
# Test Parameter Estimation Plots with Pseudo Data
# ==============================================================================
# Generates pseudo MCMC samples and plots posterior densities + group comparisons
# without needing real model fits.
#
# Usage: Rscript analysis/2_plotting/test_param_estimation_plots.R
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, patchwork)

set.seed(42)

# ==============================================================================
# Configuration
# ==============================================================================
plot_dir <- "plots/param_estimation"
if (dir.exists(plot_dir)) {
    unlink(plot_dir, recursive = TRUE)
    cat("Cleared existing param estimation plot folder.\n")
}
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Group colors
group_colors <- c(
    "HC" = "#fb9a99",
    "Amphetamine" = "#1f78b4",
    "Heroin" = "#542788"
)

# ==============================================================================
# Generate Pseudo Posterior Samples
# ==============================================================================

generate_pseudo_posterior <- function(n_samples = 4000) {
    # Simulate group-level posteriors for each model
    # These represent draws from the posterior distribution

    list(
        # PVL-Delta parameters
        pvl_delta = list(
            mu_w = rnorm(n_samples, mean = 1.2, sd = 0.3),
            mu_A = rnorm(n_samples, mean = 0.6, sd = 0.15),
            mu_theta = rnorm(n_samples, mean = 2.5, sd = 0.5),
            mu_a = rnorm(n_samples, mean = 0.4, sd = 0.1)
        ),
        # ORL parameters
        orl = list(
            mu_a_rew = rnorm(n_samples, mean = 0.35, sd = 0.08),
            mu_a_pun = rnorm(n_samples, mean = 0.25, sd = 0.06),
            mu_K = rnorm(n_samples, mean = 2.0, sd = 0.4),
            mu_theta = rnorm(n_samples, mean = 1.5, sd = 0.3),
            mu_omega_f = rnorm(n_samples, mean = 0.8, sd = 0.2),
            mu_omega_p = rnorm(n_samples, mean = -0.5, sd = 0.15)
        ),
        # EEF parameters
        eef = list(
            mu_theta = rnorm(n_samples, mean = 0.5, sd = 0.1),
            mu_lambda = rnorm(n_samples, mean = 0.3, sd = 0.08),
            mu_phi = rnorm(n_samples, mean = 0.2, sd = 0.15),
            mu_cons = rnorm(n_samples, mean = 2.0, sd = 0.4)
        )
    )
}

generate_pseudo_group_posteriors <- function(n_samples = 4000) {
    # Simulate group-level posteriors for group comparison
    # HC (control), Amphetamine, Heroin

    list(
        HC = list(
            mu_lambda = rnorm(n_samples, mean = 0.28, sd = 0.06),
            mu_theta = rnorm(n_samples, mean = 0.52, sd = 0.08),
            mu_phi = rnorm(n_samples, mean = 0.25, sd = 0.10),
            mu_cons = rnorm(n_samples, mean = 2.1, sd = 0.35)
        ),
        Amphetamine = list(
            mu_lambda = rnorm(n_samples, mean = 0.38, sd = 0.08), # Higher = faster forgetting
            mu_theta = rnorm(n_samples, mean = 0.48, sd = 0.09),
            mu_phi = rnorm(n_samples, mean = 0.18, sd = 0.12), # Lower exploration
            mu_cons = rnorm(n_samples, mean = 1.8, sd = 0.40)
        ),
        Heroin = list(
            mu_lambda = rnorm(n_samples, mean = 0.42, sd = 0.10), # Even higher forgetting
            mu_theta = rnorm(n_samples, mean = 0.45, sd = 0.10),
            mu_phi = rnorm(n_samples, mean = 0.12, sd = 0.08), # Much lower exploration
            mu_cons = rnorm(n_samples, mean = 1.5, sd = 0.45) # Lower consistency
        )
    )
}

cat("=== Generating Pseudo Parameter Estimation Data ===\n\n")

posterior_samples <- generate_pseudo_posterior()
group_posteriors <- generate_pseudo_group_posteriors()

# ==============================================================================
# PLOT 1: Posterior Densities (Single Model)
# ==============================================================================

plot_posterior_density <- function(samples, param_name, model_name) {
    df <- data.frame(value = samples)

    dens <- density(samples)
    n_samples <- length(samples)
    bw <- dens$bw

    ggplot(df, aes(x = value)) +
        geom_density(fill = NA, color = "black", linewidth = 0.8) +
        labs(
            title = param_name,
            x = NULL,
            y = "Density",
            caption = sprintf("N = %d   Bandwidth = %.4f", n_samples, bw)
        ) +
        theme_minimal(base_size = 10) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "plain", size = 12),
            plot.caption = element_text(hjust = 0.5, size = 8),
            panel.grid.major = element_line(color = "gray90"),
            panel.grid.minor = element_blank(),
            axis.title.y = element_text(size = 10)
        )
}

cat("Generating posterior density plots...\n")

for (model in names(posterior_samples)) {
    cat(sprintf("  %s...\n", model))

    params <- posterior_samples[[model]]
    plots <- list()

    for (param in names(params)) {
        plots[[param]] <- plot_posterior_density(params[[param]], param, model)
    }

    # Combine
    n_params <- length(plots)
    ncol <- min(4, n_params)

    combined <- wrap_plots(plots, ncol = ncol) +
        plot_annotation(
            title = sprintf("Posterior Distributions: %s", toupper(model)),
            theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
        )

    ggsave(
        file.path(plot_dir, sprintf("posterior_%s.png", model)),
        combined,
        width = 10, height = 6, dpi = 300
    )
}

# ==============================================================================
# PLOT 2: Group Comparison Densities (Overlapping)
# ==============================================================================

plot_group_comparison <- function(group_posteriors, param_name) {
    df <- do.call(rbind, lapply(names(group_posteriors), function(g) {
        data.frame(
            value = group_posteriors[[g]][[param_name]],
            Group = g
        )
    }))
    df$Group <- factor(df$Group, levels = c("HC", "Amphetamine", "Heroin"))

    ggplot(df, aes(x = value, color = Group)) +
        geom_density(fill = NA, linewidth = 0.8) +
        scale_color_manual(values = group_colors) +
        labs(
            title = param_name,
            x = NULL,
            y = "Density"
        ) +
        theme_minimal(base_size = 11) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "plain", size = 12),
            panel.grid.major = element_line(color = "gray90"),
            panel.grid.minor = element_blank(),
            legend.position = "none"
        )
}

cat("\nGenerating group comparison plots...\n")

param_names <- names(group_posteriors$HC)
group_plots <- lapply(param_names, function(p) plot_group_comparison(group_posteriors, p))
names(group_plots) <- param_names

# Add legend to the last plot only
group_plots[[length(group_plots)]] <- group_plots[[length(group_plots)]] +
    theme(legend.position = "bottom") +
    guides(
        fill = guide_legend(title = "Group", nrow = 1),
        color = guide_legend(title = "Group", nrow = 1)
    )

combined_group <- wrap_plots(group_plots, ncol = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(
        title = "Group Comparison: EEF Parameters",
        subtitle = "HC (Control) vs Amphetamine vs Heroin",
        theme = theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 11)
        )
    ) &
    theme(legend.position = "bottom")

ggsave(
    file.path(plot_dir, "group_comparison_densities.png"),
    combined_group,
    width = 10, height = 8, dpi = 300
)

# ==============================================================================
# PLOT 3: HDI Difference Plot (Substance - Control)
# ==============================================================================

cat("Generating HDI difference plots...\n")

compute_hdi <- function(samples, prob = 0.95) {
    sorted <- sort(samples)
    n <- length(sorted)
    ci_size <- ceiling(prob * n)

    widths <- sorted[(ci_size + 1):n] - sorted[1:(n - ci_size)]
    min_idx <- which.min(widths)

    c(lower = sorted[min_idx], upper = sorted[min_idx + ci_size])
}

diff_plots <- list()
comparisons <- list(
    "Heroin - HC" = c("Heroin", "HC"),
    "Amphetamine - HC" = c("Amphetamine", "HC")
)

for (name in names(comparisons)) {
    g1 <- comparisons[[name]][1]
    g2 <- comparisons[[name]][2]

    # Difference in mu_lambda (forgetting rate)
    diff <- group_posteriors[[g1]]$mu_lambda - group_posteriors[[g2]]$mu_lambda

    hdi <- compute_hdi(diff)
    mean_diff <- mean(diff)
    prob_pos <- mean(diff > 0)

    df <- data.frame(value = diff)

    diff_plots[[name]] <- ggplot(df, aes(x = value)) +
        geom_density(fill = NA, color = "black", linewidth = 0.8) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
        labs(
            title = name,
            x = expression(Delta * mu[lambda]),
            y = "Density",
            caption = sprintf(
                "Mean = %.3f, 95%% HDI [%.3f, %.3f], P(diff > 0) = %.2f",
                mean_diff, hdi["lower"], hdi["upper"], prob_pos
            )
        ) +
        theme_minimal(base_size = 10) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "plain", size = 12),
            plot.caption = element_text(hjust = 0.5, size = 8),
            panel.grid.major = element_line(color = "gray90"),
            panel.grid.minor = element_blank()
        )
}

combined_diff <- wrap_plots(diff_plots, ncol = 2) +
    plot_annotation(
        title = "Group Differences in Forgetting Rate",
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    )

ggsave(
    file.path(plot_dir, "hdi_differences.png"),
    combined_diff,
    width = 10, height = 4, dpi = 300
)

cat("\n=== Parameter Estimation Plots Generated ===\n")
cat("Output directory:", plot_dir, "\n")
cat("Files:\n")
cat("  - posterior_pvl_delta.png\n")
cat("  - posterior_orl.png\n")
cat("  - posterior_eef.png\n")
cat("  - group_comparison_densities.png\n")
cat("  - hdi_differences.png\n")
