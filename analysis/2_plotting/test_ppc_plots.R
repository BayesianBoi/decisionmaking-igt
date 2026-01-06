#!/usr/bin/env Rscript
# ==============================================================================
# Test PPC Plots Using Pseudo Data
# ==============================================================================
# Generates PPC plots using pseudo data to verify plotting works correctly.
#
# Usage: Rscript analysis/2_plotting/test_ppc_plots.R
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, patchwork)

# ==============================================================================
# Load Pseudo PPC Data
# ==============================================================================
ppc_file <- "analysis/outputs/ppc/pseudo_ppc_results.rds"

if (!file.exists(ppc_file)) {
    cat("Generating pseudo PPC data first...\n")
    source("analysis/utils/generate_pseudo_ppc.R")
}

ppc_results <- readRDS(ppc_file)
cat("Loaded pseudo PPC data.\n")

# ==============================================================================
# Configuration
# ==============================================================================
plot_dir <- "analysis/plots/ppc"
if (dir.exists(plot_dir)) {
    unlink(plot_dir, recursive = TRUE)
    cat("Cleared existing PPC plot folder.\n")
}
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Colors matching Example2.pdf style
point_color <- "gray30"
band_color <- "#F4A8A8" # Light pink/salmon
mean_line_color <- "black"
chance_line_color <- "black"

# ==============================================================================
# Plot Function
# ==============================================================================

plot_ppc_single <- function(ppc_data, model_name) {
    n_subj <- length(ppc_data$subject_accuracy)

    df <- data.frame(
        subject = 1:n_subj,
        accuracy = ppc_data$subject_accuracy
    )

    overall_mean <- ppc_data$overall_mean
    overall_sd <- ppc_data$overall_sd

    p <- ggplot(df, aes(x = subject, y = accuracy)) +
        # Pink ±1 SD band
        geom_rect(
            aes(
                xmin = -Inf, xmax = Inf,
                ymin = overall_mean - overall_sd,
                ymax = overall_mean + overall_sd
            ),
            fill = band_color, alpha = 0.5
        ) +
        # Mean line
        geom_hline(yintercept = overall_mean, color = mean_line_color, linewidth = 1) +
        # Chance line (dashed)
        geom_hline(yintercept = 0.25, linetype = "dashed", color = chance_line_color) +
        # Per-subject points
        geom_point(color = point_color, size = 2) +
        labs(
            title = toupper(model_name),
            x = "Subject",
            y = "Prediction Accuracy"
        ) +
        coord_cartesian(ylim = c(0.1, 0.7)) +
        theme_minimal(base_size = 11) +
        theme(
            panel.grid.major = element_line(color = "gray90"),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold")
        )

    return(p)
}

# ==============================================================================
# Generate Individual Plots
# ==============================================================================

cat("\nGenerating PPC plots...\n")

plots <- list()
for (model in names(ppc_results)) {
    cat(sprintf("  %s...\n", model))

    p <- plot_ppc_single(ppc_results[[model]], model)
    plots[[model]] <- p

    # Save individual
    ggsave(
        file.path(plot_dir, sprintf("ppc_%s.png", model)),
        p,
        width = 8, height = 4, dpi = 300
    )
    ggsave(
        file.path(plot_dir, sprintf("ppc_%s.pdf", model)),
        p,
        width = 8, height = 4
    )
}

# ==============================================================================
# Combined Comparison Plot
# ==============================================================================

cat("Generating combined comparison plot...\n")

p_combined <- plots[[1]] + plots[[2]] + plots[[3]] +
    plot_layout(ncol = 3) +
    plot_annotation(
        title = "Posterior Predictive Checks: Model Comparison",
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
    )

ggsave(
    file.path(plot_dir, "ppc_comparison.png"),
    p_combined,
    width = 14, height = 4, dpi = 300
)
ggsave(
    file.path(plot_dir, "ppc_comparison.pdf"),
    p_combined,
    width = 14, height = 4
)

cat("\n=== PPC Plots Generated ===\n")
cat("Output directory:", plot_dir, "\n")
cat("Individual plots:\n")
for (model in names(ppc_results)) {
    cat(sprintf("  - ppc_%s.png\n", model))
}
cat("Combined: ppc_comparison.png\n")
