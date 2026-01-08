#!/usr/bin/env Rscript
# ==============================================================================
# PPC SUMMARY PLOTS
# ==============================================================================
# Creates one combined plot per model showing all groups side by side
# Usage: Rscript scripts/plotting/plot_ppc_summary.R
# ==============================================================================

library(ggplot2)
library(dplyr)

cat("=== Generating PPC Summary Plots ===\n\n")

# Settings
models <- c("orl", "pvl_delta", "eef")
groups <- c("HC", "Amph", "Hero")
group_labels <- c("HC" = "Healthy Controls", "Amph" = "Amphetamine", "Hero" = "Heroin")
colors <- c("HC" = "#4CAF50", "Amph" = "#FF5722", "Hero" = "#2196F3")

output_dir <- "outputs/ppc"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Generate one plot per model
for (m in models) {
    cat("Processing model:", toupper(m), "\n")

    model_data <- data.frame()
    summary_stats <- data.frame()

    # Load data for all groups
    for (g in groups) {
        file <- paste0("outputs/ppc/ppc_individual_", m, "_", g, ".rds")
        if (file.exists(file)) {
            ppc <- readRDS(file)
            df <- ppc$df
            df$Group <- g
            df$GroupLabel <- group_labels[g]
            model_data <- rbind(model_data, df)

            # Calculate summary
            n_subjects <- nrow(df)
            summary_stats <- rbind(summary_stats, data.frame(
                Group = g,
                GroupLabel = group_labels[g],
                Mean = ppc$mean,
                SD = ppc$sd,
                N = n_subjects
            ))
        }
    }

    if (nrow(model_data) == 0) {
        cat("  No data found, skipping.\n")
        next
    }

    # Set factor order
    model_data$GroupLabel <- factor(model_data$GroupLabel,
        levels = c("Healthy Controls", "Amphetamine", "Heroin")
    )
    summary_stats$GroupLabel <- factor(summary_stats$GroupLabel,
        levels = c("Healthy Controls", "Amphetamine", "Heroin")
    )

    # Create plot
    p <- ggplot(model_data, aes(x = subject, y = accuracy)) +
        # Chance level line
        geom_hline(yintercept = 0.25, linetype = "dashed", color = "black", linewidth = 0.8) +
        # SD ribbon
        geom_ribbon(
            data = model_data %>%
                left_join(summary_stats, by = c("Group", "GroupLabel")) %>%
                group_by(GroupLabel) %>%
                mutate(xmin = 0, xmax = max(subject)),
            aes(ymin = Mean - SD, ymax = Mean + SD, fill = Group),
            alpha = 0.3
        ) +
        # Mean line
        geom_hline(
            data = summary_stats,
            aes(yintercept = Mean, color = Group),
            linewidth = 1.2
        ) +
        # Individual points
        geom_point(aes(color = Group), alpha = 0.6, size = 2.5) +
        # Facet by group
        facet_wrap(~GroupLabel, scales = "free_x", nrow = 1) +
        # Colors
        scale_color_manual(values = colors) +
        scale_fill_manual(values = colors) +
        # Axis limits
        ylim(0, 0.85) +
        # Labels
        labs(
            title = paste("Posterior Predictive Check:", toupper(m)),
            subtitle = sprintf(
                "Mean accuracy: HC = %.1f%%, Amph = %.1f%%, Hero = %.1f%% | Dashed line = chance (25%%)",
                summary_stats$Mean[summary_stats$Group == "HC"] * 100,
                summary_stats$Mean[summary_stats$Group == "Amph"] * 100,
                summary_stats$Mean[summary_stats$Group == "Hero"] * 100
            ),
            x = "Subject",
            y = "Prediction Accuracy"
        ) +
        # Theme
        theme_minimal(base_size = 12) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 10),
            strip.text = element_text(face = "bold", size = 12),
            legend.position = "none",
            panel.grid.minor = element_blank()
        )

    # Save plot
    plot_path <- file.path(output_dir, paste0("ppc_combined_", m, ".png"))
    ggsave(plot_path, p, width = 14, height = 5, dpi = 300)
    cat("  Saved:", plot_path, "\n")
}

cat("\nDone!\n")
