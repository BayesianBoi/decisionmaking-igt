#!/usr/bin/env Rscript
# ==============================================================================
# PPC BAR PLOT - ALL MODELS COMBINED
# ==============================================================================
# Creates a grouped bar plot comparing all models and groups with SD error bars
# Usage: Rscript scripts/plotting/plot_ppc_barplot.R
# ==============================================================================

library(ggplot2)
library(dplyr)

cat("=== Generating PPC Bar Plot ===\n\n")

# Settings
models <- c("pvl_delta", "orl", "eef")
model_labels <- c("pvl_delta" = "PVL-Delta", "orl" = "ORL", "eef" = "EEF")
groups <- c("HC", "Amph", "Hero")
group_labels <- c("HC" = "Healthy Controls", "Amph" = "Amphetamine", "Hero" = "Heroin")

# Color palette
colors <- c(
    "Healthy Controls" = "#4CAF50",
    "Amphetamine" = "#FF5722",
    "Heroin" = "#2196F3"
)

output_dir <- "plots/ppc"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Collect summary statistics directly
# We need to load individual files to calculate the mean and SD across subjects,
# just like in plot_ppc_summary.R

all_data <- data.frame()

for (m in models) {
    for (g in groups) {
        file <- paste0("outputs/ppc/ppc_individual_", m, "_", g, ".rds")
        if (file.exists(file)) {
            ppc <- readRDS(file)

            # The .rds file contains $mean (mean accuracy across all subjects?)
            # or is it $df which has individual subject accuracies?
            # Let's check the structure relative to plot_ppc_summary.R
            # plot_ppc_summary.R uses:
            # df <- ppc$df (which has $accuracy per subject)
            # and calculates summary from that.

            df <- ppc$df
            if (is.null(df)) {
                # Fallback if structure is different
                cat("Warning: Unexpected structure for", m, g, "\n")
                next
            }

            mean_acc <- mean(df$accuracy, na.rm = TRUE)
            sd_acc <- sd(df$accuracy, na.rm = TRUE)
            n_sub <- nrow(df)

            all_data <- rbind(all_data, data.frame(
                Model = model_labels[m],
                Group = group_labels[g],
                Mean = mean_acc,
                SD = sd_acc,
                N = n_sub
            ))
        }
    }
}

if (nrow(all_data) == 0) {
    stop("No PPC data found!")
}

# Set factor levels for ordering
all_data$Model <- factor(all_data$Model, levels = c("PVL-Delta", "ORL", "EEF"))
all_data$Group <- factor(all_data$Group, levels = c("Healthy Controls", "Amphetamine", "Heroin"))

cat("Summary Statistics:\n")
print(all_data)

# Create Grouped Bar Plot
p <- ggplot(all_data, aes(x = Model, y = Mean, fill = Group)) +
    # Chance level line
    geom_hline(yintercept = 0.25, linetype = "dashed", color = "gray40", linewidth = 0.8) +

    # Bars
    geom_bar(
        stat = "identity",
        position = position_dodge(width = 0.8),
        width = 0.7,
        alpha = 0.9
    ) +

    # Error Bars (Mean +/- SD)
    geom_errorbar(
        aes(ymin = Mean - SD, ymax = Mean + SD),
        position = position_dodge(width = 0.8),
        width = 0.25,
        color = "black",
        linewidth = 0.6
    ) +

    # Aesthetics
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.1), expand = c(0, 0)) +

    # Labels
    labs(
        title = "Posterior Predictive Accuracy",
        subtitle = "Mean prediction accuracy by model and group (Error bars = ±1 SD)",
        y = "Prediction Accuracy",
        x = "Model"
    ) +

    # Theme
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(), # Model names are clear enough
        legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
    )

# Save plot
plot_path <- file.path(output_dir, "ppc_barplot_combined.png")
ggsave(plot_path, p, width = 10, height = 7, dpi = 300)
cat("\nSaved:", plot_path, "\n")

cat("\nDone!\n")
