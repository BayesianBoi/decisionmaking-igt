#!/usr/bin/env Rscript
#
# plot_ppc_barplot.R
# Combined bar plot showing prediction accuracy for all models and groups.
# The dashed line at 0.25 marks chance level (random guessing among 4 decks).
#

library(ggplot2)
library(dplyr)

cat("=== Generating PPC Bar Plot ===\n\n")

models <- c("pvl_delta", "orl", "eef")
model_labels <- c("pvl_delta" = "PVL-Delta", "orl" = "ORL", "eef" = "EEF")
groups <- c("HC", "Amph", "Hero")
group_labels <- c("HC" = "Healthy Controls", "Amph" = "Amphetamine", "Hero" = "Heroin")

colors <- c(
    "Healthy Controls" = "#4CAF50",
    "Amphetamine" = "#FF5722",
    "Heroin" = "#2196F3"
)

output_dir <- "figures/ppc"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

all_data <- data.frame()

for (m in models) {
    for (g in groups) {
        file <- paste0("data/processed/ppc/", m, "/ppc_individual_", m, "_", g, ".rds")
        if (file.exists(file)) {
            ppc <- readRDS(file)

            df <- ppc$df
            if (is.null(df)) {
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

all_data$Model <- factor(all_data$Model, levels = c("PVL-Delta", "ORL", "EEF"))
all_data$Group <- factor(all_data$Group, levels = c("Healthy Controls", "Amphetamine", "Heroin"))

cat("Summary Statistics:\n")
print(all_data)

p <- ggplot(all_data, aes(x = Model, y = Mean, fill = Group)) +
    geom_hline(yintercept = 0.25, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_bar(
        stat = "identity",
        position = position_dodge(width = 0.8),
        width = 0.7,
        alpha = 0.9
    ) +
    geom_errorbar(
        aes(ymin = Mean - SD, ymax = Mean + SD),
        position = position_dodge(width = 0.8),
        width = 0.25,
        color = "black",
        linewidth = 0.6
    ) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.1), expand = c(0, 0)) +
    labs(
        title = "Posterior Predictive Accuracy",
        subtitle = "Mean prediction accuracy by model and group (Error bars show one SD)",
        y = "Prediction Accuracy",
        x = "Model"
    ) +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
    )

plot_path <- file.path("figures/paper", "ppc_barplot_combined.png")
ggsave(plot_path, p, width = 10, height = 7, dpi = 600)
cat("\nSaved:", plot_path, "\n")

cat("\nDone!\n")
