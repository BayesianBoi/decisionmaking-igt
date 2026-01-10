# Combined PPC comparison plot
# Run after all individual PPC runs complete
# Usage: Rscript analysis/1_PPC/plot_ppc_combined.R

library(ggplot2)

cat("=== Creating Combined PPC Plot ===\n\n")

output_dir <- "outputs/ppc"

models <- c("eef", "orl", "pvl_delta")
groups <- c("HC", "Amph", "Hero")
colors <- c(HC = "lightgreen", Amph = "salmon", Hero = "lightblue")

results_list <- list()

for (model in models) {
    for (group in groups) {
        results_path <- file.path(output_dir, paste0("ppc_individual_", model, "_", group, ".rds"))

        if (!file.exists(results_path)) {
            cat("Missing:", results_path, "\n")
            next
        }

        res <- readRDS(results_path)

        results_list[[paste0(model, "_", group)]] <- data.frame(
            accuracy = res$accuracy,
            model = toupper(model),
            group = group,
            mean_acc = res$mean,
            sd_acc = res$sd
        )

        cat(sprintf(
            "%s %s: Mean = %.3f (SD = %.3f)\n",
            toupper(model), group, res$mean, res$sd
        ))
    }
}

if (length(results_list) == 0) {
    stop("No PPC results found. Run individual PPC scripts first.")
}

combined_df <- do.call(rbind, results_list)
combined_df$group <- factor(combined_df$group, levels = c("HC", "Amph", "Hero"))
combined_df$subject <- ave(seq_len(nrow(combined_df)),
    combined_df$model, combined_df$group,
    FUN = seq_along
)

p_facet <- ggplot(combined_df, aes(x = subject, y = accuracy)) +
    geom_ribbon(aes(ymin = mean_acc - sd_acc, ymax = mean_acc + sd_acc, fill = group),
        alpha = 0.3
    ) +
    geom_hline(aes(yintercept = mean_acc), size = 0.8) +
    geom_hline(yintercept = 0.25, linetype = "dashed", color = "black", size = 0.6) +
    geom_point(aes(color = group), size = 2, alpha = 0.6) +
    facet_grid(model ~ group) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    ylim(0, 0.8) +
    labs(
        title = "Posterior Predictive Accuracy: Model Comparison",
        subtitle = "Dashed line = chance level (0.25)",
        x = "Subject",
        y = "Proportion Correctly Predicted"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none",
        strip.text = element_text(face = "bold")
    )

ggsave(file.path(output_dir, "ppc_mpd_comparison.png"), p_facet, width = 12, height = 8, dpi = 300)
ggsave(file.path(output_dir, "ppc_mpd_comparison.pdf"), p_facet, width = 12, height = 8)

p_box <- ggplot(combined_df, aes(x = group, y = accuracy, fill = group)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
    geom_hline(yintercept = 0.25, linetype = "dashed", color = "black", size = 0.8) +
    facet_wrap(~model) +
    scale_fill_manual(values = colors) +
    ylim(0, 0.8) +
    labs(
        title = "Posterior Predictive Accuracy by Group",
        x = "",
        y = "Proportion Correctly Predicted"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none"
    )

ggsave(file.path(output_dir, "ppc_mpd_boxplot.png"), p_box, width = 10, height = 5, dpi = 300)

cat("\nPlots saved to:", output_dir, "\n")

cat("\n=== Summary Statistics ===\n")
for (model in models) {
    cat("\n", toupper(model), ":\n", sep = "")
    for (group in groups) {
        key <- paste0(model, "_", group)
        if (key %in% names(results_list)) {
            res <- results_list[[key]]
            cat(sprintf("  %s: %.3f (SD = %.3f)\n", group, mean(res$accuracy), sd(res$accuracy)))
        }
    }
}
