#!/usr/bin/env Rscript
# group comparison via HDI - if 95% interval excludes zero, groups differ
# usage: Rscript compare_groups_posthoc.R <model>

if (!require("pacman")) install.packages("pacman")
pacman::p_load(HDInterval, ggplot2, dplyr, tidyr)

set.seed(69420)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    cat("Usage: Rscript compare_groups_posthoc.R <model>\n")
    cat("Models: eef, pvl_delta, orl\n")
    cat("\nExample: Rscript compare_groups_posthoc.R eef\n")
    stop("Please provide a model name.")
}

model_name <- tolower(args[1])

valid_models <- c("eef", "pvl_delta", "orl")
if (!model_name %in% valid_models) {
    stop("Invalid model. Must be one of: ", paste(valid_models, collapse = ", "))
}

# Which parameters to test for each model
param_map <- list(
    "eef" = c("mu_theta", "mu_lambda", "mu_phi", "mu_cons"),
    "pvl_delta" = c("mu_A", "mu_a", "mu_w", "mu_theta"),
    "orl" = c("mu_a_rew", "mu_a_pun", "mu_K", "mu_omega_f", "mu_omega_p")
)

# Readable labels for the plots
label_map <- list(
    "mu_theta" = "Outcome Sensitivity",
    "mu_lambda" = "Forgetting Rate",
    "mu_phi" = "Exploration Bonus",
    "mu_cons" = "Choice Consistency",
    "mu_A" = "Outcome Sensitivity",
    "mu_a" = "Learning Rate",
    "mu_w" = "Loss Aversion",
    "mu_a_rew" = "Reward Learning Rate",
    "mu_a_pun" = "Punishment Learning Rate",
    "mu_K" = "Perseverance Decay",
    "mu_omega_f" = "Frequency Weight",
    "mu_omega_p" = "Perseverance Weight"
)

params <- param_map[[model_name]]

cat("\n========================================\n")
cat("GROUP COMPARISON:", toupper(model_name), "\n")
cat("Method: Post-hoc HDI (Ahn et al., 2017)\n")
cat("========================================\n\n")

# Load fitted models
groups <- c("HC", "Amph", "Hero")
fits <- list()

fits_folder <- paste0("data/processed/fits/", model_name)

for (g in groups) {
    fit_path <- file.path(fits_folder, paste0(model_name, "_fit_", g, ".rds"))

    if (!file.exists(fit_path)) {
        stop(
            "Fit file not found: ", fit_path,
            "\nMake sure you have run the fitting scripts first."
        )
    }

    cat("Loading:", fit_path, "\n")
    fits[[g]] <- readRDS(fit_path)
}

# Pull out posterior samples for each parameter
extract_mu <- function(fit, param) {
    samples <- fit$BUGSoutput$sims.list[[param]]
    if (is.matrix(samples)) samples <- samples[, 1]
    return(as.numeric(samples))
}

posteriors <- list()
for (g in groups) {
    posteriors[[g]] <- list()
    for (p in params) {
        posteriors[[g]][[p]] <- extract_mu(fits[[g]], p)
    }
}

cat("\nExtracted", length(posteriors[[groups[1]]][[params[1]]]), "posterior samples per parameter\n\n")

# Compare all pairs of groups
group_pairs <- list(
    c("HC", "Amph"),
    c("HC", "Hero"),
    c("Amph", "Hero")
)

results <- data.frame(
    model = character(),
    parameter = character(),
    group1 = character(),
    group2 = character(),
    mean_g1 = numeric(),
    mean_g2 = numeric(),
    diff_mean = numeric(),
    diff_sd = numeric(),
    hdi_lower = numeric(),
    hdi_upper = numeric(),
    prob_g1_greater = numeric(),
    excludes_zero = logical(),
    stringsAsFactors = FALSE
)

cat("Comparing groups...\n")
cat("========================================\n\n")

for (pair in group_pairs) {
    g1 <- pair[1]
    g2 <- pair[2]

    cat("---", g1, "vs", g2, "---\n\n")

    for (p in params) {
        samples_g1 <- posteriors[[g1]][[p]]
        samples_g2 <- posteriors[[g2]][[p]]

        # Positive difference means g1 > g2
        diff_samples <- samples_g1 - samples_g2

        hdi_95 <- hdi(diff_samples, credMass = 0.95)

        # If the interval does not include zero then the groups differ
        excludes_zero <- (hdi_95[1] > 0) | (hdi_95[2] < 0)

        # What proportion of the posterior has g1 > g2?
        prob_g1_greater <- mean(diff_samples > 0)

        results <- rbind(results, data.frame(
            model = model_name,
            parameter = p,
            group1 = g1,
            group2 = g2,
            mean_g1 = mean(samples_g1),
            mean_g2 = mean(samples_g2),
            diff_mean = mean(diff_samples),
            diff_sd = sd(diff_samples),
            hdi_lower = hdi_95[1],
            hdi_upper = hdi_95[2],
            prob_g1_greater = prob_g1_greater,
            excludes_zero = excludes_zero,
            stringsAsFactors = FALSE
        ))

        param_label <- ifelse(p %in% names(label_map), label_map[[p]], p)
        cat(param_label, " (", p, ")\n", sep = "")
        cat("  ", g1, " M =", round(mean(samples_g1), 3), "\n")
        cat("  ", g2, " M =", round(mean(samples_g2), 3), "\n")
        cat("  Difference (", g1, " - ", g2, ") = ",
            round(mean(diff_samples), 3), " [",
            round(hdi_95[1], 3), ", ", round(hdi_95[2], 3), "]\n",
            sep = ""
        )
        cat("  P(", g1, " > ", g2, ") = ", round(prob_g1_greater, 3), "\n", sep = "")

        if (excludes_zero) {
            direction <- ifelse(mean(diff_samples) > 0,
                paste(g1, ">", g2),
                paste(g2, ">", g1)
            )
            cat("  >> HDI excludes zero:", direction, "\n")
        } else {
            cat("  >> HDI includes zero: no credible difference\n")
        }
        cat("\n")
    }
}

# Save CSV of results
output_dir <- "results/group_comparison"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

csv_path <- file.path(output_dir, paste0("group_comparison_", model_name, ".csv"))
write.csv(results, csv_path, row.names = FALSE)
cat("\nResults saved to:", csv_path, "\n")

# Density plots for each parameter
plot_dir <- "figures/group_comparison"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

group_colours <- c("HC" = "#4DAF4A", "Amph" = "#E41A1C", "Hero" = "#377EB8")

plot_data <- data.frame()
for (g in groups) {
    for (p in params) {
        samples <- posteriors[[g]][[p]]
        plot_data <- rbind(plot_data, data.frame(
            group = g,
            parameter = p,
            value = samples
        ))
    }
}

for (p in params) {
    param_label <- ifelse(p %in% names(label_map), label_map[[p]], p)

    p_data <- plot_data[plot_data$parameter == p, ]

    hdi_data <- data.frame()
    for (g in groups) {
        samples <- posteriors[[g]][[p]]
        hdi_95 <- hdi(samples, credMass = 0.95)
        hdi_data <- rbind(hdi_data, data.frame(
            group = g,
            hdi_lower = hdi_95[1],
            hdi_upper = hdi_95[2],
            mean_val = mean(samples),
            y = -0.05 - 0.03 * which(groups == g)
        ))
    }

    plot <- ggplot(p_data, aes(x = value, fill = group, colour = group)) +
        geom_density(alpha = 0.4, linewidth = 0.8) +
        geom_segment(
            data = hdi_data,
            aes(x = hdi_lower, xend = hdi_upper, y = y, yend = y, colour = group),
            linewidth = 2, inherit.aes = FALSE
        ) +
        geom_point(
            data = hdi_data,
            aes(x = mean_val, y = y, colour = group),
            size = 3, inherit.aes = FALSE
        ) +
        scale_fill_manual(values = group_colours) +
        scale_colour_manual(values = group_colours) +
        labs(
            title = paste(toupper(model_name), "-", param_label),
            subtitle = "Group-level posterior distributions with 95% HDI",
            x = "Parameter Value",
            y = "Density"
        ) +
        theme_minimal(base_size = 12) +
        theme(
            legend.position = "bottom",
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)
        )

    plot_path <- file.path(plot_dir, paste0(model_name, "_", p, "_posteriors.png"))
    ggsave(plot_path, plot, width = 8, height = 6, dpi = 300)
}

cat("Posterior plots saved to:", plot_dir, "\n")

# Summary of credible differences
cat("\n========================================\n")
cat("SUMMARY: Parameters with credible group differences\n")
cat("========================================\n\n")

sig_results <- results[results$excludes_zero, ]

if (nrow(sig_results) > 0) {
    for (i in 1:nrow(sig_results)) {
        row <- sig_results[i, ]
        direction <- ifelse(row$diff_mean > 0,
            paste(row$group1, ">", row$group2),
            paste(row$group2, ">", row$group1)
        )
        param_label <- ifelse(row$parameter %in% names(label_map),
            label_map[[row$parameter]], row$parameter
        )
        cat("-", row$group1, "vs", row$group2, ":", param_label, "->", direction, "\n")
    }
} else {
    cat("No parameters showed credible group differences (all 95% HDIs include zero).\n")
}

cat("\n")
