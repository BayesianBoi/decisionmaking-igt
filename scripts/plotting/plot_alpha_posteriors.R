#!/usr/bin/env Rscript
# plot_alpha_posteriors.R
# -----------------------
# Plots the posterior distributions of alpha parameters from a group comparison.
#
# Alpha represents the difference between two groups on each cognitive parameter.
# If the distribution is clearly away from zero, we have evidence of a group
# difference. These plots are useful for papers because they show the full
# posterior, not just summary statistics.
#
# The script produces:
#   - Density plots for each alpha parameter
#   - A vertical line at zero for reference
#   - Shaded 95% HDI region
#
# Usage: Rscript plot_alpha_posteriors.R <model> <group1> <group2>
# Example: Rscript plot_alpha_posteriors.R orl HC Hero

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Usage: Rscript plot_alpha_posteriors.R <model> <group1> <group2>")
}

model_name <- args[1]
g1 <- args[2]
g2 <- args[3]

# Load required packages
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("HDInterval", quietly = TRUE)) install.packages("HDInterval")

library(ggplot2)
library(HDInterval)

# Build path to the fitted comparison model
fit_file <- file.path(
    "outputs/group_comparison",
    paste0("compare_", model_name, "_", g1, "_vs_", g2, ".rds")
)

if (!file.exists(fit_file)) {
    stop("Fit file not found: ", fit_file)
}

cat("\n========================================\n")
cat("Plotting Alpha Posteriors:", toupper(model_name), "\n")
cat(g1, "vs", g2, "\n")
cat("========================================\n\n")

# Load the fitted model
fit <- readRDS(fit_file)

# Find all alpha parameters
all_params <- names(fit$BUGSoutput$sims.list)
alpha_params <- all_params[grepl("^alpha_", all_params)]

if (length(alpha_params) == 0) {
    stop("No alpha parameters found in the fitted model.")
}

# Create nicer labels for the parameters
# These will appear on the plot
label_map <- list(
    # ORL parameters
    "alpha_a_rew" = "Reward Learning Rate",
    "alpha_a_pun" = "Punishment Learning Rate",
    "alpha_K" = "Perseverance Decay",
    "alpha_theta" = "Choice Consistency",
    "alpha_omega_f" = "Frequency Weight",
    "alpha_omega_p" = "Perseverance Weight",
    # EEF parameters
    "alpha_lambda" = "Forgetting Rate",
    "alpha_phi" = "Exploration Bonus",
    "alpha_cons" = "Choice Consistency"
)

# Collect all samples into a data frame for plotting
plot_data <- data.frame()

for (param in alpha_params) {
    samples <- as.numeric(fit$BUGSoutput$sims.list[[param]])

    # Use nice label if available, otherwise use the parameter name
    label <- if (param %in% names(label_map)) label_map[[param]] else param

    df <- data.frame(
        value = samples,
        param = label,
        stringsAsFactors = FALSE
    )

    plot_data <- rbind(plot_data, df)
}

# Compute HDI for each parameter (for shading)
hdi_data <- do.call(rbind, lapply(alpha_params, function(param) {
    samples <- as.numeric(fit$BUGSoutput$sims.list[[param]])
    hdi_95 <- hdi(samples, credMass = 0.95)
    label <- if (param %in% names(label_map)) label_map[[param]] else param

    data.frame(
        param = label,
        hdi_lower = hdi_95[1],
        hdi_upper = hdi_95[2],
        mean = mean(samples),
        stringsAsFactors = FALSE
    )
}))

# Create the plot
# Each panel shows one parameter's posterior
p <- ggplot(plot_data, aes(x = value)) +
    geom_density(fill = "#3182bd", alpha = 0.6, colour = "#08519c") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.8) +
    facet_wrap(~param, scales = "free", ncol = 2) +
    labs(
        title = paste0(toupper(model_name), ": Group Difference Posteriors (", g1, " vs ", g2, ")"),
        subtitle = paste0("Dashed line = no difference. If distribution is away from zero, groups differ."),
        x = paste0("Alpha (positive = ", g2, " > ", g1, ")"),
        y = "Density"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.subtitle = element_text(colour = "grey40", size = 9)
    )

# Create output directory
output_dir <- "plots/group_comparison"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save the plot
plot_file <- file.path(
    output_dir,
    paste0("alpha_posteriors_", model_name, "_", g1, "_vs_", g2, ".png")
)

ggsave(plot_file, p, width = 10, height = 8, dpi = 300)
cat("Plot saved:", plot_file, "\n")

# Also save as PDF for publication
pdf_file <- gsub(".png$", ".pdf", plot_file)
ggsave(pdf_file, p, width = 10, height = 8)
cat("PDF saved:", pdf_file, "\n")

# Print a quick summary
cat("\n")
cat("Summary of group differences:\n")
for (i in seq_len(nrow(hdi_data))) {
    row <- hdi_data[i, ]
    excludes_zero <- (row$hdi_lower > 0) | (row$hdi_upper < 0)
    status <- if (excludes_zero) "* DIFFERS" else "  no diff"

    cat(sprintf(
        "  %s %s: mean = %.3f, 95%% HDI [%.3f, %.3f]\n",
        status, row$param, row$mean, row$hdi_lower, row$hdi_upper
    ))
}
