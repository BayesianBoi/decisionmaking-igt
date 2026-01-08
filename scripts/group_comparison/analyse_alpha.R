#!/usr/bin/env Rscript
# ANALYSE ALPHA
#
# Analyses the difference parameters (alpha) from a group comparison model.
# The alpha parameters represent the difference between two groups on each
# cognitive parameter. If the 95% HDI excludes zero then we can say there is
# a group difference worth reporting.
#
# The joint model sets up group means like this
#   Group1 mean = mu - alpha/2
#   Group2 mean = mu + alpha/2
# So alpha > 0 means Group2 has higher values than Group1.
#
# Run from project root with model name and both groups as arguments.
# Example call from terminal
#   Rscript scripts/group_comparison/analyse_alpha.R orl HC Hero

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Usage requires model and both groups. For example orl HC Hero")
}

model_name <- args[1]
g1 <- args[2]
g2 <- args[3]

# We need the HDInterval package for computing highest density intervals
if (!require("HDInterval", quietly = TRUE)) {
    install.packages("HDInterval")
}
library(HDInterval)

# Build the path to the fitted model
fit_file <- file.path(
    "outputs/group_comparison",
    paste0("compare_", model_name, "_", g1, "_vs_", g2, ".rds")
)

if (!file.exists(fit_file)) {
    stop(paste("Fit file not found at", fit_file))
}

cat("\n========================================\n")
cat("Group Comparison Analysis", toupper(model_name), "\n")
cat(g1, "vs", g2, "\n")
cat("========================================\n\n")

# Load the fitted model
fit <- readRDS(fit_file)

# Find all alpha parameters since these encode group differences
all_params <- names(fit$BUGSoutput$sims.list)
alpha_params <- all_params[grepl("^alpha_", all_params)]

if (length(alpha_params) == 0) {
    stop("No alpha parameters found in fit.")
}

cat("Parameters analysed\n")

# Set up results dataframe
results <- data.frame(
    parameter = character(),
    mean = numeric(),
    sd = numeric(),
    hdi_lower = numeric(),
    hdi_upper = numeric(),
    prob_positive = numeric(),
    prob_negative = numeric(),
    significant = character(),
    stringsAsFactors = FALSE
)

# Loop through each alpha parameter
for (param in alpha_params) {
    # Extract posterior samples for this parameter
    samples <- as.numeric(fit$BUGSoutput$sims.list[[param]])

    # The 95% Highest Density Interval is the narrowest interval containing 95% of samples
    hdi_95 <- hdi(samples, credMass = 0.95)

    # Probability of direction tells us how confident we are about the sign
    prob_pos <- mean(samples > 0)
    prob_neg <- mean(samples < 0)

    # If the HDI does not include zero we have a difference worth reporting
    excludes_zero <- (hdi_95[1] > 0) | (hdi_95[2] < 0)
    sig_label <- ifelse(excludes_zero, "Yes", "No")

    # Store results
    results <- rbind(results, data.frame(
        parameter = param,
        mean = mean(samples),
        sd = sd(samples),
        hdi_lower = hdi_95[1],
        hdi_upper = hdi_95[2],
        prob_positive = prob_pos,
        prob_negative = prob_neg,
        significant = sig_label,
        stringsAsFactors = FALSE
    ))

    # Print to console
    cat("\n", param, "\n", sep = "")
    cat("  Mean (SD)", round(mean(samples), 3), "(", round(sd(samples), 3), ")\n")
    cat("  95% HDI [", round(hdi_95[1], 3), ", ", round(hdi_95[2], 3), "]\n", sep = "")
    cat("  P(alpha > 0)", round(prob_pos, 3), "\n")
    cat("  P(alpha < 0)", round(prob_neg, 3), "\n")

    if (excludes_zero) {
        # Work out which group is higher
        direction <- ifelse(mean(samples) > 0, paste0(g2, " > ", g1), paste0(g1, " > ", g2))
        cat("  >> HDI excludes zero --", direction, "\n")
    } else {
        cat("  >> HDI includes zero -- no group difference\n")
    }
}

# Save the results to CSV
output_file <- file.path(
    "outputs/group_comparison",
    paste0("alpha_analysis_", model_name, "_", g1, "_vs_", g2, ".csv")
)
write.csv(results, output_file, row.names = FALSE)
cat("\n\nResults saved to", output_file, "\n")

# Print summary
cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")

sig_params <- results$parameter[results$significant == "Yes"]
if (length(sig_params) > 0) {
    cat("Parameters with group differences\n")
    for (p in sig_params) {
        row <- results[results$parameter == p, ]
        direction <- ifelse(row$mean > 0, paste0(g2, " > ", g1), paste0(g1, " > ", g2))
        cat("  -", p, direction, "\n")
    }
} else {
    cat("No parameters showed group differences (95% HDI includes zero for all).\n")
}

# Remind user how to interpret alpha
cat("\nHow to interpret alpha\n")
cat("  alpha > 0 means", g2, "has higher values than", g1, "\n")
cat("  alpha < 0 means", g1, "has higher values than", g2, "\n")
