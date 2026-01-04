# ===========================================================================
# Run Posterior Predictive Checks (PPC)
# ===========================================================================
#
# Usage: Rscript analysis/scripts/run_ppc.R
#
# This script performs posterior predictive checks for all fitted models
# (PVL-Delta, ORL, EEF) to validate model descriptive adequacy.
# It generates summary tables and plots comparing observed vs predicted
# choice behavior.
# ===========================================================================

library(rjags)
library(coda)
library(ggplot2)
library(gridExtra)

source("analysis/utils/load_data.R")
source("analysis/utils/ppc.R")
source("analysis/utils/prepare_jags_data.R")
source("analysis/utils/prepare_eef_data.R")

# ===========================================================================
# Configuration
# ===========================================================================

output_dir <- "results/ppc"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

models <- c("pvl_delta", "orl", "eef")
n_sim <- 500 # Number of posterior samples to use for simulation

# ===========================================================================
# Helper: Extract parameters and data for each model type
# ===========================================================================

get_model_objects <- function(model_name) {
    results_dir <- file.path("results", model_name)
    samples_path <- file.path(results_dir, "mcmc_samples.rds")
    data_path <- file.path(results_dir, "jags_data.rds")

    if (!file.exists(samples_path) || !file.exists(data_path)) {
        warning(sprintf("Fit results not found for %s", model_name))
        return(NULL)
    }

    samples <- readRDS(samples_path)
    jags_data <- readRDS(data_path)

    # Standardize list format for the PPC function
    return(list(
        samples = samples,
        data = jags_data
    ))
}

# ===========================================================================
# Main PPC Loop
# ===========================================================================

ppc_summary_list <- list()

for (model in models) {
    message(sprintf("\nRunning PPC for model: %s...", model))

    objs <- get_model_objects(model)
    if (is.null(objs)) next

    # Convert samples to matrix only once
    samples_mat <- as.matrix(objs$samples)

    # The run_ppc function in utils/ppc.R expects a fit object with $samples
    # We construct a wrapper to match that expectation
    fit_wrapper <- list(samples = objs$samples)

    # Run the PPC
    # Note: The utility function run_ppc takes (fit_result, observed_data, model_name)
    # We need to make sure observed_data has what simulate_* functions need

    # For EEF, the data structure is slightly different, but the necessary content
    # (outcome, choice) works if names match.

    ppc_out <- run_ppc(fit_wrapper, objs$data, model, n_sim = n_sim)

    if (!is.null(ppc_out)) {
        # Save raw results
        saveRDS(ppc_out, file.path(output_dir, sprintf("ppc_results_%s.rds", model)))

        # Store summary for comparison
        summ <- ppc_out$summary
        summ$model <- model
        ppc_summary_list[[model]] <- summ

        # --- Plotting ---

        # Create plot for this model
        p <- ggplot(summ, aes(x = factor(deck))) +
            # Observed bars
            geom_bar(aes(y = observed, fill = "Observed"), stat = "identity", alpha = 0.6, width = 0.6) +
            # Predicted error bars
            geom_errorbar(aes(ymin = predicted_lower, ymax = predicted_upper, color = "Predicted (95% CI)"),
                width = 0.3, size = 1
            ) +
            # Predicted points
            geom_point(aes(y = predicted_mean, color = "Predicted (95% CI)"), size = 3) +
            scale_fill_manual(values = c("Observed" = "gray30")) +
            scale_color_manual(values = c("Predicted (95% CI)" = "firebrick")) +
            labs(
                title = sprintf("Posterior Predictive Check: %s", toupper(model)),
                subtitle = "Observed vs Predicted Choice Proportions",
                x = "Deck", y = "Choice Proportion",
                fill = "", color = ""
            ) +
            theme_minimal() +
            theme(legend.position = "bottom")

        ggsave(file.path(output_dir, sprintf("ppc_plot_%s.pdf", model)), p, width = 6, height = 5)
    }
}

# ===========================================================================
# Comparison Plot
# ===========================================================================

if (length(ppc_summary_list) > 0) {
    message("\nGenerating comparison plot...")

    all_ppc <- do.call(rbind, ppc_summary_list)

    p_comp <- ggplot(all_ppc, aes(x = factor(deck))) +
        geom_bar(aes(y = observed, fill = "Observed"), stat = "identity", alpha = 0.4, position = "identity") +
        geom_point(aes(y = predicted_mean, color = model), position = position_dodge(width = 0.5), size = 2) +
        geom_errorbar(aes(ymin = predicted_lower, ymax = predicted_upper, color = model),
            position = position_dodge(width = 0.5), width = 0.2
        ) +
        facet_wrap(~model) +
        labs(
            title = "Model Fit Comparison (PPC)",
            subtitle = "Ability to recover overall choice rates",
            x = "Deck", y = "Proportion"
        ) +
        theme_bw()

    ggsave(file.path(output_dir, "ppc_comparison.pdf"), p_comp, width = 10, height = 4)
    message(sprintf("PPC complete. Results saved to %s", output_dir))
} else {
    message("No PPC results generated. Check if model fits exist.")
}
