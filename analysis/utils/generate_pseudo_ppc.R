#!/usr/bin/env Rscript
# ==============================================================================
# Generate Pseudo PPC Data for Testing Plots
# ==============================================================================
# Creates simulated per-subject prediction accuracy data to test PPC plots
# without needing real model fits.
#
# Usage: Rscript analysis/utils/generate_pseudo_ppc.R
# ==============================================================================

set.seed(42)

# ==============================================================================
# Configuration
# ==============================================================================
output_dir <- "analysis/outputs/ppc"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_subjects <- 48 # Typical IGT sample size
models <- c("pvl_delta", "orl", "eef")

# ==============================================================================
# Generate Pseudo PPC Results
# ==============================================================================

generate_pseudo_ppc <- function(n_subjects, model_quality = "good") {
    # Simulate per-subject prediction accuracy
    # Good model: mean ~0.40-0.50, moderate spread
    # Moderate model: mean ~0.35, higher spread
    # Poor model: mean ~0.30, near chance

    if (model_quality == "good") {
        mean_acc <- 0.45
        sd_acc <- 0.08
    } else if (model_quality == "moderate") {
        mean_acc <- 0.38
        sd_acc <- 0.10
    } else {
        mean_acc <- 0.32
        sd_acc <- 0.06
    }

    # Generate per-subject accuracies (bounded between 0.15 and 0.75)
    subject_accuracy <- rnorm(n_subjects, mean = mean_acc, sd = sd_acc)
    subject_accuracy <- pmax(pmin(subject_accuracy, 0.75), 0.15)

    return(list(
        subject_accuracy = subject_accuracy,
        overall_mean = mean(subject_accuracy),
        overall_sd = sd(subject_accuracy)
    ))
}

cat("=== Generating Pseudo PPC Data ===\n\n")

# Generate for each model with varying quality
model_qualities <- c("pvl_delta" = "good", "orl" = "moderate", "eef" = "good")

ppc_results <- list()

for (model in models) {
    quality <- model_qualities[model]
    ppc_results[[model]] <- generate_pseudo_ppc(n_subjects, quality)

    cat(sprintf("%s (%s quality):\n", toupper(model), quality))
    cat(sprintf("  Mean accuracy: %.3f\n", ppc_results[[model]]$overall_mean))
    cat(sprintf("  SD: %.3f\n\n", ppc_results[[model]]$overall_sd))
}

# Save for plotting
saveRDS(ppc_results, file.path(output_dir, "pseudo_ppc_results.rds"))
cat("Saved to:", file.path(output_dir, "pseudo_ppc_results.rds"), "\n")
