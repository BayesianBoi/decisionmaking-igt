#!/usr/bin/env Rscript
# Generate all publication-ready outputs
# Run AFTER fit_models.R completes
# Usage: Rscript analysis/generate_paper_outputs.R

cat("=== Generating Publication Outputs ===\n\n")
cat("This script generates all figures, tables, and summaries for manuscript preparation.\n\n")

# Check if fitted models exist
if (!file.exists("analysis/outputs/pvl_delta_fit.rds")) {
  stop("ERROR: No fitted models found. Run fit_models.R first.")
}

# Install/load required packages
required_packages <- c("coda", "ggplot2", "gridExtra")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    library(pkg, character.only = TRUE)
  }
}

# Load utilities
source("analysis/utils/diagnostics.R")
source("analysis/utils/visualization.R")
source("analysis/utils/reporting.R")

# Create output directories
dir.create("analysis/outputs/publication", recursive = TRUE, showWarnings = FALSE)
dir.create("analysis/outputs/publication/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("analysis/outputs/publication/tables", recursive = TRUE, showWarnings = FALSE)

# Step 1: Run diagnostics if not already done
cat("Step 1: Computing convergence diagnostics...\n")
if (!file.exists("analysis/outputs/all_diagnostics.rds")) {
  diagnostics <- run_full_diagnostics()
} else {
  diagnostics <- readRDS("analysis/outputs/all_diagnostics.rds")
  cat("  (Using existing diagnostics)\n")
}

# Step 2: Generate all figures
cat("\nStep 2: Generating publication figures...\n")
generate_all_figures("analysis/outputs", "analysis/outputs/publication/figures")

# Step 3: Generate all tables
cat("\nStep 3: Generating manuscript tables...\n")

# Load fitted models
fit_files <- list.files("analysis/outputs", pattern = "_fit\\.rds$", full.names = TRUE)
fit_results <- list()
for (fit_file in fit_files) {
  model_name <- gsub("_fit\\.rds$", "", basename(fit_file))
  fit_results[[model_name]] <- readRDS(fit_file)
  cat(sprintf("  Loaded: %s\n", model_name))
}

tables <- generate_manuscript_tables(fit_results, diagnostics,
                                     "analysis/outputs/publication/tables")

# Step 4: Generate results summary
cat("\nStep 4: Generating results summary...\n")
write_results_summary(fit_results, "analysis/outputs/publication/RESULTS_SUMMARY.md")

# Step 5: Create formatted output for copying to manuscript
cat("\nStep 5: Creating formatted text outputs...\n")

# Parameter estimates in APA format
apa_file <- "analysis/outputs/publication/APA_FORMATTED_RESULTS.txt"
sink(apa_file)

cat("=== FORMATTED RESULTS FOR MANUSCRIPT ===\n\n")

for (model_name in names(fit_results)) {
  samples <- as.matrix(fit_results[[model_name]]$samples)
  mu_params <- grep("^mu_", colnames(samples), value = TRUE)

  cat(sprintf("\n%s Model\n", toupper(model_name)))
  cat(paste(rep("-", 50), collapse = ""), "\n\n")

  for (param in mu_params) {
    values <- samples[, param]
    mean_val <- mean(values)
    sd_val <- sd(values)
    ci_lower <- quantile(values, 0.025)
    ci_upper <- quantile(values, 0.975)

    # APA format: M = X.XX, SD = X.XX, 95% CI [X.XX, X.XX]
    cat(sprintf("%s: M = %.3f, SD = %.3f, 95%% CI [%.3f, %.3f]\n",
                param, mean_val, sd_val, ci_lower, ci_upper))
  }
}

cat("\n\n=== CONVERGENCE SUMMARY ===\n\n")

for (model_name in names(diagnostics)) {
  conv <- diagnostics[[model_name]]$convergence

  cat(sprintf("\n%s:\n", toupper(model_name)))
  cat(sprintf("  Convergence: %s\n", ifelse(conv$converged, "YES", "NO")))
  cat(sprintf("  Max R-hat: %.3f (threshold: %.2f)\n",
              max(diagnostics[[model_name]]$diagnostics$rhat), conv$threshold))
  cat(sprintf("  Problematic parameters: %d/%d\n",
              conv$n_problematic, conv$n_total))
}

sink()

cat(sprintf("  Formatted results saved to: %s\n", apa_file))

# Step 6: Generate model comparison summary
cat("\nStep 6: Creating model comparison summary...\n")

comparison_file <- "analysis/outputs/publication/MODEL_COMPARISON.md"
writeLines(c(
  "# Model Comparison Summary",
  "",
  "## Parameter Complexity",
  ""
), comparison_file)

# Add parameter counts
for (model_name in names(fit_results)) {
  samples <- as.matrix(fit_results[[model_name]]$samples)
  n_group <- length(grep("^(mu_|sigma_)", colnames(samples)))
  n_total <- ncol(samples)

  cat(sprintf("- **%s**: %d group-level parameters, %d total parameters\n",
              toupper(model_name), n_group, n_total),
      file = comparison_file, append = TRUE)
}

writeLines(c(
  "",
  "## Key Findings",
  "",
  "*(To be filled based on specific research questions)*",
  "",
  "## Model Performance",
  "",
  "See `tables/table3_convergence.csv` for detailed convergence metrics.",
  ""
), comparison_file, append = TRUE)

cat(sprintf("  Model comparison saved to: %s\n", comparison_file))

# Summary of outputs
cat("\n=== Publication Package Complete ===\n\n")
cat("Generated outputs:\n\n")
cat("FIGURES (analysis/outputs/publication/figures/):\n")
cat("  - *_traces.pdf         : MCMC trace plots for convergence checking\n")
cat("  - *_posteriors.pdf     : Posterior density plots\n")
cat("  - model_comparison_forest.pdf : Parameter estimates across models\n\n")

cat("TABLES (analysis/outputs/publication/tables/):\n")
cat("  - table1_parameter_estimates.csv : All parameter estimates with CIs\n")
cat("  - table2_model_comparison.csv    : Model complexity comparison\n")
cat("  - table3_convergence.csv         : Convergence diagnostics\n\n")

cat("SUMMARIES (analysis/outputs/publication/):\n")
cat("  - RESULTS_SUMMARY.md            : Quick overview of results\n")
cat("  - APA_FORMATTED_RESULTS.txt     : Copy-paste ready APA format\n")
cat("  - MODEL_COMPARISON.md           : Model comparison notes\n\n")

cat("Next steps:\n")
cat("  1. Review trace plots for convergence\n")
cat("  2. Check table3_convergence.csv (all R-hat < 1.1?)\n")
cat("  3. Review posterior densities\n")
cat("  4. Copy APA_FORMATTED_RESULTS.txt to manuscript\n")
cat("  5. Insert figures and tables into manuscript\n\n")

cat("Done! ✅\n")
