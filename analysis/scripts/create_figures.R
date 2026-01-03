# Generate publication-ready figures for the paper
# Creates all figures from model fits and comparisons
#
# Run: Rscript analysis/scripts/create_figures.R

library(ggplot2)
library(coda)

#==============================================================================
# CONFIGURATION
#==============================================================================

# Input directories
model_comparison_dir <- "results/model_comparison"
group_comparison_dir <- "results/group_comparison"
parameter_recovery_dir <- "results/parameter_recovery"

# Output directory
figures_dir <- "results/figures"
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Figure settings
fig_width <- 7
fig_height <- 5
fig_dpi <- 300

#==============================================================================
# FIGURE 1: MODEL COMPARISON
#==============================================================================

message("=== CREATING FIGURE 1: MODEL COMPARISON ===\n")

comparison_file <- file.path(model_comparison_dir, "model_comparison_table.csv")
if (file.exists(comparison_file)) {

  comparison_df <- read.csv(comparison_file)

  pdf(file.path(figures_dir, "figure1_model_comparison.pdf"),
      width = fig_width, height = fig_height)

  par(mfrow = c(1, 2), mar = c(5, 5, 3, 2))

  # Panel A: Convergence (R-hat)
  barplot(comparison_df$mean_rhat,
          names.arg = toupper(comparison_df$model),
          main = "A. Model Convergence",
          ylab = "Mean R-hat",
          ylim = c(0.9, 1.2),
          col = c("steelblue", "orange", "darkgreen"),
          border = "white")
  abline(h = 1.1, col = "red", lty = 2, lwd = 2)
  text(x = 0.7, y = 1.15, labels = "Threshold", col = "red", cex = 0.8)

  # Panel B: Effective sample size
  barplot(comparison_df$median_ess,
          names.arg = toupper(comparison_df$model),
          main = "B. Effective Sample Size",
          ylab = "Median ESS",
          col = c("steelblue", "orange", "darkgreen"),
          border = "white")
  abline(h = 1000, col = "red", lty = 2, lwd = 2)

  dev.off()

  message("Figure 1 saved.\n")

} else {
  warning("Model comparison file not found. Skipping Figure 1.")
}

#==============================================================================
# FIGURE 2: GROUP DIFFERENCES IN FORGETTING
#==============================================================================

message("=== CREATING FIGURE 2: GROUP DIFFERENCES ===\n")

group_comp_file <- file.path(group_comparison_dir, "group_comparisons.rds")
if (file.exists(group_comp_file)) {

  comparisons <- readRDS(group_comp_file)

  # Load EEF samples to get full distributions
  eef_samples_file <- "results/eef_clinical/mcmc_samples.rds"
  eef_jags_file <- "results/eef_clinical/jags_data.rds"

  if (file.exists(eef_samples_file) && file.exists(eef_jags_file)) {

    samples <- readRDS(eef_samples_file)
    jags_data <- readRDS(eef_jags_file)

    samples_mat <- as.matrix(samples)
    lambda_cols <- grep("^lambda_forget\\[", colnames(samples_mat), value = TRUE)
    lambda_samples <- samples_mat[, lambda_cols]

    group_assignments <- jags_data$group
    group_names <- c("HC", "Amphetamine", "Heroin", "Cannabis")
    n_groups <- max(group_assignments)

    pdf(file.path(figures_dir, "figure2_group_differences.pdf"),
        width = 10, height = 6)

    par(mfrow = c(1, 2), mar = c(5, 5, 3, 2))

    # Panel A: Violin plot
    group_colors <- c("steelblue", "orange", "darkred", "darkgreen")

    plot(NULL, xlim = c(0.5, n_groups + 0.5), ylim = c(0, 1),
         xlab = "Group", ylab = "Forgetting Rate (lambda)",
         main = "A. Posterior Distributions by Group",
         xaxt = "n", cex.lab = 1.2)
    axis(1, at = 1:n_groups, labels = group_names[1:n_groups])

    for (g in 1:n_groups) {
      group_subjects <- which(group_assignments == g)
      group_mean_per_iter <- rowMeans(lambda_samples[, group_subjects, drop = FALSE])

      overall_mean <- mean(group_mean_per_iter)
      ci_lower <- quantile(group_mean_per_iter, 0.025)
      ci_upper <- quantile(group_mean_per_iter, 0.975)

      # Violin (density)
      dens <- density(group_mean_per_iter)
      dens_scaled <- dens$y / max(dens$y) * 0.35
      polygon(c(g - dens_scaled, rev(g + dens_scaled)),
              c(dens$x, rev(dens$x)),
              col = adjustcolor(group_colors[g], alpha = 0.4),
              border = group_colors[g], lwd = 2)

      # Mean and CI
      points(g, overall_mean, pch = 19, cex = 2, col = group_colors[g])
      arrows(g, ci_lower, g, ci_upper, angle = 90, code = 3,
             length = 0.1, lwd = 3, col = group_colors[g])
    }

    legend("topright", legend = group_names[1:n_groups],
           col = group_colors[1:n_groups], pch = 19, cex = 1.1, bty = "n")

    # Panel B: Effect sizes (differences from HC)
    comp_names <- c("amph_vs_hc", "heroin_vs_hc", "cannabis_vs_hc")
    comp_labels <- c("Amphetamine", "Heroin", "Cannabis")
    n_comp <- length(comp_names)

    mean_diffs <- sapply(comp_names, function(x) {
      if (x %in% names(comparisons)) comparisons[[x]]$mean_diff else NA
    })
    ci_lower <- sapply(comp_names, function(x) {
      if (x %in% names(comparisons)) comparisons[[x]]$ci_lower else NA
    })
    ci_upper <- sapply(comp_names, function(x) {
      if (x %in% names(comparisons)) comparisons[[x]]$ci_upper else NA
    })
    prob_greater <- sapply(comp_names, function(x) {
      if (x %in% names(comparisons)) comparisons[[x]]$prob_greater else NA
    })

    # Remove NAs
    valid_idx <- !is.na(mean_diffs)
    mean_diffs <- mean_diffs[valid_idx]
    ci_lower <- ci_lower[valid_idx]
    ci_upper <- ci_upper[valid_idx]
    prob_greater <- prob_greater[valid_idx]
    comp_labels <- comp_labels[valid_idx]

    if (length(mean_diffs) > 0) {
      bp <- barplot(mean_diffs,
                    names.arg = comp_labels,
                    main = "B. Group Differences vs HC",
                    ylab = "Difference in Forgetting Rate",
                    ylim = range(c(0, ci_lower, ci_upper)) * 1.2,
                    col = c("orange", "darkred", "darkgreen")[valid_idx],
                    border = "white",
                    cex.lab = 1.2)

      arrows(x0 = bp, y0 = ci_lower, x1 = bp, y1 = ci_upper,
             angle = 90, code = 3, length = 0.1, lwd = 2)

      abline(h = 0, lty = 2, col = "gray50", lwd = 2)

      # Add probability labels
      for (i in seq_along(bp)) {
        text(bp[i], ci_upper[i] + 0.02,
             sprintf("P=%.2f", prob_greater[i]),
             cex = 0.9)
      }
    }

    dev.off()

    message("Figure 2 saved.\n")

  } else {
    warning("EEF samples not found. Skipping Figure 2.")
  }

} else {
  warning("Group comparison file not found. Skipping Figure 2.")
}

#==============================================================================
# FIGURE 3: PARAMETER RECOVERY
#==============================================================================

message("=== CREATING FIGURE 3: PARAMETER RECOVERY ===\n")

recovery_file <- file.path(parameter_recovery_dir, "recovery_results.rds")
if (file.exists(recovery_file)) {

  recovery_results <- readRDS(recovery_file)

  pdf(file.path(figures_dir, "figure3_parameter_recovery.pdf"),
      width = 10, height = 10)

  par(mfrow = c(2, 2), mar = c(5, 5, 3, 2))

  param_names <- c("theta", "lambda_forget", "phi", "cons")
  param_labels <- c("Outcome Sensitivity (theta)",
                    "Forgetting Rate (lambda)",
                    "Exploration Incentive (phi)",
                    "Choice Consistency (cons)")

  for (i in seq_along(param_names)) {
    param <- param_names[i]
    res <- recovery_results[[param]]

    plot(res$true, res$recovered,
         xlab = sprintf("True %s", param_labels[i]),
         ylab = sprintf("Recovered %s", param_labels[i]),
         main = sprintf("%s", param_labels[i]),
         pch = 19, col = adjustcolor("steelblue", alpha = 0.5),
         cex = 1.2,
         xlim = range(c(res$true, res$recovered)),
         ylim = range(c(res$true, res$recovered)),
         cex.lab = 1.2)

    abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)

    fit <- lm(res$recovered ~ res$true)
    abline(fit, col = "blue", lwd = 2)

    legend("topleft",
           legend = c(sprintf("r = %.3f", res$correlation),
                      sprintf("RMSE = %.3f", res$rmse)),
           bty = "n", cex = 1.1)
  }

  dev.off()

  message("Figure 3 saved.\n")

} else {
  warning("Parameter recovery file not found. Skipping Figure 3.")
}

#==============================================================================
# FIGURE 4: POSTERIOR PREDICTIVE CHECK
#==============================================================================

message("=== CREATING FIGURE 4: POSTERIOR PREDICTIVE CHECK ===\n")

# Load EEF data and samples
eef_samples_file <- "results/eef_clinical/mcmc_samples.rds"
eef_jags_file <- "results/eef_clinical/jags_data.rds"

if (file.exists(eef_samples_file) && file.exists(eef_jags_file)) {

  samples <- readRDS(eef_samples_file)
  jags_data <- readRDS(eef_jags_file)

  # Compute observed choice frequencies per deck
  observed_choices <- jags_data$choice
  n_subjects <- nrow(observed_choices)
  n_trials <- ncol(observed_choices)

  # Count choices for each deck (ignoring NAs)
  deck_counts <- matrix(0, nrow = n_subjects, ncol = 4)
  for (s in 1:n_subjects) {
    valid_trials <- !is.na(observed_choices[s, ])
    choices_s <- observed_choices[s, valid_trials]
    for (d in 1:4) {
      deck_counts[s, d] <- sum(choices_s == d)
    }
  }

  # Compute proportions
  deck_props <- deck_counts / rowSums(deck_counts)

  # Aggregate across subjects
  mean_props <- colMeans(deck_props)

  pdf(file.path(figures_dir, "figure4_posterior_predictive.pdf"),
      width = fig_width, height = fig_height)

  par(mar = c(5, 5, 3, 2))

  barplot(mean_props,
          names.arg = c("Deck A", "Deck B", "Deck C", "Deck D"),
          main = "Choice Proportions by Deck",
          ylab = "Proportion of Choices",
          ylim = c(0, 0.4),
          col = c("darkred", "red", "lightblue", "steelblue"),
          border = "white",
          cex.lab = 1.2)

  # Add error bars (SD across subjects)
  sd_props <- apply(deck_props, 2, sd)
  se_props <- sd_props / sqrt(n_subjects)

  bp <- barplot(mean_props,
                names.arg = c("Deck A", "Deck B", "Deck C", "Deck D"),
                main = "Observed Choice Proportions by Deck",
                ylab = "Proportion of Choices",
                ylim = c(0, max(mean_props + se_props) * 1.2),
                col = c("darkred", "red", "lightblue", "steelblue"),
                border = "white",
                cex.lab = 1.2)

  arrows(x0 = bp, y0 = mean_props - se_props,
         x1 = bp, y1 = mean_props + se_props,
         angle = 90, code = 3, length = 0.1, lwd = 2)

  dev.off()

  message("Figure 4 saved.\n")

} else {
  warning("EEF data not found. Skipping Figure 4.")
}

#==============================================================================
# SUMMARY
#==============================================================================

message("\n=== FIGURE GENERATION COMPLETE ===\n")
message(sprintf("Output directory: %s", figures_dir))
message("Figures created:")
message("  - figure1_model_comparison.pdf")
message("  - figure2_group_differences.pdf")
message("  - figure3_parameter_recovery.pdf")
message("  - figure4_posterior_predictive.pdf")
message("\nThese figures are publication-ready PDFs suitable for manuscript submission.")
