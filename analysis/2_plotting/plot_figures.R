# ===========================================================================
# Publication-Quality Figures
# ===========================================================================
#
# What this script does:
#   Generates publication-ready figures for our IGT modeling paper.
#   Style matches Example2.pdf (salmon points, R² annotations, clean theme).
#
# Output figures:
#   1. Parameter Recovery - Can we recover known parameters from simulated data?
#   2. Learning Curves - How do groups differ in behavioral learning?
#   3. Posterior Densities - What are the estimated group-level parameters?
#   4. Model Comparison - Which model fits best? (R-hat, ESS)
#
# Usage:
#   Rscript analysis/2_plotting/plot_figures.R
#
# ===========================================================================

library(ggplot2)
library(gridExtra)
library(coda)

source("analysis/utils/load_data.R")
source("analysis/utils/plotting_utils.R")

# ===========================================================================
# Configuration
# ===========================================================================

figures_dir <- "results/figures"
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# FIGURE 1: PARAMETER RECOVERY
# ===========================================================================
# Purpose: Demonstrate that our models can accurately estimate parameters.
#
# Method: We simulate data from known "true" parameters, fit the model,
# and compare recovered parameters to the true values.
#
# Good recovery = points lie close to diagonal line, R² close to 1.0
# Poor recovery = scattered points, low R²
#
# Style: Salmon points with gray linear fit + SE ribbon (matches Example2.pdf)

message("=== CREATING FIGURE 1: PARAMETER RECOVERY ===\n")

recovery_file <- "results/parameter_recovery/recovery_results.rds"

if (file.exists(recovery_file)) {
  recovery_data <- readRDS(recovery_file)

  # Combine into one dataframe (handles both list and single df formats)
  if (inherits(recovery_data, "list") && is.data.frame(recovery_data[[1]])) {
    plot_df <- do.call(rbind, recovery_data)
  } else if (is.data.frame(recovery_data)) {
    plot_df <- recovery_data
  } else {
    warning("Recovery results format unrecognized.")
    plot_df <- NULL
  }

  if (!is.null(plot_df)) {
    # Compute R² for each Model × Parameter combination
    # R² tells us what proportion of variance in recovered values is explained by true values
    r2_df <- do.call(rbind, by(plot_df, list(plot_df$Model, plot_df$Parameter), function(d) {
      if (nrow(d) > 2 && var(d$True) > 0) {
        r2 <- compute_r2(d$True, d$Recovered)
        return(data.frame(
          Model = d$Model[1],
          Parameter = d$Parameter[1],
          r2 = r2,
          label = sprintf("italic(R)^2 == %.2f", r2) # For plotmath parsing
        ))
      }
      return(NULL)
    }))

    # Create the recovery plot
    p_recovery <- ggplot(plot_df, aes(x = True, y = Recovered)) +
      # Salmon-colored points (our signature style)
      geom_point(color = salmon_color, alpha = 0.6, size = 1.5) +
      # Linear fit with standard error ribbon
      geom_smooth(
        method = "lm", color = "gray50", fill = "gray80",
        se = TRUE, linewidth = 0.8
      ) +
      # Separate panel for each Model × Parameter combination
      facet_wrap(Model ~ Parameter, scales = "free", ncol = 4) +
      # R² annotation in top-left corner of each panel
      geom_text(
        data = r2_df,
        aes(x = -Inf, y = Inf, label = label),
        hjust = -0.1, vjust = 1.5,
        parse = TRUE, color = salmon_color, size = 3
      ) +
      labs(
        title = "Parameter Recovery (Hierarchical)",
        x = "True Value",
        y = "Inferred Value"
      ) +
      theme_publication()

    # Save to PDF (vector format for publication)
    pdf(file.path(figures_dir, "figure1_parameter_recovery.pdf"), width = 12, height = 10)
    print(p_recovery)
    dev.off()

    message("Figure 1 saved.\n")
  }
} else {
  warning("Recovery results not found. Run parameter_recovery.R first.")
}

# ===========================================================================
# FIGURE 2: LEARNING CURVES
# ===========================================================================
# Purpose: Show behavioral differences between clinical groups.
#
# Method: Compute proportion of "advantageous" choices (Decks C & D)
# in blocks of 20 trials. This is the classic IGT learning metric.
#
# Interpretation:
#   - Healthy controls typically learn to prefer good decks over time
#   - Substance users may show flatter curves (poor learning)
#   - The 0.5 dashed line is chance level

message("=== CREATING FIGURE 2: LEARNING CURVES ===\n")

dat_all <- load_all_igt_data()

# Focus on clinical groups from Ahn (2014)
clinical_studies <- c("Ahn2014_HC", "Ahn2014_Amph", "Ahn2014_Hero")
dat <- dat_all[dat_all$study %in% clinical_studies, ]

# Create readable group labels
dat$group_label <- NA
dat$group_label[dat$study == "Ahn2014_HC"] <- "Control"
dat$group_label[dat$study == "Ahn2014_Amph"] <- "Amphetamine"
dat$group_label[dat$study == "Ahn2014_Hero"] <- "Heroin"

# Decks C (3) and D (4) are the "good" decks with positive expected value
dat$is_good <- ifelse(dat$choice %in% c(3, 4), 1, 0)

# Create blocks of 20 trials (standard in IGT literature)
dat$block <- ceiling(dat$trial / 20)

# First, average within each subject (to avoid pseudoreplication)
dat_agg <- aggregate(is_good ~ group_label + block + subj_unique, data = dat, mean)

# Then, compute group means and standard errors
dat_summary <- aggregate(is_good ~ group_label + block,
  data = dat_agg,
  FUN = function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x)))
)
dat_summary <- do.call(data.frame, dat_summary)
colnames(dat_summary) <- c("group", "block", "mean", "se")

# Create the learning curve plot
p_learn <- ggplot(dat_summary, aes(x = block, y = mean, color = group, group = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  # Shaded ribbon for standard error
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = group), alpha = 0.2, color = NA) +
  # Dashed line at chance level (0.5)
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  labs(
    title = "Learning Curves by Group",
    subtitle = "Proportion of Advantageous Choices (Decks C & D)",
    x = "Block (20 Trials)",
    y = "P(Advantageous)",
    color = "Group",
    fill = "Group"
  ) +
  theme_publication()

pdf(file.path(figures_dir, "figure2_learning_curves.pdf"), width = 8, height = 6)
print(p_learn)
dev.off()

message("Figure 2 saved.\n")

# ===========================================================================
# FIGURE 3: POSTERIOR DENSITIES
# ===========================================================================
# Purpose: Visualize the posterior distributions of group-level parameters.
#
# Style: Black line density (no fill) matching Example2.pdf Figure 2a.
# Each panel shows one parameter with N and bandwidth annotations.
#
# We also create a summary table with MAP (mode) and 95% CI for each parameter.

message("=== CREATING FIGURE 3: POSTERIOR DENSITIES ===\n")

eef_samples_file <- "results/eef/mcmc_samples.rds"

if (file.exists(eef_samples_file)) {
  samples <- readRDS(eef_samples_file)
  mat <- as.matrix(samples)

  # Group-level parameters to plot
  params <- c("mu_theta", "mu_lambda_forget", "mu_phi", "mu_cons")
  param_labels <- c(
    "mu_theta" = expression(mu[theta]),
    "mu_lambda_forget" = expression(mu[lambda]),
    "mu_phi" = expression(mu[phi]),
    "mu_cons" = expression(mu[c])
  )

  plot_list <- list()

  for (param in params) {
    if (param %in% colnames(mat)) {
      vals <- mat[, param]
      df <- data.frame(value = vals)

      # Compute summary statistics
      dens <- density(vals)
      map_val <- dens$x[which.max(dens$y)] # Maximum a posteriori (mode)
      ci <- quantile(vals, c(0.025, 0.975)) # 95% credible interval
      n_samples <- length(vals)
      bw <- dens$bw # Bandwidth used for smoothing

      p <- ggplot(df, aes(x = value)) +
        # Black line density (matching Example2.pdf style)
        geom_density(color = "black", linewidth = 0.8, fill = NA) +
        labs(
          title = param_labels[[param]],
          x = NULL,
          y = "Density"
        ) +
        # N and bandwidth annotation (matching Example2.pdf)
        annotate(
          "text",
          x = Inf, y = -Inf,
          label = sprintf("N = %d  Bandwidth = %.4f", n_samples, bw),
          hjust = 1.1, vjust = -0.5, size = 2.5
        ) +
        theme_minimal(base_size = 10) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 12),
          axis.title.y = element_text(size = 9)
        )

      plot_list[[param]] <- p
    }
  }

  if (length(plot_list) > 0) {
    # Create summary table (like Example2.pdf Figure 2b)
    summary_table <- data.frame(
      Parameter = character(),
      MAP = numeric(),
      CI_Lower = numeric(),
      CI_Upper = numeric(),
      stringsAsFactors = FALSE
    )

    for (param in params) {
      if (param %in% colnames(mat)) {
        vals <- mat[, param]
        dens <- density(vals)
        map_val <- dens$x[which.max(dens$y)]
        ci <- quantile(vals, c(0.025, 0.975))
        summary_table <- rbind(summary_table, data.frame(
          Parameter = param,
          MAP = round(map_val, 4),
          CI_Lower = round(ci[1], 4),
          CI_Upper = round(ci[2], 4)
        ))
      }
    }

    # Save plots
    pdf(file.path(figures_dir, "figure3_posterior_densities.pdf"), width = 8, height = 6)
    do.call(grid.arrange, c(plot_list, ncol = 2))
    dev.off()

    # Save summary table as CSV
    write.csv(summary_table, file.path(figures_dir, "posterior_summary_table.csv"), row.names = FALSE)

    message("Figure 3 saved.\n")
  }
} else {
  warning("EEF samples not found. Run fit_eef.R first.")
}

# ===========================================================================
# FIGURE 4: MODEL COMPARISON
# ===========================================================================
# Purpose: Compare model quality across PVL-Delta, ORL, and EEF.
#
# Metrics shown:
#   - Mean R-hat: convergence diagnostic (should be < 1.1)
#   - Median ESS: effective sample size (should be > 1000)
#
# Better models have lower R-hat and higher ESS.

message("=== CREATING FIGURE 4: MODEL COMPARISON ===\n")

comparison_file <- "results/model_comparison/model_comparison_table.csv"

if (file.exists(comparison_file)) {
  comp_df <- read.csv(comparison_file)

  # Panel A: R-hat bar plot (lower is better, threshold at 1.1)
  p_rhat <- ggplot(comp_df, aes(x = model, y = mean_rhat)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.6) +
    geom_hline(yintercept = 1.1, linetype = "dashed", color = "red") +
    geom_text(aes(label = sprintf("%.3f", mean_rhat)), vjust = -0.5, size = 3) +
    coord_cartesian(ylim = c(0.95, 1.2)) +
    labs(title = "Convergence (Mean R-hat)", x = NULL, y = "Mean R-hat") +
    theme_publication()

  # Panel B: ESS bar plot (higher is better, target at 1000)
  p_ess <- ggplot(comp_df, aes(x = model, y = median_ess)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.6) +
    geom_hline(yintercept = 1000, linetype = "dashed", color = "red") +
    geom_text(aes(label = sprintf("%.0f", median_ess)), vjust = -0.5, size = 3) +
    labs(title = "Median Effective Sample Size", x = NULL, y = "Median ESS") +
    theme_publication()

  pdf(file.path(figures_dir, "figure4_model_comparison.pdf"), width = 10, height = 5)
  grid.arrange(p_rhat, p_ess, ncol = 2)
  dev.off()

  message("Figure 4 saved.\n")
} else {
  warning("Model comparison results not found. Run compare_models.R first.")
}

message("\n=== FIGURE GENERATION COMPLETE ===\n")
message("Output directory: ", figures_dir, "\n")
