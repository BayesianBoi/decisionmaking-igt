# Generate publication-ready figures for the paper
# Creates all figures from model fits and comparisons
# Standardized to IGT field norms (Learning Curves, Forest Plots, etc.)
#
# Run: Rscript analysis/scripts/create_figures.R

library(ggplot2)
library(coda)
library(gridExtra)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Input directories
model_comparison_dir <- "results/model_comparison"
group_comparison_dir <- "results/group_comparison"
parameter_recovery_dir <- "results/parameter_recovery"

# Output directory
figures_dir <- "results/figures"
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Theme for publication
theme_pub <- function() {
  theme_classic() +
    theme(
      text = element_text(size = 12, color = "black"),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "top",
      panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
    )
}

# ==============================================================================
# FIGURE 1: MODEL COMPARISON (WAIC & PARAMS)
# ==============================================================================

message("=== CREATING FIGURE 1: MODEL COMPARISON ===\n")

comparison_file <- file.path(model_comparison_dir, "model_comparison_table.csv")
param_est_file <- file.path(model_comparison_dir, "parameter_estimates.csv")

if (file.exists(comparison_file) && file.exists(param_est_file)) {
  comp_df <- read.csv(comparison_file)
  param_df <- read.csv(param_est_file)

  # Panel A: WAIC Bar Plot (if available, otherwise R-hat)
  # Note: Assuming WAIC will be populated after cloud run
  # Fallback to R-hat for now if WAIC is NA

  if (all(is.na(comp_df$dic))) {
    p_conv <- ggplot(comp_df, aes(x = toupper(model), y = mean_rhat, fill = model)) +
      geom_bar(stat = "identity", width = 0.6, alpha = 0.8) +
      geom_hline(yintercept = 1.1, linetype = "dashed", color = "red") +
      scale_fill_brewer(palette = "Set2") +
      coord_cartesian(ylim = c(0.95, 1.2)) +
      labs(title = "A. Convergence (R-hat)", x = NULL, y = "Mean R-hat") +
      theme_pub() +
      theme(legend.position = "none")
  } else {
    # If WAIC/DIC exists, plot that (lower is better)
    p_conv <- ggplot(comp_df, aes(x = toupper(model), y = dic, fill = model)) +
      geom_bar(stat = "identity", width = 0.6, alpha = 0.8) +
      scale_fill_brewer(palette = "Set2") +
      labs(title = "A. Model Fit (DIC)", x = NULL, y = "DIC (Lower is Better)") +
      theme_pub() +
      theme(legend.position = "none")
  }

  # Panel B: Key Parameters (Focus on EEF)
  # Filter for EEF parameters to show prior/posterior
  eef_params <- subset(param_df, model == "eef" & parameter %in% c("mu_theta", "mu_lambda_forget", "mu_phi"))

  p_param <- ggplot(eef_params, aes(x = parameter, y = mean)) +
    geom_point(size = 4, color = "darkblue") +
    geom_errorbar(aes(ymin = mean - 2 * sd, ymax = mean + 2 * sd), width = 0.2, color = "darkblue") +
    labs(title = "B. EEF Parameter Estimates (Mean ± 2SD)", x = NULL, y = "Estimate") +
    scale_x_discrete(labels = c("lambda_forget" = "Forgetting", "phi" = "Exploration", "theta" = "Sensitivity")) +
    theme_pub()

  pdf(file.path(figures_dir, "figure1_model_comparison.pdf"), width = 10, height = 5)
  grid.arrange(p_conv, p_param, ncol = 2)
  dev.off()

  message("Figure 1 saved.\n")
}

# ==============================================================================
# FIGURE 2: LEARNING CURVES (BEHAVIORAL CHECK)
# ==============================================================================
# Standard IGT Figure: 5 blocks of 20 trials, P(Good Deck)

message("=== CREATING FIGURE 2: LEARNING CURVES ===\n")

# Load Raw Data
source("analysis/utils/load_data.R")
dat_all <- load_all_igt_data()
clinical_studies <- c("Ahn2014_HC", "Ahn2014_Amph", "Ahn2014_Hero", "Fridberg2010_HC", "Fridberg2010_Cbis")
dat <- dat_all[dat_all$study %in% clinical_studies, ]

# Map groups
dat$group_label <- NA
dat$group_label[dat$study %in% c("Ahn2014_HC", "Fridberg2010_HC")] <- "Healthy Control"
dat$group_label[dat$study == "Ahn2014_Amph"] <- "Amphetamine"
dat$group_label[dat$study == "Ahn2014_Hero"] <- "Heroin"
dat$group_label[dat$study == "Fridberg2010_Cbis"] <- "Cannabis"

# Define Good Decks (C & D) vs Bad Decks (A & B)
# Note: Usually C&D are advantageous. A&B are disadvantageous.
dat$is_good <- ifelse(dat$deck %in% c("C", "D"), 1, 0)

# Create Blocks (20 trials each)
dat$block <- ceiling(dat$trial / 20)
dat_agg <- aggregate(is_good ~ group_label + block + subj_unique, data = dat, mean)
dat_summary <- aggregate(is_good ~ group_label + block,
  data = dat_agg,
  FUN = function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x)))
)
dat_summary <- do.call(data.frame, dat_summary)
colnames(dat_summary) <- c("group", "block", "mean", "se")

# Plot
p_learn <- ggplot(dat_summary, aes(x = block, y = mean, color = group, group = group)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Learning Curves by Group",
    subtitle = "Proportion of Advantageous Choices (Decks C & D)",
    x = "Block (20 Trials)",
    y = "P(Advantageous)"
  ) +
  theme_pub()

pdf(file.path(figures_dir, "figure2_learning_curves.pdf"), width = 8, height = 6)
print(p_learn)
dev.off()

message("Figure 2 saved.\n")

# ==============================================================================
# FIGURE 3: GROUP DIFFERENCES (FOREST PLOT STYLE)
# ==============================================================================

message("=== CREATING FIGURE 3: GROUP DIFFERENCES ===\n")

group_comp_file <- file.path(group_comparison_dir, "group_comparisons.rds")

if (file.exists(group_comp_file)) {
  comparisons <- readRDS(group_comp_file)

  # Extract data for plotting
  comps <- c(
    "Heroin vs HC" = "heroin_vs_hc",
    "Amphetamine vs HC" = "amph_vs_hc",
    "Cannabis vs HC" = "cannabis_vs_hc"
  )

  res_list <- list()
  for (n in names(comps)) {
    k <- comps[[n]]
    if (k %in% names(comparisons)) {
      res_list[[n]] <- data.frame(
        Comparison = n,
        Diff = comparisons[[k]]$mean_diff,
        Lower = comparisons[[k]]$ci_lower,
        Upper = comparisons[[k]]$ci_upper,
        Prob = comparisons[[k]]$prob_greater
      )
    }
  }

  df_forest <- do.call(rbind, res_list)

  p_forest <- ggplot(df_forest, aes(x = Diff, y = Comparison)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = 4, color = "black") +
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, color = "black") +
    labs(
      title = "Difference in Forgetting Rate (Lambda)",
      subtitle = "Posterior Difference (Substance - Comp)",
      x = "Difference in Lambda",
      y = NULL
    ) +
    theme_pub()

  pdf(file.path(figures_dir, "figure3_group_differences.pdf"), width = 8, height = 4)
  print(p_forest)
  dev.off()

  message("Figure 3 saved.\n")
}

# ==============================================================================
# FIGURE 4: PARAMETER RECOVERY (SCATTER)
# ==============================================================================

message("=== CREATING FIGURE 4: PARAMETER RECOVERY ===\n")

# Load recovery results data
rec_data_file <- file.path(parameter_recovery_dir, "recovery_results.rds")

if (file.exists(rec_data_file)) {
  recovery_data <- readRDS(rec_data_file)

  # Check structure: Is it the full list (pvl, orl, eef) or old format?
  # We expect a list where each element is a dataframe or a list of dataframes

  # Combine into one dataframe
  if (inherits(recovery_data, "list") && is.data.frame(recovery_data[[1]])) {
    # New format: list of dataframes
    plot_df <- do.call(rbind, recovery_data)
  } else {
    # Old format handling or fallback
    warning("Recovery results format unrecognized. Skipping Figure 4 regeneration.")
    plot_df <- NULL
  }

  if (!is.null(plot_df)) {
    p_rec <- ggplot(plot_df, aes(x = True, y = Recovered)) +
      geom_point(alpha = 0.5, color = "blue") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      facet_wrap(~ Model + Parameter, scales = "free") +
      labs(
        title = "Parameter Recovery (All Models)",
        x = "True Parameter Value",
        y = "Recovered Posterior Mean"
      ) +
      theme_pub()

    pdf(file.path(figures_dir, "figure4_parameter_recovery.pdf"), width = 12, height = 10)
    print(p_rec)
    dev.off()

    message("Figure 4 saved (regenerated from data).\n")
  }
} else {
  warning("recovery_results.rds not found. Run parameter_recovery.R first.")
}

message("\n=== FIGURE GENERATION COMPLETE ===\n")
