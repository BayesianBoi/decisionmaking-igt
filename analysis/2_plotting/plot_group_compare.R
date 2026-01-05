# ===========================================================================
# Group Comparison Plots (Reference r/ Style)
# ===========================================================================
#
# Usage: Rscript analysis/2_plotting/plot_group_compare.R
#
# Generates:
#   - Overlapping posterior density plots by group (HC vs Amph vs Hero)
#   - HDI difference plots (Substance - Control)
#
# ===========================================================================

library(ggplot2)
library(gridExtra)
library(coda)
library(reshape2)

source("analysis/utils/plotting_utils.R")

# ===========================================================================
# Configuration
# ===========================================================================

output_dir <- "results/group_comparison"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Group colors matching reference r/ scripts
group_colors <- c(
  "Control" = "#fb9a99",
  "HC" = "#fb9a99",
  "Amphetamine" = "#1f78b4",
  "Heroin" = "#542788"
)

# ===========================================================================
# Load EEF Samples and Metadata
# ===========================================================================

message("=== GENERATING GROUP COMPARISON PLOTS ===\n")

samples_file <- "results/eef/mcmc_samples.rds"
metadata_file <- "results/subject_metadata.csv"

if (!file.exists(samples_file)) {
  stop("EEF samples not found. Run fit_eef.R first.")
}

if (!file.exists(metadata_file)) {
  # Create metadata if missing
  message("Creating subject metadata...")
  source("analysis/utils/load_data.R")

  dat <- load_all_igt_data()
  clinical <- c("Ahn2014_HC", "Ahn2014_Amph", "Ahn2014_Hero")
  dat <- dat[dat$study %in% clinical, ]

  dat$group <- NA
  dat$group[dat$study == "Ahn2014_HC"] <- "HC"
  dat$group[dat$study == "Ahn2014_Amph"] <- "Amphetamine"
  dat$group[dat$study == "Ahn2014_Hero"] <- "Heroin"

  # Get unique subjects
  subj_df <- unique(dat[, c("subj_unique", "group")])
  subj_df$subj_idx <- as.integer(factor(subj_df$subj_unique, levels = unique(subj_df$subj_unique)))

  dir.create("results", showWarnings = FALSE)
  write.csv(subj_df, metadata_file, row.names = FALSE)
}

samples <- readRDS(samples_file)
mat <- as.matrix(samples)
metadata <- read.csv(metadata_file)

# ===========================================================================
# Helper: Get group-level posterior for a parameter
# ===========================================================================

get_group_posterior <- function(param_name, group_name, mat, metadata) {
  subj_indices <- metadata$subj_idx[metadata$group == group_name]
  col_names <- sprintf("%s[%d]", param_name, subj_indices)
  valid_cols <- col_names[col_names %in% colnames(mat)]

  if (length(valid_cols) == 0) {
    return(NULL)
  }

  # Average across subjects for each posterior draw
  group_samples <- mat[, valid_cols, drop = FALSE]
  return(rowMeans(group_samples))
}

# ===========================================================================
# FIGURE 1: Overlapping Posterior Densities by Group
# ===========================================================================

message("Creating overlapping density plots...")

# Parameters to compare (EEF model)
params <- c("lambda_forget", "theta", "phi", "cons")
param_labels <- list(
  "lambda_forget" = expression(lambda[forget]),
  "theta" = expression(theta),
  "phi" = expression(phi),
  "cons" = expression(c)
)

groups <- c("HC", "Amphetamine", "Heroin")
density_plots <- list()

for (param in params) {
  # Collect posteriors for each group
  group_data <- lapply(groups, function(g) {
    post <- get_group_posterior(param, g, mat, metadata)
    if (!is.null(post)) {
      data.frame(value = post, Group = g)
    } else {
      NULL
    }
  })

  df <- do.call(rbind, group_data)

  if (!is.null(df) && nrow(df) > 0) {
    df$Group <- factor(df$Group, levels = c("HC", "Amphetamine", "Heroin"))

    p <- ggplot(df, aes(x = value, fill = Group, color = Group)) +
      geom_density(alpha = 0.5, linewidth = 0.8) +
      scale_fill_manual(values = group_colors) +
      scale_color_manual(values = group_colors) +
      labs(
        title = param_labels[[param]],
        x = "",
        y = ""
      ) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5)
      )

    density_plots[[param]] <- p
  }
}

if (length(density_plots) > 0) {
  # Create legend from one plot
  legend_plot <- density_plots[[1]] +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(title = "Group"), color = guide_legend(title = "Group"))

  # Extract legend
  g <- ggplotGrob(legend_plot)
  legend_idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  legend <- g$grobs[[legend_idx]]

  # Combine plots
  pdf(file.path(output_dir, "group_densities.pdf"), width = 10, height = 8)
  grid.arrange(
    arrangeGrob(grobs = density_plots, ncol = 2),
    legend,
    heights = c(0.9, 0.1)
  )
  dev.off()

  message("Density plots saved.\n")
}

# ===========================================================================
# FIGURE 2: HDI Difference Plots (Substance - Control)
# ===========================================================================

message("Creating HDI difference plots...")

comparisons <- list(
  "Heroin vs HC" = c("Heroin", "HC"),
  "Amphetamine vs HC" = c("Amphetamine", "HC")
)

# Focus on lambda_forget (forgetting rate)
target_param <- "lambda_forget"

diff_plots <- list()
diff_results <- list()

for (name in names(comparisons)) {
  groups_pair <- comparisons[[name]]
  g1_post <- get_group_posterior(target_param, groups_pair[1], mat, metadata)
  g2_post <- get_group_posterior(target_param, groups_pair[2], mat, metadata)

  if (!is.null(g1_post) && !is.null(g2_post)) {
    diff_post <- g1_post - g2_post

    # Statistics
    mean_diff <- mean(diff_post)
    hdi <- compute_hdi(diff_post, 0.95)
    prob_greater <- mean(diff_post > 0)

    diff_results[[name]] <- list(
      mean_diff = mean_diff,
      ci_lower = hdi["lower"],
      ci_upper = hdi["upper"],
      prob_greater = prob_greater
    )

    message(sprintf(
      "  %s: Mean Diff = %.3f, 95%% HDI [%.3f, %.3f], P(>0) = %.2f",
      name, mean_diff, hdi["lower"], hdi["upper"], prob_greater
    ))

    # Create HDI difference plot
    diff_plots[[name]] <- plot_hdi_diff(diff_post, title = name)
  }
}

if (length(diff_plots) > 0) {
  pdf(file.path(output_dir, "hdi_differences.pdf"), width = 10, height = 4)
  do.call(grid.arrange, c(diff_plots, ncol = length(diff_plots)))
  dev.off()

  message("HDI difference plots saved.\n")
}

# Save results
saveRDS(diff_results, file.path(output_dir, "group_comparisons.rds"))

# ===========================================================================
# FIGURE 3: Forest Plot (Summary)
# ===========================================================================

if (length(diff_results) > 0) {
  message("Creating forest plot...")

  forest_df <- do.call(rbind, lapply(names(diff_results), function(n) {
    data.frame(
      Comparison = n,
      Diff = diff_results[[n]]$mean_diff,
      Lower = diff_results[[n]]$ci_lower,
      Upper = diff_results[[n]]$ci_upper,
      Prob = diff_results[[n]]$prob_greater
    )
  }))

  p_forest <- ggplot(forest_df, aes(x = Diff, y = Comparison)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, color = "black", linewidth = 1) +
    geom_point(size = 4, color = "black") +
    labs(
      title = expression(paste("Difference in Forgetting Rate (", lambda[forget], ")")),
      subtitle = "Posterior Difference (Substance - Control)",
      x = expression(paste(Delta, lambda[forget])),
      y = NULL
    ) +
    theme_publication()

  ggsave(
    file.path(output_dir, "forest_plot.pdf"),
    p_forest,
    width = 8, height = 4
  )

  message("Forest plot saved.\n")
}

message("\n=== GROUP COMPARISON COMPLETE ===\n")
message("Output directory: ", output_dir, "\n")
