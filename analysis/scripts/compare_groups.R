# Compare forgetting rates across clinical groups
# Tests hypothesis: Substance users have higher lambda_forget than HC
#
# Run: Rscript analysis/scripts/compare_groups.R

library(coda)
library(ggplot2)

# Source utility functions
source("analysis/utils/load_data.R")
source("analysis/utils/prepare_eef_data.R")

#==============================================================================
# CONFIGURATION
#==============================================================================

# Model directory
eef_dir <- "results/eef_clinical"

# Output directory
output_dir <- "results/group_comparison"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Group names (in order used during fitting)
group_names <- c("HC", "Amphetamine", "Heroin", "Cannabis")

#==============================================================================
# LOAD DATA AND SAMPLES
#==============================================================================

message("=== LOADING EEF MODEL RESULTS ===\n")

# Load MCMC samples
samples_file <- file.path(eef_dir, "mcmc_samples.rds")
if (!file.exists(samples_file)) {
  stop("EEF model results not found. Run analysis/scripts/fit_eef_clinical.R first.")
}

samples <- readRDS(samples_file)
message("MCMC samples loaded.")

# Load JAGS data to get group assignments
jags_data_file <- file.path(eef_dir, "jags_data.rds")
if (!file.exists(jags_data_file)) {
  stop("JAGS data not found.")
}

jags_data <- readRDS(jags_data_file)
message("JAGS data loaded.")

# Get group assignments
if (is.null(jags_data$group)) {
  stop("No group structure found in JAGS data. Was the model fit with group_var specified?")
}

group_assignments <- jags_data$group
n_groups <- jags_data$n_groups

message(sprintf("Found %d groups: %s\n",
                n_groups, paste(group_names[1:n_groups], collapse = ", ")))

#==============================================================================
# EXTRACT SUBJECT-LEVEL FORGETTING RATES
#==============================================================================

message("=== EXTRACTING FORGETTING RATES ===\n")

# Convert samples to matrix
samples_mat <- as.matrix(samples)

# Extract lambda_forget parameters
lambda_cols <- grep("^lambda_forget\\[", colnames(samples_mat), value = TRUE)
lambda_samples <- samples_mat[, lambda_cols]

n_subjects <- ncol(lambda_samples)
n_iterations <- nrow(lambda_samples)

message(sprintf("Extracted %d subjects x %d iterations", n_subjects, n_iterations))

# Create data frame with posterior means and group labels
lambda_means <- colMeans(lambda_samples)
lambda_sds <- apply(lambda_samples, 2, sd)

subject_df <- data.frame(
  subject = 1:n_subjects,
  group_code = group_assignments,
  group = group_names[group_assignments],
  lambda_mean = lambda_means,
  lambda_sd = lambda_sds
)

# Print group summaries
message("\nGroup summaries (posterior means):")
for (g in 1:n_groups) {
  group_data <- subject_df[subject_df$group_code == g, ]
  message(sprintf("  %s (n=%d): %.3f ± %.3f",
                  group_names[g],
                  nrow(group_data),
                  mean(group_data$lambda_mean),
                  sd(group_data$lambda_mean)))
}

#==============================================================================
# BAYESIAN GROUP COMPARISONS
#==============================================================================

message("\n=== BAYESIAN GROUP COMPARISONS ===\n")

# For each substance group, compute P(lambda_substance > lambda_HC)
# This is done by comparing posterior distributions at the iteration level

# Identify subjects in each group
hc_subjects <- which(group_assignments == 1)  # HC is group 1
amph_subjects <- which(group_assignments == 2)
heroin_subjects <- which(group_assignments == 3)
cannabis_subjects <- which(group_assignments == 4)

# Function to compute group difference probability
compute_group_diff_prob <- function(lambda_samples, group1_idx, group2_idx, group1_name, group2_name) {

  # For each iteration, compute mean difference
  n_iter <- nrow(lambda_samples)
  diff_samples <- numeric(n_iter)

  for (i in 1:n_iter) {
    mean_group1 <- mean(lambda_samples[i, group1_idx])
    mean_group2 <- mean(lambda_samples[i, group2_idx])
    diff_samples[i] <- mean_group1 - mean_group2
  }

  # Compute probability that group1 > group2
  prob_greater <- mean(diff_samples > 0)

  # Compute credible interval for difference
  ci_lower <- quantile(diff_samples, 0.025)
  ci_upper <- quantile(diff_samples, 0.975)

  result <- list(
    comparison = sprintf("%s vs %s", group1_name, group2_name),
    mean_diff = mean(diff_samples),
    sd_diff = sd(diff_samples),
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    prob_greater = prob_greater,
    diff_samples = diff_samples
  )

  return(result)
}

# Compare each substance group to HC
comparisons <- list()

if (length(amph_subjects) > 0) {
  comparisons$amph_vs_hc <- compute_group_diff_prob(
    lambda_samples, amph_subjects, hc_subjects, "Amphetamine", "HC"
  )
}

if (length(heroin_subjects) > 0) {
  comparisons$heroin_vs_hc <- compute_group_diff_prob(
    lambda_samples, heroin_subjects, hc_subjects, "Heroin", "HC"
  )
}

if (length(cannabis_subjects) > 0) {
  comparisons$cannabis_vs_hc <- compute_group_diff_prob(
    lambda_samples, cannabis_subjects, hc_subjects, "Cannabis", "HC"
  )
}

# Print results
message("Group comparisons (substance vs HC):\n")
for (comp_name in names(comparisons)) {
  comp <- comparisons[[comp_name]]

  message(sprintf("%s:", comp$comparison))
  message(sprintf("  Mean difference: %.3f (95%% CI: [%.3f, %.3f])",
                  comp$mean_diff, comp$ci_lower, comp$ci_upper))
  message(sprintf("  P(substance > HC): %.3f", comp$prob_greater))

  if (comp$prob_greater > 0.95) {
    message("  => STRONG evidence for higher forgetting in substance group")
  } else if (comp$prob_greater > 0.80) {
    message("  => MODERATE evidence for higher forgetting in substance group")
  } else if (comp$prob_greater < 0.20) {
    message("  => Evidence for LOWER forgetting in substance group")
  } else {
    message("  => No clear difference")
  }
  message("")
}

# Save comparison results
saveRDS(comparisons, file.path(output_dir, "group_comparisons.rds"))

# Create summary table
comparison_df <- data.frame(
  comparison = sapply(comparisons, function(x) x$comparison),
  mean_diff = sapply(comparisons, function(x) x$mean_diff),
  sd_diff = sapply(comparisons, function(x) x$sd_diff),
  ci_lower = sapply(comparisons, function(x) x$ci_lower),
  ci_upper = sapply(comparisons, function(x) x$ci_upper),
  prob_greater = sapply(comparisons, function(x) x$prob_greater)
)

write.csv(comparison_df, file.path(output_dir, "group_comparison_table.csv"),
          row.names = FALSE)

#==============================================================================
# VISUALIZATIONS
#==============================================================================

message("=== CREATING VISUALIZATIONS ===\n")

# Violin plot of posterior distributions
pdf(file.path(output_dir, "forgetting_by_group.pdf"), width = 10, height = 6)

par(mfrow = c(1, 1))
group_colors <- c("steelblue", "orange", "darkred", "darkgreen")

# Create violin-like density plots
plot(NULL, xlim = c(0.5, n_groups + 0.5), ylim = c(0, 1),
     xlab = "Group", ylab = "Forgetting Rate (lambda)",
     main = "Posterior Distributions of Forgetting Rate by Group",
     xaxt = "n")
axis(1, at = 1:n_groups, labels = group_names[1:n_groups])

for (g in 1:n_groups) {
  group_subjects <- which(group_assignments == g)
  group_lambda <- lambda_samples[, group_subjects]

  # Plot density for each subject (semi-transparent)
  for (s in group_subjects) {
    dens <- density(lambda_samples[, s])
    # Scale density to fit in narrow strip
    dens_scaled <- dens$y / max(dens$y) * 0.3
    polygon(c(g - dens_scaled, rev(g + dens_scaled)),
            c(dens$x, rev(dens$x)),
            col = adjustcolor(group_colors[g], alpha = 0.1),
            border = NA)
  }

  # Add group mean and CI
  group_mean_per_iter <- rowMeans(group_lambda)
  overall_mean <- mean(group_mean_per_iter)
  ci_lower <- quantile(group_mean_per_iter, 0.025)
  ci_upper <- quantile(group_mean_per_iter, 0.975)

  points(g, overall_mean, pch = 19, cex = 2, col = group_colors[g])
  arrows(g, ci_lower, g, ci_upper, angle = 90, code = 3,
         length = 0.1, lwd = 3, col = group_colors[g])
}

legend("topright", legend = group_names[1:n_groups],
       col = group_colors[1:n_groups], pch = 19, cex = 1.2)

dev.off()

message(sprintf("Violin plot saved to: %s",
                file.path(output_dir, "forgetting_by_group.pdf")))

# Difference distributions
pdf(file.path(output_dir, "group_differences.pdf"), width = 10, height = 8)

n_comparisons <- length(comparisons)
par(mfrow = c(ceiling(n_comparisons / 2), 2))

for (comp_name in names(comparisons)) {
  comp <- comparisons[[comp_name]]

  hist(comp$diff_samples, breaks = 50,
       main = comp$comparison,
       xlab = "Difference in Forgetting Rate",
       col = "steelblue",
       border = "white")

  abline(v = 0, col = "red", lty = 2, lwd = 2)
  abline(v = comp$mean_diff, col = "darkblue", lwd = 2)
  abline(v = comp$ci_lower, col = "darkblue", lty = 2)
  abline(v = comp$ci_upper, col = "darkblue", lty = 2)

  legend("topright",
         legend = c(sprintf("P(diff > 0) = %.3f", comp$prob_greater),
                    sprintf("Mean = %.3f", comp$mean_diff)),
         bty = "n")
}

dev.off()

message(sprintf("Difference plots saved to: %s",
                file.path(output_dir, "group_differences.pdf")))

#==============================================================================
# SUMMARY
#==============================================================================

message("\n=== GROUP COMPARISON COMPLETE ===\n")
message(sprintf("Output directory: %s", output_dir))
message("Files created:")
message("  - group_comparisons.rds (full comparison results)")
message("  - group_comparison_table.csv (summary table)")
message("  - forgetting_by_group.pdf (violin plots)")
message("  - group_differences.pdf (posterior difference distributions)")
