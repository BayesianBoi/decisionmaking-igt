#!/usr/bin/env Rscript
# ppc_blockwise.R
# ---------------
# Block-wise posterior predictive check for IGT models.
#
# The IGT has 100 trials which we split into 5 blocks of 20 trials each.
# This lets us see whether participants (and the model) learn to prefer
# the advantageous decks (C and D) over time. A key signature of healthy
# IGT performance is shifting from A/B to C/D across blocks.
#
# If the model captures learning dynamics well, predicted block preferences
# should match observed preferences.
#
# Usage: Rscript ppc_blockwise.R <model> <group>
# Example: Rscript ppc_blockwise.R orl HC

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: Rscript ppc_blockwise.R <model> <group>")
}

model_name <- args[1]
group <- args[2]

library(ggplot2)

# Map short group names to full dataset identifiers
groups_map <- c(
    "HC" = "Ahn2014_HC",
    "Amph" = "Ahn2014_Amph",
    "Hero" = "Ahn2014_Hero"
)

if (!group %in% names(groups_map)) {
    stop("Invalid group. Use: HC, Amph, or Hero")
}

# Load behavioural data
source("utils/load_data.R")
all_data <- load_all_igt_data()
group_data <- all_data[all_data$study == groups_map[[group]], ]

# Load the fitted model
fit_file <- file.path("outputs", model_name, paste0(model_name, "_fit_", group, ".rds"))
if (!file.exists(fit_file)) {
    stop("Fit file not found: ", fit_file)
}

fit <- readRDS(fit_file)

cat("\n========================================\n")
cat("Block-wise PPC:", toupper(model_name), "-", group, "\n")
cat("========================================\n\n")

# Standard IGT blocking: 5 blocks of 20 trials
n_blocks <- 5
trials_per_block <- 20

# Get subject IDs
subIDs <- unique(group_data$subj)
n_subs <- length(subIDs)

# Count deck selections per block per subject
# observed_block[subject, block, deck] = count
observed_block <- array(0, dim = c(n_subs, n_blocks, 4))

for (s in seq_len(n_subs)) {
    subj_data <- group_data[group_data$subj == subIDs[s], ]

    for (b in seq_len(n_blocks)) {
        # Trial indices for this block
        start_trial <- (b - 1) * trials_per_block + 1
        end_trial <- min(b * trials_per_block, nrow(subj_data))

        if (start_trial <= nrow(subj_data)) {
            block_choices <- subj_data$choice[start_trial:end_trial]

            # Count how many times each deck was chosen
            for (d in 1:4) {
                observed_block[s, b, d] <- sum(block_choices == d, na.rm = TRUE)
            }
        }
    }
}

# Check if we have the p parameter (trial-by-trial choice probabilities)
if (!"p" %in% names(fit$BUGSoutput$sims.list)) {
    cat("Note: The 'p' parameter was not saved during fitting.\n")
    cat("Showing observed data only. To get predictions, re-fit with 'p' in parameters.to.save.\n\n")

    # Report observed block preferences
    cat("Observed block-wise deck selection (mean across subjects):\n\n")

    obs_means <- apply(observed_block, c(2, 3), mean)
    colnames(obs_means) <- c("A", "B", "C", "D")
    rownames(obs_means) <- paste0("Block", 1:n_blocks)
    print(round(obs_means, 2))

    # The key metric: advantageous (C+D) vs disadvantageous (A+B)
    cat("\nAdvantageous deck preference (C+D) by block:\n")
    adv_pref <- obs_means[, 3] + obs_means[, 4]
    disadv_pref <- obs_means[, 1] + obs_means[, 2]

    for (b in seq_len(n_blocks)) {
        cat(
            "  Block", b, ": C+D =", round(adv_pref[b], 2),
            "  A+B =", round(disadv_pref[b], 2),
            "  Diff =", round(adv_pref[b] - disadv_pref[b], 2), "\n"
        )
    }
} else {
    # We have the p parameter, so we can compare predicted vs observed
    p_array <- fit$BUGSoutput$sims.list$p
    n_samples <- dim(p_array)[1]
    n_trials_p <- dim(p_array)[3]

    cat("Computing predicted block preferences from", n_samples, "posterior samples...\n")
    cat("p array has", n_trials_p, "trials\n\n")

    # Average choice probability across posterior samples
    p_mean <- apply(p_array, c(2, 3, 4), mean)

    # Sum up predicted probabilities within each block
    predicted_block <- array(0, dim = c(n_subs, n_blocks, 4))

    for (s in seq_len(n_subs)) {
        for (b in seq_len(n_blocks)) {
            start_trial <- (b - 1) * trials_per_block + 1
            end_trial <- min(b * trials_per_block, n_trials_p)

            for (d in 1:4) {
                # Handle case where p array might have fewer trials
                if (end_trial >= start_trial && end_trial <= n_trials_p) {
                    predicted_block[s, b, d] <- sum(p_mean[s, start_trial:end_trial, d], na.rm = TRUE)
                }
            }
        }
    }

    pred_means <- apply(predicted_block, c(2, 3), mean)
    obs_means <- apply(observed_block, c(2, 3), mean)

    colnames(pred_means) <- c("A", "B", "C", "D")
    colnames(obs_means) <- c("A", "B", "C", "D")
    rownames(pred_means) <- paste0("Block", 1:n_blocks)
    rownames(obs_means) <- paste0("Block", 1:n_blocks)

    cat("Observed:\n")
    print(round(obs_means, 2))
    cat("\nPredicted:\n")
    print(round(pred_means, 2))

    # How well do predictions match observations?
    cor_val <- cor(as.vector(obs_means), as.vector(pred_means))
    cat("\nCorrelation (observed vs predicted):", round(cor_val, 3), "\n")
}

# Generate a plot showing learning curve
obs_means <- apply(observed_block, c(2, 3), mean)

# Calculate proportion choosing advantageous vs disadvantageous
plot_data <- data.frame(
    block = rep(1:n_blocks, 2),
    proportion = c(
        rowSums(obs_means[, 3:4]) / rowSums(obs_means),
        rowSums(obs_means[, 1:2]) / rowSums(obs_means)
    ),
    type = rep(c("Advantageous (C+D)", "Disadvantageous (A+B)"), each = n_blocks)
)

output_dir <- file.path("plots/ppc")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

p <- ggplot(plot_data, aes(x = block, y = proportion, colour = type, linetype = type)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(breaks = 1:5) +
    scale_colour_manual(values = c("Advantageous (C+D)" = "#2E7D32", "Disadvantageous (A+B)" = "#C62828")) +
    labs(
        title = paste0(toupper(model_name), " - ", group, ": Observed Deck Preference by Block"),
        x = "Block (20 trials each)",
        y = "Proportion of choices",
        colour = "", linetype = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

plot_file <- file.path(output_dir, paste0("blockwise_", model_name, "_", group, ".png"))
ggsave(plot_file, p, width = 8, height = 5, dpi = 300)
cat("\nPlot saved:", plot_file, "\n")
