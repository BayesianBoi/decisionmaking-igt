#!/usr/bin/env Rscript
# ppc_blockwise.R
# ---------------
# Block-wise analysis of IGT deck preferences.
#
# The IGT has 100 trials which we split into 5 blocks of 20 trials each.
# This lets us see whether participants learn to prefer the advantageous
# decks (C and D) over time. A key signature of healthy IGT performance
# is shifting from A/B to C/D across blocks.
#
# Usage: Rscript ppc_blockwise.R <group>
# Example: Rscript ppc_blockwise.R HC

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("Usage: Rscript ppc_blockwise.R <group>")
}

group <- args[1]

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

cat("\n========================================\n")
cat("Block-wise Analysis:", group, "\n")
cat("========================================\n\n")

# Standard IGT blocking: 5 blocks of 20 trials
n_blocks <- 5
trials_per_block <- 20

# Get subject IDs
subIDs <- unique(group_data$subj)
n_subs <- length(subIDs)

cat("Subjects:", n_subs, "\n\n")

# Count deck selections per block per subject
# observed_block[subject, block, deck] = count
observed_block <- array(0, dim = c(n_subs, n_blocks, 4))

for (s in seq_len(n_subs)) {
    subj_data <- group_data[group_data$subj == subIDs[s], ]

    for (b in seq_len(n_blocks)) {
        start_trial <- (b - 1) * trials_per_block + 1
        end_trial <- min(b * trials_per_block, nrow(subj_data))

        if (start_trial <= nrow(subj_data)) {
            block_choices <- subj_data$choice[start_trial:end_trial]

            for (d in 1:4) {
                observed_block[s, b, d] <- sum(block_choices == d, na.rm = TRUE)
            }
        }
    }
}

# Compute means across subjects
obs_means <- apply(observed_block, c(2, 3), mean)
colnames(obs_means) <- c("A", "B", "C", "D")
rownames(obs_means) <- paste0("Block", 1:n_blocks)

cat("Mean deck selections by block:\n\n")
print(round(obs_means, 2))

# The key metric: advantageous (C+D) vs disadvantageous (A+B)
cat("\nAdvantageous vs Disadvantageous by block:\n")
adv_pref <- obs_means[, 3] + obs_means[, 4]
disadv_pref <- obs_means[, 1] + obs_means[, 2]

for (b in seq_len(n_blocks)) {
    cat(
        "  Block", b, ": C+D =", round(adv_pref[b], 2),
        "  A+B =", round(disadv_pref[b], 2),
        "  Diff =", round(adv_pref[b] - disadv_pref[b], 2), "\n"
    )
}

# Generate a plot showing learning curve
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
        title = paste0(group, ": Observed Deck Preference by Block"),
        x = "Block (20 trials each)",
        y = "Proportion of choices",
        colour = "", linetype = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

plot_file <- file.path(output_dir, paste0("blockwise_", group, ".png"))
ggsave(plot_file, p, width = 8, height = 5, dpi = 300)
cat("\nPlot saved:", plot_file, "\n")
