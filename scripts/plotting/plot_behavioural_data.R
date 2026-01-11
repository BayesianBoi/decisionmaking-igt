# Block-wise Deck Proportions from Raw Behavioural Data
# Creates Figure 1 style plot showing learning curves by group

library(ggplot2)
library(dplyr)

# Load the raw data
data_dir <- "data/raw/Ahn_2014"
output_dir <- "figures/paper"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Group colours
group_colors <- c("HC" = "#4DAF4A", "Amph" = "#E41A1C", "Hero" = "#377EB8")

# Read all three groups
hc <- read.table(file.path(data_dir, "IGTdata_HC.txt"), header = TRUE, sep = "\t")
hc$group <- "HC"

amph <- read.table(file.path(data_dir, "IGTdata_Amph.txt"), header = TRUE, sep = "\t")
amph$group <- "Amph"

hero <- read.table(file.path(data_dir, "IGTdata_Hero.txt"), header = TRUE, sep = "\t")
hero$group <- "Hero"

# Combine
all_data <- rbind(hc, amph, hero)

# Decks 1-2 are disadvantageous (A, B), decks 3-4 are advantageous (C, D)
all_data$advantageous <- ifelse(all_data$deck %in% c(3, 4), 1, 0)

# Assign trials to blocks (10 trials per block)
all_data$block <- ceiling(all_data$trial / 10)

# Compute proportion of advantageous choices per block per subject
subject_block <- all_data %>%
    group_by(group, subjID, block) %>%
    summarise(
        prop_advantageous = mean(advantageous),
        n_trials = n(),
        .groups = "drop"
    )

# Compute group means and SE per block
block_summary <- subject_block %>%
    group_by(group, block) %>%
    summarise(
        mean_prop = mean(prop_advantageous),
        se = sd(prop_advantageous) / sqrt(n()),
        n_subj = n(),
        .groups = "drop"
    )

block_summary$group <- factor(block_summary$group, levels = c("HC", "Amph", "Hero"))

# Create the plot
p_behaviour <- ggplot(block_summary, aes(x = block, y = mean_prop, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = mean_prop - se, ymax = mean_prop + se), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    scale_x_continuous(breaks = 1:10) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
        title = "Learning Curves on the Iowa Gambling Task",
        x = "Block (10 trials per block)",
        y = "Proportion of Advantageous Deck Choices",
        colour = "Group",
        fill = "Group"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank()
    )

ggsave(file.path(output_dir, "behavioural_learning_curves.png"), p_behaviour, width = 10, height = 6, dpi = 600)

cat("Behavioural learning curves saved to:", output_dir, "\n")

# Print sample sizes
cat("\nSample sizes per group:\n")
all_data %>%
    group_by(group) %>%
    summarise(n_subjects = n_distinct(subjID)) %>%
    print()
