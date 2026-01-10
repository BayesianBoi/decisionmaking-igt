#!/usr/bin/env Rscript
# ==============================================================================
# PPC: Deck Proportions Over Blocks (One-Step-Ahead Prediction)
# ==============================================================================
# Uses observed choice history to compute predicted choice probabilities
# at each trial, then averages to get expected deck proportions per block.
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(42)

# ==============================================================================
# Settings
# ==============================================================================
n_posterior_samples <- 100
block_size <- 20
n_trials <- 100
n_blocks <- n_trials / block_size

groups <- c("HC", "Amph", "Hero")
models <- c("orl", "pvl_delta", "eef")
deck_colors <- c("A" = "#E41A1C", "B" = "#FF7F00", "C" = "#377EB8", "D" = "#4DAF4A")

output_dir <- "outputs/figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Load Data
# ==============================================================================
source("utils/load_data.R")
all_data <- load_all_igt_data()

# ==============================================================================
# One-Step-Ahead Prediction Functions
# ==============================================================================
# These compute p(choice_t | history_1:t-1, params) using observed outcomes

softmax <- function(x) {
    exp_x <- exp(x - max(x))
    return(exp_x / sum(exp_x))
}

# ORL: Compute choice probabilities for all trials given observed history
predict_orl <- function(params, choices, outcomes, n_trials) {
    a_rew <- params$a_rew
    a_pun <- params$a_pun
    K <- params$K
    omega_f <- params$omega_f
    omega_p <- params$omega_p
    theta <- 1

    # Initialize
    Ev <- rep(0, 4)
    Ef <- rep(0, 4)
    PS <- rep(1, 4)

    # Store probabilities for each trial
    probs <- matrix(0.25, nrow = n_trials, ncol = 4) # Trial 1 = uniform

    for (t in 2:n_trials) {
        prev_choice <- choices[t - 1]
        X <- outcomes[t - 1]
        signX <- ifelse(X < 0, -1, 1)

        # Update for chosen deck
        if (X >= 0) {
            Ev[prev_choice] <- Ev[prev_choice] + a_rew * (X - Ev[prev_choice])
            Ef[prev_choice] <- Ef[prev_choice] + a_rew * (signX - Ef[prev_choice])
        } else {
            Ev[prev_choice] <- Ev[prev_choice] + a_pun * (X - Ev[prev_choice])
            Ef[prev_choice] <- Ef[prev_choice] + a_pun * (signX - Ef[prev_choice])
        }

        # Update Ef for unchosen decks
        for (d in 1:4) {
            if (d != prev_choice) {
                if (X >= 0) {
                    Ef[d] <- Ef[d] + a_pun * ((-signX / 3) - Ef[d])
                } else {
                    Ef[d] <- Ef[d] + a_rew * ((-signX / 3) - Ef[d])
                }
            }
        }

        # Update perseveration
        PS <- PS / (1 + K)
        PS[prev_choice] <- 1 / (1 + K)

        # Compute choice probabilities for trial t
        V <- Ev + Ef * omega_f + PS * omega_p
        probs[t, ] <- softmax(theta * V)
    }

    return(probs)
}

# PVL-Delta: Compute choice probabilities
predict_pvl_delta <- function(params, choices, outcomes, n_trials) {
    A <- params$A
    w <- params$w
    a <- params$a
    theta <- params$theta

    Ev <- rep(0, 4)
    probs <- matrix(0.25, nrow = n_trials, ncol = 4)

    for (t in 2:n_trials) {
        prev_choice <- choices[t - 1]
        X <- outcomes[t - 1]

        # Prospect utility
        if (X >= 0) {
            u <- X^A
        } else {
            u <- -w * abs(X)^A
        }

        # Delta rule update
        Ev[prev_choice] <- Ev[prev_choice] + a * (u - Ev[prev_choice])

        # Choice probabilities
        probs[t, ] <- softmax(theta * Ev)
    }

    return(probs)
}

# EEF: Compute choice probabilities
predict_eef <- function(params, choices, outcomes, n_trials) {
    theta_param <- params$theta
    lambda <- params$lambda
    phi <- params$phi
    cons <- params$cons

    Ev <- rep(0, 4)
    Explore <- rep(0, 4)

    probs <- matrix(0.25, nrow = n_trials, ncol = 4)

    for (t in 2:n_trials) {
        prev_choice <- choices[t - 1]
        X <- outcomes[t - 1]

        # Apply forgetting
        Ev <- Ev * lambda
        Explore <- Explore * lambda

        # Update chosen deck
        Ev[prev_choice] <- Ev[prev_choice] + X
        Explore[prev_choice] <- 0

        # Update exploration for unchosen
        for (d in 1:4) {
            if (d != prev_choice) {
                Explore[d] <- Explore[d] + 1
            }
        }

        # Choice probabilities
        V <- theta_param * Ev + phi * Explore + cons
        probs[t, ] <- softmax(V)
    }

    return(probs)
}

# ==============================================================================
# Compute Observed Proportions
# ==============================================================================
cat("Computing observed proportions...\n")

observed_list <- list()

for (group in groups) {
    study_label <- paste0("Ahn2014_", group)
    group_data <- all_data[all_data$study == study_label, ]
    subjects <- unique(group_data$subj)

    all_block_props <- list()

    for (i in seq_along(subjects)) {
        subj_data <- group_data[group_data$subj == subjects[i], ]
        n_t <- nrow(subj_data)
        n_b <- floor(n_t / block_size)

        if (n_b < 1) next

        for (b in 1:min(n_b, n_blocks)) {
            start_idx <- (b - 1) * block_size + 1
            end_idx <- b * block_size
            block_choices <- subj_data$choice[start_idx:end_idx]

            for (d in 1:4) {
                key <- paste(b, d, sep = "_")
                if (is.null(all_block_props[[key]])) {
                    all_block_props[[key]] <- numeric(0)
                }
                all_block_props[[key]] <- c(
                    all_block_props[[key]],
                    sum(block_choices == d) / block_size
                )
            }
        }
    }

    for (b in 1:n_blocks) {
        for (d in 1:4) {
            key <- paste(b, d, sep = "_")
            props <- all_block_props[[key]]
            if (!is.null(props) && length(props) > 0) {
                observed_list[[length(observed_list) + 1]] <- data.frame(
                    group = group,
                    source = "Behavior",
                    block = b,
                    deck = LETTERS[d],
                    mean = mean(props, na.rm = TRUE),
                    lower = mean(props, na.rm = TRUE),
                    upper = mean(props, na.rm = TRUE)
                )
            }
        }
    }
}

observed_df <- do.call(rbind, observed_list)

# ==============================================================================
# Compute Model Predictions (One-Step-Ahead)
# ==============================================================================
cat("Computing model predictions (one-step-ahead)...\n")

model_list <- list()

for (model in models) {
    cat(sprintf("  Processing %s...\n", model))

    for (group in groups) {
        cat(sprintf("    Group %s...\n", group))

        fit_path <- sprintf("outputs/parameter_estimation/%s_fit_%s.rds", model, group)
        if (!file.exists(fit_path)) {
            cat(sprintf("    Skipping: %s\n", fit_path))
            next
        }

        fit <- readRDS(fit_path)
        sims <- fit$BUGSoutput$sims.list

        study_label <- paste0("Ahn2014_", group)
        group_data <- all_data[all_data$study == study_label, ]
        subjects <- unique(group_data$subj)
        n_subjects <- length(subjects)

        # Sample posterior indices
        n_iter <- nrow(sims[[1]])
        if (is.null(dim(sims[[1]]))) n_iter <- length(sims[[1]])
        sample_idx <- sample(1:n_iter, min(n_posterior_samples, n_iter))

        # Store: [sample, block, deck] - average predicted probabilities
        pred_props <- array(NA, dim = c(length(sample_idx), n_blocks, 4))

        for (s_idx in seq_along(sample_idx)) {
            idx <- sample_idx[s_idx]

            # Collect predicted probabilities across all subjects
            all_probs <- list()

            for (subj_i in 1:n_subjects) {
                subj_data <- group_data[group_data$subj == subjects[subj_i], ]
                choices <- subj_data$choice
                outcomes <- (subj_data$gain + subj_data$loss) / 100
                n_t <- nrow(subj_data)

                # Extract parameters
                if (model == "orl") {
                    params <- list(
                        a_rew = sims$a_rew[idx, subj_i],
                        a_pun = sims$a_pun[idx, subj_i],
                        K = sims$K[idx, subj_i],
                        omega_f = sims$omega_f[idx, subj_i],
                        omega_p = sims$omega_p[idx, subj_i]
                    )
                    probs <- predict_orl(params, choices, outcomes, n_t)
                } else if (model == "pvl_delta") {
                    params <- list(
                        A = sims$A[idx, subj_i],
                        w = sims$w[idx, subj_i],
                        a = sims$a[idx, subj_i],
                        theta = sims$theta[idx, subj_i]
                    )
                    probs <- predict_pvl_delta(params, choices, outcomes, n_t)
                } else if (model == "eef") {
                    params <- list(
                        theta = sims$theta[idx, subj_i],
                        lambda = sims$lambda[idx, subj_i],
                        phi = sims$phi[idx, subj_i],
                        cons = sims$cons[idx, subj_i]
                    )
                    probs <- predict_eef(params, choices, outcomes, n_t)
                }

                all_probs[[subj_i]] <- probs
            }

            # Average across subjects per block
            for (b in 1:n_blocks) {
                start_idx <- (b - 1) * block_size + 1
                end_idx <- min(b * block_size, n_trials)

                # Collect block probabilities from all subjects
                block_probs <- matrix(0, nrow = 0, ncol = 4)
                for (subj_i in 1:n_subjects) {
                    subj_probs <- all_probs[[subj_i]]
                    if (nrow(subj_probs) >= end_idx) {
                        block_probs <- rbind(block_probs, subj_probs[start_idx:end_idx, ])
                    }
                }

                # Mean probability per deck in this block
                if (nrow(block_probs) > 0) {
                    pred_props[s_idx, b, ] <- colMeans(block_probs)
                }
            }
        }

        # Compute mean and 95% CI across posterior samples
        mean_props <- apply(pred_props, c(2, 3), mean, na.rm = TRUE)
        lower_props <- apply(pred_props, c(2, 3), quantile, probs = 0.025, na.rm = TRUE)
        upper_props <- apply(pred_props, c(2, 3), quantile, probs = 0.975, na.rm = TRUE)

        model_label <- toupper(model)
        if (model == "pvl_delta") model_label <- "PVL-Delta"

        for (b in 1:n_blocks) {
            for (d in 1:4) {
                model_list[[length(model_list) + 1]] <- data.frame(
                    group = group,
                    source = model_label,
                    block = b,
                    deck = LETTERS[d],
                    mean = mean_props[b, d],
                    lower = lower_props[b, d],
                    upper = upper_props[b, d]
                )
            }
        }
    }
}

model_df <- do.call(rbind, model_list)

# ==============================================================================
# Combine and Plot
# ==============================================================================
cat("Creating plot...\n")

plot_df <- rbind(observed_df, model_df)
plot_df$source <- factor(plot_df$source, levels = c("Behavior", "ORL", "PVL-Delta", "EEF"))
plot_df$group <- factor(plot_df$group,
    levels = c("HC", "Amph", "Hero"),
    labels = c("Healthy Controls", "Amphetamine", "Heroin")
)
plot_df$deck <- factor(plot_df$deck, levels = c("A", "B", "C", "D"))

p <- ggplot(plot_df, aes(x = block, y = mean, color = deck, fill = deck)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    facet_grid(group ~ source) +
    scale_color_manual(values = deck_colors, name = "Deck") +
    scale_fill_manual(values = deck_colors, name = "Deck") +
    scale_x_continuous(breaks = 1:5, labels = 1:5) +
    labs(x = "Block", y = "Choice Proportion") +
    ylim(0, 0.5) +
    theme_minimal(base_size = 11) +
    theme(
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.5, "lines")
    )

ggsave(file.path(output_dir, "ppc_deck_proportions.pdf"), p, width = 12, height = 8)
ggsave(file.path(output_dir, "ppc_deck_proportions.png"), p, width = 12, height = 8, dpi = 300)

cat(sprintf("\nPlot saved to: %s\n", file.path(output_dir, "ppc_deck_proportions.pdf")))
cat("Done!\n")
