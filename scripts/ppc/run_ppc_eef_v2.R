#!/usr/bin/env Rscript
# =============================================================================
# INDIVIDUAL PPC - EEF v2 Model
# =============================================================================
# Fits EEF v2 model individually for each subject and computes PPC
# This version uses separate Gain/Loss inputs per Yang et al. 2025
# Usage: Rscript scripts/ppc/run_ppc_eef_v2.R <group>
# Example: Rscript scripts/ppc/run_ppc_eef_v2.R HC
# =============================================================================

# Setup
required_packages <- c("R2jags", "parallel", "ggplot2")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")
invisible(lapply(required_packages, library, character.only = TRUE))

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript run_ppc_eef_v2.R <group>\nGroup: HC, Amph, Hero")
}

target_group <- args[1]

groups_map <- c(
    "HC" = "Ahn2014_HC",
    "Amph" = "Ahn2014_Amph",
    "Hero" = "Ahn2014_Hero"
)

if (!target_group %in% names(groups_map)) {
    stop("Invalid group. Must be one of: ", paste(names(groups_map), collapse = ", "))
}
study_label <- groups_map[[target_group]]

cat("\n================================================\n")
cat("INDIVIDUAL PPC: EEF v2 -", target_group, "\n")
cat("================================================\n")

# Load Data
source("utils/load_data.R")
all_data <- load_all_igt_data()
raw_data <- all_data[all_data$study == study_label, ]

subIDs <- unique(raw_data$subj)
nsubs <- length(subIDs)

cat("Found", nsubs, "subjects\n")

# MPD Function
MPD <- function(x) {
    x <- as.numeric(x)
    if (length(x) < 2) {
        return(NA)
    }
    dens <- density(x)
    dens$x[which.max(dens$y)]
}

model_file <- "models/eef_v2_individual.txt"
if (!file.exists(model_file)) {
    stop("Model file not found: ", model_file)
}

# Parallel processing
n_cores <- min(detectCores() - 1, 10)
cat("Using", n_cores, "cores for parallel processing\n")

start_time <- Sys.time()

process_subject <- function(s) {
    subj_df <- raw_data[raw_data$subj == subIDs[s], ]
    ntrials <- nrow(subj_df)

    # Prepare data with separate Gain and Loss
    x <- subj_df$choice
    Gain <- subj_df$gain / 100
    Loss <- subj_df$loss / 100

    jags_data <- list(
        "x" = x,
        "Gain" = Gain,
        "Loss" = Loss,
        "ntrials" = ntrials
    )

    params <- c("p")

    result <- tryCatch(
        {
            fit <- jags(
                data = jags_data,
                inits = NULL,
                parameters.to.save = params,
                model.file = model_file,
                n.chains = 3,
                n.iter = 6000,
                n.burnin = 1000,
                n.thin = 1,
                quiet = TRUE
            )

            # Extract p and compute one-step-ahead accuracy
            p_post <- fit$BUGSoutput$sims.list$p

            x_predict <- numeric(ntrials)
            x_predict[1] <- NA

            for (t in 2:ntrials) {
                p_predict <- c(
                    MPD(p_post[, t, 1]),
                    MPD(p_post[, t, 2]),
                    MPD(p_post[, t, 3]),
                    MPD(p_post[, t, 4])
                )
                x_predict[t] <- which.max(p_predict)
            }

            acc <- sum(x_predict[2:ntrials] == x[2:ntrials], na.rm = TRUE) / (ntrials - 1)
            cat(sprintf("Subject %d/%d: accuracy = %.3f\n", s, nsubs, acc))
            acc
        },
        error = function(e) {
            cat(sprintf("Subject %d: ERROR - %s\n", s, e$message))
            NA
        }
    )
    return(result)
}

pred_success <- unlist(mclapply(1:nsubs, process_subject, mc.cores = n_cores))

end_time <- Sys.time()
cat("\nTotal time:", difftime(end_time, start_time, units = "mins"), "minutes\n")

# Summary
mean_acc <- mean(pred_success, na.rm = TRUE)
sd_acc <- sd(pred_success, na.rm = TRUE)

cat("\n=== EEF v2 INDIVIDUAL PPC RESULTS ===\n")
cat(sprintf("Group: %s\n", target_group))
cat(sprintf("Mean accuracy: %.3f (SD: %.3f)\n", mean_acc, sd_acc))
cat(sprintf("Range: %.3f - %.3f\n", min(pred_success, na.rm = TRUE), max(pred_success, na.rm = TRUE)))

# Save Results
output_dir <- "outputs/ppc"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pred_df <- data.frame(
    subject = 1:nsubs,
    accuracy = pred_success,
    mean_acc = mean_acc,
    sd_acc = sd_acc
)

results_path <- file.path(output_dir, paste0("ppc_individual_eef_v2_", target_group, ".rds"))
saveRDS(list(accuracy = pred_success, mean = mean_acc, sd = sd_acc, df = pred_df), results_path)
cat("Results saved to:", results_path, "\n")

# Plot
colors <- list(HC = "lightgreen", Amph = "salmon", Hero = "lightblue")

p <- ggplot(pred_df, aes(x = subject, y = accuracy)) +
    geom_ribbon(aes(ymin = mean_acc - sd_acc, ymax = mean_acc + sd_acc),
        fill = colors[[target_group]], alpha = 0.3
    ) +
    geom_hline(aes(yintercept = mean_acc), size = 1) +
    geom_hline(yintercept = 0.25, linetype = "dashed", color = "black", size = 0.8) +
    geom_point(color = colors[[target_group]], size = 2.5, alpha = 0.6) +
    ylim(0, 0.8) +
    labs(
        title = paste("Individual PPC: EEF v2 -", target_group),
        subtitle = sprintf("Mean Accuracy: %.3f (SD: %.3f)", mean_acc, sd_acc),
        x = "Subject",
        y = "Proportion Correctly Predicted"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
    )

plot_path <- file.path(output_dir, paste0("ppc_individual_eef_v2_", target_group, ".png"))
ggsave(plot_path, plot = p, width = 12, height = 6, dpi = 300)
cat("Plot saved to:", plot_path, "\n")
