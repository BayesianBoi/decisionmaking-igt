#!/usr/bin/env Rscript
# compare_orl.R
# -------------
# Joint hierarchical ORL model for comparing two groups.
#
# This script fits a single model to both groups simultaneously, with
# shared grand means (mu) and group difference parameters (alpha).
# The parameterisation is:
#   Group1 mean = mu - alpha/2
#   Group2 mean = mu + alpha/2
#
# This approach lets us directly test whether alpha differs from zero,
# which tells us if the groups differ on each parameter.
#
# Usage: Rscript compare_orl.R [Group1] [Group2]
# Example: Rscript compare_orl.R HC Hero

if (!require("pacman")) install.packages("pacman")
pacman::p_load(R2jags, parallel)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Must specify TWO groups. Example: HC Hero")
}

g1_name <- args[1]
g2_name <- args[2]

# Map short names to dataset labels
groups_map <- c(
    "HC" = "Ahn2014_HC",
    "Amph" = "Ahn2014_Amph",
    "Hero" = "Ahn2014_Hero"
)

if (!g1_name %in% names(groups_map) || !g2_name %in% names(groups_map)) {
    stop("Invalid groups. Valid options: HC, Amph, Hero")
}

g1_label <- groups_map[[g1_name]]
g2_label <- groups_map[[g2_name]]

cat("========================================\n")
cat("Joint Comparison: ORL\n")
cat(g1_name, "vs", g2_name, "\n")
cat("========================================\n\n")

# Load behavioural data
source("utils/load_data.R")
all_data <- load_all_igt_data()

# Function to prepare data for one group
prep_group_data <- function(data, study_label) {
    raw <- data[data$study == study_label, ]
    subIDs <- unique(raw$subj)
    nsubs <- length(subIDs)
    ntrials_max <- 100

    # Pre-allocate arrays
    ntrials_all <- array(0, c(nsubs))
    x_all <- array(0, c(nsubs, ntrials_max))
    X_all <- array(0, c(nsubs, ntrials_max))

    for (s in seq_len(nsubs)) {
        subj_df <- raw[raw$subj == subIDs[s], ]
        ntrials_all[s] <- nrow(subj_df)

        # Choices (1-4 for decks A-D)
        x_sub <- subj_df$choice
        length(x_sub) <- ntrials_max
        x_all[s, ] <- x_sub

        # Net outcomes (gain + loss, where loss is negative)
        X_raw <- subj_df$gain + subj_df$loss

        # Scale by dividing by 100 (matches hBayesDM convention)
        X_sub <- X_raw / 100
        length(X_sub) <- ntrials_max
        X_all[s, ] <- X_sub
    }

    return(list(x = x_all, X = X_all, ntrials = ntrials_all, nsubs = nsubs))
}

# Prepare data for both groups
g1_data <- prep_group_data(all_data, g1_label)
g2_data <- prep_group_data(all_data, g2_label)

cat("Group 1 (", g1_name, "): n =", g1_data$nsubs, "subjects\n")
cat("Group 2 (", g2_name, "): n =", g2_data$nsubs, "subjects\n\n")

# Bundle data for JAGS
jags_data <- list(
    "x_g1" = g1_data$x, "X_g1" = g1_data$X,
    "ntrials_g1" = g1_data$ntrials, "nsubs_g1" = g1_data$nsubs,
    "x_g2" = g2_data$x, "X_g2" = g2_data$X,
    "ntrials_g2" = g2_data$ntrials, "nsubs_g2" = g2_data$nsubs
)

# Parameters to monitor
# mu = grand mean across both groups
# alpha = difference (Group2 - Group1)
params <- c(
    "mu_a_rew", "alpha_a_rew",
    "mu_a_pun", "alpha_a_pun",
    "mu_K", "alpha_K",
    "mu_theta", "alpha_theta",
    "mu_omega_f", "alpha_omega_f",
    "mu_omega_p", "alpha_omega_p"
)

output_dir <- "outputs/group_comparison"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
model_file <- "models/orl_compare.txt"

cat("Fitting joint model...\n")
cat("This may take a while.\n\n")

start_time <- Sys.time()

fit <- jags.parallel(
    data = jags_data,
    inits = NULL,
    parameters.to.save = params,
    model.file = model_file,
    n.chains = 3,
    n.iter = 10000,
    n.burnin = 2000,
    n.thin = 2
)

end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

cat("Fitting complete. Duration:", round(duration, 1), "minutes\n")

# Save results
save_name <- paste0("compare_orl_", g1_name, "_vs_", g2_name, ".rds")
save_path <- file.path(output_dir, save_name)
saveRDS(fit, file = save_path)
cat("Results saved to:", save_path, "\n")

# Quick summary of alpha parameters
cat("\nQuick look at group differences (alpha):\n")
alpha_summary <- fit$BUGSoutput$summary[grepl("^alpha_", rownames(fit$BUGSoutput$summary)), ]
print(round(alpha_summary[, c("mean", "sd", "2.5%", "97.5%", "Rhat")], 3))

cat("\nRun analyse_alpha.R for a full HDI analysis.\n")
