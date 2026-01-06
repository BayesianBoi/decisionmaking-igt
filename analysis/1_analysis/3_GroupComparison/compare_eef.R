#!/usr/bin/env Rscript
# ==============================================================================
# COMPARE EEF: Joint Hierarchical Model
# ==============================================================================
# Usage: Rscript compare_eef.R [Group1] [Group2]

if (!require("pacman")) install.packages("pacman")
pacman::p_load(R2jags, parallel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Must specify TWO groups.")
}

g1_name <- args[1]
g2_name <- args[2]

groups_map <- c(
    "HC" = "Ahn2014_HC",
    "Amph" = "Ahn2014_Amph",
    "Hero" = "Ahn2014_Hero"
)

if (!g1_name %in% names(groups_map) || !g2_name %in% names(groups_map)) {
    stop("Invalid groups.")
}

g1_label <- groups_map[[g1_name]]
g2_label <- groups_map[[g2_name]]

cat("JOINT COMPARISON: EEF\n")
cat(g1_name, "vs", g2_name, "\n")

source("analysis/utils/load_data.R")
all_data <- load_all_igt_data()

prep_group_data <- function(data, study_label) {
    raw <- data[data$study == study_label, ]
    subIDs <- unique(raw$subj)
    nsubs <- length(subIDs)
    ntrials_max <- 100

    ntrials_all <- array(0, c(nsubs))
    x_all <- array(0, c(nsubs, ntrials_max))
    X_all <- array(0, c(nsubs, ntrials_max))

    for (s in 1:nsubs) {
        subj_df <- raw[raw$subj == subIDs[s], ]
        ntrials_all[s] <- nrow(subj_df)

        x_sub <- subj_df$choice
        length(x_sub) <- ntrials_max
        x_all[s, ] <- x_sub

        # RAW Outcomes
        X_raw <- subj_df$gain + subj_df$loss
        # SCALE by /100 (EEF Requirement)
        X_sub <- X_raw / 100

        length(X_sub) <- ntrials_max
        X_all[s, ] <- X_sub
    }

    return(list(x = x_all, X = X_all, ntrials = ntrials_all, nsubs = nsubs))
}

g1_data <- prep_group_data(all_data, g1_label)
g2_data <- prep_group_data(all_data, g2_label)

jags_data <- list(
    "x_g1" = g1_data$x, "X_g1" = g1_data$X, "ntrials_g1" = g1_data$ntrials, "nsubs_g1" = g1_data$nsubs,
    "x_g2" = g2_data$x, "X_g2" = g2_data$X, "ntrials_g2" = g2_data$ntrials, "nsubs_g2" = g2_data$nsubs
)

params <- c(
    "mu_theta", "alpha_theta",
    "mu_lambda", "alpha_lambda",
    "mu_phi", "alpha_phi",
    "mu_cons", "alpha_cons"
)

output_dir <- "analysis/outputs/group_comparison"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
model_file <- "analysis/models/eef_compare.txt"

cat("Fitting...\n")
fit <- jags.parallel(
    data = jags_data,
    inits = NULL,
    parameters.to.save = params,
    model.file = model_file,
    n.chains = 3,
    n.iter = 5000, n.burnin = 1000, n.thin = 1
)

save_name <- paste0("compare_eef_", g1_name, "_vs_", g2_name, ".rds")
saveRDS(fit, file = file.path(output_dir, save_name))
cat("Done.\n")
