#!/usr/bin/env Rscript
# ==============================================================================
# COMPARE PVL-DELTA: Joint Hierarchical Model
# ==============================================================================
# Usage: Rscript compare_pvl_delta.R [Group1] [Group2]
# Example: Rscript compare_pvl_delta.R HC Hero

if (!require("pacman")) install.packages("pacman")
pacman::p_load(R2jags, parallel)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Must specify TWO groups to compare. Example: HC Hero")
}

g1_name <- args[1]
g2_name <- args[2]

# Map aliases
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

cat("\n================================================\n")
cat("JOINT COMPARISON: PVL-DELTA\n")
cat("Group 1:", g1_name, "(", g1_label, ")\n")
cat("Group 2:", g2_name, "(", g2_label, ")\n")
cat("================================================\n")

# 1. Load Data
source("analysis/utils/load_data.R")
all_data <- load_all_igt_data()

# 2. Prep Function
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

        # Pad
        x_sub <- subj_df$choice
        length(x_sub) <- ntrials_max
        x_all[s, ] <- x_sub

        # RAW Outcomes -> SCALED
        X_raw <- subj_df$gain + subj_df$loss
        X_sub <- X_raw / 100
        length(X_sub) <- ntrials_max
        X_all[s, ] <- X_sub
    }

    return(list(
        x = x_all,
        X = X_all,
        ntrials = ntrials_all,
        nsubs = nsubs
    ))
}

g1_data <- prep_group_data(all_data, g1_label)
g2_data <- prep_group_data(all_data, g2_label)

# 3. JAGS Prep
jags_data <- list(
    "x_g1" = g1_data$x,
    "X_g1" = g1_data$X,
    "ntrials_g1" = g1_data$ntrials,
    "nsubs_g1" = g1_data$nsubs,
    "x_g2" = g2_data$x,
    "X_g2" = g2_data$X,
    "ntrials_g2" = g2_data$ntrials,
    "nsubs_g2" = g2_data$nsubs
)

params <- c(
    "mu_w", "alpha_w",
    "mu_A", "alpha_A",
    "mu_theta", "alpha_theta",
    "mu_a", "alpha_a"
)

# 4. Fit
output_dir <- "analysis/outputs/group_comparison"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

model_file <- "analysis/models/pvl_delta_compare.txt"

cat("Fitting Joint Model...\n")
start_time <- Sys.time()

fit <- jags.parallel(
    data = jags_data,
    inits = NULL,
    parameters.to.save = params,
    model.file = model_file,
    n.chains = 3,
    n.iter = 5000,
    n.burnin = 1000,
    n.thin = 1
)

end_time <- Sys.time()
cat("Duration:", end_time - start_time, "\n")

save_name <- paste0("compare_pvl_delta_", g1_name, "_vs_", g2_name, ".rds")
save_path <- file.path(output_dir, save_name)
saveRDS(fit, file = save_path)
cat("Saved to:", save_path, "\n")
