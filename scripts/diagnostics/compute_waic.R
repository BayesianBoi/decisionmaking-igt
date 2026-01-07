#!/usr/bin/env Rscript
# compute_waic.R
# --------------
# Computes WAIC (Widely Applicable Information Criterion) for model comparison.
#
# WAIC is a Bayesian measure of out-of-sample predictive accuracy. Lower is better.
# It uses the pointwise log-likelihood computed from the posterior predictive
# distribution. The formula is: WAIC = -2 * (LPPD - p_waic)
#   LPPD = log pointwise predictive density (how well the model predicts each trial)
#   p_waic = effective number of parameters (complexity penalty)
#
# This script requires the p parameter (trial-by-trial choice probabilities)
# to be saved during model fitting.
#
# Usage: Rscript compute_waic.R <group>
# Example: Rscript compute_waic.R HC

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("Usage: Rscript compute_waic.R <group>")
}

group <- args[1]

groups_map <- c(
    "HC" = "Ahn2014_HC",
    "Amph" = "Ahn2014_Amph",
    "Hero" = "Ahn2014_Hero"
)

if (!group %in% names(groups_map)) {
    stop("Invalid group. Use: HC, Amph, or Hero")
}

cat("\n========================================\n")
cat("WAIC Computation:", group, "\n")
cat("========================================\n\n")

# Load behavioural data
source("utils/load_data.R")
all_data <- load_all_igt_data()
group_data <- all_data[all_data$study == groups_map[[group]], ]

subIDs <- unique(group_data$subj)
n_subs <- length(subIDs)

# Build choice matrix to match against model predictions
ntrials_max <- 100
x_all <- array(0, c(n_subs, ntrials_max))
ntrials_all <- array(0, c(n_subs))

for (s in seq_len(n_subs)) {
    subj_data <- group_data[group_data$subj == subIDs[s], ]
    ntrials_all[s] <- nrow(subj_data)

    x_sub <- subj_data$choice
    length(x_sub) <- ntrials_max
    x_all[s, ] <- x_sub
}

# Function to compute WAIC from the p parameter
compute_waic_from_p <- function(p_array, x_obs, ntrials) {
    # p_array has shape [samples, subjects, trials, decks]
    # x_obs has shape [subjects, trials]
    # ntrials is a vector of trial counts per subject

    n_samples <- dim(p_array)[1]
    n_subs <- dim(p_array)[2]

    # We need the log-likelihood for each observation across all posterior samples
    ll_matrix <- list()
    obs_idx <- 1

    for (s in seq_len(n_subs)) {
        # Models start at t=2 because t=1 has no prior to condition on
        for (t in 2:ntrials[s]) {
            choice <- x_obs[s, t]

            # Skip if missing or invalid
            if (!is.na(choice) && choice >= 1 && choice <= 4) {
                # Log probability of the observed choice, across all samples
                # We add a small constant to avoid log(0)
                ll_samples <- log(p_array[, s, t, choice] + 1e-10)
                ll_matrix[[obs_idx]] <- ll_samples
                obs_idx <- obs_idx + 1
            }
        }
    }

    # Stack into a matrix: rows = samples, columns = observations
    ll_mat <- do.call(cbind, ll_matrix)

    # LPPD: for each observation, average likelihood across samples, then take log
    lppd <- sum(log(colMeans(exp(ll_mat))))

    # Effective parameters: variance of log-likelihood across samples
    p_waic <- sum(apply(ll_mat, 2, var))

    # WAIC formula
    waic <- -2 * (lppd - p_waic)

    return(list(
        waic = waic,
        lppd = lppd,
        p_waic = p_waic,
        n_obs = ncol(ll_mat)
    ))
}

# Results table
results <- data.frame(
    model = character(),
    group = character(),
    waic = numeric(),
    lppd = numeric(),
    p_waic = numeric(),
    n_obs = numeric(),
    stringsAsFactors = FALSE
)

# Compute WAIC for each model
for (model_name in c("orl", "eef")) {
    fit_file <- file.path("outputs", model_name, paste0(model_name, "_fit_", group, ".rds"))

    if (!file.exists(fit_file)) {
        cat("Skipping", model_name, "- file not found\n")
        next
    }

    cat("Processing", toupper(model_name), "...\n")
    fit <- readRDS(fit_file)

    # Check if p was saved
    if (!"p" %in% names(fit$BUGSoutput$sims.list)) {
        cat("  The 'p' parameter was not saved during fitting.\n")
        cat("  Re-run with 'p' in parameters.to.save to compute WAIC.\n")
        next
    }

    p_array <- fit$BUGSoutput$sims.list$p

    waic_result <- compute_waic_from_p(p_array, x_all, ntrials_all)

    results <- rbind(results, data.frame(
        model = model_name,
        group = group,
        waic = waic_result$waic,
        lppd = waic_result$lppd,
        p_waic = waic_result$p_waic,
        n_obs = waic_result$n_obs,
        stringsAsFactors = FALSE
    ))

    cat("  WAIC:", round(waic_result$waic, 2), "\n")
    cat("  LPPD:", round(waic_result$lppd, 2), "\n")
    cat("  p_waic:", round(waic_result$p_waic, 2), "\n\n")
}

# Compare models if we have both
if (nrow(results) >= 2) {
    cat("\n========================================\n")
    cat("Model Comparison:", group, "\n")
    cat("========================================\n")

    # Sort by WAIC (lower is better)
    results <- results[order(results$waic), ]
    print(results)

    delta_waic <- diff(results$waic)
    cat("\nDelta WAIC:", round(delta_waic, 2), "\n")
    cat("Preferred model:", results$model[1], "\n")

    # Rough interpretation guidelines
    if (abs(delta_waic) < 5) {
        cat("Interpretation: Models are roughly equivalent.\n")
    } else if (abs(delta_waic) < 10) {
        cat("Interpretation: Some evidence for preferred model.\n")
    } else {
        cat("Interpretation: Strong evidence for preferred model.\n")
    }
} else {
    cat("Could not compare models (need both ORL and EEF fits).\n")
}

# Save results
output_file <- file.path("outputs", paste0("waic_comparison_", group, ".csv"))
write.csv(results, output_file, row.names = FALSE)
cat("\nResults saved:", output_file, "\n")
