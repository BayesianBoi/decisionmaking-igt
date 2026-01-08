#!/usr/bin/env Rscript
# CHECK CONVERGENCE
#
# Reports MCMC diagnostics for a fitted model. Two main things to look at.
#
# Rhat (Gelman-Rubin statistic) tells us whether the chains agree. Values
# above 1.1 mean they explored different regions of the posterior. Below
# 1.05 is ideal.
#
# n.eff (effective sample size) reflects how many independent samples we
# have after accounting for autocorrelation. Below 100 is a worry since
# it means we lack precision.
#
# Run from project root with model name and group as arguments.
# Example call from terminal
#   Rscript scripts/diagnostics/check_convergence.R orl HC

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage requires model and group. For example orl HC")
}

model_name <- args[1]
target_group <- args[2]

cat("=== Convergence Diagnostics", model_name, "-", target_group, "===\n\n")

library(coda)

# Build path to the fitted model
fit_path <- file.path(
    "outputs", model_name,
    paste0(model_name, "_fit_", target_group, ".rds")
)

if (!file.exists(fit_path)) {
    stop(paste("Fit file not found at", fit_path))
}

fit <- readRDS(fit_path)

# The summary table has all the diagnostics we need
summ <- fit$BUGSoutput$summary

# Get parameter names but exclude the p array if it was saved
param_names <- rownames(summ)
non_p_params <- param_names[!grepl("^p\\[", param_names)]

cat("Total parameters", length(param_names), "\n")
cat("Parameters checked (excluding p)", length(non_p_params), "\n\n")

# Work only with non-p parameters
summ_filtered <- summ[non_p_params, , drop = FALSE]

# --- Rhat diagnostics ---
cat("=== RHAT DIAGNOSTICS ===\n")
cat("Target is below 1.1 (below 1.05 is ideal)\n\n")

rhat_vals <- summ_filtered[, "Rhat"]
rhat_vals <- rhat_vals[!is.na(rhat_vals)]

cat("  Min Rhat", round(min(rhat_vals), 3), "\n")
cat("  Max Rhat", round(max(rhat_vals), 3), "\n")
cat("  Mean Rhat", round(mean(rhat_vals), 3), "\n")
cat("  Median Rhat", round(median(rhat_vals), 3), "\n\n")

# Count how many are above thresholds
n_high_rhat <- sum(rhat_vals > 1.1)
n_medium_rhat <- sum(rhat_vals > 1.05 & rhat_vals <= 1.1)

cat("  Parameters with Rhat > 1.1", n_high_rhat, "/", length(rhat_vals), "\n")
cat("  Parameters with 1.05 < Rhat <= 1.1", n_medium_rhat, "/", length(rhat_vals), "\n\n")

# Show problematic parameters if any
if (n_high_rhat > 0) {
    cat("  Problematic parameters (Rhat > 1.1)\n")
    high_rhat <- rhat_vals[rhat_vals > 1.1]
    high_rhat <- sort(high_rhat, decreasing = TRUE)

    n_to_show <- min(10, length(high_rhat))
    for (i in seq_len(n_to_show)) {
        cat("    ", names(high_rhat)[i], " ", round(high_rhat[i], 3), "\n")
    }
    if (length(high_rhat) > 10) {
        cat("    ... and", length(high_rhat) - 10, "more\n")
    }
    cat("\n")
}

# --- Effective sample size ---
cat("=== EFFECTIVE SAMPLE SIZE (n.eff) ===\n")
cat("Target is above 100 (above 400 is ideal)\n\n")

neff_vals <- summ_filtered[, "n.eff"]
neff_vals <- neff_vals[!is.na(neff_vals)]

cat("  Min n.eff", round(min(neff_vals)), "\n")
cat("  Max n.eff", round(max(neff_vals)), "\n")
cat("  Mean n.eff", round(mean(neff_vals)), "\n")
cat("  Median n.eff", round(median(neff_vals)), "\n\n")

n_low_neff <- sum(neff_vals < 100)
n_medium_neff <- sum(neff_vals >= 100 & neff_vals < 400)

cat("  Parameters with n.eff < 100", n_low_neff, "/", length(neff_vals), "\n")
cat("  Parameters with 100 <= n.eff < 400", n_medium_neff, "/", length(neff_vals), "\n\n")

if (n_low_neff > 0) {
    cat("  Low n.eff parameters (< 100)\n")
    low_neff <- neff_vals[neff_vals < 100]
    low_neff <- sort(low_neff)

    n_to_show <- min(10, length(low_neff))
    for (i in seq_len(n_to_show)) {
        cat("    ", names(low_neff)[i], " ", round(low_neff[i]), "\n")
    }
    if (length(low_neff) > 10) {
        cat("    ... and", length(low_neff) - 10, "more\n")
    }
    cat("\n")
}

# --- Group-level parameters ---
cat("=== GROUP-LEVEL PARAMETERS ===\n")

# These are the ones we care most about
group_params <- non_p_params[grepl("^(mu_|lambda_)", non_p_params)]

if (length(group_params) > 0) {
    group_summ <- summ_filtered[group_params, c("mean", "sd", "2.5%", "97.5%", "Rhat", "n.eff"), drop = FALSE]
    print(round(group_summ, 3))
}

# --- Overall verdict ---
cat("\n=== CONVERGENCE VERDICT ===\n")

if (n_high_rhat == 0 && n_low_neff == 0) {
    cat("EXCELLENT. All parameters converged well\n")
} else if (n_high_rhat == 0 && n_low_neff <= 5) {
    cat("GOOD. Minor n.eff issues but Rhat is fine\n")
} else if (n_high_rhat <= 5) {
    cat("ACCEPTABLE. A few parameters may need more iterations\n")
} else {
    cat("POOR. Consider running longer chains or adjusting priors\n")
}

# Suggestions if there are problems
if (n_high_rhat > 0 || n_low_neff > 10) {
    cat("\nSuggestions\n")
    cat("  - Try increasing n.iter (e.g. 20000 instead of 10000)\n")
    cat("  - Increase n.burnin proportionally\n")
    cat("  - Check for prior-likelihood conflicts\n")
}

cat("\n=== DONE ===\n")
