# Script to summarize parameter estimates from fitted .rds files
library(dplyr)
library(tidyr)
library(R2jags)

# Function to extract summary stats
get_param_summary <- function(fit_file, model_name, group_name) {
    if (!file.exists(fit_file)) {
        return(NULL)
    }

    fit <- readRDS(fit_file)

    # Check if fit is valid
    if (is.null(fit$BUGSoutput)) {
        return(NULL)
    }

    summary_stats <- fit$BUGSoutput$summary

    # Parameter names depend on model
    if (model_name == "eef") {
        params <- c("mu_theta", "mu_lambda", "mu_phi", "mu_cons")
    } else if (model_name == "orl") {
        params <- c("mu_a_rew", "mu_a_pun", "mu_K", "mu_theta", "mu_omega_f", "mu_omega_p")
    } else if (model_name == "pvl_delta") {
        params <- c("mu_w", "mu_A", "mu_a", "mu_theta")
    }

    results <- data.frame()

    for (p in params) {
        if (p %in% rownames(summary_stats)) {
            mean_val <- summary_stats[p, "mean"]
            lower <- summary_stats[p, "2.5%"]
            upper <- summary_stats[p, "97.5%"]

            results <- rbind(results, data.frame(
                Model = model_name,
                Group = group_name,
                Parameter = p,
                Mean = round(mean_val, 3),
                CI_2.5 = round(lower, 3),
                CI_97.5 = round(upper, 3)
            ))
        }
    }
    return(results)
}

# Process all files
all_results <- data.frame()

# Groups
groups <- c("HC", "Amph", "Hero")

# EEF
for (g in groups) {
    res <- get_param_summary(paste0("analysis/outputs/eef/eef_fit_", g, ".rds"), "eef", g)
    if (!is.null(res)) all_results <- rbind(all_results, res)
}

# ORL
for (g in groups) {
    res <- get_param_summary(paste0("analysis/outputs/orl/orl_fit_", g, ".rds"), "orl", g)
    if (!is.null(res)) all_results <- rbind(all_results, res)
}

# PVL-Delta
for (g in groups) {
    res <- get_param_summary(paste0("analysis/outputs/pvl_delta/pvl_delta_fit_", g, ".rds"), "pvl_delta", g)
    if (!is.null(res)) all_results <- rbind(all_results, res)
}

# Print results nicely
print(all_results)
