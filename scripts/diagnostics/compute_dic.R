#!/usr/bin/env Rscript
# COMPUTE DIC
#
# Compares models using DIC (Deviance Information Criterion). JAGS computes
# this automatically during fitting so we just need to read it out.
#
# DIC = D_bar + pD
#   D_bar is the mean deviance across posterior samples (goodness of fit)
#   pD is the effective number of parameters (penalty for complexity)
#
# Lower DIC means better. A difference of 5 to 10 is worth noting.
# Differences under 2 are not meaningful.
#
# Run from project root with group as argument.
# Example call from terminal
#   Rscript scripts/diagnostics/compute_dic.R HC

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("Usage requires group argument. For example HC")
}

group <- args[1]

valid_groups <- c("HC", "Amph", "Hero")
if (!group %in% valid_groups) {
    stop("Invalid group. Must be one of HC, Amph, or Hero")
}

cat("\n========================================\n")
cat("DIC Model Comparison", group, "\n")
cat("========================================\n\n")

# Models we want to compare
models <- c("orl", "eef", "pvl_delta")

results <- data.frame(
    model = character(),
    group = character(),
    DIC = numeric(),
    pD = numeric(),
    Dbar = numeric(),
    stringsAsFactors = FALSE
)

for (model_name in models) {
    fit_file <- file.path("outputs", model_name, paste0(model_name, "_fit_", group, ".rds"))

    # Handle alternative naming/path for ORL if needed (e.g. outputs/ORL vs outputs/orl)
    if (!file.exists(fit_file) && model_name == "orl") {
        fit_file <- file.path("outputs", "ORL", paste0("orl_fit_", group, ".rds"))
    }

    if (!file.exists(fit_file)) {
        cat("Skipping", model_name, "- file not found\n")
        next
    }

    tryCatch(
        {
            fit <- readRDS(fit_file)

            # DIC and pD are stored in the BUGSoutput slot
            dic <- fit$BUGSoutput$DIC

            # R2jags sometimes stores pD as pV or pD depending on the version/method
            if ("pD" %in% names(fit$BUGSoutput)) {
                pd <- fit$BUGSoutput$pD
            } else if ("pV" %in% names(fit$BUGSoutput)) {
                pd <- fit$BUGSoutput$pV
            } else {
                pd <- NA
            }

            # D_bar is the fit component before the complexity penalty
            if (!is.na(dic) && !is.na(pd)) {
                dbar <- dic - pd
            } else {
                dbar <- NA
            }

            results <- rbind(results, data.frame(
                model = model_name,
                group = group,
                DIC = ifelse(is.null(dic), NA, dic),
                pD = ifelse(is.null(pd), NA, pd),
                Dbar = ifelse(is.null(dbar), NA, dbar),
                stringsAsFactors = FALSE
            ))

            cat(model_name, "\n")
            cat("  DIC", round(dic, 2), "\n")
            if (!is.na(pd)) cat("  pD (effective params)", round(pd, 2), "\n")
            if (!is.na(dbar)) cat("  D_bar (mean deviance)", round(dbar, 2), "\n\n")
        },
        error = function(e) {
            cat("Error processing", model_name, ":", e$message, "\n")
        }
    )
}

if (nrow(results) < 2) {
    cat("Need at least 2 models to compare. Found:", nrow(results), "\n")
    if (nrow(results) == 1) print(results)
    quit()
}

# Sort by DIC where lower is better
results <- results[order(results$DIC), ]

cat("========================================\n")
cat("Model Ranking (lower DIC = better fit)\n")
cat("========================================\n\n")

for (i in seq_len(nrow(results))) {
    row <- results[i, ]
    if (i == 1) {
        cat("  1.", row$model, "- DIC", round(row$DIC, 2), "(BEST)\n")
    } else {
        delta <- row$DIC - results$DIC[1]
        cat(
            " ", paste0(i, "."), row$model, "- DIC", round(row$DIC, 2),
            "(+", round(delta, 2), ")\n"
        )
    }
}

# Interpretation
best_model <- results$model[1]
if (nrow(results) >= 2) {
    delta_dic <- results$DIC[2] - results$DIC[1]

    cat("\nDelta DIC (2nd - 1st)", round(delta_dic, 2), "\n")

    if (delta_dic < 2) {
        cat("The models are too close to distinguish\n")
    } else if (delta_dic < 5) {
        cat("Weak evidence favouring", best_model, "\n")
    } else if (delta_dic < 10) {
        cat("Moderate evidence favouring", best_model, "\n")
    } else {
        cat("Strong evidence favouring", best_model, "\n")
    }
}

# Save results to disk
output_file <- file.path("outputs", paste0("dic_comparison_", group, ".csv"))
write.csv(results, output_file, row.names = FALSE)
cat("\nResults saved to", output_file, "\n")
