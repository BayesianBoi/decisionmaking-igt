#!/usr/bin/env Rscript
# compute_dic.R
# -------------
# Compares models using DIC (Deviance Information Criterion).
#
# DIC is computed automatically by JAGS during fitting and stored in the
# fit object. It balances goodness of fit against model complexity:
#   DIC = D_bar + pD
#   D_bar = mean deviance (how well the model fits)
#   pD = effective number of parameters (complexity penalty)
#
# Lower DIC = better model. A difference of 5-10 is meaningful.
#
# Usage: Rscript compute_dic.R <group>
# Example: Rscript compute_dic.R HC

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("Usage: Rscript compute_dic.R <group>")
}

group <- args[1]

valid_groups <- c("HC", "Amph", "Hero")
if (!group %in% valid_groups) {
    stop("Invalid group. Use: HC, Amph, or Hero")
}

cat("\n========================================\n")
cat("DIC Model Comparison:", group, "\n")
cat("========================================\n\n")

# Models to compare
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

    if (!file.exists(fit_file)) {
        cat("Skipping", model_name, "- file not found\n")
        next
    }

    fit <- readRDS(fit_file)

    # DIC is stored in the BUGSoutput
    dic <- fit$BUGSoutput$DIC
    pd <- fit$BUGSoutput$pD

    # D_bar = DIC - pD
    dbar <- dic - pd

    results <- rbind(results, data.frame(
        model = model_name,
        group = group,
        DIC = dic,
        pD = pd,
        Dbar = dbar,
        stringsAsFactors = FALSE
    ))

    cat(model_name, ":\n")
    cat("  DIC:", round(dic, 2), "\n")
    cat("  pD (effective params):", round(pd, 2), "\n")
    cat("  D_bar (mean deviance):", round(dbar, 2), "\n\n")
}

if (nrow(results) < 2) {
    cat("Need at least 2 models to compare.\n")
    quit()
}

# Sort by DIC (lower is better)
results <- results[order(results$DIC), ]

cat("========================================\n")
cat("Model Ranking (lower DIC = better fit)\n")
cat("========================================\n\n")

for (i in seq_len(nrow(results))) {
    row <- results[i, ]
    if (i == 1) {
        cat("  1.", row$model, "- DIC:", round(row$DIC, 2), "(BEST)\n")
    } else {
        delta <- row$DIC - results$DIC[1]
        cat(
            " ", paste0(i, "."), row$model, "- DIC:", round(row$DIC, 2),
            "(+", round(delta, 2), ")\n"
        )
    }
}

# Interpretation
best_model <- results$model[1]
if (nrow(results) >= 2) {
    delta_dic <- results$DIC[2] - results$DIC[1]

    cat("\nDelta DIC (2nd - 1st):", round(delta_dic, 2), "\n")

    if (delta_dic < 2) {
        cat("Interpretation: Models are essentially equivalent.\n")
    } else if (delta_dic < 5) {
        cat("Interpretation: Weak evidence for", best_model, "\n")
    } else if (delta_dic < 10) {
        cat("Interpretation: Moderate evidence for", best_model, "\n")
    } else {
        cat("Interpretation: Strong evidence for", best_model, "\n")
    }
}

# Save results
output_file <- file.path("outputs", paste0("dic_comparison_", group, ".csv"))
write.csv(results, output_file, row.names = FALSE)
cat("\nResults saved:", output_file, "\n")
