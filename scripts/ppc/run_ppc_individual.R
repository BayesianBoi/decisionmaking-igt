#!/usr/bin/env Rscript
# ==============================================================================
# INDIVIDUAL PPC - All Models
# ==============================================================================
# Fits model INDIVIDUALLY for each subject (non-hierarchical)
# Then computes PPC to compare with hierarchical results
# Usage: Rscript scripts/ppc/run_ppc_individual.R <model> <group>
# Example: Rscript scripts/ppc/run_ppc_individual.R orl HC
# ==============================================================================

# 1. Setup
# ------------------------------------------------------------------------------
required_packages <- c("R2jags", "parallel", "ggplot2")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")
invisible(lapply(required_packages, library, character.only = TRUE))

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript run_ppc_individual.R <model> <group>\nModel: orl, pvl_delta, eef\nGroup: HC, Amph, Hero")
}

model_name <- args[1]
target_group <- args[2]

# Validate model
valid_models <- c("orl", "pvl_delta", "eef")
if (!model_name %in% valid_models) {
    stop("Invalid model. Must be one of: ", paste(valid_models, collapse = ", "))
}

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
cat("INDIVIDUAL PPC:", toupper(model_name), "-", target_group, "\n")
cat("================================================\n")

# 2. Load Data
# ------------------------------------------------------------------------------
source("utils/load_data.R")
all_data <- load_all_igt_data()
raw_data <- all_data[all_data$study == study_label, ]

subIDs <- unique(raw_data$subj)
nsubs <- length(subIDs)

cat("Found", nsubs, "subjects\n")

# 3. MPD Function
# ------------------------------------------------------------------------------
MPD <- function(x) {
    x <- as.numeric(x)
    if (length(x) < 2) {
        return(NA)
    }
    dens <- density(x)
    dens$x[which.max(dens$y)]
}

# 4. Model-specific settings
# ------------------------------------------------------------------------------
model_file <- paste0("models/", model_name, "_individual.txt")

if (!file.exists(model_file)) {
    stop("Model file not found: ", model_file)
}

# 5. Fit Each Subject Individually and Compute PPC
# ------------------------------------------------------------------------------
pred_success <- numeric(nsubs)

start_time <- Sys.time()

for (s in 1:nsubs) {
    subj_df <- raw_data[raw_data$subj == subIDs[s], ]
    ntrials <- nrow(subj_df)

    # Prepare data for single subject
    x <- subj_df$choice
    X <- (subj_df$gain + subj_df$loss) / 100 # Scale outcomes

    jags_data <- list(
        "x" = x,
        "X" = X,
        "ntrials" = ntrials
    )

    params <- c("p")

    # Fit this subject individually
    tryCatch(
        {
            fit <- jags.parallel(
                data = jags_data,
                inits = NULL,
                parameters.to.save = params,
                model.file = model_file,
                n.chains = 4,
                n.iter = 8000,
                n.burnin = 1600,
                n.thin = 1,
                n.cluster = 4
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

            pred_success[s] <- sum(x_predict[2:ntrials] == x[2:ntrials], na.rm = TRUE) / (ntrials - 1)
        },
        error = function(e) {
            cat(sprintf("  Subject %d: ERROR - %s\n", s, e$message))
            pred_success[s] <<- NA
        }
    )

    cat(sprintf("Subject %d/%d: accuracy = %.3f\n", s, nsubs, pred_success[s]))
}

end_time <- Sys.time()
cat("\nTotal time:", difftime(end_time, start_time, units = "mins"), "minutes\n")

# 6. Summary
# ------------------------------------------------------------------------------
mean_acc <- mean(pred_success, na.rm = TRUE)
sd_acc <- sd(pred_success, na.rm = TRUE)

cat("\n=== INDIVIDUAL PPC RESULTS ===\n")
cat(sprintf("Model: %s | Group: %s\n", toupper(model_name), target_group))
cat(sprintf("Mean accuracy: %.3f (SD: %.3f)\n", mean_acc, sd_acc))
cat(sprintf("Range: %.3f - %.3f\n", min(pred_success, na.rm = TRUE), max(pred_success, na.rm = TRUE)))

# 7. Save Results
# ------------------------------------------------------------------------------
output_dir <- "outputs/ppc"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pred_df <- data.frame(
    subject = 1:nsubs,
    accuracy = pred_success,
    mean_acc = mean_acc,
    sd_acc = sd_acc
)

results_path <- file.path(output_dir, paste0("ppc_individual_", model_name, "_", target_group, ".rds"))
saveRDS(list(accuracy = pred_success, mean = mean_acc, sd = sd_acc, df = pred_df), results_path)
cat("Results saved to:", results_path, "\n")

# 8. Plot
# ------------------------------------------------------------------------------
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
        title = paste("Individual PPC:", toupper(model_name), "-", target_group),
        subtitle = sprintf("Mean Accuracy: %.3f (SD: %.3f)", mean_acc, sd_acc),
        x = "Subject",
        y = "Proportion Correctly Predicted"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
    )

plot_path <- file.path(output_dir, paste0("ppc_individual_", model_name, "_", target_group, ".png"))
ggsave(plot_path, plot = p, width = 12, height = 6, dpi = 300)
cat("Plot saved to:", plot_path, "\n")
