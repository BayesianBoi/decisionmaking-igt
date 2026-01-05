# ===========================================================================
# Posterior Predictive Checks (Example2.pdf Style)
# ===========================================================================
#
# Usage: Rscript analysis/2_plotting/plot_ppc.R
#
# Generates per-subject prediction accuracy plots:
#   - Each dot = one subject's prediction accuracy
#   - Solid line = mean across subjects
#   - Pink band = ±1 SD
#   - Dashed line = chance level (0.25)
#
# ===========================================================================

library(ggplot2)
library(gridExtra)
library(coda)

source("analysis/utils/load_data.R")
source("analysis/utils/ppc.R")
source("analysis/utils/prepare_jags_data.R")
source("analysis/utils/prepare_eef_data.R")
source("analysis/utils/plotting_utils.R")

# ===========================================================================
# Configuration
# ===========================================================================

output_dir <- "results/ppc"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

models <- c("pvl_delta", "orl", "eef")
n_sim <- 100 # Posterior samples for simulation

# Colors matching Example2.pdf
point_color <- "gray30"
band_color <- "#F4A8A8" # Light pink/salmon
mean_line_color <- "black"
chance_line_color <- "black"

# ===========================================================================
# Helper: Compute per-subject prediction accuracy
# ===========================================================================

compute_subject_accuracy <- function(observed_choices, simulated_choices, Tsubj) {
    n_subj <- nrow(observed_choices)
    accuracy <- rep(NA, n_subj)

    for (s in 1:n_subj) {
        n_trials <- Tsubj[s]
        obs <- observed_choices[s, 1:n_trials]
        sim <- simulated_choices[s, 1:n_trials]

        # Remove NA
        valid <- !is.na(obs) & !is.na(sim)
        if (sum(valid) > 0) {
            accuracy[s] <- mean(obs[valid] == sim[valid])
        }
    }

    return(accuracy)
}

# ===========================================================================
# Main PPC Loop
# ===========================================================================

message("=== RUNNING POSTERIOR PREDICTIVE CHECKS ===\n")

ppc_results <- list()

for (model in models) {
    message(sprintf("Processing model: %s", model))

    # Load samples and data
    samples_file <- file.path("results", model, "mcmc_samples.rds")
    data_file <- file.path("results", model, "jags_data.rds")

    if (!file.exists(samples_file) || !file.exists(data_file)) {
        warning(sprintf("Fit results not found for %s", model))
        next
    }

    samples <- readRDS(samples_file)
    jags_data <- readRDS(data_file)
    mat <- as.matrix(samples)

    # Sample posterior indices
    n_posterior <- nrow(mat)
    posterior_idx <- sample(1:n_posterior, size = min(n_sim, n_posterior), replace = FALSE)

    # Store per-subject accuracies across simulations
    n_subj <- jags_data$N
    Tsubj <- jags_data$Tsubj
    all_accuracies <- matrix(NA, nrow = length(posterior_idx), ncol = n_subj)

    for (i in seq_along(posterior_idx)) {
        idx <- posterior_idx[i]

        # Extract parameters and simulate
        if (model == "pvl_delta") {
            params <- list(
                a = mat[idx, grep("^a\\[", colnames(mat))],
                A = mat[idx, grep("^A\\[", colnames(mat))],
                theta = mat[idx, grep("^theta\\[", colnames(mat))],
                w = mat[idx, grep("^w\\[", colnames(mat))]
            )
            sim_choices <- simulate_pvl_delta(params, jags_data$outcome, Tsubj)
        } else if (model == "orl") {
            params <- list(
                Arew = mat[idx, grep("^Arew\\[", colnames(mat))],
                Apun = mat[idx, grep("^Apun\\[", colnames(mat))],
                K = mat[idx, grep("^K\\[", colnames(mat))],
                betaF = mat[idx, grep("^betaF\\[", colnames(mat))],
                betaP = mat[idx, grep("^betaP\\[", colnames(mat))]
            )
            sim_choices <- simulate_orl(params, jags_data$outcome, Tsubj)
        } else if (model == "eef") {
            params <- list(
                theta = mat[idx, grep("^theta\\[", colnames(mat))],
                lambda = mat[idx, grep("^lambda_forget\\[", colnames(mat))],
                phi = mat[idx, grep("^phi\\[", colnames(mat))],
                cons = mat[idx, grep("^cons\\[", colnames(mat))]
            )
            w_ini <- if (!is.null(jags_data$w_ini)) jags_data$w_ini else rep(0, 4)
            sim_choices <- simulate_eef(params, jags_data$outcome, Tsubj, w_ini)
        } else {
            next
        }

        # Compute per-subject accuracy
        all_accuracies[i, ] <- compute_subject_accuracy(jags_data$choice, sim_choices, Tsubj)
    }

    # Average across simulations
    mean_accuracy <- colMeans(all_accuracies, na.rm = TRUE)
    overall_mean <- mean(mean_accuracy, na.rm = TRUE)
    overall_sd <- sd(mean_accuracy, na.rm = TRUE)

    ppc_results[[model]] <- list(
        subject_accuracy = mean_accuracy,
        overall_mean = overall_mean,
        overall_sd = overall_sd
    )

    # Create plot (Example2.pdf style)
    df <- data.frame(
        subject = 1:n_subj,
        accuracy = mean_accuracy
    )

    p <- ggplot(df, aes(x = subject, y = accuracy)) +
        # Pink ±1 SD band
        geom_rect(
            aes(
                xmin = -Inf, xmax = Inf,
                ymin = overall_mean - overall_sd,
                ymax = overall_mean + overall_sd
            ),
            fill = band_color, alpha = 0.5
        ) +
        # Mean line
        geom_hline(yintercept = overall_mean, color = mean_line_color, linewidth = 1) +
        # Chance line (dashed)
        geom_hline(yintercept = 0.25, linetype = "dashed", color = chance_line_color) +
        # Per-subject points
        geom_point(color = point_color, size = 1.5) +
        labs(
            title = sprintf("Posterior Predictive Checks: %s", toupper(model)),
            x = "Subject",
            y = "Predicted Success Percentage"
        ) +
        coord_cartesian(ylim = c(0.1, 1.0)) +
        theme_gray(base_size = 11) +
        theme(
            panel.grid.major = element_line(color = "white"),
            panel.grid.minor = element_blank()
        )

    # Save individual plot
    ggsave(
        file.path(output_dir, sprintf("ppc_%s.pdf", model)),
        p,
        width = 8, height = 4
    )

    message(sprintf(
        "  %s: Mean accuracy = %.2f (SD = %.2f)",
        model, overall_mean, overall_sd
    ))
}

# ===========================================================================
# Comparison Plot (all models)
# ===========================================================================

if (length(ppc_results) > 0) {
    message("\nGenerating comparison plot...")

    # Combine data
    all_df <- do.call(rbind, lapply(names(ppc_results), function(m) {
        data.frame(
            model = toupper(m),
            subject = seq_along(ppc_results[[m]]$subject_accuracy),
            accuracy = ppc_results[[m]]$subject_accuracy,
            mean = ppc_results[[m]]$overall_mean,
            sd = ppc_results[[m]]$overall_sd
        )
    }))

    p_comp <- ggplot(all_df, aes(x = subject, y = accuracy)) +
        # Band and mean per model
        geom_rect(
            aes(
                xmin = -Inf, xmax = Inf,
                ymin = mean - sd, ymax = mean + sd
            ),
            fill = band_color, alpha = 0.4
        ) +
        geom_hline(aes(yintercept = mean), color = mean_line_color, linewidth = 0.8) +
        geom_hline(yintercept = 0.25, linetype = "dashed", color = chance_line_color) +
        geom_point(color = point_color, size = 1) +
        facet_wrap(~model, ncol = 3) +
        labs(
            title = "Posterior Predictive Checks: Model Comparison",
            x = "Subject",
            y = "Predicted Success Percentage"
        ) +
        coord_cartesian(ylim = c(0.1, 1.0)) +
        theme_gray(base_size = 10) +
        theme(
            panel.grid.major = element_line(color = "white"),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "grey90"),
            strip.text = element_text(face = "bold")
        )

    ggsave(
        file.path(output_dir, "ppc_comparison.pdf"),
        p_comp,
        width = 12, height = 4
    )

    # Save results
    saveRDS(ppc_results, file.path(output_dir, "ppc_results.rds"))

    message(sprintf("\nPPC complete. Results saved to %s", output_dir))
}
