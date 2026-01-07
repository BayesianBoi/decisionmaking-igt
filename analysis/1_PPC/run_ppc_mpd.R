# PPC with MPD (Maximum Posterior Density) approach
# Usage: Rscript analysis/1_PPC/run_ppc_mpd.R <model> <group>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript run_ppc_mpd.R <model> <group>\nModel: eef, orl\nGroup: HC, Amph, Hero")
}

model_name <- args[1]
target_group <- args[2]

cat("=== PPC (MPD) for", model_name, "-", target_group, "===\n\n")

library(ggplot2)

MPD <- function(x) {
    x <- as.numeric(x)
    if (length(x) < 2) {
        return(NA)
    }
    dens <- density(x)
    dens$x[which.max(dens$y)]
}

source("analysis/utils/load_data.R")

fit_path <- file.path(
    "analysis/outputs", model_name,
    paste0(model_name, "_fit_", target_group, ".rds")
)

if (!file.exists(fit_path)) {
    stop("Fit file not found: ", fit_path)
}

cat("Loading fit from:", fit_path, "\n")
fit <- readRDS(fit_path)

p_array <- fit$BUGSoutput$sims.list$p

if (is.null(p_array)) {
    stop("No choice probabilities (p) found in fit. Re-run fitting with p monitored.")
}

cat("p array dimensions:", dim(p_array), "\n")
n_iterations <- dim(p_array)[1]
n_subjects <- dim(p_array)[2]
n_trials <- dim(p_array)[3]
n_decks <- dim(p_array)[4]

cat("Iterations:", n_iterations, "\n")
cat("Subjects:", n_subjects, "\n")
cat("Trials:", n_trials, "\n")

all_data <- load_all_igt_data()
all_data$group <- gsub("Ahn2014_", "", all_data$study)
group_data <- all_data[all_data$group == target_group, ]

subj_list <- unique(group_data$subj_unique)
nsubs <- length(subj_list)

x_all <- matrix(NA, nrow = nsubs, ncol = 100)
ntrials_all <- numeric(nsubs)

for (s in 1:nsubs) {
    subj_data <- group_data[group_data$subj_unique == subj_list[s], ]
    subj_data <- subj_data[order(subj_data$trial), ]
    n_trials_s <- nrow(subj_data)
    ntrials_all[s] <- n_trials_s
    x_all[s, 1:n_trials_s] <- subj_data$choice
}

pred_success <- numeric(nsubs)

for (s in 1:nsubs) {
    ntrials <- ntrials_all[s]
    x <- x_all[s, 1:ntrials]

    p_post <- p_array[, s, , ]
    n_posterior_trials <- dim(p_post)[2]

    x_predict <- numeric(ntrials)
    x_predict[1] <- NA

    for (t in 2:min(ntrials, n_posterior_trials)) {
        p_predict <- c(
            MPD(p_post[, t, 1]),
            MPD(p_post[, t, 2]),
            MPD(p_post[, t, 3]),
            MPD(p_post[, t, 4])
        )
        p_predict <- p_predict / sum(p_predict)
        x_predict[t] <- which.max(p_predict)
    }

    n_compare <- min(ntrials, n_posterior_trials)
    pred_success[s] <- sum(x_predict[2:n_compare] == x[2:n_compare], na.rm = TRUE) / (n_compare - 1)

    if (s %% 10 == 0) {
        cat(sprintf("Subject %d/%d: accuracy = %.3f\n", s, nsubs, pred_success[s]))
    }
}

mean_acc <- mean(pred_success, na.rm = TRUE)
sd_acc <- sd(pred_success, na.rm = TRUE)

cat("\n=== RESULTS ===\n")
cat(sprintf("Mean accuracy: %.3f (SD: %.3f)\n", mean_acc, sd_acc))
cat(sprintf("Range: %.3f - %.3f\n", min(pred_success), max(pred_success)))

output_dir <- "analysis/outputs/ppc"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pred_df <- data.frame(
    subject = 1:nsubs,
    accuracy = pred_success,
    mean_acc = mean_acc,
    sd_acc = sd_acc
)

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
        title = paste("Posterior Predictive Accuracy:", toupper(model_name), "-", target_group),
        subtitle = sprintf("Mean Accuracy: %.3f (SD: %.3f)", mean_acc, sd_acc),
        x = "Subject",
        y = "Proportion Correctly Predicted"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
    )

plot_path <- file.path(output_dir, paste0("ppc_mpd_", model_name, "_", target_group, ".png"))
ggsave(plot_path, plot = p, width = 12, height = 6, dpi = 300)
cat("\nPlot saved to:", plot_path, "\n")

results_path <- file.path(output_dir, paste0("ppc_mpd_", model_name, "_", target_group, ".rds"))
saveRDS(list(accuracy = pred_success, mean = mean_acc, sd = sd_acc, df = pred_df), results_path)
cat("Results saved to:", results_path, "\n")
