#!/usr/bin/env Rscript
#
# plot_trace.R
# Generates combined trace plots for MCMC convergence diagnostics.
# Each figure shows all three participant groups as rows with group-level
# parameters as columns. Useful for visual inspection of chain mixing.
#
# Run with: Rscript scripts/plotting/plot_trace.R <model>
# Example: Rscript scripts/plotting/plot_trace.R eef
#

library(ggplot2)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript plot_trace.R <model>\nModel: eef, orl, pvl_delta")
}

model <- tolower(args[1])
groups <- c("HC", "Amph", "Hero")
group_labels <- c("HC" = "Healthy Controls", "Amph" = "Amphetamine", "Hero" = "Heroin")

output_dir <- "figures/paper"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# parameter names and display labels for each model
if (model == "eef") {
    mu_params <- c("mu_theta", "mu_lambda", "mu_phi", "mu_cons")
    labels <- c(
        "mu_theta" = expression(theta),
        "mu_lambda" = expression(lambda),
        "mu_phi" = expression(phi),
        "mu_cons" = expression(beta)
    )
} else if (model == "orl") {
    mu_params <- c("mu_a_rew", "mu_a_pun", "mu_K", "mu_omega_f", "mu_omega_p")
    labels <- c(
        "mu_a_rew" = expression(A[rew]),
        "mu_a_pun" = expression(A[pun]),
        "mu_K" = expression(K),
        "mu_omega_f" = expression(omega[F]),
        "mu_omega_p" = expression(omega[P])
    )
} else if (model == "pvl_delta") {
    mu_params <- c("mu_A", "mu_a", "mu_w", "mu_theta")
    labels <- c(
        "mu_A" = expression(A),
        "mu_a" = expression(a),
        "mu_w" = expression(w),
        "mu_theta" = expression(theta)
    )
} else {
    stop("Unknown model: ", model)
}

# single trace plot for one parameter
plot_trace <- function(chains, param, label) {
    n_iter <- dim(chains)[1]
    n_chains <- dim(chains)[2]

    df <- data.frame()
    for (c in 1:n_chains) {
        chain_df <- data.frame(
            iteration = 1:n_iter,
            value = chains[, c, param],
            chain = factor(c)
        )
        df <- rbind(df, chain_df)
    }

    ggplot(df, aes(x = iteration, y = value, color = chain)) +
        geom_line(alpha = 0.6, linewidth = 0.2) +
        labs(title = label, x = NULL, y = NULL) +
        scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
        theme_minimal(base_size = 8) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
            legend.position = "none",
            axis.text = element_text(size = 6),
            plot.margin = margin(2, 4, 2, 4)
        )
}

# build grid with rows = groups and columns = parameters
all_plots <- list()

for (group in groups) {
    fit_path <- file.path("data/processed/fits", model, paste0(model, "_fit_", group, ".rds"))
    if (!file.exists(fit_path)) {
        stop("Fit file not found: ", fit_path)
    }

    cat("Loading:", fit_path, "\n")
    fit <- readRDS(fit_path)
    chains <- fit$BUGSoutput$sims.array

    for (param in mu_params) {
        label <- labels[param]
        p <- plot_trace(chains, param, label)
        all_plots[[paste0(group, "_", param)]] <- p
    }
}

# arrange into final figure
n_params <- length(mu_params)

row_plots <- list()
for (i in seq_along(groups)) {
    group <- groups[i]
    group_plots <- list()
    for (j in seq_along(mu_params)) {
        param <- mu_params[j]
        group_plots[[j]] <- all_plots[[paste0(group, "_", param)]]
    }

    row <- ggarrange(plotlist = group_plots, ncol = n_params, nrow = 1)
    row <- annotate_figure(row, left = text_grob(group_labels[group], rot = 90, size = 10, face = "bold"))
    row_plots[[i]] <- row
}

combined <- ggarrange(plotlist = row_plots, ncol = 1, nrow = 3)

model_title <- switch(model,
    "eef" = "EEF Model",
    "orl" = "ORL Model",
    "pvl_delta" = "PVL-Delta Model"
)
combined <- annotate_figure(combined, top = text_grob(model_title, face = "bold", size = 14))

out_file <- file.path(output_dir, paste0("trace_", model, ".png"))
ggsave(out_file, combined, width = 12, height = 7, dpi = 600)
cat("Saved:", out_file, "\n")

cat("Done!\n")
