# Plot Group Posteriors (Violin Plots)
# Generates violin plots comparing parameter posteriors across groups
# Usage: Rscript analysis/2_plotting/plot_group_posteriors.R

library(ggplot2)
library(gridExtra)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

output_dir <- "plots/group_comparison"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Colors
group_colors <- c("HC" = "#4DAF4A", "Amph" = "#E41A1C", "Hero" = "#377EB8")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

extract_mu <- function(fit, param) {
    # Handles both R2jags output structures
    if ("sims.list" %in% names(fit$BUGSoutput)) {
        samples <- fit$BUGSoutput$sims.list[[param]]
    } else {
        samples <- fit$BUGSoutput$sims.matrix[, param]
    }

    if (is.matrix(samples)) samples <- samples[, 1]
    return(as.numeric(samples))
}

plot_violins <- function(data, title, scales = "free_y", ncol = 3) {
    ggplot(data, aes(x = group, y = value, fill = group)) +
        geom_violin(alpha = 0.6, trim = TRUE) +
        geom_boxplot(width = 0.2, alpha = 0.8, outlier.size = 0.5) +
        facet_wrap(~param, scales = scales, ncol = ncol) +
        scale_fill_manual(values = group_colors) +
        labs(
            title = title,
            x = "Group",
            y = "Posterior Value"
        ) +
        theme_minimal(base_size = 11) +
        theme(
            legend.position = "none",
            strip.text = element_text(face = "bold"),
            panel.grid.minor = element_blank()
        )
}

# ==============================================================================
# ORL PLOTTING
# ==============================================================================

cat("Processing ORL...\n")

orl_files <- list(
    HC = "outputs/orl/orl_fit_HC.rds",
    Amph = "outputs/orl/orl_fit_Amph.rds",
    Hero = "outputs/orl/orl_fit_Hero.rds"
)

# Check if files exist
if (all(file.exists(unlist(orl_files)))) {
    orl_params <- c("mu_a_rew", "mu_a_pun", "mu_K", "mu_omega_f", "mu_omega_p")
    orl_labels <- c("Reward Learning (A+)", "Punishment Learning (A-)", "Decay (K)", "Frequency Weight (omega_f)", "Perseverance (omega_p)")

    all_posteriors <- data.frame()

    for (grp in names(orl_files)) {
        fit <- readRDS(orl_files[[grp]])

        for (i in seq_along(orl_params)) {
            vals <- extract_mu(fit, orl_params[i])

            # Handle potentially missing params (e.g. if naming changed)
            if (!is.null(vals)) {
                df <- data.frame(
                    value = vals,
                    group = grp,
                    param = orl_labels[i]
                )
                all_posteriors <- rbind(all_posteriors, df)
            }
        }
    }

    if (nrow(all_posteriors) > 0) {
        all_posteriors$group <- factor(all_posteriors$group, levels = c("HC", "Amph", "Hero"))

        p_orl <- plot_violins(all_posteriors, "ORL Model: Group-Level Parameter Posteriors", ncol = 3)

        ggsave(file.path(output_dir, "orl_group_posteriors.png"), p_orl, width = 12, height = 8, dpi = 300)
        ggsave(file.path(output_dir, "orl_group_posteriors.pdf"), p_orl, width = 12, height = 8)
        cat("ORL plots saved.\n")
    }
} else {
    cat("Skipping ORL (files missing)\n")
}

# ==============================================================================
# EEF PLOTTING
# ==============================================================================

cat("Processing EEF...\n")

eef_files <- list(
    HC = "outputs/eef/eef_fit_HC.rds",
    Amph = "outputs/eef/eef_fit_Amph.rds",
    Hero = "outputs/eef/eef_fit_Hero.rds"
)

if (all(file.exists(unlist(eef_files)))) {
    eef_params <- c("mu_theta", "mu_lambda", "mu_phi", "mu_cons")
    eef_labels <- c("Outcome Sensitivity (theta)", "Forgetting Rate (lambda)", "Exploration Bonus (phi)", "Consistency (c)")

    eef_posteriors <- data.frame()

    for (grp in names(eef_files)) {
        fit <- readRDS(eef_files[[grp]])

        for (i in seq_along(eef_params)) {
            vals <- extract_mu(fit, eef_params[i])

            if (!is.null(vals)) {
                df <- data.frame(
                    value = vals,
                    group = grp,
                    param = eef_labels[i]
                )
                eef_posteriors <- rbind(eef_posteriors, df)
            }
        }
    }

    if (nrow(eef_posteriors) > 0) {
        eef_posteriors$group <- factor(eef_posteriors$group, levels = c("HC", "Amph", "Hero"))

        p_eef <- plot_violins(eef_posteriors, "EEF Model: Group-Level Parameter Posteriors", ncol = 2)

        ggsave(file.path(output_dir, "eef_group_posteriors.png"), p_eef, width = 10, height = 8, dpi = 300)
        ggsave(file.path(output_dir, "eef_group_posteriors.pdf"), p_eef, width = 10, height = 8)
        cat("EEF plots saved.\n")
    }
} else {
    cat("Skipping EEF (files missing)\n")
}

# ==============================================================================
# PVL-DELTA PLOTTING
# ==============================================================================

cat("Processing PVL-Delta...\n")

pvl_files <- list(
    HC = "outputs/pvl_delta/pvl_delta_fit_HC.rds",
    Amph = "outputs/pvl_delta/pvl_delta_fit_Amph.rds",
    Hero = "outputs/pvl_delta/pvl_delta_fit_Hero.rds"
)

if (all(file.exists(unlist(pvl_files)))) {
    pvl_params <- c("mu_A", "mu_a", "mu_w", "mu_theta")
    pvl_labels <- c(
        "Outcome Sensitivity (A)",
        "Learning Rate (a)",
        "Loss Aversion (w)",
        "Response Consistency (theta)"
    )

    pvl_posteriors <- data.frame()

    for (grp in names(pvl_files)) {
        fit <- readRDS(pvl_files[[grp]])

        for (i in seq_along(pvl_params)) {
            vals <- extract_mu(fit, pvl_params[i])

            if (!is.null(vals)) {
                df <- data.frame(
                    value = vals,
                    group = grp,
                    param = pvl_labels[i]
                )
                pvl_posteriors <- rbind(pvl_posteriors, df)
            }
        }
    }

    if (nrow(pvl_posteriors) > 0) {
        pvl_posteriors$group <- factor(pvl_posteriors$group, levels = c("HC", "Amph", "Hero"))

        p_pvl <- plot_violins(pvl_posteriors, "PVL-Delta Model: Group-Level Parameter Posteriors", ncol = 2)

        ggsave(file.path(output_dir, "pvl_delta_group_posteriors.png"), p_pvl, width = 10, height = 8, dpi = 300)
        ggsave(file.path(output_dir, "pvl_delta_group_posteriors.pdf"), p_pvl, width = 10, height = 8)
        cat("PVL-Delta plots saved.\n")
    }
} else {
    cat("Skipping PVL-Delta (files missing)\n")
}

cat("Done.\n")
