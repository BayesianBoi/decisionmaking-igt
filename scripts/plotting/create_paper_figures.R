#!/usr/bin/env Rscript
#
# create_paper_figures.R
# Generates the main figures for the paper. Includes DIC model comparison
# and group-level posterior distributions for all three models.
#

library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(patchwork)

output_dir <- "figures/paper"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Consistent colours across all plots
group_colors <- c("HC" = "#4DAF4A", "Amph" = "#E41A1C", "Hero" = "#377EB8")

# Grabs the group-level mean parameter samples from a fitted model
extract_mu <- function(fit, param) {
    if ("sims.list" %in% names(fit$BUGSoutput)) {
        samples <- fit$BUGSoutput$sims.list[[param]]
    } else {
        samples <- fit$BUGSoutput$sims.matrix[, param]
    }
    if (is.matrix(samples)) samples <- samples[, 1]
    return(as.numeric(samples))
}

# DIC comparison bar chart
# Shows which model fits best for each group. Lower is better.
cat("Creating DIC comparison figure...\n")

# Updated with actual extracted values
dic_data <- data.frame(
    Group = rep(c("HC", "Amph", "Hero"), each = 3),
    Model = rep(c("ORL", "EEF", "PVL-Delta"), 3),
    DIC = c(
        11718, 11766, 12587, # HC: ORL best
        9129, 9111, 9741, # Amph: EEF best
        10276, 10211, 11135 # Hero: EEF best
    )
)

dic_data$Group <- factor(dic_data$Group, levels = c("HC", "Amph", "Hero"))
dic_data$Model <- factor(dic_data$Model, levels = c("ORL", "EEF", "PVL-Delta"))

dic_data <- dic_data %>%
    group_by(Group) %>%
    mutate(
        Best_DIC = min(DIC),
        Delta_DIC = DIC - Best_DIC
    ) %>%
    ungroup()

p_dic <- ggplot(dic_data, aes(x = Group, y = DIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = round(DIC, 0)),
        position = position_dodge(width = 0.8), vjust = -0.5, size = 3
    ) +
    scale_fill_manual(values = c("ORL" = "#2166AC", "EEF" = "#B2182B", "PVL-Delta" = "#4DAF4A")) +
    labs(
        title = "Model Comparison by DIC",
        x = "Group",
        y = "Deviance Information Criterion (DIC)",
        fill = "Model"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank()
    ) +
    coord_cartesian(ylim = c(8500, 13000))

ggsave(file.path(output_dir, "dic_model_comparison.png"), p_dic, width = 8, height = 6, dpi = 600)
cat("DIC figure saved.\n")

# EEF posteriors
# The forgetting rate (lambda) is our main hypothesis parameter
cat("Creating EEF group posteriors figure...\n")

eef_files <- list(
    HC = "data/processed/fits/eef/eef_fit_HC.rds",
    Amph = "data/processed/fits/eef/eef_fit_Amph.rds",
    Hero = "data/processed/fits/eef/eef_fit_Hero.rds"
)

eef_params <- c("mu_lambda", "mu_theta", "mu_phi", "mu_cons")
eef_labels <- c(
    "mu_lambda" = "Forgetting Rate (lambda)",
    "mu_theta" = "Outcome Sensitivity (theta)",
    "mu_phi" = "Exploration Bonus (phi)",
    "mu_cons" = "Choice Consistency (c)"
)

eef_posteriors <- data.frame()

for (grp in names(eef_files)) {
    fit <- readRDS(eef_files[[grp]])
    for (p in eef_params) {
        vals <- extract_mu(fit, p)
        if (!is.null(vals)) {
            df <- data.frame(
                value = vals,
                group = grp,
                parameter = eef_labels[p]
            )
            eef_posteriors <- rbind(eef_posteriors, df)
        }
    }
}

eef_posteriors$group <- factor(eef_posteriors$group, levels = c("HC", "Amph", "Hero"))
eef_posteriors$parameter <- factor(eef_posteriors$parameter,
    levels = c(
        "Forgetting Rate (lambda)", "Outcome Sensitivity (theta)",
        "Exploration Bonus (phi)", "Choice Consistency (c)"
    )
)

p_eef <- ggplot(eef_posteriors, aes(x = value, fill = group, colour = group)) +
    geom_density(alpha = 0.4, linewidth = 0.8) +
    facet_wrap(~parameter, scales = "free", ncol = 2) +
    scale_fill_manual(values = group_colors) +
    scale_colour_manual(values = group_colors) +
    labs(
        title = "EEF Model Group-Level Parameter Posteriors",
        x = "Parameter Value",
        y = "Density",
        fill = "Group",
        colour = "Group"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")
    )

ggsave(file.path(output_dir, "eef_group_posteriors.png"), p_eef, width = 10, height = 8, dpi = 600)
cat("EEF posteriors saved.\n")

# ORL posteriors
cat("Creating ORL group posteriors figure...\n")

orl_files <- list(
    HC = "data/processed/fits/orl/orl_fit_HC.rds",
    Amph = "data/processed/fits/orl/orl_fit_Amph.rds",
    Hero = "data/processed/fits/orl/orl_fit_Hero.rds"
)

orl_params <- c("mu_a_rew", "mu_a_pun", "mu_K", "mu_omega_f", "mu_omega_p")
orl_labels <- c(
    "mu_a_rew" = "Reward Learning Rate",
    "mu_a_pun" = "Punishment Learning Rate",
    "mu_K" = "Perseverance Decay (K)",
    "mu_omega_f" = "Frequency Weight",
    "mu_omega_p" = "Perseverance Weight"
)

orl_posteriors <- data.frame()

for (grp in names(orl_files)) {
    fit <- readRDS(orl_files[[grp]])
    for (p in orl_params) {
        vals <- extract_mu(fit, p)
        if (!is.null(vals)) {
            df <- data.frame(
                value = vals,
                group = grp,
                parameter = orl_labels[p]
            )
            orl_posteriors <- rbind(orl_posteriors, df)
        }
    }
}

orl_posteriors$group <- factor(orl_posteriors$group, levels = c("HC", "Amph", "Hero"))
orl_posteriors$parameter <- factor(orl_posteriors$parameter,
    levels = c(
        "Reward Learning Rate", "Punishment Learning Rate",
        "Perseverance Decay (K)", "Frequency Weight", "Perseverance Weight"
    )
)

p_orl <- ggplot(orl_posteriors, aes(x = value, fill = group, colour = group)) +
    geom_density(alpha = 0.4, linewidth = 0.8) +
    facet_wrap(~parameter, scales = "free", ncol = 3) +
    scale_fill_manual(values = group_colors) +
    scale_colour_manual(values = group_colors) +
    labs(
        title = "ORL Model Group-Level Parameter Posteriors",
        x = "Parameter Value",
        y = "Density",
        fill = "Group",
        colour = "Group"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")
    )

ggsave(file.path(output_dir, "orl_group_posteriors.png"), p_orl, width = 12, height = 8, dpi = 600)
cat("ORL posteriors saved.\n")

# PVL-Delta posteriors
cat("Creating PVL-Delta group posteriors figure...\n")

pvl_files <- list(
    HC = "data/processed/fits/pvl_delta/pvl_delta_fit_HC.rds",
    Amph = "data/processed/fits/pvl_delta/pvl_delta_fit_Amph.rds",
    Hero = "data/processed/fits/pvl_delta/pvl_delta_fit_Hero.rds"
)

pvl_params <- c("mu_A", "mu_a", "mu_w", "mu_theta")
pvl_labels <- c(
    "mu_A" = "Outcome Sensitivity (A)",
    "mu_a" = "Learning Rate (a)",
    "mu_w" = "Loss Aversion (w)",
    "mu_theta" = "Response Consistency (theta)"
)

pvl_posteriors <- data.frame()

for (grp in names(pvl_files)) {
    fit <- readRDS(pvl_files[[grp]])
    for (p in pvl_params) {
        vals <- extract_mu(fit, p)
        if (!is.null(vals)) {
            df <- data.frame(
                value = vals,
                group = grp,
                parameter = pvl_labels[p]
            )
            pvl_posteriors <- rbind(pvl_posteriors, df)
        }
    }
}

pvl_posteriors$group <- factor(pvl_posteriors$group, levels = c("HC", "Amph", "Hero"))
pvl_posteriors$parameter <- factor(pvl_posteriors$parameter,
    levels = c(
        "Outcome Sensitivity (A)", "Learning Rate (a)",
        "Loss Aversion (w)", "Response Consistency (theta)"
    )
)

p_pvl <- ggplot(pvl_posteriors, aes(x = value, fill = group, colour = group)) +
    geom_density(alpha = 0.4, linewidth = 0.8) +
    facet_wrap(~parameter, scales = "free", ncol = 2) +
    scale_fill_manual(values = group_colors) +
    scale_colour_manual(values = group_colors) +
    labs(
        title = "PVL-Delta Model Group-Level Parameter Posteriors",
        x = "Parameter Value",
        y = "Density",
        fill = "Group",
        colour = "Group"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")
    )

ggsave(file.path(output_dir, "pvl_delta_group_posteriors.png"), p_pvl, width = 10, height = 8, dpi = 600)
cat("PVL-Delta posteriors saved.\n")

cat("\n========================================\n")
cat("All figures saved to:", output_dir, "\n")
cat("========================================\n")
list.files(output_dir)
