# Plotting Utilities for Decision Making Analysis
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, bayestestR, dplyr, tidyr, ggpubr)

#' Plot Parameter Recovery (Scatter Plot)
#'
#' @param true_vals Numeric vector of ground-truth values
#' @param infer_vals Numeric vector of recovered parameters (means)
#' @param param_name String, name of the parameter (e.g., "Learning Rate")
#' @return A ggplot object
plot_recovery <- function(true_vals, infer_vals, param_name) {
    df <- data.frame(True = true_vals, Inferred = infer_vals)

    # Calculate Correlation
    corr <- cor(true_vals, infer_vals)
    rmse <- sqrt(mean((true_vals - infer_vals)^2))

    ggplot(df, aes(x = True, y = Inferred)) +
        geom_point(alpha = 0.6, color = "#2C3E50") +
        geom_abline(intercept = 0, slope = 1, color = "#E74C3C", linetype = "dashed") +
        geom_smooth(method = "lm", color = "#3498DB", fill = "lightblue", alpha = 0.2) +
        labs(
            title = paste0(param_name),
            subtitle = paste0("r = ", round(corr, 2), ", RMSE = ", round(rmse, 2)),
            x = "True Parameter",
            y = "Recovered (Mean)"
        ) +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold"))
}

#' Plot Posterior Difference (Group Comparison)
#'
#' @param alpha_samples Numeric vector of posterior samples for the difference parameter (alpha)
#' @param param_name String, name of parameter
#' @param group1_name Name of Group 1 (Ref)
#' @param group2_name Name of Group 2
#' @return A ggplot object
plot_difference <- function(alpha_samples, param_name, group1_name, group2_name) {
    df <- data.frame(Difference = alpha_samples)

    # HDI and ROPE
    hdi_res <- hdi(alpha_samples, ci = 0.95)
    pd <- p_direction(alpha_samples)

    fill_color <- ifelse(hdi_res$CI_low > 0 | hdi_res$CI_high < 0, "#27AE60", "#95A5A6")

    ggplot(df, aes(x = Difference)) +
        geom_density(fill = fill_color, alpha = 0.5, color = "black") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
        labs(
            title = paste0("Difference: ", param_name),
            subtitle = paste0(
                group2_name, " - ", group1_name, "\n",
                "PD: ", round(pd * 100, 1), "% | 95% HDI: [", round(hdi_res$CI_low, 2), ", ", round(hdi_res$CI_high, 2), "]"
            ),
            x = "Estimated Difference (alpha)"
        ) +
        theme_minimal()
}

#' Combine Recovery Plots
#'
#' @param plot_list List of ggplots
#' @param title Main title
combine_plots <- function(plot_list, title) {
    ggarrange(plotlist = plot_list, common.legend = TRUE, legend = "bottom") %>%
        annotate_figure(top = text_grob(title, face = "bold", size = 14))
}
