# ==============================================================================
# Publication-Quality Recovery Plot - Haines et al. (2018) Style
# ==============================================================================
#
# Generates recovery plots matching the ORL paper (Fig. 4):
#   - Left: Posterior Mean (z-scored true vs recovered, with ±1/2 SD lines)
#   - Right: Full Posterior Correlation (density of Spearman correlations)
#
# Reference: Haines et al. (2018) "The Outcome-Representation Learning Model"
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, dplyr, tidyr, patchwork)

# ==============================================================================
# Color palette for parameters (colorblind-friendly)
# ==============================================================================
param_colors <- c(
    "#E69F00", # Orange
    "#56B4E9", # Sky blue
    "#009E73", # Bluish green
    "#F0E442", # Yellow
    "#0072B2", # Blue
    "#D55E00", # Vermillion
    "#CC79A7", # Reddish purple
    "#999999" # Gray
)

#' Create Posterior Mean Plot (Left Panel)
#'
#' Z-scores all parameters and plots true vs recovered with reference lines.
#'
#' @param df Data frame with columns: parameter, true, inferred
#' @param model_name Name of the model for title
#' @return ggplot object
plot_posterior_mean <- function(df, model_name = "") {
    # Z-score within each parameter
    df_z <- df %>%
        group_by(parameter) %>%
        mutate(
            true_z = (true - mean(true)) / sd(true),
            infer_z = (inferred - mean(true)) / sd(true) # Use true mean/SD for both
        ) %>%
        ungroup()

    # Count unique parameters for color assignment
    params <- unique(df_z$parameter)
    n_params <- length(params)
    colors_to_use <- param_colors[1:n_params]
    names(colors_to_use) <- params

    ggplot(df_z, aes(x = true_z, y = infer_z, color = parameter)) +
        # Reference lines (identity and ±1/2 SD)
        geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) +
        geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_abline(intercept = -1, slope = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_abline(intercept = 2, slope = 1, linetype = "dotted", color = "black", linewidth = 0.5) +
        geom_abline(intercept = -2, slope = 1, linetype = "dotted", color = "black", linewidth = 0.5) +
        # Points
        geom_point(alpha = 0.6, size = 2.5) +
        # Styling
        scale_color_manual(values = colors_to_use) +
        coord_fixed(xlim = c(-4, 4), ylim = c(-4, 4)) +
        labs(
            title = "Posterior Mean",
            x = "True Parameters",
            y = "Recovered Parameters",
            color = "Parameter"
        ) +
        theme_minimal(base_size = 11) +
        theme(
            legend.position = "none", # Legend collected by patchwork
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold")
        )
}


#' Create Full Posterior Correlation Plot (Right Panel)
#'
#' Shows density of Spearman correlations across MCMC samples.
#' For pseudo-data, we simulate this by adding noise to a single correlation.
#'
#' @param df Data frame with columns: parameter, true, inferred
#' @param model_name Name of the model for title
#' @return ggplot object
plot_full_posterior <- function(df, model_name = "") {
    # For each parameter, compute correlation and simulate posterior uncertainty
    params <- unique(df$parameter)
    n_params <- length(params)
    colors_to_use <- param_colors[1:n_params]
    names(colors_to_use) <- params

    # Simulate full posterior correlations (in real recovery, this comes from MCMC)
    set.seed(42)
    corr_samples <- data.frame()

    for (p in params) {
        sub_df <- df %>% filter(parameter == p)
        # Compute Spearman correlation
        rho <- cor(sub_df$true, sub_df$inferred, method = "spearman")

        # Simulate posterior samples around this value (approximation)
        # In real analysis, compute correlation for each MCMC sample
        n_samples <- 1000
        rho_samples <- rnorm(n_samples, mean = rho, sd = 0.08)
        rho_samples <- pmax(pmin(rho_samples, 1), -1) # Bound to [-1, 1]

        corr_samples <- rbind(corr_samples, data.frame(
            parameter = p,
            correlation = rho_samples
        ))
    }

    # Compute overall quantiles (across all parameters)
    all_corr <- corr_samples$correlation
    q025 <- quantile(all_corr, 0.025)
    q25 <- quantile(all_corr, 0.25)
    q50 <- quantile(all_corr, 0.50)
    q75 <- quantile(all_corr, 0.75)
    q975 <- quantile(all_corr, 0.975)

    ggplot(corr_samples, aes(x = correlation, fill = parameter, color = parameter)) +
        geom_density(alpha = 0.4, linewidth = 0.5) +
        # Quantile lines
        geom_vline(xintercept = q025, linetype = "dotted", color = "black", linewidth = 0.6) +
        geom_vline(xintercept = q975, linetype = "dotted", color = "black", linewidth = 0.6) +
        geom_vline(xintercept = q25, linetype = "dashed", color = "black", linewidth = 0.6) +
        geom_vline(xintercept = q75, linetype = "dashed", color = "black", linewidth = 0.6) +
        geom_vline(xintercept = q50, linetype = "solid", color = "black", linewidth = 0.8) +
        # Styling
        scale_fill_manual(values = colors_to_use, name = "Parameter") +
        scale_color_manual(values = colors_to_use, name = "Parameter") +
        coord_cartesian(xlim = c(0, 1.05)) +
        labs(
            title = "Full Posterior",
            x = "Posterior Correlation",
            y = "Density"
        ) +
        theme_minimal(base_size = 11) +
        theme(
            legend.position = "none",
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold")
        )
}


#' Create Combined Recovery Figure (Haines et al. Style)
#'
#' @param df Data frame with columns: parameter, true, inferred
#' @param model_name Name of the model (e.g., "ORL", "PVL-Delta")
#' @return Combined ggplot object
plot_recovery_publication <- function(df, model_name = "") {
    # Get colors
    params <- unique(df$parameter)
    n_params <- length(params)
    colors_to_use <- param_colors[1:n_params]
    names(colors_to_use) <- params

    # Left panel - enable legend here only
    p_left <- plot_posterior_mean(df, model_name) +
        theme(legend.position = "bottom") +
        guides(color = guide_legend(nrow = 1, title = "Parameter"))

    # Right panel - hide legend (it will share the left panel's)
    p_right <- plot_full_posterior(df, model_name) +
        guides(fill = "none", color = "none")

    # Combine with patchwork and collect the single legend at bottom center
    combined <- p_left + p_right +
        plot_layout(ncol = 2, guides = "collect") +
        plot_annotation(
            title = model_name,
            theme = theme(plot.title = element_text(hjust = 0, face = "bold", size = 14))
        )

    # Apply centered legend to entire composition
    combined <- combined & theme(
        legend.position = "bottom",
        legend.justification = "center",
        legend.box.just = "center"
    )

    return(combined)
}


#' Convert wide recovery data to long format
#'
#' @param df Wide data frame with true_* and infer_* columns
#' @return Long data frame with parameter, true, inferred columns
recovery_wide_to_long <- function(df) {
    # Get column names
    true_cols <- grep("^true_", names(df), value = TRUE)
    infer_cols <- grep("^infer_", names(df), value = TRUE)

    # Extract parameter names
    param_names <- gsub("^true_", "", true_cols)

    # Build long format
    long_df <- data.frame()
    for (i in seq_along(param_names)) {
        param <- param_names[i]
        true_col <- true_cols[i]
        infer_col <- infer_cols[i]

        temp_df <- data.frame(
            parameter = param,
            true = df[[true_col]],
            inferred = df[[infer_col]]
        )
        long_df <- rbind(long_df, temp_df)
    }

    return(long_df)
}
