# Effect Size Reporting for Group Comparisons
# Run: source("utils/effect_size.R")
#
# Computes Bayesian effect sizes from posterior distributions for
# group comparison analysis in IGT hierarchical models.

if (!require("pacman")) install.packages("pacman")
pacman::p_load(coda, ggplot2)

# ==============================================================================
# EFFECT SIZE FUNCTIONS
# ==============================================================================

#' Compute Cohen's d from two posterior distributions
#' @param posterior_g1 Vector of posterior samples for group 1
#' @param posterior_g2 Vector of posterior samples for group 2
#' @param pooled_sd Optional: provide pooled SD, otherwise estimated from posteriors
#' @return Distribution of Cohen's d values
compute_cohens_d <- function(posterior_g1, posterior_g2, pooled_sd = NULL) {
    # Mean difference for each posterior sample
    mean_diff <- posterior_g2 - posterior_g1

    if (is.null(pooled_sd)) {
        # Estimate pooled SD from posterior variance
        # This is approximate - better to use within-group SDs if available
        combined <- c(posterior_g1, posterior_g2)
        pooled_sd <- sd(combined)
    }

    # Cohen's d for each sample
    d <- mean_diff / pooled_sd

    return(d)
}

#' Compute Cohen's d from alpha (difference) parameter
#' Uses the joint model parameterization: alpha = group2 - group1
#' @param alpha_samples Vector of posterior samples for alpha (difference)
#' @param sigma_samples Optional: posterior samples for group-level SD
#' @return Distribution of Cohen's d values
compute_cohens_d_from_alpha <- function(alpha_samples, sigma_samples = NULL) {
    if (!is.null(sigma_samples)) {
        # Use posterior samples of sigma for proper uncertainty propagation
        d <- alpha_samples / sigma_samples
    } else {
        # Use posterior SD of alpha as approximate denominator
        d <- alpha_samples / sd(alpha_samples)
    }

    return(d)
}

#' Compute probability of superiority (common language effect size)
#' P(X_g2 > X_g1) - probability that a random draw from g2 exceeds g1
#' @param posterior_g1 Vector of posterior samples for group 1
#' @param posterior_g2 Vector of posterior samples for group 2
#' @return Probability of superiority
compute_prob_superiority <- function(posterior_g1, posterior_g2) {
    # For each posterior sample, compute whether g2 > g1
    # This naturally integrates over parameter uncertainty
    prob <- mean(posterior_g2 > posterior_g1)
    return(prob)
}

#' Compute probability of superiority from alpha
#' @param alpha_samples Vector of posterior samples for alpha
#' @return Probability that alpha > 0 (group 2 > group 1)
compute_prob_superiority_from_alpha <- function(alpha_samples) {
    prob <- mean(alpha_samples > 0)
    return(prob)
}

#' Compute credible interval
#' @param samples Vector of posterior samples
#' @param prob Coverage probability (default 0.95)
#' @return Named vector with lower and upper bounds
compute_credible_interval <- function(samples, prob = 0.95) {
    alpha <- (1 - prob) / 2
    ci <- quantile(samples, probs = c(alpha, 1 - alpha))
    names(ci) <- c("lower", "upper")
    return(ci)
}

#' Compute HDI (Highest Density Interval)
#' @param samples Vector of posterior samples
#' @param prob Coverage probability (default 0.95)
#' @return Named vector with lower and upper bounds
compute_hdi <- function(samples, prob = 0.95) {
    sorted <- sort(samples)
    n <- length(sorted)
    ci_width <- floor(prob * n)

    # Find narrowest interval
    widths <- sorted[(ci_width + 1):n] - sorted[1:(n - ci_width)]
    min_idx <- which.min(widths)

    hdi <- c(lower = sorted[min_idx], upper = sorted[min_idx + ci_width])
    return(hdi)
}

# ==============================================================================
# COMPREHENSIVE EFFECT SIZE REPORT
# ==============================================================================

#' Generate effect size report for group comparison
#' @param fit JAGS fit object from compare_*.R script
#' @param model_name Model name for labeling
#' @param group1_name Name of group 1 (reference)
#' @param group2_name Name of group 2 (comparison)
#' @return Data frame with effect size statistics
generate_effect_size_report <- function(fit, model_name,
                                        group1_name = "HC",
                                        group2_name = "Clinical") {
    cat(sprintf("\n=== Effect Size Report: %s ===\n", model_name))
    cat(sprintf("Comparison: %s vs %s\n\n", group1_name, group2_name))

    samples <- fit$BUGSoutput$sims.list

    # Find alpha parameters (group differences)
    alpha_params <- names(samples)[grep("^alpha_", names(samples))]

    if (length(alpha_params) == 0) {
        # Try alternative naming
        alpha_params <- names(samples)[grep("^delta_", names(samples))]
    }

    results_list <- list()

    for (param in alpha_params) {
        alpha_samples <- samples[[param]]

        # Get corresponding mu parameter for sigma estimation
        param_base <- gsub("^alpha_", "", param)
        mu_param <- paste0("mu_", param_base)

        sigma_samples <- NULL
        if (mu_param %in% names(samples)) {
            # Approximate sigma from posterior variance of mu
            sigma_samples <- rep(sd(samples[[mu_param]]), length(alpha_samples))
        }

        # Compute effect sizes
        d <- compute_cohens_d_from_alpha(alpha_samples, sigma_samples)
        prob_sup <- compute_prob_superiority_from_alpha(alpha_samples)
        ci_alpha <- compute_credible_interval(alpha_samples)
        ci_d <- compute_credible_interval(d)
        hdi_alpha <- compute_hdi(alpha_samples)

        results_list[[param]] <- data.frame(
            model = model_name,
            parameter = param_base,
            group1 = group1_name,
            group2 = group2_name,

            # Raw difference (alpha)
            alpha_mean = mean(alpha_samples),
            alpha_sd = sd(alpha_samples),
            alpha_ci_lower = ci_alpha["lower"],
            alpha_ci_upper = ci_alpha["upper"],
            alpha_hdi_lower = hdi_alpha["lower"],
            alpha_hdi_upper = hdi_alpha["upper"],

            # Cohen's d
            d_mean = mean(d),
            d_sd = sd(d),
            d_ci_lower = ci_d["lower"],
            d_ci_upper = ci_d["upper"],

            # Probability of superiority
            prob_g2_gt_g1 = prob_sup,

            # Evidence interpretation
            evidence = interpret_bayes_factor(prob_sup),
            stringsAsFactors = FALSE
        )
    }

    results_df <- do.call(rbind, results_list)
    rownames(results_df) <- NULL

    # Print summary
    cat("Parameter Effect Sizes:\n")
    cat("------------------------\n")
    for (i in 1:nrow(results_df)) {
        row <- results_df[i, ]
        cat(sprintf("\n%s:\n", row$parameter))
        cat(sprintf(
            "  Difference (α): %.3f [%.3f, %.3f]\n",
            row$alpha_mean, row$alpha_ci_lower, row$alpha_ci_upper
        ))
        cat(sprintf(
            "  Cohen's d: %.3f [%.3f, %.3f]\n",
            row$d_mean, row$d_ci_lower, row$d_ci_upper
        ))
        cat(sprintf(
            "  P(%s > %s): %.3f (%s)\n",
            row$group2, row$group1, row$prob_g2_gt_g1, row$evidence
        ))
    }

    return(results_df)
}

#' Interpret probability as Bayes Factor evidence
#' @param prob Probability of superiority
#' @return Character string interpretation
interpret_bayes_factor <- function(prob) {
    # Convert probability to odds ratio (approximates Bayes Factor)
    if (prob > 0.5) {
        bf <- prob / (1 - prob)
    } else {
        bf <- (1 - prob) / prob
    }

    direction <- ifelse(prob > 0.5, "+", "-")

    if (bf < 1) {
        return("No evidence")
    } else if (bf < 3) {
        return(paste0("Anecdotal", direction))
    } else if (bf < 10) {
        return(paste0("Moderate", direction))
    } else if (bf < 30) {
        return(paste0("Strong", direction))
    } else if (bf < 100) {
        return(paste0("Very strong", direction))
    } else {
        return(paste0("Extreme", direction))
    }
}

#' Interpret Cohen's d magnitude
#' @param d Cohen's d value
#' @return Character string interpretation
interpret_cohens_d <- function(d) {
    abs_d <- abs(d)

    if (abs_d < 0.2) {
        return("negligible")
    } else if (abs_d < 0.5) {
        return("small")
    } else if (abs_d < 0.8) {
        return("medium")
    } else {
        return("large")
    }
}

# ==============================================================================
# PLOTTING FUNCTIONS
# ==============================================================================

#' Plot effect size forest plot
#' @param results_df Data frame from generate_effect_size_report
#' @param output_file Optional file to save plot
#' @return ggplot object
plot_effect_sizes <- function(results_df, output_file = NULL) {
    p <- ggplot(results_df, aes(x = d_mean, y = parameter)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        geom_vline(xintercept = c(-0.2, 0.2), linetype = "dotted", color = "gray70") +
        geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "gray70") +
        geom_vline(xintercept = c(-0.8, 0.8), linetype = "dotted", color = "gray70") +
        geom_errorbarh(aes(xmin = d_ci_lower, xmax = d_ci_upper), height = 0.2) +
        geom_point(size = 3) +
        labs(
            x = "Cohen's d",
            y = "Parameter",
            title = paste0("Effect Sizes: ", unique(results_df$group1), " vs ", unique(results_df$group2)),
            subtitle = "95% Credible Intervals"
        ) +
        theme_minimal() +
        theme(
            panel.grid.minor = element_blank()
        )

    if (!is.null(output_file)) {
        ggsave(output_file, p, width = 8, height = 6, dpi = 150)
        cat("\nPlot saved to:", output_file, "\n")
    }

    return(p)
}

#' Plot probability of superiority
#' @param results_df Data frame from generate_effect_size_report
#' @param output_file Optional file to save plot
#' @return ggplot object
plot_prob_superiority <- function(results_df, output_file = NULL) {
    p <- ggplot(results_df, aes(x = prob_g2_gt_g1, y = parameter)) +
        geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50") +
        geom_vline(xintercept = c(0.75, 0.90, 0.95), linetype = "dotted", color = "gray70", alpha = 0.7) +
        geom_vline(xintercept = c(0.25, 0.10, 0.05), linetype = "dotted", color = "gray70", alpha = 0.7) +
        geom_point(size = 4, aes(color = evidence)) +
        geom_text(aes(label = sprintf("%.2f", prob_g2_gt_g1)), hjust = -0.3, size = 3) +
        scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
        scale_color_manual(
            values = c(
                "No evidence" = "gray50",
                "Anecdotal+" = "lightblue",
                "Anecdotal-" = "lightblue",
                "Moderate+" = "steelblue",
                "Moderate-" = "steelblue",
                "Strong+" = "darkblue",
                "Strong-" = "darkblue",
                "Very strong+" = "purple",
                "Very strong-" = "purple",
                "Extreme+" = "red",
                "Extreme-" = "red"
            )
        ) +
        labs(
            x = paste0("P(", unique(results_df$group2), " > ", unique(results_df$group1), ")"),
            y = "Parameter",
            title = "Probability of Superiority",
            color = "Evidence"
        ) +
        theme_minimal()

    if (!is.null(output_file)) {
        ggsave(output_file, p, width = 10, height = 6, dpi = 150)
        cat("\nPlot saved to:", output_file, "\n")
    }

    return(p)
}

# ==============================================================================
# RUN ALL COMPARISONS
# ==============================================================================

#' Generate effect size reports for all group comparisons
#' @param comparison_dir Directory with comparison fit files
#' @param output_dir Directory for output files
#' @return List of effect size results
run_all_effect_sizes <- function(comparison_dir = "analysis/outputs/group_comparison",
                                 output_dir = "analysis/outputs/effect_sizes") {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    # Find comparison files
    fit_files <- list.files(comparison_dir, pattern = "^compare_.*\\.rds$", full.names = TRUE)

    if (length(fit_files) == 0) {
        cat("No comparison files found in", comparison_dir, "\n")
        return(NULL)
    }

    all_results <- list()

    for (fit_file in fit_files) {
        # Parse filename
        basename <- gsub("\\.rds$", "", basename(fit_file))
        parts <- strsplit(basename, "_vs_")[[1]]

        # Extract model and groups
        pre_vs <- parts[1]
        group2 <- parts[2]

        model_parts <- strsplit(pre_vs, "_")[[1]]
        model_name <- paste(model_parts[2:(length(model_parts) - 1)], collapse = "_")
        group1 <- model_parts[length(model_parts)]

        cat(sprintf("\n\nProcessing: %s (%s vs %s)\n", model_name, group1, group2))

        # Load fit
        fit <- readRDS(fit_file)

        # Generate report
        results <- generate_effect_size_report(fit, model_name, group1, group2)

        # Save individual results
        output_name <- sprintf("effect_sizes_%s_%s_vs_%s", model_name, group1, group2)
        write.csv(results, file.path(output_dir, paste0(output_name, ".csv")), row.names = FALSE)

        # Generate plots
        plot_effect_sizes(results, file.path(output_dir, paste0(output_name, "_forest.png")))
        plot_prob_superiority(results, file.path(output_dir, paste0(output_name, "_prob.png")))

        all_results[[basename]] <- results
    }

    # Combine all results
    combined <- do.call(rbind, all_results)
    write.csv(combined, file.path(output_dir, "all_effect_sizes.csv"), row.names = FALSE)

    cat("\n\n=== All Effect Size Reports Complete ===\n")
    cat("Results saved to:", output_dir, "\n")

    return(all_results)
}

# ==============================================================================
# LATEX TABLE GENERATION
# ==============================================================================

#' Generate LaTeX table for effect sizes
#' @param results_df Data frame from generate_effect_size_report
#' @param caption Table caption
#' @param label Table label
#' @return Character string with LaTeX code
effect_size_to_latex <- function(results_df,
                                 caption = "Effect sizes for group comparison",
                                 label = "tab:effect_sizes") {
    # Format for display
    df_fmt <- data.frame(
        Parameter = results_df$parameter,
        `$\\alpha$` = sprintf(
            "%.3f [%.3f, %.3f]",
            results_df$alpha_mean,
            results_df$alpha_ci_lower,
            results_df$alpha_ci_upper
        ),
        `Cohen's $d$` = sprintf(
            "%.3f [%.3f, %.3f]",
            results_df$d_mean,
            results_df$d_ci_lower,
            results_df$d_ci_upper
        ),
        `$P(G_2 > G_1)$` = sprintf("%.3f", results_df$prob_g2_gt_g1),
        Evidence = results_df$evidence,
        check.names = FALSE
    )

    # Build LaTeX
    header <- paste(names(df_fmt), collapse = " & ")

    rows <- apply(df_fmt, 1, function(row) {
        paste(row, collapse = " & ")
    })

    latex <- c(
        "\\begin{table}[htbp]",
        "\\centering",
        sprintf("\\caption{%s}", caption),
        sprintf("\\label{%s}", label),
        paste0("\\begin{tabular}{l", paste(rep("c", ncol(df_fmt) - 1), collapse = ""), "}"),
        "\\hline",
        paste0(header, " \\\\"),
        "\\hline",
        paste0(rows, " \\\\"),
        "\\hline",
        "\\end{tabular}",
        "\\end{table}"
    )

    return(paste(latex, collapse = "\n"))
}

# ==============================================================================
# USAGE EXAMPLE
# ==============================================================================
#
# source("utils/effect_size.R")
#
# # For single comparison
# fit <- readRDS("analysis/outputs/group_comparison/compare_pvl_delta_HC_vs_Hero.rds")
# results <- generate_effect_size_report(fit, "pvl_delta", "HC", "Hero")
#
# # Plot
# plot_effect_sizes(results, "effect_sizes_pvl_delta.png")
#
# # LaTeX
# latex_table <- effect_size_to_latex(results)
# cat(latex_table)
#
# # Run all comparisons
# all_results <- run_all_effect_sizes()
