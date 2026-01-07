# Visualization utilities for publication-quality figures
# Run: source("analysis/utils/visualization.R")

library(ggplot2)
library(gridExtra)
library(coda)

#' Create publication-quality trace plots
#' @param mcmc_samples MCMC samples
#' @param params Parameters to plot
#' @param output_file Output file path
#' @param ncol Number of columns in layout
create_trace_plots <- function(mcmc_samples, params = NULL,
                               output_file = "trace_plots.pdf",
                               ncol = 2) {

  if (is.null(params)) {
    # Default to group-level parameters
    all_params <- colnames(mcmc_samples[[1]])
    params <- grep("^(mu_|sigma_)", all_params, value = TRUE)
  }

  # Convert to data frame
  samples_df <- as.data.frame(as.matrix(mcmc_samples))
  samples_df$iteration <- 1:nrow(samples_df)
  samples_df$chain <- rep(1:length(mcmc_samples),
                          each = nrow(samples_df) / length(mcmc_samples))

  # Create plots
  plot_list <- list()

  for (param in params) {
    p <- ggplot(samples_df, aes(x = iteration, y = .data[[param]], color = factor(chain))) +
      geom_line(alpha = 0.7) +
      labs(title = param, x = "Iteration", y = "Value") +
      theme_minimal() +
      theme(legend.position = "none")

    plot_list[[param]] <- p
  }

  # Save to PDF
  pdf(output_file, width = 12, height = 3 * ceiling(length(params) / ncol))
  grid.arrange(grobs = plot_list, ncol = ncol)
  dev.off()

  cat(sprintf("Trace plots saved to: %s\n", output_file))
}

#' Create posterior density plots
#' @param mcmc_samples MCMC samples
#' @param params Parameters to plot
#' @param output_file Output file path
create_density_plots <- function(mcmc_samples, params = NULL,
                                 output_file = "density_plots.pdf") {

  if (is.null(params)) {
    all_params <- colnames(mcmc_samples[[1]])
    params <- grep("^mu_", all_params, value = TRUE)
  }

  samples_df <- as.data.frame(as.matrix(mcmc_samples))

  plot_list <- list()

  for (param in params) {
    # Get summary stats
    mean_val <- mean(samples_df[[param]])
    hdi_lower <- quantile(samples_df[[param]], 0.025)
    hdi_upper <- quantile(samples_df[[param]], 0.975)

    p <- ggplot(samples_df, aes(x = .data[[param]])) +
      geom_density(fill = "steelblue", alpha = 0.5) +
      geom_vline(xintercept = mean_val, linetype = "dashed", color = "red") +
      geom_vline(xintercept = c(hdi_lower, hdi_upper),
                linetype = "dotted", color = "darkgray") +
      labs(title = param,
           subtitle = sprintf("Mean: %.3f, 95%% HDI: [%.3f, %.3f]",
                            mean_val, hdi_lower, hdi_upper),
           x = "Value", y = "Density") +
      theme_minimal()

    plot_list[[param]] <- p
  }

  pdf(output_file, width = 10, height = 3 * ceiling(length(params) / 2))
  grid.arrange(grobs = plot_list, ncol = 2)
  dev.off()

  cat(sprintf("Density plots saved to: %s\n", output_file))
}

#' Create parameter comparison plot across models
#' @param fit_results List of fitted models
#' @param param Parameter name (e.g., "mu_A")
#' @param output_file Output file path
compare_parameter_across_models <- function(fit_results, param,
                                           output_file = NULL) {

  # Extract parameter from each model
  param_data <- list()

  for (model_name in names(fit_results)) {
    samples <- as.matrix(fit_results[[model_name]]$samples)

    if (param %in% colnames(samples)) {
      param_data[[model_name]] <- data.frame(
        value = samples[, param],
        model = model_name
      )
    }
  }

  if (length(param_data) == 0) {
    warning(sprintf("Parameter %s not found in any model", param))
    return(NULL)
  }

  # Combine data
  combined_data <- do.call(rbind, param_data)

  # Create plot
  p <- ggplot(combined_data, aes(x = value, fill = model)) +
    geom_density(alpha = 0.5) +
    labs(title = sprintf("Posterior Distribution: %s", param),
         x = "Value", y = "Density") +
    theme_minimal() +
    theme(legend.position = "bottom")

  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 8, height = 6)
    cat(sprintf("Comparison plot saved to: %s\n", output_file))
  }

  return(p)
}

#' Create correlation plot for parameter recovery
#' @param recovery_results Parameter recovery results
#' @param output_file Output file path
plot_parameter_recovery <- function(recovery_results, output_file = "recovery_plots.pdf") {

  plot_list <- list()

  for (param_name in names(recovery_results$recovery_results)) {
    recovery_data <- recovery_results$recovery_results[[param_name]]

    # Correlation
    cor_val <- cor(recovery_data$true, recovery_data$recovered)

    p <- ggplot(recovery_data, aes(x = true, y = recovered)) +
      geom_point(alpha = 0.6) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      geom_smooth(method = "lm", se = TRUE, color = "blue") +
      labs(title = param_name,
           subtitle = sprintf("r = %.3f", cor_val),
           x = "True Value", y = "Recovered Value") +
      theme_minimal()

    plot_list[[param_name]] <- p
  }

  pdf(output_file, width = 10, height = 3 * ceiling(length(plot_list) / 2))
  grid.arrange(grobs = plot_list, ncol = 2)
  dev.off()

  cat(sprintf("Parameter recovery plots saved to: %s\n", output_file))
}

#' Create forest plot for parameter estimates
#' @param fit_results List of fitted models
#' @param params Parameters to plot
#' @param output_file Output file path
create_forest_plot <- function(fit_results, params = NULL,
                               output_file = "forest_plot.pdf") {

  # Extract parameter estimates
  param_estimates <- list()

  for (model_name in names(fit_results)) {
    samples <- as.matrix(fit_results[[model_name]]$samples)

    if (is.null(params)) {
      params <- grep("^mu_", colnames(samples), value = TRUE)
    }

    for (param in params) {
      if (param %in% colnames(samples)) {
        values <- samples[, param]
        param_estimates[[length(param_estimates) + 1]] <- data.frame(
          model = model_name,
          parameter = param,
          mean = mean(values),
          lower = quantile(values, 0.025),
          upper = quantile(values, 0.975)
        )
      }
    }
  }

  # Combine data
  estimates_df <- do.call(rbind, param_estimates)

  # Create plot
  p <- ggplot(estimates_df, aes(x = mean, y = parameter, color = model)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbarh(aes(xmin = lower, xmax = upper),
                   position = position_dodge(width = 0.5),
                   height = 0.2) +
    labs(title = "Parameter Estimates (95% Credible Intervals)",
         x = "Estimate", y = "Parameter") +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(output_file, p, width = 10, height = 6)
  cat(sprintf("Forest plot saved to: %s\n", output_file))

  return(p)
}

#' Create choice proportion plots
#' @param observed_data Observed choice data
#' @param predicted_data Predicted choice data (from PPC)
#' @param output_file Output file path
plot_choice_proportions <- function(observed_data, predicted_data = NULL,
                                   output_file = "choice_proportions.pdf") {

  # Compute observed proportions
  obs_props <- table(observed_data$choice) / nrow(observed_data)

  plot_data <- data.frame(
    deck = names(obs_props),
    observed = as.numeric(obs_props),
    type = "Observed"
  )

  if (!is.null(predicted_data)) {
    pred_props <- colMeans(predicted_data)
    pred_lower <- apply(predicted_data, 2, quantile, probs = 0.025)
    pred_upper <- apply(predicted_data, 2, quantile, probs = 0.975)

    pred_df <- data.frame(
      deck = 1:4,
      predicted = pred_props,
      lower = pred_lower,
      upper = pred_upper
    )

    plot_data <- merge(plot_data, pred_df, by = "deck")

    p <- ggplot(plot_data, aes(x = deck)) +
      geom_bar(aes(y = observed), stat = "identity", fill = "steelblue", alpha = 0.5) +
      geom_point(aes(y = predicted), color = "red", size = 3) +
      geom_errorbar(aes(ymin = lower, ymax = upper), color = "red", width = 0.2) +
      labs(title = "Observed vs Predicted Choice Proportions",
           subtitle = "Red: Model predictions with 95% CI | Blue: Observed",
           x = "Deck", y = "Proportion") +
      theme_minimal()
  } else {
    p <- ggplot(plot_data, aes(x = deck, y = observed)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      labs(title = "Observed Choice Proportions",
           x = "Deck", y = "Proportion") +
      theme_minimal()
  }

  ggsave(output_file, p, width = 8, height = 6)
  cat(sprintf("Choice proportion plot saved to: %s\n", output_file))

  return(p)
}

#' Generate all publication figures
#' @param models_dir Directory with fitted models
#' @param output_dir Output directory for figures
generate_all_figures <- function(models_dir = "analysis/outputs",
                                output_dir = "analysis/outputs/figures") {

  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cat("=== Generating Publication Figures ===\n\n")

  # Load fitted models
  fit_files <- list.files(models_dir, pattern = "_fit\\.rds$", full.names = TRUE)

  if (length(fit_files) == 0) {
    stop("No fitted models found in: ", models_dir)
  }

  fit_results <- list()
  for (fit_file in fit_files) {
    model_name <- gsub("_fit\\.rds$", "", basename(fit_file))
    fit_results[[model_name]] <- readRDS(fit_file)
  }

  # 1. Trace plots for each model
  cat("1. Creating trace plots...\n")
  for (model_name in names(fit_results)) {
    trace_file <- file.path(output_dir, sprintf("%s_traces.pdf", model_name))
    create_trace_plots(fit_results[[model_name]]$samples, output_file = trace_file)
  }

  # 2. Density plots for each model
  cat("\n2. Creating posterior density plots...\n")
  for (model_name in names(fit_results)) {
    density_file <- file.path(output_dir, sprintf("%s_posteriors.pdf", model_name))
    create_density_plots(fit_results[[model_name]]$samples, output_file = density_file)
  }

  # 3. Forest plot comparing models
  cat("\n3. Creating forest plot...\n")
  forest_file <- file.path(output_dir, "model_comparison_forest.pdf")
  create_forest_plot(fit_results, output_file = forest_file)

  cat("\n=== All figures generated ===\n")
  cat(sprintf("Saved to: %s\n", output_dir))

  invisible(fit_results)
}
