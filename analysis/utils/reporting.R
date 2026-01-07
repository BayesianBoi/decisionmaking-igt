# Reporting utilities for manuscript preparation
# Run: source("analysis/utils/reporting.R")

library(coda)

#' Generate parameter summary table
#' @param fit_result Fitted model
#' @param model_name Model name
#' @return Data frame with parameter summaries
generate_parameter_table <- function(fit_result, model_name) {

  samples <- as.matrix(fit_result$samples)

  # Get group-level parameters
  mu_params <- grep("^mu_", colnames(samples), value = TRUE)

  # Create summary table
  summary_list <- list()

  for (param in mu_params) {
    values <- samples[, param]

    summary_list[[param]] <- data.frame(
      Model = model_name,
      Parameter = param,
      Mean = mean(values),
      SD = sd(values),
      CI_lower = quantile(values, 0.025),
      CI_upper = quantile(values, 0.975),
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL

  return(summary_df)
}

#' Generate model comparison table
#' @param fit_results List of fitted models
#' @param output_file Optional CSV file to save
#' @return Data frame with model comparison
generate_model_comparison_table <- function(fit_results, output_file = NULL) {

  comparison_list <- list()

  for (model_name in names(fit_results)) {
    fit <- fit_results[[model_name]]
    samples <- as.matrix(fit$samples)

    # Count parameters
    n_params <- ncol(samples)
    n_group_params <- length(grep("^(mu_|sigma_)", colnames(samples)))
    n_subj_params <- n_params - n_group_params

    comparison_list[[model_name]] <- data.frame(
      Model = model_name,
      N_group_params = n_group_params,
      N_subject_params = n_subj_params / length(unique(gsub("\\[.*\\]", "", colnames(samples)[grep("\\[", colnames(samples))]))),
      Total_params = n_params,
      stringsAsFactors = FALSE
    )
  }

  comparison_df <- do.call(rbind, comparison_list)
  rownames(comparison_df) <- NULL

  if (!is.null(output_file)) {
    write.csv(comparison_df, output_file, row.names = FALSE)
    cat(sprintf("Model comparison table saved to: %s\n", output_file))
  }

  return(comparison_df)
}

#' Generate convergence diagnostics table
#' @param diagnostics_results Diagnostics from run_full_diagnostics
#' @param output_file Optional CSV file to save
#' @return Data frame with convergence summaries
generate_convergence_table <- function(diagnostics_results, output_file = NULL) {

  convergence_list <- list()

  for (model_name in names(diagnostics_results)) {
    diag <- diagnostics_results[[model_name]]

    # Get group-level parameters only
    group_params <- grep("^(mu_|sigma_)", diag$diagnostics$parameter, value = TRUE)
    group_diag <- diag$diagnostics[diag$diagnostics$parameter %in% group_params, ]

    convergence_list[[model_name]] <- data.frame(
      Model = model_name,
      Max_Rhat = max(group_diag$rhat),
      Mean_Rhat = mean(group_diag$rhat),
      Min_ESS = min(group_diag$ess),
      Mean_ESS = mean(group_diag$ess),
      Converged = diag$convergence$converged,
      stringsAsFactors = FALSE
    )
  }

  convergence_df <- do.call(rbind, convergence_list)
  rownames(convergence_df) <- NULL

  if (!is.null(output_file)) {
    write.csv(convergence_df, output_file, row.names = FALSE)
    cat(sprintf("Convergence table saved to: %s\n", output_file))
  }

  return(convergence_df)
}

#' Generate formatted results for manuscript
#' @param fit_results List of fitted models
#' @param diagnostics_results Diagnostics results
#' @param output_dir Output directory
generate_manuscript_tables <- function(fit_results,
                                      diagnostics_results = NULL,
                                      output_dir = "analysis/outputs/tables") {

  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cat("=== Generating Manuscript Tables ===\n\n")

  # Table 1: Parameter estimates
  cat("1. Generating parameter estimates table...\n")
  param_tables <- list()
  for (model_name in names(fit_results)) {
    param_tables[[model_name]] <- generate_parameter_table(fit_results[[model_name]], model_name)
  }
  all_params <- do.call(rbind, param_tables)

  param_file <- file.path(output_dir, "table1_parameter_estimates.csv")
  write.csv(all_params, param_file, row.names = FALSE)
  cat(sprintf("  Saved to: %s\n", param_file))

  # Table 2: Model comparison
  cat("\n2. Generating model comparison table...\n")
  comparison_file <- file.path(output_dir, "table2_model_comparison.csv")
  comparison_df <- generate_model_comparison_table(fit_results, comparison_file)

  # Table 3: Convergence diagnostics
  if (!is.null(diagnostics_results)) {
    cat("\n3. Generating convergence diagnostics table...\n")
    convergence_file <- file.path(output_dir, "table3_convergence.csv")
    convergence_df <- generate_convergence_table(diagnostics_results, convergence_file)
  }

  cat("\n=== All tables generated ===\n")
  cat(sprintf("Tables saved to: %s\n", output_dir))

  # Return combined results
  results <- list(
    parameters = all_params,
    comparison = comparison_df
  )

  if (!is.null(diagnostics_results)) {
    results$convergence <- convergence_df
  }

  return(results)
}

#' Generate LaTeX table from data frame
#' @param df Data frame
#' @param caption Table caption
#' @param label Table label
#' @return LaTeX string
generate_latex_table <- function(df, caption = "", label = "") {

  # Round numeric columns
  df_formatted <- df
  for (col in names(df)) {
    if (is.numeric(df[[col]])) {
      df_formatted[[col]] <- sprintf("%.3f", df[[col]])
    }
  }

  # Start table
  latex <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\begin{tabular}{",
    paste(rep("l", ncol(df)), collapse = ""),
    "}",
    "\\hline"
  )

  # Header
  header <- paste(names(df), collapse = " & ")
  latex <- c(latex, paste(header, "\\\\"), "\\hline")

  # Rows
  for (i in 1:nrow(df)) {
    row <- paste(df_formatted[i, ], collapse = " & ")
    latex <- c(latex, paste(row, "\\\\"))
  }

  # End table
  latex <- c(latex,
             "\\hline",
             "\\end{tabular}",
             "\\end{table}")

  return(paste(latex, collapse = "\n"))
}

#' Write results summary for manuscript
#' @param fit_results List of fitted models
#' @param output_file Output markdown file
write_results_summary <- function(fit_results, output_file = "analysis/outputs/RESULTS_SUMMARY.md") {

  cat("=== Writing Results Summary ===\n")

  lines <- c(
    "# Results Summary for Manuscript",
    "",
    sprintf("Generated: %s", Sys.time()),
    "",
    "## Model Fits",
    ""
  )

  for (model_name in names(fit_results)) {
    samples <- as.matrix(fit_results[[model_name]]$samples)

    # Get group-level parameters
    mu_params <- grep("^mu_", colnames(samples), value = TRUE)

    lines <- c(lines,
               sprintf("### %s Model", toupper(model_name)),
               "")

    for (param in mu_params) {
      values <- samples[, param]
      mean_val <- mean(values)
      ci_lower <- quantile(values, 0.025)
      ci_upper <- quantile(values, 0.975)

      lines <- c(lines,
                 sprintf("- **%s**: M = %.3f, 95%% CI [%.3f, %.3f]",
                        param, mean_val, ci_lower, ci_upper))
    }

    lines <- c(lines, "")
  }

  # Write to file
  writeLines(lines, output_file)
  cat(sprintf("Results summary saved to: %s\n", output_file))
}

#' Generate complete publication package
#' @param models_dir Directory with fitted models
#' @param output_dir Output directory
generate_publication_package <- function(models_dir = "analysis/outputs",
                                        output_dir = "analysis/outputs/publication") {

  cat("=== Generating Complete Publication Package ===\n\n")

  # Create subdirectories
  fig_dir <- file.path(output_dir, "figures")
  table_dir <- file.path(output_dir, "tables")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  # Load fitted models
  fit_files <- list.files(models_dir, pattern = "_fit\\.rds$", full.names = TRUE)
  fit_results <- list()
  for (fit_file in fit_files) {
    model_name <- gsub("_fit\\.rds$", "", basename(fit_file))
    fit_results[[model_name]] <- readRDS(fit_file)
  }

  # Load diagnostics if available
  diag_file <- file.path(models_dir, "all_diagnostics.rds")
  diagnostics_results <- NULL
  if (file.exists(diag_file)) {
    diagnostics_results <- readRDS(diag_file)
  }

  # Generate figures
  cat("Generating figures...\n")
  source("analysis/2_plotting/visualization.R")
  generate_all_figures(models_dir, fig_dir)

  # Generate tables
  cat("\nGenerating tables...\n")
  tables <- generate_manuscript_tables(fit_results, diagnostics_results, table_dir)

  # Generate results summary
  cat("\nGenerating results summary...\n")
  summary_file <- file.path(output_dir, "RESULTS_SUMMARY.md")
  write_results_summary(fit_results, summary_file)

  cat("\n=== Publication Package Complete ===\n")
  cat(sprintf("All outputs saved to: %s\n", output_dir))
  cat("\nContents:\n")
  cat(sprintf("  - Figures: %s/\n", fig_dir))
  cat(sprintf("  - Tables: %s/\n", table_dir))
  cat(sprintf("  - Summary: %s\n", summary_file))

  return(invisible(tables))
}
