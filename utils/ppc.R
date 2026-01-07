# Posterior Predictive Checks for IGT models
# Run: source("utils/ppc.R")

library(coda)

#' Simulate choice data from PVL-Delta model
#' nomenclature: a (learning rate), A (shape), theta (consistency), w (loss aversion)
#' Reference: Ahn et al. (2008)
#' @param params List with a, A, theta, w
#' @param outcomes Matrix of outcomes (subj x trial)
#' @param n_trials Number of trials per subject
#' @return Matrix of simulated choices
simulate_pvl_delta <- function(params, outcomes, n_trials) {
  n_subj <- length(params$a)
  choices_sim <- matrix(NA, nrow = n_subj, ncol = ncol(outcomes))

  # Handle n_trials being scalar or vector
  if (length(n_trials) == 1) n_trials <- rep(n_trials, n_subj)

  for (s in 1:n_subj) {
    # Initialize expected values
    ev <- rep(0, 4)

    for (t in 1:n_trials[s]) {
      # Compute choice probabilities (matching JAGS model formulation)
      v <- params$theta[s] * ev
      p_choice <- exp(v - max(v)) / sum(exp(v - max(v)))

      # Sample choice
      choices_sim[s, t] <- sample(1:4, size = 1, prob = p_choice)

      # Get outcome for chosen deck
      outcome <- outcomes[s, t]

      # Compute utility using A (sensitivity) and w (loss aversion)
      if (outcome >= 0) {
        util <- outcome^params$A[s]
      } else {
        util <- -params$w[s] * abs(outcome)^params$A[s]
      }

      # Update expected value using a (learning rate)
      chosen <- choices_sim[s, t]
      ev[chosen] <- ev[chosen] + params$a[s] * (util - ev[chosen])
    }
  }

  return(choices_sim)
}

#' Simulate choice data from ORL model
#' @param params List with Arew, Apun, K, betaF, betaP
#' @param outcomes Matrix of outcomes (subj x trial)
#' @param n_trials Number of trials per subject
#' @return Matrix of simulated choices
simulate_orl <- function(params, outcomes, n_trials) {
  n_subj <- length(params$Arew)
  choices_sim <- matrix(NA, nrow = n_subj, ncol = ncol(outcomes))

  # Handle n_trials being scalar or vector
  if (length(n_trials) == 1) n_trials <- rep(n_trials, n_subj)

  for (s in 1:n_subj) {
    # Initialize
    ev <- rep(0, 4)
    ef <- rep(0, 4)
    pers <- rep(0, 4)

    for (t in 1:n_trials[s]) {
      # Compute choice probabilities
      util_combined <- ev + ef * params$betaF[s] + pers * params$betaP[s]
      p_choice <- exp(util_combined - max(util_combined)) / sum(exp(util_combined - max(util_combined)))

      # Sample choice
      choices_sim[s, t] <- sample(1:4, size = 1, prob = p_choice)

      # Get outcome
      outcome <- outcomes[s, t]
      sign_out <- ifelse(outcome > 0, 1, ifelse(outcome < 0, -1, 0))

      # Prediction errors
      chosen <- choices_sim[s, t]
      PE_val <- outcome - ev[chosen]
      PE_freq <- sign_out - ef[chosen]

      # Learning rates
      lr_val <- ifelse(outcome >= 0, params$Arew[s], params$Apun[s])
      lr_freq <- ifelse(outcome >= 0, params$Arew[s], params$Apun[s])
      lr_fic <- ifelse(outcome >= 0, params$Apun[s], params$Arew[s])

      # Update chosen deck
      ev[chosen] <- ev[chosen] + lr_val * PE_val
      ef[chosen] <- ef[chosen] + lr_freq * PE_freq

      # Update unchosen decks (fictive)
      for (d in setdiff(1:4, chosen)) {
        PE_freq_fic <- -sign_out / 3 - ef[d]
        ef[d] <- ef[d] + lr_fic * PE_freq_fic
      }

      # Update perseverance
      pers[chosen] <- 1 / (1 + params$K[s])
      for (d in setdiff(1:4, chosen)) {
        pers[d] <- pers[d] / (1 + params$K[s])
      }
    }
  }

  return(choices_sim)
}

#' Simulate choice data from EEF model
#' @param params List with theta, lambda (forget), phi, cons
#' @param outcomes Matrix of outcomes (subj x trial)
#' @param n_trials Vector of number of trials per subject
#' @param w_ini Vector of initial weights (default rep(0,4) if missing)
#' @return Matrix of simulated choices
simulate_eef <- function(params, outcomes, n_trials, w_ini = NULL) {
  n_subj <- length(params$theta)
  # Handle n_trials being scalar or vector
  if (length(n_trials) == 1) n_trials <- rep(n_trials, n_subj)
  max_trials <- ncol(outcomes)

  choices_sim <- matrix(NA, nrow = n_subj, ncol = max_trials)

  # Default w_ini if not provided
  if (is.null(w_ini)) {
    w_ini <- rep(0, 4)
  }

  for (s in 1:n_subj) {
    # Initialize
    exploit <- w_ini # Start with initial priors
    explore <- rep(0, 4)

    # Individual parameters
    p_theta <- params$theta[s]
    p_lambda <- params$lambda[s] # learning/forgetting rate
    p_phi <- params$phi[s]
    p_cons <- params$cons[s]

    for (t in 1:n_trials[s]) {
      # Compute choice probabilities
      # Composite weight: Exploitation + Exploration
      w <- exploit + explore
      v <- p_cons * w
      p_choice <- exp(v - max(v)) / sum(exp(v - max(v)))

      # Sample choice
      choices_sim[s, t] <- sample(1:4, size = 1, prob = p_choice)
      chosen <- choices_sim[s, t]

      # Get outcome
      outcome <- outcomes[s, t]
      gain <- max(outcome, 0)
      loss <- max(-outcome, 0)

      # Utility (symmetric power function)
      util <- (gain + 0.001)^p_theta - (loss + 0.001)^p_theta

      # Update Exploitation
      # Chosen: (1-lambda)*old + util
      # Unchosen: (1-lambda)*old
      for (d in 1:4) {
        if (d == chosen) {
          exploit[d] <- (1 - p_lambda) * exploit[d] + util
        } else {
          exploit[d] <- (1 - p_lambda) * exploit[d]
        }
      }

      # Update Exploration
      # Chosen: reset to 0
      # Unchosen: lambda*old + (1-lambda)*phi
      for (d in 1:4) {
        if (d == chosen) {
          explore[d] <- 0
        } else {
          explore[d] <- (1 - p_lambda) * explore[d] + p_lambda * p_phi
        }
      }
    }
  }

  return(choices_sim)
}

#' Compute choice proportions for each deck
#' @param choices Matrix of choices (subj x trial)
#' @return Vector of proportions for decks 1-4
compute_choice_proportions <- function(choices) {
  all_choices <- as.vector(choices)
  all_choices <- all_choices[!is.na(all_choices)]

  props <- table(all_choices) / length(all_choices)

  # Ensure all 4 decks are represented
  result <- rep(0, 4)
  names(result) <- 1:4
  result[names(props)] <- props

  return(result)
}

#' Compute block-wise choice patterns
#' @param choices Matrix of choices (subj x trial)
#' @param block_size Number of trials per block
#' @return Data frame with block-wise proportions
compute_block_proportions <- function(choices, block_size = 20) {
  n_trials <- ncol(choices)
  n_blocks <- floor(n_trials / block_size)

  block_props <- list()

  for (b in 1:n_blocks) {
    block_start <- (b - 1) * block_size + 1
    block_end <- b * block_size
    block_choices <- choices[, block_start:block_end]

    props <- compute_choice_proportions(block_choices)

    block_props[[b]] <- data.frame(
      block = b,
      deck = 1:4,
      proportion = as.numeric(props)
    )
  }

  return(do.call(rbind, block_props))
}

#' Run posterior predictive check for a model
#' @param fit_result Fitted model object
#' @param observed_data JAGS data list with observed choices and outcomes
#' @param model_name Model name ("pvl_delta", "orl", "eef")
#' @param n_sim Number of posterior samples to use for simulation
#' @return List with PPC results
run_ppc <- function(fit_result, observed_data, model_name, n_sim = 100) {
  cat(sprintf("\n=== Posterior Predictive Check: %s ===\n", model_name))

  # Extract posterior samples (handle both jags and rjags output formats)
  if (!is.null(fit_result$BUGSoutput$sims.matrix)) {
    samples_matrix <- fit_result$BUGSoutput$sims.matrix
  } else if (!is.null(fit_result$samples)) {
    samples_matrix <- as.matrix(fit_result$samples)
  } else {
    stop("Cannot extract posterior samples from fit object")
  }
  n_posterior <- nrow(samples_matrix)

  # Sample from posterior
  posterior_idx <- sample(1:n_posterior, size = n_sim, replace = FALSE)

  # Observed choices
  obs_choices <- observed_data$choice

  # Store simulated data
  sim_props_list <- list()

  for (i in 1:n_sim) {
    idx <- posterior_idx[i]

    # Extract parameters for this posterior sample
    if (model_name == "pvl_delta") {
      params <- list(
        a = samples_matrix[idx, grep("^a\\[", colnames(samples_matrix))],
        A = samples_matrix[idx, grep("^A\\[", colnames(samples_matrix))],
        theta = samples_matrix[idx, grep("^theta\\[", colnames(samples_matrix))],
        w = samples_matrix[idx, grep("^w\\[", colnames(samples_matrix))]
      )

      # Simulate
      sim_choices <- simulate_pvl_delta(params, observed_data$outcome, observed_data$Tsubj)
    } else if (model_name == "orl") {
      params <- list(
        Arew = samples_matrix[idx, grep("^a_rew\\[", colnames(samples_matrix))],
        Apun = samples_matrix[idx, grep("^a_pun\\[", colnames(samples_matrix))],
        K = samples_matrix[idx, grep("^K\\[", colnames(samples_matrix))],
        betaF = samples_matrix[idx, grep("^omega_f\\[", colnames(samples_matrix))],
        betaP = samples_matrix[idx, grep("^omega_p\\[", colnames(samples_matrix))]
      )

      # Simulate
      sim_choices <- simulate_orl(params, observed_data$outcome, observed_data$Tsubj)
    } else if (model_name == "eef") {
      params <- list(
        theta = samples_matrix[idx, grep("^theta\\[", colnames(samples_matrix))],
        lambda = samples_matrix[idx, grep("^lambda\\[", colnames(samples_matrix))],
        phi = samples_matrix[idx, grep("^phi\\[", colnames(samples_matrix))],
        cons = samples_matrix[idx, grep("^cons\\[", colnames(samples_matrix))]
      )

      # Handle w_ini
      w_ini <- if (!is.null(observed_data$w_ini)) observed_data$w_ini else rep(0, 4)

      # Handle trials
      n_trials_vec <- if (!is.null(observed_data$Tsubj)) observed_data$Tsubj else observed_data$T

      # Simulate
      sim_choices <- simulate_eef(params, observed_data$outcome, n_trials_vec, w_ini)
    } else {
      warning(sprintf("PPC not implemented for model: %s", model_name))
      return(NULL)
    }

    # Compute choice proportions
    sim_props_list[[i]] <- compute_choice_proportions(sim_choices)
  }

  # Aggregate simulated proportions
  sim_props_matrix <- do.call(rbind, sim_props_list)
  sim_props_mean <- colMeans(sim_props_matrix)
  sim_props_lower <- apply(sim_props_matrix, 2, quantile, probs = 0.025)
  sim_props_upper <- apply(sim_props_matrix, 2, quantile, probs = 0.975)

  # Observed proportions
  obs_props <- compute_choice_proportions(obs_choices)

  # Summary
  ppc_summary <- data.frame(
    deck = 1:4,
    observed = as.numeric(obs_props),
    predicted_mean = sim_props_mean,
    predicted_lower = sim_props_lower,
    predicted_upper = sim_props_upper
  )

  cat("\nPosterior Predictive Check Summary:\n")
  print(round(ppc_summary, 3))

  # Check if observed falls within credible interval
  within_ci <- ppc_summary$observed >= ppc_summary$predicted_lower &
    ppc_summary$observed <= ppc_summary$predicted_upper

  cat(sprintf("\nObserved within 95%% CI: %d/4 decks\n", sum(within_ci)))

  return(list(
    summary = ppc_summary,
    simulated_props = sim_props_matrix,
    within_ci = within_ci
  ))
}

#' Run PPC for all fitted models
#' @param models_dir Directory containing fitted models
#' @return List of PPC results
run_all_ppc <- function(models_dir = "analysis/outputs") {
  # Load data
  source("utils/load_data.R")
  source("utils/prepare_jags_data.R")
  source("utils/prepare_eef_data.R")

  all_data <- load_all_igt_data()

  # Find fitted models
  fit_files <- list.files(models_dir, pattern = "_fit\\.rds$", full.names = TRUE)

  ppc_results <- list()

  for (fit_file in fit_files) {
    model_name <- gsub("_fit\\.rds$", "", basename(fit_file))

    cat(sprintf("\nRunning PPC for %s...\n", model_name))
    fit <- readRDS(fit_file)

    # Prepare model-specific data
    if (model_name == "eef") {
      model_data <- prepare_eef_jags_data_pooled(all_data)
    } else {
      model_data <- prepare_jags_data_for_model(all_data, model_name)
    }

    # Run PPC (only for models with simulation functions)
    if (model_name %in% c("pvl_delta", "orl", "eef")) {
      ppc_results[[model_name]] <- run_ppc(fit, model_data, model_name)
    }
  }

  # Save results
  output_file <- file.path(models_dir, "ppc_results.rds")
  saveRDS(ppc_results, output_file)
  cat(sprintf("\nPPC results saved to: %s\n", output_file))

  return(ppc_results)
}
