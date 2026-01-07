# Parameter Recovery Analysis for IGT models
# Run: source("utils/parameter_recovery.R")

library(rjags)
library(coda)

#' Sample parameters from prior distributions
#' @param model_name Model name ("pvl_delta", "orl")
#' @param n_subj Number of subjects to simulate
#' @return List of sampled parameters
sample_from_prior <- function(model_name, n_subj = 20) {
  if (model_name == "pvl_delta") {
    # Sample group-level parameters (using reference naming: a, A, theta, w)
    mu_a <- runif(1, 0, 1) # Learning rate
    mu_A <- runif(1, 0, 2) # Outcome sensitivity
    mu_theta <- runif(1, 0, 5) # Choice consistency
    mu_w <- runif(1, 0, 10) # Loss aversion

    sigma_a <- runif(1, 0, 0.3)
    sigma_A <- runif(1, 0, 0.5)
    sigma_theta <- runif(1, 0, 1)
    sigma_w <- runif(1, 0, 2)

    # Sample subject-level parameters
    a <- pmin(pmax(rnorm(n_subj, mu_a, sigma_a), 0), 1)
    A <- pmax(rnorm(n_subj, mu_A, sigma_A), 0.01)
    theta <- pmax(rnorm(n_subj, mu_theta, sigma_theta), 0.01)
    w <- pmax(rnorm(n_subj, mu_w, sigma_w), 0.01)

    return(list(
      a = a, A = A, theta = theta, w = w,
      mu_a = mu_a, mu_A = mu_A, mu_theta = mu_theta, mu_w = mu_w
    ))
  } else if (model_name == "orl") {
    mu_Arew <- runif(1, 0, 1)
    mu_Apun <- runif(1, 0, 1)
    mu_K <- runif(1, 0, 5)
    mu_betaF <- runif(1, -5, 5)
    mu_betaP <- runif(1, -5, 5)

    sigma_Arew <- runif(1, 0, 0.3)
    sigma_Apun <- runif(1, 0, 0.3)
    sigma_K <- runif(1, 0, 1)
    sigma_betaF <- runif(1, 0, 2)
    sigma_betaP <- runif(1, 0, 2)

    Arew <- pmin(pmax(rnorm(n_subj, mu_Arew, sigma_Arew), 0), 1)
    Apun <- pmin(pmax(rnorm(n_subj, mu_Apun, sigma_Apun), 0), 1)
    K <- pmax(rnorm(n_subj, mu_K, sigma_K), 0.01)
    betaF <- rnorm(n_subj, mu_betaF, sigma_betaF)
    betaP <- rnorm(n_subj, mu_betaP, sigma_betaP)

    return(list(
      Arew = Arew, Apun = Apun, K = K, betaF = betaF, betaP = betaP,
      mu_Arew = mu_Arew, mu_Apun = mu_Apun, mu_K = mu_K,
      mu_betaF = mu_betaF, mu_betaP = mu_betaP
    ))
  }
}

#' Generate standard IGT reward structure
#' @param n_trials Number of trials
#' @return Array of deck outcomes (trial x outcome_type x deck)
generate_igt_outcomes <- function(n_trials = 100) {
  # Standard IGT deck structure
  # Deck A: Immediate reward 100, occasional large loss
  # Deck B: Immediate reward 100, occasional very large loss
  # Deck C: Immediate reward 50, occasional small loss
  # Deck D: Immediate reward 50, very occasional small loss

  deck_outcomes <- array(0, dim = c(n_trials, 2, 4))

  for (t in 1:n_trials) {
    # Deck A: +100, loss every 10 trials
    deck_outcomes[t, 1, 1] <- 100
    deck_outcomes[t, 2, 1] <- ifelse(t %% 10 == 0, -250, 0)

    # Deck B: +100, loss every 10 trials (larger)
    deck_outcomes[t, 1, 2] <- 100
    deck_outcomes[t, 2, 2] <- ifelse(t %% 10 == 0, -1250, 0)

    # Deck C: +50, loss every 10 trials
    deck_outcomes[t, 1, 3] <- 50
    deck_outcomes[t, 2, 3] <- ifelse(t %% 10 == 0, -50, 0)

    # Deck D: +50, loss every 10 trials (smaller)
    deck_outcomes[t, 1, 4] <- 50
    deck_outcomes[t, 2, 4] <- ifelse(t %% 10 == 0, -250, 0)
  }

  return(deck_outcomes)
}

#' Simulate data from PVL-Delta model (using reference naming: a, A, theta, w)
#' @param params Parameter list with a, A, theta, w
#' @param deck_outcomes IGT deck structure
#' @param n_trials Number of trials
#' @return List with simulated choices and outcomes
simulate_data_pvl_delta <- function(params, deck_outcomes, n_trials) {
  n_subj <- length(params$a)

  choices <- matrix(NA, nrow = n_subj, ncol = n_trials)
  outcomes <- matrix(NA, nrow = n_subj, ncol = n_trials)

  for (s in 1:n_subj) {
    ev <- rep(0, 4)

    for (t in 1:n_trials) {
      # Choice probabilities (matching JAGS model formulation)
      v <- params$theta[s] * ev
      exp_util <- exp(v)
      p_choice <- exp_util / sum(exp_util)

      # Sample choice
      choice <- sample(1:4, size = 1, prob = p_choice)
      choices[s, t] <- choice

      # Get outcome
      gain <- deck_outcomes[t, 1, choice]
      loss <- deck_outcomes[t, 2, choice]
      outcome <- (gain + loss) / 100 # Scale
      outcomes[s, t] <- outcome

      # Compute utility using A (sensitivity) and w (loss aversion)
      if (outcome >= 0) {
        util <- outcome^params$A[s]
      } else {
        util <- -params$w[s] * abs(outcome)^params$A[s]
      }

      # Update EV using a (learning rate)
      ev[choice] <- ev[choice] + params$a[s] * (util - ev[choice])
    }
  }

  return(list(choices = choices, outcomes = outcomes))
}

#' Simulate data from ORL model
#' @param params Parameter list
#' @param deck_outcomes IGT deck structure
#' @param n_trials Number of trials
#' @return List with simulated choices and outcomes
simulate_data_orl <- function(params, deck_outcomes, n_trials) {
  n_subj <- length(params$Arew)

  choices <- matrix(NA, nrow = n_subj, ncol = n_trials)
  outcomes <- matrix(NA, nrow = n_subj, ncol = n_trials)

  for (s in 1:n_subj) {
    ev <- rep(0, 4)
    ef <- rep(0, 4)
    pers <- rep(0, 4)

    for (t in 1:n_trials) {
      # Combined utility (matching JAGS formulation)
      util_combined <- ev + ef * params$betaF[s] + pers * params$betaP[s]
      exp_util <- exp(util_combined)
      p_choice <- exp_util / sum(exp_util)

      # Sample choice
      choice <- sample(1:4, size = 1, prob = p_choice)
      choices[s, t] <- choice

      # Get outcome
      gain <- deck_outcomes[t, 1, choice]
      loss <- deck_outcomes[t, 2, choice]
      outcome <- (gain + loss) / 100 # Scale
      outcomes[s, t] <- outcome

      # Sign of outcome
      sign_out <- ifelse(outcome > 0, 1, ifelse(outcome < 0, -1, 0))

      # Prediction errors
      PE_val <- outcome - ev[choice]
      PE_freq <- sign_out - ef[choice]

      # Learning rates
      lr_val <- ifelse(outcome >= 0, params$Arew[s], params$Apun[s])
      lr_freq <- ifelse(outcome >= 0, params$Arew[s], params$Apun[s])
      lr_fic <- ifelse(outcome >= 0, params$Apun[s], params$Arew[s])

      # Update chosen deck
      ev[choice] <- ev[choice] + lr_val * PE_val
      ef[choice] <- ef[choice] + lr_freq * PE_freq

      # Update unchosen decks (fictive updating)
      for (d in setdiff(1:4, choice)) {
        PE_freq_fic <- -sign_out / 3 - ef[d]
        ef[d] <- ef[d] + lr_fic * PE_freq_fic
      }

      # Update perseverance (matching JAGS formulation)
      pers[choice] <- 1 / (1 + params$K[s])
      for (d in setdiff(1:4, choice)) {
        pers[d] <- pers[d] / (1 + params$K[s])
      }
    }
  }

  return(list(choices = choices, outcomes = outcomes))
}

#' Run parameter recovery analysis
#' @param model_name Model name
#' @param n_subj Number of subjects to simulate
#' @param n_trials Number of trials per subject
#' @param output_dir Output directory
#' @return Parameter recovery results
run_parameter_recovery <- function(model_name = "pvl_delta",
                                   n_subj = 20,
                                   n_trials = 100,
                                   output_dir = "analysis/outputs") {
  cat(sprintf("=== Parameter Recovery: %s ===\n\n", model_name))

  # Step 1: Sample true parameters
  cat("Step 1: Sampling true parameters from prior...\n")
  true_params <- sample_from_prior(model_name, n_subj)

  # Step 2: Simulate data
  cat("Step 2: Simulating choice data...\n")
  deck_outcomes <- generate_igt_outcomes(n_trials)

  if (model_name == "pvl_delta") {
    sim_data <- simulate_data_pvl_delta(true_params, deck_outcomes, n_trials)
  } else if (model_name == "orl") {
    sim_data <- simulate_data_orl(true_params, deck_outcomes, n_trials)
  } else {
    warning(sprintf("Simulation not implemented for %s", model_name))
    return(NULL)
  }

  # Step 3: Prepare JAGS data
  cat("Step 3: Preparing data for JAGS...\n")
  jags_data <- list(
    N = n_subj,
    T = n_trials,
    Tsubj = rep(n_trials, n_subj),
    choice = sim_data$choices,
    outcome = sim_data$outcomes
  )

  # Step 4: Fit model to simulated data
  cat("Step 4: Fitting model to recover parameters...\n")
  model_file <- sprintf("analysis/models/%s_v2.jags", model_name)

  jags_model <- jags.model(
    file = model_file,
    data = jags_data,
    n.chains = 2,
    n.adapt = 1000,
    quiet = FALSE
  )

  update(jags_model, n.iter = 1000)

  # Monitor parameters
  if (model_name == "pvl_delta") {
    params_to_monitor <- c("a", "A", "theta", "w")
  } else if (model_name == "orl") {
    params_to_monitor <- c("Arew", "Apun", "K", "betaF", "betaP")
  }

  samples <- coda.samples(
    model = jags_model,
    variable.names = params_to_monitor,
    n.iter = 2000
  )

  # Step 5: Compare true vs recovered
  cat("Step 5: Comparing true vs. recovered parameters...\n")
  samples_matrix <- as.matrix(samples)

  recovery_results <- list()

  for (param in params_to_monitor) {
    # Extract recovered parameters (posterior means)
    param_cols <- grep(sprintf("^%s\\[", param), colnames(samples_matrix))
    recovered <- colMeans(samples_matrix[, param_cols])

    # True parameters
    true <- true_params[[param]]

    # Correlation
    correlation <- cor(true, recovered)

    recovery_results[[param]] <- data.frame(
      subject = 1:n_subj,
      true = true,
      recovered = recovered,
      error = recovered - true
    )

    cat(sprintf("  %s: r = %.3f\n", param, correlation))
  }

  # Save results
  output_file <- file.path(output_dir, sprintf("%s_parameter_recovery.rds", model_name))
  recovery_output <- list(
    true_params = true_params,
    simulated_data = sim_data,
    recovered_samples = samples,
    recovery_results = recovery_results
  )
  saveRDS(recovery_output, output_file)
  cat(sprintf("\nParameter recovery results saved to: %s\n", output_file))

  return(recovery_output)
}

# Run if executed as script
if (!interactive()) {
  recovery <- run_parameter_recovery("pvl_delta", n_subj = 10, n_trials = 100)
}
