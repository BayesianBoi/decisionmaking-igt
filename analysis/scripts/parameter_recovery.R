# ===========================================================================
# Parameter Recovery Analysis (All Models)
# ===========================================================================
#
# Simulates data from known parameters for PVL-Delta, ORL, and EEF models,
# then re-fits the models to verify parameters can be reliably recovered.
#
# Run: Rscript analysis/scripts/parameter_recovery.R
#
# ===========================================================================

library(rjags)
library(coda)
library(parallel)
library(ggplot2)

source("analysis/utils/load_data.R")
source("analysis/utils/prepare_jags_data.R")
source("analysis/utils/prepare_eef_data.R")

# ===========================================================================
# CONFIGURATION
# ===========================================================================

n_subjects <- 30 # Subjects per model
n_trials <- 100 # Standard IGT length
n_chains <- 3
n_adapt <- 1000
n_burnin <- 2000
n_iter <- 5000
thin <- 2

output_dir <- "results/parameter_recovery"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(123)

# Payoff scheme (Standard IGT A,B,C,D)
# A: -250 net (freq), B: -250 net (infreq), C: +250 net (freq), D: +250 net (infreq)
# Simplified approximations for simulation:
payoffs <- list(
  A = list(wins = 100, losses = -250, freq = 0.5), # Bad
  B = list(wins = 100, losses = -1250, freq = 0.1), # Bad
  C = list(wins = 50, losses = -50, freq = 0.5), # Good
  D = list(wins = 50, losses = -250, freq = 0.1) # Good
)

# ===========================================================================
# SIMULATION FUNCTIONS
# ===========================================================================

# --- PVL-Delta Simulation ---
simulate_pvl <- function(n_subj, n_trials) {
  # True params: A (learning), alpha (shape), cons (consistency), lambda (loss av)
  params <- data.frame(
    A = runif(n_subj, 0.1, 0.5), # Learning rate
    alpha = runif(n_subj, 0.3, 0.7), # Curve shape
    cons = runif(n_subj, 1.0, 3.0), # Consistency
    lambda = runif(n_subj, 1.0, 3.0) # Loss aversion
  )

  choice_mat <- matrix(NA, n_subj, n_trials)
  outcome_mat <- matrix(NA, n_subj, n_trials)

  for (s in 1:n_subj) {
    Ev <- rep(0, 4) # Expected values
    p <- params[s, ]

    for (t in 1:n_trials) {
      # Choice prob (Softmax)
      theta <- 3^p$cons - 1
      prob <- exp(theta * Ev) / sum(exp(theta * Ev))
      ch <- sample(1:4, 1, prob = prob)
      choice_mat[s, t] <- ch

      # Outcome
      deck <- LETTERS[ch]
      win <- payoffs[[deck]]$wins
      is_loss <- runif(1) < payoffs[[deck]]$freq
      loss <- if (is_loss) abs(payoffs[[deck]]$losses) else 0
      net <- win - loss
      outcome_mat[s, t] <- net

      # Update
      util <- ifelse(net >= 0, net^p$alpha, -p$lambda * (abs(net)^p$alpha))
      Ev[ch] <- (1 - p$A) * Ev[ch] + p$A * util
    }
  }
  return(list(choice = choice_mat, outcome = outcome_mat, params = params))
}

# --- ORL Simulation ---
# --- ORL Simulation ---
simulate_orl <- function(n_subj, n_trials) {
  # A_rew, A_pun, K (decay), betaF (freq weight), betaP (pers)
  params <- data.frame(
    A_rew = runif(n_subj, 0.1, 0.5),
    A_pun = runif(n_subj, 0.05, 0.3),
    K = runif(n_subj, 0.1, 0.5),
    betaF = runif(n_subj, 1.0, 3.0),
    betaP = runif(n_subj, 0, 1.0)
  )

  choice_mat <- matrix(NA, n_subj, n_trials)
  outcome_mat <- matrix(NA, n_subj, n_trials)
  sign_mat <- matrix(NA, n_subj, n_trials)

  for (s in 1:n_subj) {
    Ev <- rep(0, 4)
    Ef <- rep(0, 4)
    Pers <- rep(0, 4)
    p <- params[s, ]

    for (t in 1:n_trials) {
      # Total Utility (Theta fixed at 1 per Haines 2018)
      V <- Ev + p$betaF * Ef + p$betaP * Pers
      prob <- exp(V) / sum(exp(V))
      ch <- sample(1:4, 1, prob = prob)
      choice_mat[s, t] <- ch

      # Outcome
      deck <- LETTERS[ch]
      win <- payoffs[[deck]]$wins
      is_loss <- runif(1) < payoffs[[deck]]$freq
      loss <- if (is_loss) abs(payoffs[[deck]]$losses) else 0
      net <- win - loss
      outcome_mat[s, t] <- net
      sign_val <- ifelse(net > 0, 1, ifelse(net < 0, -1, 0))
      sign_mat[s, t] <- sign_val

      # --- Update EV ---
      # Asymmetric learning
      if (net >= 0) {
        Ev[ch] <- Ev[ch] + p$A_rew * (net - Ev[ch])
      } else {
        Ev[ch] <- Ev[ch] + p$A_pun * (net - Ev[ch])
      }

      # --- Update EF (Expected Frequency) ---
      # Chosen deck: update toward sign
      if (net >= 0) {
        Ef[ch] <- Ef[ch] + p$A_rew * (sign_val - Ef[ch])
      } else {
        Ef[ch] <- Ef[ch] + p$A_pun * (sign_val - Ef[ch])
      }

      # Unchosen decks: Counterfactual updating (-sign/3)
      # If win: unchosen updated with Apun towards negative
      # If loss: unchosen updated with Arew towards positive (Opposite)
      cf_sign <- -sign_val / 3
      for (d in 1:4) {
        if (d != ch) {
          if (net >= 0) {
            Ef[d] <- Ef[d] + p$A_pun * (cf_sign - Ef[d])
          } else {
            Ef[d] <- Ef[d] + p$A_rew * (cf_sign - Ef[d])
          }
        }
      }

      # --- Update Pers (Perseverance) ---
      # Decay rule: 1/(1+K)
      # Chosen: Set to 1 then decay. Unchosen: Keep then decay.
      for (d in 1:4) {
        if (d == ch) {
          Pers[d] <- 1 / (1 + p$K)
        } else {
          Pers[d] <- Pers[d] / (1 + p$K)
        }
      }
    }
  }
  return(list(choice = choice_mat, outcome = outcome_mat, sign = sign_mat, params = params))
}


# --- EEF Simulation ---
simulate_eef <- function(n_subj, n_trials) {
  params <- data.frame(
    id = 1:n_subj,
    theta = runif(n_subj, 0.1, 0.5),
    lambda = runif(n_subj, 0.1, 0.5), # lambda_forget
    phi = runif(n_subj, -1, 1),
    cons = runif(n_subj, 1.0, 3.0)
  )

  choice_mat <- matrix(NA, n_subj, n_trials)
  outcome_mat <- matrix(NA, n_subj, n_trials)

  for (s in 1:n_subj) {
    p <- params[s, ]
    exploit <- rep(0, 4)
    explore <- rep(0, 4)

    for (t in 1:n_trials) {
      w <- exploit + explore
      prob <- exp(p$cons * w) / sum(exp(p$cons * w))
      ch <- sample(1:4, 1, prob = prob)
      choice_mat[s, t] <- ch

      deck <- LETTERS[ch]
      win <- payoffs[[deck]]$wins
      is_loss <- runif(1) < payoffs[[deck]]$freq
      loss <- if (is_loss) abs(payoffs[[deck]]$losses) else 0
      net <- win - loss
      outcome_mat[s, t] <- net

      # Update
      util <- (max(net, 0) + 0.001)^p$theta - (max(-net, 0) + 0.001)^p$theta
      for (d in 1:4) {
        if (d == ch) {
          exploit[d] <- (1 - p$lambda) * exploit[d] + util
          explore[d] <- 0
        } else {
          exploit[d] <- (1 - p$lambda) * exploit[d]
          explore[d] <- p$lambda * explore[d] + (1 - p$lambda) * p$phi
        }
      }
    }
  }
  return(list(choice = choice_mat, outcome = outcome_mat, params = params))
}

# ===========================================================================
# RECOVERY LOOP
# ===========================================================================

models_to_test <- list(
  list(
    name = "pvl_delta", sim_func = simulate_pvl,
    params = c("mu_A", "mu_alpha", "mu_cons", "mu_lambda"),
    indiv_params = c("A", "alpha", "cons", "lambda"),
    file = "analysis/models/pvl_delta.jags"
  ),
  list(
    name = "orl", sim_func = simulate_orl,
    params = c("mu_Arew", "mu_Apun", "mu_K", "mu_betaF", "mu_betaP"),
    indiv_params = c("Arew", "Apun", "K", "betaF", "betaP"),
    file = "analysis/models/orl.jags"
  ),
  list(
    name = "eef", sim_func = simulate_eef,
    params = c("mu_theta", "mu_lambda_forget", "mu_phi", "mu_cons"),
    indiv_params = c("theta", "lambda_forget", "phi", "cons"),
    file = "analysis/models/eef.jags"
  )
)

all_results <- list()

for (m in models_to_test) {
  message(sprintf("\n=== RECOVERING %s ===", toupper(m$name)))

  # 1. Simulate
  sim <- m$sim_func(n_subjects, n_trials)

  # 2. Prepare JAGS Data
  # Use the standard prep function but mock the group
  mock_df <- data.frame(
    subj_unique = rep(1:n_subjects, each = n_trials),
    trial = rep(1:n_trials, n_subjects),
    choice = as.vector(t(sim$choice)),
    deck = LETTERS[as.vector(t(sim$choice))],
    outcome = as.vector(t(sim$outcome)),
    net = as.vector(t(sim$outcome)),
    gain = pmax(as.vector(t(sim$outcome)), 0),
    loss = pmax(-as.vector(t(sim$outcome)), 0),
    sign = if (!is.null(sim$sign)) as.vector(t(sim$sign)) else NA,
    study = "Sim",
    group = "Sim"
  )

  if (m$name == "eef") {
    jags_data <- prepare_eef_jags_data(mock_df, group_var = "group", reference_group = "Sim")
  } else {
    jags_data <- prepare_jags_data(mock_df)
  }

  # 3. Fit
  message("Fitting...")
  mod <- jags.model(m$file, jags_data, n.chains = n_chains, n.adapt = n_adapt, quiet = TRUE)
  update(mod, n_burnin, progress.bar = "none")

  # Monitor individual parameters too
  monitors <- c(m$params, m$indiv_params)
  samples <- coda.samples(mod, monitors, n_iter, thin = thin)

  # 4. Evaluate
  sum_stats <- summary(samples)$statistics

  res_df <- data.frame()

  # Check individual parameter recovery (Correlation)
  for (p in m$indiv_params) {
    # Extract posterior means for all subjects: e.g. theta[1], theta[2]...
    # Note: Regex needs to be precise
    # For params like "A" in PVL, need to avoid matching "mu_A"

    # Construct regex: e.g. "^A\["
    p_regex <- paste0("^", p, "\\[")
    est_rows <- grep(p_regex, rownames(sum_stats))
    estimated <- sum_stats[est_rows, "Mean"]

    # Handle parameter name mismatches for EEF
    true_p_name <- p
    if (m$name == "eef" && p == "lambda_forget") true_p_name <- "lambda"
    if (m$name == "orl" && p == "Arew") true_p_name <- "A_rew"
    if (m$name == "orl" && p == "Apun") true_p_name <- "A_pun"

    true_vals <- sim$params[[true_p_name]]

    # Check alignment
    if (length(estimated) == n_subjects) {
      cor_val <- cor(true_vals, estimated)
      rmse_val <- sqrt(mean((true_vals - estimated)^2))
      message(sprintf("  %s: r = %.2f, RMSE = %.2f", p, cor_val, rmse_val))

      res_df <- rbind(res_df, data.frame(
        Model = m$name,
        Parameter = p,
        True = true_vals,
        Recovered = estimated,
        Correlation = cor_val,
        RMSE = rmse_val
      ))
    }
  }
  all_results[[m$name]] <- res_df
}

# ===========================================================================
# PLOT
# ===========================================================================

message("\nGenerating Recovery Plots...")
combined_df <- do.call(rbind, all_results)

pdf(file.path(output_dir, "parameter_recovery.pdf"), width = 12, height = 10)

ggplot(combined_df, aes(x = True, y = Recovered)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  facet_wrap(~ Model + Parameter, scales = "free") +
  theme_minimal() +
  labs(
    title = "Parameter Recovery (All Models)",
    subtitle = "Simulated vs Recovered Parameters",
    x = "True Parameter Value",
    y = "Recovered Posterior Mean"
  )

dev.off()

saveRDS(all_results, file.path(output_dir, "recovery_results.rds"))
message("Recovery complete. Results saved.")
