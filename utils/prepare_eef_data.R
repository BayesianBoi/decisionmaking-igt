# Prepare data specifically for EEF model with first-choice priors
# Run: source("utils/prepare_eef_data.R")

#' Calculate first-choice priors from data (Yang et al., 2025 Section 2.4)
#'
#' Yang's EEF model uses empirical first-choice frequencies to initialize
#' deck weights (w_ini). This reflects participants' initial biases before
#' any learning has occurred.
#'
#' @param dat Data frame with harmonized IGT data (from load_data.R)
#' @param beta Inverse temperature for reverse softmax (default = 3, following Yang)
#' @return Vector of 4 initial weights for decks A, B, C, D
#' @details
#' Process:
#' 1. Count first-choice frequencies for each deck
#' 2. Convert to probabilities
#' 3. Reverse softmax: w_ini = log(p) / beta
#'
#' This ensures model starts with realistic deck preferences
#' rather than assuming uniform initial values (0,0,0,0)
calculate_first_choice_priors <- function(dat, beta = 3) {

  # Get first choice for each unique subject
  subj_list <- unique(dat$subj_unique)
  first_choices <- sapply(subj_list, function(s) {
    subj_data <- dat[dat$subj_unique == s, ]
    subj_data <- subj_data[order(subj_data$trial), ]
    return(subj_data$choice[1])
  })

  # Count frequencies
  freq_table <- table(factor(first_choices, levels = 1:4))
  freq_props <- as.numeric(freq_table) / sum(freq_table)

  # Add small constant to avoid log(0)
  freq_props <- pmax(freq_props, 1e-6)

  # Reverse softmax transformation
  w_ini <- log(freq_props) / beta

  # Report frequencies for transparency
  message("First-choice frequencies:")
  message(sprintf("  Deck A: %.1f%% (n=%d)", freq_props[1]*100, freq_table[1]))
  message(sprintf("  Deck B: %.1f%% (n=%d)", freq_props[2]*100, freq_table[2]))
  message(sprintf("  Deck C: %.1f%% (n=%d)", freq_props[3]*100, freq_table[3]))
  message(sprintf("  Deck D: %.1f%% (n=%d)", freq_props[4]*100, freq_table[4]))
  message(sprintf("Initial weights: [%.4f, %.4f, %.4f, %.4f]",
                  w_ini[1], w_ini[2], w_ini[3], w_ini[4]))

  # Statistical test: Are first choices non-uniform?
  chisq_test <- chisq.test(freq_table)
  message(sprintf("Chi-square test for uniformity: X2=%.2f, p=%.2e",
                  chisq_test$statistic, chisq_test$p.value))

  return(w_ini)
}


#' Prepare JAGS data for EEF model with group-level structure
#'
#' This function prepares data for hierarchical Bayesian modeling
#' with group comparisons (HC vs. substance users)
#'
#' @param dat Data frame with harmonized IGT data
#' @param group_var Name of column indicating group membership (e.g., "group")
#' @param reference_group Name of reference group (e.g., "HC")
#' @param study_filter Optional vector of study names to include
#' @return List formatted for JAGS input with group structure
prepare_eef_jags_data <- function(dat,
                                   group_var = NULL,
                                   reference_group = NULL,
                                   study_filter = NULL) {

  # Filter by study if requested
  if (!is.null(study_filter)) {
    dat <- dat[dat$study %in% study_filter, ]
  }

  # Calculate first-choice priors from full dataset
  w_ini <- calculate_first_choice_priors(dat)

  # Get unique subjects
  subj_list <- unique(dat$subj_unique)
  n_subj <- length(subj_list)

  # Get trial counts
  n_trials_per_subj <- sapply(subj_list, function(s) {
    sum(dat$subj_unique == s)
  })
  max_trials <- max(n_trials_per_subj)

  # Initialize arrays
  choice_mat <- matrix(NA, nrow = n_subj, ncol = max_trials)
  outcome_mat <- matrix(NA, nrow = n_subj, ncol = max_trials)

  # Fill arrays
  for (i in seq_along(subj_list)) {
    subj_id <- subj_list[i]
    subj_data <- dat[dat$subj_unique == subj_id, ]
    subj_data <- subj_data[order(subj_data$trial), ]

    n_trials <- nrow(subj_data)
    choice_mat[i, 1:n_trials] <- subj_data$choice
    # Scale outcomes for numerical stability (cents -> dollars)
    outcome_mat[i, 1:n_trials] <- subj_data$net / 100
  }

  # Create base JAGS data list
  jags_data <- list(
    N = n_subj,
    Tsubj = n_trials_per_subj,
    choice = choice_mat,
    outcome = outcome_mat,
    w_ini = w_ini  # First-choice priors
  )

  # Add group structure if requested
  if (!is.null(group_var)) {

    if (!group_var %in% names(dat)) {
      stop(sprintf("Group variable '%s' not found in data", group_var))
    }

    # Get group membership for each subject
    group_membership <- sapply(subj_list, function(s) {
      unique(dat[[group_var]][dat$subj_unique == s])
    })

    # Convert to numeric codes
    group_levels <- unique(group_membership)
    if (!is.null(reference_group)) {
      # Put reference group first
      group_levels <- c(reference_group, setdiff(group_levels, reference_group))
    }
    group_codes <- as.numeric(factor(group_membership, levels = group_levels))

    # Add to JAGS data
    jags_data$n_groups <- length(group_levels)
    jags_data$group <- group_codes
    jags_data$group_names <- group_levels

    # Report group sizes
    message("\nGroup structure:")
    for (g in 1:length(group_levels)) {
      n_in_group <- sum(group_codes == g)
      message(sprintf("  %s (code %d): n=%d", group_levels[g], g, n_in_group))
    }
  }

  return(jags_data)
}


#' Prepare JAGS data for pooled analysis (all groups combined)
#'
#' Use this for initial model validation before running group comparisons
#'
#' @param dat Data frame with harmonized IGT data
#' @param study_filter Optional vector of study names to include
#' @return List formatted for JAGS input (no group structure)
prepare_eef_jags_data_pooled <- function(dat, study_filter = NULL) {
  prepare_eef_jags_data(dat,
                        group_var = NULL,
                        study_filter = study_filter)
}


#' Compare first-choice patterns between groups
#'
#' This is useful for checking whether different populations
#' (HC vs. substance users) show different initial deck preferences
#'
#' @param dat Data frame with harmonized IGT data
#' @param group_var Name of column indicating group membership
#' @return Data frame with first-choice proportions by group
compare_first_choices_by_group <- function(dat, group_var = "group") {

  if (!group_var %in% names(dat)) {
    stop(sprintf("Group variable '%s' not found in data", group_var))
  }

  # Get first choice for each subject with group label
  subj_list <- unique(dat$subj_unique)
  first_choice_data <- data.frame(
    subj = subj_list,
    group = sapply(subj_list, function(s) {
      unique(dat[[group_var]][dat$subj_unique == s])
    }),
    first_choice = sapply(subj_list, function(s) {
      subj_data <- dat[dat$subj_unique == s, ]
      subj_data <- subj_data[order(subj_data$trial), ]
      return(subj_data$choice[1])
    })
  )

  # Create contingency table
  contingency_table <- table(first_choice_data$group, first_choice_data$first_choice)

  # Convert to proportions
  prop_table <- prop.table(contingency_table, margin = 1)

  # Statistical test
  if (nrow(contingency_table) > 1) {
    chisq_test <- chisq.test(contingency_table)
    message("\nChi-square test for group differences in first choice:")
    message(sprintf("  X2(%.0f) = %.2f, p = %.4f",
                    chisq_test$parameter,
                    chisq_test$statistic,
                    chisq_test$p.value))

    if (chisq_test$p.value < 0.05) {
      message("  => Groups show DIFFERENT first-choice patterns")
      message("  => Consider using GROUP-SPECIFIC w_ini priors")
    } else {
      message("  => Groups show SIMILAR first-choice patterns")
      message("  => Can use pooled w_ini priors")
    }
  }

  # Print formatted table
  message("\nFirst-choice proportions by group:")
  print(round(prop_table * 100, 1))

  return(as.data.frame.matrix(prop_table))
}


#' Calculate group-specific first-choice priors
#'
#' Use this if groups show significantly different first-choice patterns
#' Returns separate w_ini vectors for each group
#'
#' @param dat Data frame with harmonized IGT data
#' @param group_var Name of column indicating group membership
#' @param beta Inverse temperature for reverse softmax
#' @return Named list of w_ini vectors, one per group
calculate_first_choice_priors_by_group <- function(dat,
                                                     group_var = "group",
                                                     beta = 3) {

  groups <- unique(dat[[group_var]])
  w_ini_list <- list()

  for (g in groups) {
    message(sprintf("\n--- Group: %s ---", g))
    dat_group <- dat[dat[[group_var]] == g, ]
    w_ini_list[[g]] <- calculate_first_choice_priors(dat_group, beta)
  }

  return(w_ini_list)
}


#' Check EEF JAGS data integrity
#'
#' Validates that data is correctly formatted for EEF model
#'
#' @param jags_data List with JAGS data
#' @return TRUE if valid, stops with error otherwise
check_eef_jags_data <- function(jags_data) {

  # Basic dimension checks
  if (jags_data$N != length(jags_data$Tsubj)) {
    stop("Number of subjects doesn't match Tsubj length")
  }

  if (jags_data$N != nrow(jags_data$choice)) {
    stop("Number of subjects doesn't match choice matrix rows")
  }

  # Check w_ini exists and has 4 elements
  if (is.null(jags_data$w_ini)) {
    stop("w_ini (first-choice priors) missing from JAGS data")
  }

  if (length(jags_data$w_ini) != 4) {
    stop("w_ini must have exactly 4 elements (one per deck)")
  }

  # Check choice values
  valid_choices <- jags_data$choice[!is.na(jags_data$choice)]
  if (any(valid_choices < 1 | valid_choices > 4)) {
    stop("Invalid choice values (must be 1-4)")
  }

  # Check outcome scaling (should be roughly -2 to +2 after /100 scaling)
  valid_outcomes <- jags_data$outcome[!is.na(jags_data$outcome)]
  outcome_range <- range(valid_outcomes)
  if (outcome_range[1] < -10 | outcome_range[2] > 10) {
    warning("Outcomes may not be properly scaled (expected range: -2 to +2)")
    message(sprintf("  Actual range: [%.2f, %.2f]", outcome_range[1], outcome_range[2]))
  }

  # Check group structure if present
  if (!is.null(jags_data$group)) {
    if (length(jags_data$group) != jags_data$N) {
      stop("Group vector length doesn't match number of subjects")
    }

    if (max(jags_data$group) != jags_data$n_groups) {
      stop("Group codes don't match n_groups")
    }

    message(sprintf("Group structure validated: %d groups, %d subjects",
                    jags_data$n_groups, jags_data$N))
  }

  message(sprintf("\nEEF JAGS data validated successfully:"))
  message(sprintf("  N = %d subjects", jags_data$N))
  message(sprintf("  Max trials = %d", ncol(jags_data$choice)))
  message(sprintf("  Total observations = %d", sum(jags_data$Tsubj)))
  message(sprintf("  w_ini = [%.4f, %.4f, %.4f, %.4f]",
                  jags_data$w_ini[1], jags_data$w_ini[2],
                  jags_data$w_ini[3], jags_data$w_ini[4]))

  return(TRUE)
}
