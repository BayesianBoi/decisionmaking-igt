# Prepare data for JAGS model fitting
# Run: source("analysis/utils/prepare_jags_data.R")

#' Prepare JAGS data list from harmonized IGT data
#' @param dat Data frame with harmonized IGT data (from load_data.R)
#' @param study_filter Optional vector of study names to include
#' @return List formatted for JAGS input
prepare_jags_data <- function(dat, study_filter = NULL) {
  # Filter by study if requested
  if (!is.null(study_filter)) {
    dat <- dat[dat$study %in% study_filter, ]
  }

  # Get unique subjects (Respect subj_idx if present for consistent ordering)
  if ("subj_idx" %in% names(dat)) {
    subj_list <- unique(dat$subj_unique[order(dat$subj_idx)])
  } else {
    subj_list <- unique(dat$subj_unique)
  }
  n_subj <- length(subj_list)

  # Get maximum number of trials
  n_trials_per_subj <- sapply(subj_list, function(s) {
    sum(dat$subj_unique == s)
  })
  max_trials <- max(n_trials_per_subj)

  # Initialize arrays
  choice_mat <- matrix(NA, nrow = n_subj, ncol = max_trials)
  reward_mat <- matrix(NA, nrow = n_subj, ncol = max_trials)
  loss_mat <- matrix(NA, nrow = n_subj, ncol = max_trials)
  net_mat <- matrix(NA, nrow = n_subj, ncol = max_trials)

  # Fill arrays
  for (i in seq_along(subj_list)) {
    subj_id <- subj_list[i]
    subj_data <- dat[dat$subj_unique == subj_id, ]
    subj_data <- subj_data[order(subj_data$trial), ]

    n_trials <- nrow(subj_data)
    choice_mat[i, 1:n_trials] <- subj_data$choice
    reward_mat[i, 1:n_trials] <- subj_data$gain
    loss_mat[i, 1:n_trials] <- subj_data$loss
    net_mat[i, 1:n_trials] <- subj_data$net
  }

  # Create JAGS data list
  # Note: Only including variables used by all models
  # T, reward, and loss are computed but only needed for specific models
  jags_data <- list(
    N = n_subj,
    Tsubj = n_trials_per_subj,
    choice = choice_mat,
    outcome = net_mat / 100 # Scaled outcomes for numerical stability
  )

  return(jags_data)
}

#' Prepare JAGS data by individual study
#' @param dat Data frame with harmonized IGT data
#' @return Named list of JAGS data lists, one per study
prepare_jags_data_by_study <- function(dat) {
  studies <- unique(dat$study)
  jags_data_list <- list()

  for (study in studies) {
    jags_data_list[[study]] <- prepare_jags_data(dat, study_filter = study)
  }

  return(jags_data_list)
}

#' Create sign-coded outcome matrix for ORL model
#' @param outcome_mat Matrix of net outcomes
#' @return Matrix with values coded as 1 (positive), 0 (zero), -1 (negative)
create_sign_outcome <- function(outcome_mat) {
  sign_mat <- matrix(0, nrow = nrow(outcome_mat), ncol = ncol(outcome_mat))
  sign_mat[outcome_mat > 0] <- 1
  sign_mat[outcome_mat < 0] <- -1
  return(sign_mat)
}

#' Prepare JAGS data with additional variables for specific models
#' @param dat Data frame with harmonized IGT data
#' @param model One of \"pvl_delta\", \"orl\", \"eef\"
#' @param study_filter Optional vector of study names to include
#' @return List formatted for JAGS input with model-specific variables
prepare_jags_data_for_model <- function(dat, model, study_filter = NULL) {
  # Get base JAGS data
  jags_data <- prepare_jags_data(dat, study_filter)

  # Add model-specific variables
  if (model == "orl") {
    # ORL needs sign-coded outcomes
    jags_data$sign_outcome <- create_sign_outcome(jags_data$outcome)
  }

  return(jags_data)
}

#' Check JAGS data integrity
#' @param jags_data List with JAGS data
#' @return TRUE if valid, stops with error otherwise
check_jags_data <- function(jags_data) {
  # Get max trials from choice matrix
  max_trials <- ncol(jags_data$choice)

  # Check dimensions
  if (jags_data$N != length(jags_data$Tsubj)) {
    stop("Number of subjects doesn't match Tsubj length")
  }

  if (jags_data$N != nrow(jags_data$choice)) {
    stop("Number of subjects doesn't match choice matrix rows")
  }

  # Check choice values
  valid_choices <- jags_data$choice[!is.na(jags_data$choice)]
  if (any(valid_choices < 1 | valid_choices > 4)) {
    stop("Invalid choice values (must be 1-4)")
  }

  # Check for NAs in appropriate places
  for (i in 1:jags_data$N) {
    n_trials <- jags_data$Tsubj[i]

    # Should have data for trials 1:n_trials
    if (any(is.na(jags_data$choice[i, 1:n_trials]))) {
      warning(sprintf("Missing choice data for subject %d within Tsubj range", i))
    }

    # Should be NA after n_trials
    if (n_trials < max_trials) {
      if (!all(is.na(jags_data$choice[i, (n_trials + 1):max_trials]))) {
        warning(sprintf("Non-NA data for subject %d beyond Tsubj", i))
      }
    }
  }

  message(sprintf(
    "JAGS data check passed: N=%d, max_T=%d, total observations=%d",
    jags_data$N, max_trials, sum(jags_data$Tsubj)
  ))

  return(TRUE)
}
