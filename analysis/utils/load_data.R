# Load and harmonize IGT data from multiple sources
# Run: source("analysis/utils/load_data.R")
#
# Datasets: Ahn et al. (2014)
# Note: Fridberg (2010) and other datasets removed from analysis.

#' Load Ahn et al. (2014) data
#' @param group One of "HC", "Amph", "Hero"
#' @return Data frame with standardized columns
load_ahn_2014 <- function(group = "HC") {
  file_path <- sprintf("data/raw/Ahn_2014/IGTdata_%s.txt", group)
  dat <- read.table(file_path, header = TRUE)

  # Standardize column names
  dat_std <- data.frame(
    subj = dat$subjID,
    trial = dat$trial,
    choice = dat$deck,
    gain = dat$gain,
    loss = dat$loss,
    net = dat$gain - abs(dat$loss),
    study = paste0("Ahn2014_", group)
  )

  return(dat_std)
}

#' Load all IGT datasets
#' @return Data frame with all studies combined
load_all_igt_data <- function() {
  # Load all datasets (Ahn 2014)
  datasets <- list(
    load_ahn_2014("HC"),
    load_ahn_2014("Amph"),
    load_ahn_2014("Hero")
  )

  # Remove NULL entries (missing datasets)
  datasets <- datasets[!sapply(datasets, is.null)]

  if (length(datasets) == 0) {
    stop("No data could be loaded! Check data directory.")
  }

  # Combine all datasets
  all_data <- do.call(rbind, datasets)

  # Create unique subject IDs across studies
  all_data$subj_unique <- paste(all_data$study, all_data$subj, sep = "_")

  return(all_data)
}

#' Add contiguous integer subject index (Safe for JAGS)
#' @param dat Data frame with subj_unique
#' @return Data frame with 'subj_idx' column
add_subject_index <- function(dat) {
  dat$subj_idx <- as.integer(factor(dat$subj_unique, levels = unique(dat$subj_unique)))
  return(dat)
}

#' Validate IGT data
#' @param dat Data frame with IGT data
#' @return TRUE if valid, stops with error otherwise
validate_igt_data <- function(dat) {
  # Check required columns
  required_cols <- c("subj", "trial", "choice", "gain", "loss", "net", "study")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Check choice range (should be 1-4 for decks)
  if (any(dat$choice < 1 | dat$choice > 4, na.rm = TRUE)) {
    stop("Choice values must be between 1 and 4")
  }

  # Check for missing values
  if (any(is.na(dat$choice))) {
    warning("Missing values detected in choice column")
  }

  # Check trial numbering
  for (subj_id in unique(dat$subj)) {
    subj_data <- dat[dat$subj == subj_id, ]
    if (!all(subj_data$trial == seq_along(subj_data$trial))) {
      warning(sprintf("Trial numbering issue for subject %s", subj_id))
    }
  }

  message(sprintf(
    "Data validation passed: %d subjects, %d total trials",
    length(unique(dat$subj_unique)), nrow(dat)
  ))

  return(TRUE)
}
