#!/usr/bin/env Rscript
# ==============================================================================
# ==============================================================================
# FIT ORL MODEL (JAGS)
# ==============================================================================
# Fits the Outcome-Representation Learning (ORL) model to IGT data using JAGS.
# - Hierarchical Bayesian estimation
# - Group-level parameters
# - Outcomes scaled by 1/100
# Usage: Rscript fit_orl.R [HC|Amph|Hero]
# ==============================================================================

# 1. Setup
# ------------------------------------------------------------------------------
# Dependencies
required_packages <- c("R2jags", "parallel")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")
invisible(lapply(required_packages, library, character.only = TRUE))

set.seed(69420)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Must specify group argument: HC, Amph, or Hero")
}
target_group <- args[1]

groups_map <- c(
  "HC" = "Ahn2014_HC",
  "Amph" = "Ahn2014_Amph",
  "Hero" = "Ahn2014_Hero"
)

if (!target_group %in% names(groups_map)) {
  stop("Invalid group. Must be one of: ", paste(names(groups_map), collapse = ", "))
}
study_label <- groups_map[[target_group]]

cat("\n================================================\n")
cat("PROCESSING GROUP:", target_group, "(", study_label, ")\n")
cat("================================================\n")

# 2. Load Data
# ------------------------------------------------------------------------------
source("utils/load_data.R")
all_data <- load_all_igt_data()

# Filter data for this group
raw_data <- all_data[all_data$study == study_label, ]

# 3. Data Preparation
# ------------------------------------------------------------------------------
subIDs <- unique(raw_data$subj)
nsubs <- length(subIDs)
ntrials_max <- 100

ntrials_all <- array(0, c(nsubs))
x_all <- array(0, c(nsubs, ntrials_max))
X_all <- array(0, c(nsubs, ntrials_max))

cat("Preparing data for", nsubs, "subjects...\n")

for (s in 1:nsubs) {
  subj_df <- raw_data[raw_data$subj == subIDs[s], ]
  ntrials_all[s] <- nrow(subj_df)

  # Pad choices (x)
  x_sub <- subj_df$choice
  length(x_sub) <- ntrials_max

  # Pad outcomes (X)
  X_sub <- subj_df$gain + subj_df$loss
  length(X_sub) <- ntrials_max

  # Assign to arrays
  x_all[s, ] <- x_sub
  X_all[s, ] <- X_sub
}

# SCALING: Divide outcomes by 100
cat("Scaling outcomes (X / 100)...\n")
X_all <- X_all / 100

# 4. JAGS Fitting
# ------------------------------------------------------------------------------
output_dir <- "outputs/orl"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

model_file <- "models/orl.txt"

params <- c(
  "mu_a_rew", "mu_a_pun", "mu_K", "mu_omega_f", "mu_omega_p",
  "lambda_a_rew", "lambda_a_pun", "lambda_K", "lambda_omega_f", "lambda_omega_p",
  "a_rew", "a_pun", "K", "omega_f", "omega_p"
)

cat("Fitting hierarchical ORL model for group:", target_group, "...\n")

jags_data <- list(
  "x" = x_all,
  "X" = X_all,
  "ntrials" = ntrials_all,
  "nsubs" = nsubs
)

start_time <- Sys.time()

fit <- jags.parallel(
  data = jags_data,
  inits = NULL,
  parameters.to.save = params,
  model.file = model_file,
  n.chains = 3,
  n.iter = 50000,
  n.burnin = 10000,
  n.thin = 5
)

end_time <- Sys.time()
cat("Fitting complete. Duration:", end_time - start_time, "\n")

# 5. Save Results
# ------------------------------------------------------------------------------
save_path <- file.path(output_dir, paste0("orl_fit_", target_group, ".rds"))
saveRDS(fit, file = save_path)
cat("Results saved to:", save_path, "\n")

# 6. Basic Diagnostics
# ------------------------------------------------------------------------------
print(fit)
