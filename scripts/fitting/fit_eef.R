#!/usr/bin/env Rscript
# =============================================================================
# FIT EEF MODEL (JAGS)
# =============================================================================
# Fits the Exploitation-Exploration with Forgetting model (Yang et al. 2025)
#   - V(t) = Gain(t)^theta - Loss(t)^theta (separate gain/loss)
#   - Exploitation chosen: (1-lambda)*Exploit + V(t)
#   - Exploration unchosen: lambda*Explore + (1-lambda)*phi
#
# Usage: Rscript fit_eef.R [HC|Amph|Hero]
# =============================================================================

# Setup
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

# Load Data
source("utils/load_data.R")
all_data <- load_all_igt_data()
raw_data <- all_data[all_data$study == study_label, ]

# Data Preparation
subIDs <- unique(raw_data$subj)
nsubs <- length(subIDs)
ntrials_max <- 100

ntrials_all <- array(0, c(nsubs))
x_all <- array(0, c(nsubs, ntrials_max))
Gain_all <- array(0, c(nsubs, ntrials_max))
Loss_all <- array(0, c(nsubs, ntrials_max))

cat("Preparing data for", nsubs, "subjects...\n")

for (s in 1:nsubs) {
  subj_df <- raw_data[raw_data$subj == subIDs[s], ]
  ntrials_all[s] <- nrow(subj_df)

  # Choices
  x_sub <- subj_df$choice
  length(x_sub) <- ntrials_max
  x_all[s, ] <- x_sub

  # Separate gain and loss
  gain_sub <- subj_df$gain
  loss_sub <- subj_df$loss
  length(gain_sub) <- ntrials_max
  length(loss_sub) <- ntrials_max

  Gain_all[s, ] <- gain_sub
  Loss_all[s, ] <- loss_sub
}

# Scale outcomes by 100 to prevent numerical issues
cat("Scaling outcomes (divide by 100)...\n")
Gain_all <- Gain_all / 100
Loss_all <- Loss_all / 100

# JAGS Fitting
output_dir <- "outputs/parameter_estimation/eef"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

model_file <- "models/eef.txt"

params <- c(
  "mu_theta", "mu_lambda", "mu_phi", "mu_cons",
  "lambda_theta", "lambda_lambda", "lambda_phi", "lambda_cons",
  "theta", "lam", "phi", "cons"
)

cat("Fitting hierarchical EEF model for group:", target_group, "...\n")

jags_data <- list(
  "x" = x_all,
  "Gain" = Gain_all,
  "Loss" = Loss_all,
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
cat("Fitting complete. Duration:", difftime(end_time, start_time, units = "mins"), "minutes\n")

# Save Results
save_path <- file.path(output_dir, paste0("eef_fit_", target_group, ".rds"))
saveRDS(fit, file = save_path)
cat("Results saved to:", save_path, "\n")

# Basic Diagnostics
print(fit)
