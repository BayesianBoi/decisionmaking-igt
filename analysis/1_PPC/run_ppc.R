# Run Posterior Predictive Checks for EEF and ORL models
# Usage: Rscript analysis/1_analysis/1_PPC/run_ppc.R <model> <group>
# Example: Rscript analysis/1_analysis/1_PPC/run_ppc.R eef HC

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: Rscript run_ppc.R <model> <group>\nModel: eef, orl\nGroup: HC, Amph, Hero")
}

model_name <- args[1]
target_group <- args[2]

cat("=== Running PPC for", model_name, "-", target_group, "===\n\n")

# Load libraries
library(R2jags)

# Source utilities
source("analysis/utils/load_data.R")
source("analysis/utils/ppc.R")

# Load fitted model
fit_path <- file.path(
    "analysis/outputs", model_name,
    paste0(model_name, "_fit_", target_group, ".rds")
)

if (!file.exists(fit_path)) {
    stop("Fit file not found: ", fit_path)
}

cat("Loading fit from:", fit_path, "\n")
fit <- readRDS(fit_path)

# Load and prepare data for this group
all_data <- load_all_igt_data()
group_data <- all_data[all_data$group == target_group, ]

cat("Subjects in group:", length(unique(group_data$subj_unique)), "\n")

# Prepare data in PPC format
subj_list <- unique(group_data$subj_unique)
n_subj <- length(subj_list)
max_trials <- 100

choice_mat <- matrix(NA, nrow = n_subj, ncol = max_trials)
outcome_mat <- matrix(NA, nrow = n_subj, ncol = max_trials)
n_trials_vec <- numeric(n_subj)

for (i in seq_along(subj_list)) {
    subj_id <- subj_list[i]
    subj_data <- group_data[group_data$subj_unique == subj_id, ]
    subj_data <- subj_data[order(subj_data$trial), ]

    n_trials <- nrow(subj_data)
    n_trials_vec[i] <- n_trials
    choice_mat[i, 1:n_trials] <- subj_data$choice
    outcome_mat[i, 1:n_trials] <- subj_data$net / 100 # Scale outcomes
}

observed_data <- list(
    choice = choice_mat,
    outcome = outcome_mat,
    Tsubj = n_trials_vec,
    N = n_subj
)

# Run PPC
cat("\nRunning posterior predictive simulation (100 draws)...\n")
ppc_result <- run_ppc(fit, observed_data, model_name, n_sim = 100)

# Save results
output_dir <- "analysis/outputs/ppc"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_file <- file.path(output_dir, paste0("ppc_", model_name, "_", target_group, ".rds"))
saveRDS(ppc_result, output_file)
cat("\nResults saved to:", output_file, "\n")

# Print summary
cat("\n=== PPC SUMMARY ===\n")
print(ppc_result$summary)
cat("\nDecks within 95% CI:", sum(ppc_result$within_ci), "/ 4\n")
