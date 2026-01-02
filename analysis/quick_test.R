#!/usr/bin/env Rscript
# Quick test of model fitting with small dataset
# Run: Rscript analysis/quick_test.R

library(rjags)
library(coda)

# Load utility functions
source("analysis/utils/load_data.R")
source("analysis/utils/prepare_jags_data.R")

set.seed(12345)

cat("=== Quick Model Test ===\n\n")

# Load only one study with few subjects
cat("Loading Ahn2014_HC data...\n")
all_data <- load_all_igt_data()
jags_data <- prepare_jags_data(all_data, study_filter = "Ahn2014_HC")

cat(sprintf("Testing with N=%d subjects, max T=%d trials\n\n",
            jags_data$N, jags_data$T))

# Test PVL-Delta model
cat("--- Testing PVL-Delta Model ---\n")
model_file <- "analysis/models/pvl_delta_v2.jags"

cat("Initializing JAGS model...\n")
jags_model <- jags.model(
  file = model_file,
  data = jags_data,
  n.chains = 2,
  n.adapt = 500,
  quiet = FALSE
)

cat("Running burn-in (500 iterations)...\n")
update(jags_model, n.iter = 500)

cat("Sampling (1000 iterations)...\n")
samples <- coda.samples(
  model = jags_model,
  variable.names = c("mu_A", "mu_alpha", "mu_cons", "mu_lambda"),
  n.iter = 1000
)

cat("\n--- Results ---\n")
print(summary(samples))

cat("\nGelman-Rubin diagnostic:\n")
print(gelman.diag(samples))

cat("\n=== Test Complete! ===\n")
cat("If you see reasonable parameter estimates above, the model is working.\n")
cat("Now you can run the full pipeline with: Rscript analysis/fit_models.R\n")
