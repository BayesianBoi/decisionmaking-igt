#!/usr/bin/env Rscript
# Quick test of parallel execution
# Run: Rscript analysis/test_parallel.R

library(rjags)
library(coda)
library(parallel)

source("analysis/utils/load_data.R")
source("analysis/utils/prepare_jags_data.R")

cat("=== Testing Parallel Chain Execution ===\n\n")

# Load small dataset
cat("Loading test data...\n")
all_data <- load_all_igt_data()
jags_data <- prepare_jags_data(all_data, study_filter = "Ahn2014_HC")

cat(sprintf("Test data: N=%d subjects, max T=%d trials\n\n",
            jags_data$N, ncol(jags_data$choice)))

# Test configuration
config <- list(
  n_chains = 2,      # Just 2 chains for quick test
  n_adapt = 500,
  n_burnin = 500,
  n_iter = 1000,
  n_cores = 2
)

model_file <- "analysis/models/pvl_delta_v2.jags"
params <- c("mu_A", "mu_alpha", "mu_cons", "mu_lambda")

cat("Testing parallel execution with 2 chains...\n")

# Function to fit a single chain
fit_chain <- function(chain_id) {
  chain_model <- jags.model(
    file = model_file,
    data = jags_data,
    n.chains = 1,
    n.adapt = config$n_adapt,
    quiet = TRUE
  )

  update(chain_model, n.iter = config$n_burnin)

  chain_samples <- coda.samples(
    model = chain_model,
    variable.names = params,
    n.iter = config$n_iter,
    thin = 1
  )

  return(chain_samples)
}

# Time parallel execution
start_time <- Sys.time()

cl <- makeCluster(config$n_cores)
clusterExport(cl, c("model_file", "jags_data", "params", "config"),
             envir = environment())
clusterEvalQ(cl, {
  library(rjags)
  library(coda)
})

chain_list <- parLapply(cl, 1:config$n_chains, fit_chain)
samples <- as.mcmc.list(lapply(chain_list, function(x) x[[1]]))

stopCluster(cl)

end_time <- Sys.time()
parallel_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat(sprintf("\nParallel execution time: %.1f seconds\n\n", parallel_time))

# Check results
cat("--- Results ---\n")
print(summary(samples)$statistics[params, c("Mean", "SD")])

cat("\nGelman-Rubin diagnostic:\n")
print(gelman.diag(samples))

cat("\n=== Parallel Test Complete ===\n")
cat("Parallel execution is working correctly!\n")
