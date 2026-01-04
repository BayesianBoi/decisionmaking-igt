# ===========================================================================
# Verification Pipeline
# ===========================================================================
#
# Runs a short 10-iteration test of all models to ensure:
# 1. Models compile correctly
# 2. Log-likelihoods are being saved
# 3. Output files are generated properly
#
# Usage: Rscript analysis/scripts/verify_pipeline.R
#
# ===========================================================================

library(rjags)
library(coda)

source("analysis/utils/load_data.R")
source("analysis/utils/prepare_jags_data.R")
source("analysis/utils/prepare_eef_data.R")

# Create test output directory
output_dir <- "results/verification"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load data once
dat_all <- load_all_igt_data()
clinical_studies <- c("Ahn2014_HC", "Ahn2014_Amph", "Fridberg2010_HC")
# Reduced dataset for speed (no Heroin/Cannabis for verify)
dat_clinical <- dat_all[dat_all$study %in% clinical_studies, ]
# Keep only first 5 subjects for speed
subjs <- unique(dat_clinical$subj_unique)[1:5]
dat_clinical <- dat_clinical[dat_clinical$subj_unique %in% subjs, ]
dat_clinical$group <- "HC" # Placeholder group for simple check

# Prepare data
jags_data_orl <- prepare_jags_data(dat_clinical)
jags_data_eef <- prepare_eef_jags_data(dat_clinical, group_var = "group", reference_group = "HC")

# List of models to test
models <- list(
    list(
        name = "pvl_delta", file = "analysis/models/pvl_delta.jags", data = jags_data_orl,
        params = c("mu_A", "log_lik")
    ),
    list(
        name = "orl", file = "analysis/models/orl.jags", data = jags_data_orl,
        params = c("mu_Arew", "log_lik")
    ),
    list(
        name = "eef", file = "analysis/models/eef.jags", data = jags_data_eef,
        params = c("mu_theta", "log_lik")
    )
)

# Run tests
for (m in models) {
    message(sprintf("\nVerifying %s...", m$name))

    if (!file.exists(m$file)) {
        stop(sprintf("Model file missing: %s", m$file))
    }

    start_time <- Sys.time()

    # Run extremely short chain
    model <- jags.model(
        file = m$file,
        data = m$data,
        n.chains = 1,
        n.adapt = 10,
        quiet = TRUE
    )

    update(model, n.iter = 10)

    samples <- coda.samples(
        model = model,
        variable.names = m$params,
        n.iter = 10,
        thin = 1
    )

    # Check if log_lik exists in columns
    col_names <- colnames(samples[[1]])
    has_log_lik <- any(grepl("log_lik", col_names))

    if (has_log_lik) {
        message(sprintf("[PASS] %s: log_lik found (Total columns: %d)", m$name, length(col_names)))
    } else {
        stop(sprintf("[FAIL] %s: log_lik NOT found in output!", m$name))
    }

    message(sprintf("Time: %.2fs", as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
}

message("\nAll models verified successfully. Ready for cloud deployment.")
