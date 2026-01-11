# Format Raw IGT Data
# This script loads and formats raw IGT data from all studies.
# It outputs a cleaned dataset ready for model fitting.
#
# Usage: Rscript analysis/0_preprocess/0_format_raw_data.R
#

source("analysis/utils/load_data.R")

message("Loading and formatting IGT data...\n")

# Load all datasets (Ahn 2014)
dat <- load_all_igt_data()

# Validate data
validate_igt_data(dat)

# Add contiguous subject index (required for JAGS)
dat <- add_subject_index(dat)

# Summary statistics
message("\n=== Data Summary ===")
message(sprintf("Total subjects: %d", length(unique(dat$subj_unique))))
message(sprintf("Total trials: %d", nrow(dat)))
message(sprintf("Studies: %s", paste(unique(dat$study), collapse = ", ")))

# Save formatted data
output_dir <- "data/processed"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(dat, file.path(output_dir, "igt_data_formatted.rds"))

message(sprintf("\nFormatted data saved to: %s/igt_data_formatted.rds", output_dir))
