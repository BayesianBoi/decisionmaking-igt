# Quick test script to verify data loading
# Run: Rscript analysis/test_data_loading.R

cat("=== Testing Data Loading ===\n\n")

# Load utilities
source("analysis/utils/load_data.R")
source("analysis/utils/prepare_jags_data.R")

# Test 1: Load individual datasets
cat("Test 1: Loading Ahn 2014 HC data...\n")
ahn_hc <- load_ahn_2014("HC")
cat(sprintf("  Loaded %d rows, %d subjects\n",
            nrow(ahn_hc),
            length(unique(ahn_hc$subj))))

cat("\nTest 2: Loading Fridberg 2010 HC data...\n")
frid_hc <- load_fridberg_2010("HC")
cat(sprintf("  Loaded %d rows, %d subjects\n",
            nrow(frid_hc),
            length(unique(frid_hc$subj))))

cat("\nTest 3: Loading Steingroever 2014 data (95 trials)...\n")
stein_95 <- load_steingroever_2014(95)
cat(sprintf("  Loaded %d rows, %d subjects\n",
            nrow(stein_95),
            length(unique(stein_95$subj))))

# Test 4: Load all data
cat("\nTest 4: Loading all datasets combined...\n")
all_data <- load_all_igt_data()
cat(sprintf("  Total: %d rows, %d unique subjects\n",
            nrow(all_data),
            length(unique(all_data$subj_unique))))

cat("\nStudy breakdown:\n")
study_summary <- table(all_data$study)
print(study_summary)

# Test 5: Validate data
cat("\nTest 5: Validating data...\n")
validate_igt_data(all_data)

# Test 6: Prepare JAGS data
cat("\nTest 6: Preparing JAGS data...\n")
jags_data <- prepare_jags_data(all_data)
cat(sprintf("  N subjects: %d\n", jags_data$N))
cat(sprintf("  Max trials: %d\n", jags_data$T))
cat(sprintf("  Total observations: %d\n", sum(jags_data$Tsubj)))

# Test 7: Check JAGS data integrity
cat("\nTest 7: Checking JAGS data integrity...\n")
check_jags_data(jags_data)

cat("\n=== All Tests Passed! ===\n")
cat("\nYou can now run the main fitting pipeline with:\n")
cat("  Rscript analysis/fit_models.R\n")
