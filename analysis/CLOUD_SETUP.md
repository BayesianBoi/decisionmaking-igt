# Cloud Setup Guide

## Running the IGT Model Fitting Pipeline on Cloud Compute

This guide explains how to run the computationally intensive model fitting on a cloud machine.

## Prerequisites

### Cloud Instance Recommendations

**Minimum Specifications:**
- 8 CPU cores (for parallel chains)
- 16 GB RAM
- 20 GB storage
- Ubuntu 20.04/22.04 or similar Linux distribution

**Recommended Providers:**
- AWS EC2: `c6i.2xlarge` or similar compute-optimized instance
- Google Cloud: `c2-standard-8` or similar
- Digital Ocean: CPU-Optimized Droplet with 8 vCPUs

**Expected Runtime:**
- With 173 subjects, 4 chains, 2000 burn-in, 5000 samples:
  - PVL-Delta: ~30-60 minutes
  - VSE: ~60-120 minutes
  - ORL: ~45-90 minutes
  - **Total: ~3-5 hours**

## Setup Instructions

### 1. Transfer Code to Cloud

```bash
# From your local machine
cd /Users/nielsvaerbak/Desktop/decision_making
tar -czf decision_making.tar.gz analysis/ data/ R/

# Upload to cloud (replace with your instance details)
scp decision_making.tar.gz user@your-cloud-ip:~/

# On cloud machine
ssh user@your-cloud-ip
tar -xzf decision_making.tar.gz
cd decision_making
```

### 2. Install Dependencies

```bash
# Update package list
sudo apt-get update

# Install R
sudo apt-get install -y r-base r-base-dev

# Install JAGS
sudo apt-get install -y jags

# Install required R packages
sudo R -e "install.packages(c('rjags', 'coda', 'dplyr'), repos='https://cloud.r-project.org/')"
```

### 3. Verify Installation

```bash
# Quick test (should complete in ~2 minutes)
Rscript analysis/quick_test.R
```

If this completes successfully with R-hat values around 1.0, you're ready to run the full pipeline.

### 4. Run Full Pipeline

```bash
# Run in background with output logging
nohup Rscript analysis/fit_models.R > fitting_output.log 2>&1 &

# Monitor progress
tail -f fitting_output.log

# Check if still running
ps aux | grep fit_models.R
```

## Configuration Options

Before running, you can edit `analysis/fit_models.R` to adjust:

```R
config <- list(
  n_chains = 4,           # Number of MCMC chains (use number of cores)
  n_adapt = 2000,         # Adaptation iterations
  n_burnin = 2000,        # Burn-in iterations
  n_iter = 5000,          # Sampling iterations
  thin = 1,               # Thinning (1 = no thinning)
  models = c("pvl_delta", "vse", "orl"),  # Models to fit
  fit_all_studies = TRUE  # TRUE = all 173 subjects, FALSE = 48 subjects
)
```

**For faster testing:**
```R
config <- list(
  n_chains = 2,
  n_adapt = 500,
  n_burnin = 500,
  n_iter = 1000,
  thin = 1,
  models = c("pvl_delta"),  # Fit only one model first
  fit_all_studies = FALSE    # Use smaller dataset
)
```

## Monitoring Progress

The pipeline outputs progress messages:

```
=== IGT Model Fitting Pipeline ===

Step 1: Loading and validating data...
Data validation passed: 173 subjects, 81831 total trials

=== Step 2: Fitting Models ===

--- Fitting PVL_DELTA model ---
Initializing JAGS model...
Compiling model graph...
Running burn-in (2000 iterations)...
Sampling from posterior (5000 iterations)...
Saved results to: analysis/outputs/pvl_delta_fit.rds

--- Fitting VSE model ---
...
```

## Retrieving Results

After completion, download the results:

```bash
# On cloud machine, compress outputs
cd ~/decision_making/analysis/outputs
tar -czf model_fits.tar.gz *.rds

# From local machine
scp user@your-cloud-ip:~/decision_making/analysis/outputs/model_fits.tar.gz .

# Extract locally
tar -xzf model_fits.tar.gz -C /Users/nielsvaerbak/Desktop/decision_making/analysis/outputs/
```

## Post-Processing Locally

Once you have the `.rds` files locally, run diagnostics and analysis:

```bash
cd /Users/nielsvaerbak/Desktop/decision_making

# Run diagnostics
Rscript -e "source('analysis/utils/diagnostics.R'); run_full_diagnostics()"

# Run posterior predictive checks
Rscript -e "source('analysis/utils/ppc.R'); run_all_ppc()"

# Compare models
Rscript -e "source('analysis/utils/model_comparison.R'); compare_models()"
```

## Troubleshooting

### Issue: Out of Memory

**Solution 1:** Reduce subjects
```R
jags_data <- prepare_jags_data(all_data, study_filter = "Ahn2014_HC")
```

**Solution 2:** Reduce iterations
```R
config$n_iter <- 2000  # instead of 5000
```

**Solution 3:** Use larger instance with more RAM

### Issue: Takes Too Long

**Solution 1:** Fit models sequentially instead of all at once
```R
config$models <- c("pvl_delta")  # Run once for each model
```

**Solution 2:** Use multiple instances (one per model)

**Solution 3:** Reduce thinning interval (saves fewer samples)
```R
config$thin <- 2  # Save every 2nd sample
```

### Issue: JAGS Not Found

```bash
# Check JAGS installation
which jags
jags -h

# If not found, reinstall
sudo apt-get install --reinstall jags
```

### Issue: Model Won't Initialize

Check the log for specific errors. Common fixes:
- Ensure using `*_v2.jags` models (not the original versions)
- Verify data loaded correctly with `quick_test.R` first

## Cost Optimization

### AWS EC2 Example

Using `c6i.2xlarge` (8 vCPU, 16GB RAM):
- On-demand: ~$0.34/hour × 5 hours = ~$1.70
- Spot instance: ~$0.10/hour × 5 hours = ~$0.50

**Recommendation:** Use spot instances for cost savings (but be prepared for interruption)

### Automated Shutdown

Add to the end of your script to auto-shutdown after completion:

```bash
# After running Rscript
sudo shutdown -h now
```

Or use cloud provider's auto-stop features.

## Alternative: Use Screen/Tmux

Instead of `nohup`, use screen for better session management:

```bash
# Start screen session
screen -S model_fitting

# Run pipeline
Rscript analysis/fit_models.R

# Detach: Ctrl+A, then D
# Reattach later
screen -r model_fitting
```

## Parallel Processing (Advanced)

To fit models in parallel, modify the pipeline to use multiple R sessions:

```bash
# Terminal 1
Rscript -e "source('analysis/fit_models.R'); config\$models <- c('pvl_delta'); ..."

# Terminal 2
Rscript -e "source('analysis/fit_models.R'); config\$models <- c('vse'); ..."

# Terminal 3
Rscript -e "source('analysis/fit_models.R'); config\$models <- c('orl'); ..."
```

Or use GNU parallel (not demonstrated here).

## Summary Checklist

- [ ] Cloud instance with 8+ cores, 16GB RAM
- [ ] R, JAGS, and packages installed
- [ ] Code and data transferred
- [ ] `quick_test.R` runs successfully
- [ ] Configuration adjusted if needed
- [ ] Pipeline running in background
- [ ] Monitor progress with `tail -f`
- [ ] Download results when complete
- [ ] Run diagnostics locally

## Support

If you encounter issues:
1. Check `fitting_output.log` for error messages
2. Verify with `quick_test.R` on smaller data
3. Review `TROUBLESHOOTING.md` for common fixes
4. Check JAGS and R package versions match local environment
