# Cloud Deployment Checklist

Use this checklist when deploying to your cloud machine.

## Pre-Deployment (Local)

- [x] ✅ All models fixed (pvl_delta_v2, vse_v2, orl_v2)
- [x] ✅ Quick test passes locally
- [x] ✅ fit_models.R updated to use v2 models
- [ ] Package code for transfer
  ```bash
  cd /Users/nielsvaerbak/Desktop/decision_making
  tar -czf decision_making.tar.gz analysis/ data/ R/
  ```

## Cloud Setup

- [ ] Provision cloud instance (8+ cores, 16GB RAM)
- [ ] Transfer code
  ```bash
  scp decision_making.tar.gz user@your-cloud-ip:~/
  ```
- [ ] SSH into instance
  ```bash
  ssh user@your-cloud-ip
  ```
- [ ] Extract code
  ```bash
  tar -xzf decision_making.tar.gz
  cd decision_making
  ```
- [ ] Install R
  ```bash
  sudo apt-get update
  sudo apt-get install -y r-base r-base-dev
  ```
- [ ] Install JAGS
  ```bash
  sudo apt-get install -y jags
  ```
- [ ] Install R packages
  ```bash
  sudo R -e "install.packages(c('rjags', 'coda', 'dplyr'), repos='https://cloud.r-project.org/')"
  ```
- [ ] Verify installation
  ```bash
  Rscript analysis/quick_test.R
  ```
  **Expected:** Completes in ~2 minutes with R-hat values around 1.0

## Running Pipeline

- [ ] Review configuration (optional)
  ```bash
  nano analysis/fit_models.R
  # Check config section - defaults are good for full run
  ```
- [ ] Start pipeline in background
  ```bash
  nohup Rscript analysis/fit_models.R > fitting_output.log 2>&1 &
  ```
- [ ] Note the process ID
  ```bash
  echo $! > pipeline.pid
  ```
- [ ] Monitor progress
  ```bash
  tail -f fitting_output.log
  ```

## Monitoring

Expected progress messages:

```
=== IGT Model Fitting Pipeline ===

Step 1: Loading and validating data...
✓ Data validation passed: 173 subjects, 81831 total trials

=== Step 2: Fitting Models ===

--- Fitting PVL_DELTA model ---
✓ Initializing JAGS model...
✓ Running burn-in (2000 iterations)...
✓ Sampling from posterior (5000 iterations)...
✓ Saved results to: analysis/outputs/pvl_delta_fit.rds

--- Fitting VSE model ---
...
```

- [ ] Check for errors
  ```bash
  grep -i "error\|failed" fitting_output.log
  ```

- [ ] Verify process still running
  ```bash
  ps aux | grep fit_models.R
  ```

## Expected Timeline

| Time    | Event                                |
|---------|--------------------------------------|
| 0:00    | Start pipeline                       |
| 0:02    | Data loaded and validated            |
| 0:05    | PVL-Delta model compilation complete |
| 0:30    | PVL-Delta fitting complete           |
| 0:35    | VSE model compilation complete       |
| 1:45    | VSE fitting complete                 |
| 1:50    | ORL model compilation complete       |
| 2:45    | ORL fitting complete                 |
| 2:45    | Pipeline summary generated           |
| **~3-5 hours total**                       |

## Completion

- [ ] Verify all models completed
  ```bash
  ls -lh analysis/outputs/*.rds
  ```
  **Expected:** 3 files (pvl_delta_fit.rds, vse_fit.rds, orl_fit.rds)

- [ ] Check final summary in log
  ```bash
  tail -50 fitting_output.log
  ```

- [ ] Compress results
  ```bash
  cd analysis/outputs
  tar -czf model_fits.tar.gz *.rds
  ls -lh model_fits.tar.gz
  ```

## Download Results

- [ ] From local machine
  ```bash
  scp user@your-cloud-ip:~/decision_making/analysis/outputs/model_fits.tar.gz .
  ```

- [ ] Extract locally
  ```bash
  tar -xzf model_fits.tar.gz -C /Users/nielsvaerbak/Desktop/decision_making/analysis/outputs/
  ```

- [ ] Verify files
  ```bash
  ls -lh /Users/nielsvaerbak/Desktop/decision_making/analysis/outputs/*.rds
  ```

## Post-Processing (Local)

- [ ] Run diagnostics
  ```bash
  cd /Users/nielsvaerbak/Desktop/decision_making
  Rscript -e "source('analysis/utils/diagnostics.R'); run_full_diagnostics()"
  ```

- [ ] Check convergence
  ```bash
  # Look for R-hat < 1.1 for all parameters
  ```

- [ ] Run posterior predictive checks
  ```bash
  Rscript -e "source('analysis/utils/ppc.R'); run_all_ppc()"
  ```

- [ ] Compare models
  ```bash
  Rscript -e "source('analysis/utils/model_comparison.R'); compare_models()"
  ```

- [ ] Review results
  ```bash
  cat analysis/outputs/model_comparison.md
  ```

## Cleanup

- [ ] Shutdown cloud instance (if using pay-per-use)
  ```bash
  sudo shutdown -h now
  ```

- [ ] Or stop instance from cloud provider console

- [ ] Archive local results
  ```bash
  cd /Users/nielsvaerbak/Desktop/decision_making/analysis/outputs
  tar -czf all_results_$(date +%Y%m%d).tar.gz *.rds *.pdf *.md
  ```

## Troubleshooting

### Pipeline stuck/no progress for >30 minutes

```bash
# Check if process is still running
ps aux | grep fit_models.R

# Check system resources
top
free -h

# Check last log entries
tail -20 fitting_output.log
```

### Out of memory errors

```bash
# Edit config to use smaller dataset
nano analysis/fit_models.R
# Set: fit_all_studies = FALSE

# Restart pipeline
pkill -f fit_models.R
nohup Rscript analysis/fit_models.R > fitting_output.log 2>&1 &
```

### Model won't initialize

```bash
# Verify using v2 models
grep "v2.jags" analysis/fit_models.R

# Test specific model
Rscript analysis/quick_test.R
```

## Success Criteria

✅ Pipeline completes without errors
✅ 3 .rds files created (pvl_delta_fit, vse_fit, orl_fit)
✅ R-hat values < 1.1 for all parameters
✅ Effective sample sizes > 400 for key parameters
✅ Posterior predictive checks show reasonable fit
✅ Model comparison report generated

## Notes

- Expected cost: $0.50-$1.70 (depending on spot vs on-demand)
- Expected time: 3-5 hours
- File sizes: ~100-500 MB per model
- Keep instance type at 8 cores for optimal MCMC chain parallelization

## Quick Reference

**Check status:**
```bash
tail -f fitting_output.log
```

**Kill pipeline:**
```bash
pkill -f fit_models.R
```

**Restart pipeline:**
```bash
nohup Rscript analysis/fit_models.R > fitting_output.log 2>&1 &
```

---

Last updated: 2026-01-02
