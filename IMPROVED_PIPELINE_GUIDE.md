# Improved Pipeline Guide

## What's New

The pipeline now has:
1. **Data caching**: Data loaded once and reused
2. **Resume capability**: Already-fitted models are skipped
3. **Detailed logging**: Step-by-step progress tracking
4. **Single model mode**: Fit models one at a time
5. **Immediate saving**: Results saved after each model completes

## Cloud Deployment

### 1. Pull Latest Code

```bash
cd /work/SorNie/decisionmaking-igt/
git pull origin main
```

### 2. Option A: Fit All Models (Recommended)

```bash
# Fit all 3 models (pvl_delta, vse, orl)
nohup Rscript analysis/fit_models.R > fitting.log 2>&1 &
echo $! > pipeline.pid
tail -f fitting.log
```

**Key features:**
- If interrupted, rerun the same command - it will skip completed models
- Data is cached after first load
- Each model saved immediately after completion

### 3. Option B: Fit One Model at a Time

```bash
# Fit PVL-Delta only
nohup Rscript analysis/fit_single_model.R pvl_delta > pvl_delta.log 2>&1 &

# Once complete, fit VSE
nohup Rscript analysis/fit_single_model.R vse > vse.log 2>&1 &

# Finally, fit ORL
nohup Rscript analysis/fit_single_model.R orl > orl.log 2>&1 &
```

**When to use:**
- Debugging a specific model
- Testing before running all models
- If you want to monitor each model separately

## New Log Output Format

You'll see clear progress indicators:

```
=================================================
MODEL: PVL_DELTA
=================================================
Start time: 2026-01-02 14:30:00

Step 1/4: Preparing model-specific data...
  ✓ Data prepared for pvl_delta model

Step 2/4: Initializing parallel execution...
  • Chains: 4
  • Cores: 4
  • Adaptation: 2000 iterations
  • Burn-in: 2000 iterations
  • Sampling: 5000 iterations

Step 3/4: Running MCMC chains in parallel...
  • Sampling started at: 2026-01-02 14:30:15
  • This may take 15-30 minutes depending on model complexity...
  ✓ Sampling complete! Duration: 18.3 minutes

Step 4/4: Combining chains and saving results...
  • Saving results to: analysis/outputs/pvl_delta_fit.rds
  ✓ Saved successfully (145.2 MB)

  Quick convergence check (group-level parameters):
  ✓ Convergence looks good! Max R-hat = 1.023

  Model complete at: 2026-01-02 14:48:18
=================================================
```

## Monitoring Progress

```bash
# Watch live output
tail -f fitting.log

# Check which models are complete
ls -lh analysis/outputs/*_fit.rds

# Check elapsed time
ps -p $(cat pipeline.pid) -o etime=

# Search for errors
grep -i "error\|warning" fitting.log
```

## Resume After Interruption

If the pipeline stops for any reason:

```bash
# Just rerun - it will skip completed models
nohup Rscript analysis/fit_models.R > fitting.log 2>&1 &
```

You'll see:
```
MODEL: PVL_DELTA
⚠ Model already fitted. Loading existing results from:
  analysis/outputs/pvl_delta_fit.rds
To refit, delete this file and re-run the pipeline.
```

## Refit a Specific Model

```bash
# Delete the old fit
rm analysis/outputs/pvl_delta_fit.rds

# Rerun pipeline (will only fit the deleted model)
nohup Rscript analysis/fit_models.R > fitting.log 2>&1 &
```

Or use single model mode:
```bash
rm analysis/outputs/pvl_delta_fit.rds
nohup Rscript analysis/fit_single_model.R pvl_delta > pvl_delta.log 2>&1 &
```

## Expected Timeline (64 cores, parallel mode)

| Model | Expected Duration |
|-------|------------------|
| PVL-Delta | 15-20 min |
| VSE | 25-35 min |
| ORL | 20-25 min |
| **Total** | **~60-80 min** |

## Cached Files

After first run, you'll have:

```
analysis/outputs/
├── cached_data.rds          # Data cache (reused on reruns)
├── pvl_delta_fit.rds        # Model 1 results
├── vse_fit.rds              # Model 2 results
└── orl_fit.rds              # Model 3 results
```

## Troubleshooting

**Problem: "Model already fitted" but you want to refit**

Solution: Delete the specific model file
```bash
rm analysis/outputs/pvl_delta_fit.rds
```

**Problem: Want to start completely fresh**

Solution: Delete all outputs
```bash
rm analysis/outputs/*.rds
```

**Problem: Stuck at "Sampling started" for hours**

Solution:
1. Check if processes are actually running: `top` or `htop`
2. Kill and restart: `kill $(cat pipeline.pid)`
3. Try single model mode to isolate issue

**Problem: Out of memory**

Unlikely with 96GB RAM, but if it happens:
```R
# Edit analysis/fit_models.R or fit_single_model.R
config <- list(
  n_chains = 2,    # Reduce from 4 to 2
  n_iter = 3000,   # Reduce from 5000 to 3000
  ...
)
```

## Next Steps After Completion

```bash
# Compress results for download
cd /work/SorNie/decisionmaking-igt/analysis/outputs
tar -czf model_fits.tar.gz pvl_delta_fit.rds vse_fit.rds orl_fit.rds
ls -lh model_fits.tar.gz

# Download to local machine
# (On your local machine)
scp user@cloud:/work/SorNie/decisionmaking-igt/analysis/outputs/model_fits.tar.gz .

# Extract and generate publication outputs
tar -xzf model_fits.tar.gz
Rscript analysis/generate_paper_outputs.R
```

---

**The improved pipeline is production-ready!** 🚀
