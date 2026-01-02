# Parallel Execution Guide

## What Changed

The pipeline now runs MCMC chains in **parallel** for ~4x faster execution.

**Old approach (sequential):**
- 4 chains run one after another
- Uses 1 core at a time
- Full pipeline: ~2-3 hours on 16-core machine

**New approach (parallel):**
- 4 chains run simultaneously
- Uses 4 cores concurrently
- Full pipeline: **~45-60 minutes** on 16-core machine

## Cloud Deployment

### 1. Pull Latest Code

```bash
cd /work/SorNie/decisionmaking-igt/
git pull origin main
```

### 2. Verify Parallel Package

The `parallel` package is included in base R, but let's verify:

```bash
Rscript -e "if (require('parallel', quietly=TRUE)) cat('✓ parallel available\n') else cat('✗ parallel NOT available\n')"
```

### 3. Run Pipeline (Parallel Mode)

```bash
# Full pipeline with parallel chains
nohup Rscript analysis/fit_models.R > fitting.log 2>&1 &

# Save process ID
echo $! > pipeline.pid

# Monitor progress
tail -f fitting.log
```

You should see:
```
Parallel mode enabled: using 4 cores (out of 16 available)
=== IGT Model Fitting Pipeline ===
...
Running 4 chains in parallel on 4 cores...
```

### 4. Optional: Test Parallel Mode First

```bash
# Quick 2-minute test
Rscript analysis/test_parallel.R
```

Expected output: ~2 minutes, parameters around:
- mu_A: 0.04-0.05
- mu_alpha: 0.25-0.30
- mu_cons: 2.3-2.5
- mu_lambda: 0.30-0.35

## Configuration Options

Edit `analysis/fit_models.R` if needed:

```R
config <- list(
  n_chains = 4,        # Number of MCMC chains
  n_adapt = 2000,      # Adaptation iterations
  n_burnin = 2000,     # Burn-in iterations
  n_iter = 5000,       # Sampling iterations
  thin = 1,
  models = c("pvl_delta", "vse", "orl"),
  fit_all_studies = TRUE,
  parallel = TRUE,     # ← Set to FALSE to disable parallelization
  n_cores = 4          # ← Cores to use (max: n_chains)
)
```

## Expected Timeline (16-core machine)

**With parallel = TRUE:**
- Data loading: ~2 min
- PVL-Delta: ~10-15 min
- VSE: ~20-30 min
- ORL: ~15-20 min
- **Total: ~45-60 minutes**

**With parallel = FALSE (old):**
- Total: ~2-3 hours

## Monitoring

```bash
# Check if using multiple cores
top  # Look for multiple Rscript processes

# Check progress
tail -50 fitting.log

# Check time elapsed
ps -p $(cat pipeline.pid) -o etime=
```

## Troubleshooting

**Issue: "parallel package not available"**

The pipeline will automatically fall back to sequential mode. This is rare since `parallel` is in base R.

**Issue: Slower than expected**

Check system load:
```bash
top
htop  # If available
```

If other processes are using CPU, consider:
- Reducing `n_cores` to 2
- Setting `parallel = FALSE` temporarily

**Issue: Out of memory**

Unlikely with 96GB RAM, but if it happens:
- Set `fit_all_studies = FALSE` to test on smaller subset
- Reduce `n_chains` to 2

## Verification

After pipeline completes, check outputs:

```bash
ls -lh analysis/outputs/*.rds

# Should see:
# pvl_delta_fit.rds (~100-200 MB)
# vse_fit.rds       (~200-300 MB)
# orl_fit.rds       (~150-250 MB)
```

## Performance Comparison

| Setup | Cores Used | Total Time | Speedup |
|-------|------------|------------|---------|
| Sequential | 1 | ~2-3 hours | 1x |
| Parallel (2 chains) | 2 | ~1-1.5 hours | ~2x |
| Parallel (4 chains) | 4 | ~45-60 min | ~3-4x |

---

**Ready to deploy with parallel execution!** 🚀
