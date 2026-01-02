# Parallel Model Fitting for High-Core Systems

For systems with many cores (32+), you can dramatically speed up fitting by running models in parallel across separate processes.

## Option 1: Run Models in Parallel (Recommended for 64 cores)

Instead of fitting all models sequentially, launch each model in a separate terminal:

```bash
# Terminal 1
Rscript analysis/fit_single_model.R pvl_delta

# Terminal 2
Rscript analysis/fit_single_model.R vse

# Terminal 3
Rscript analysis/fit_single_model.R orl
```

Each model will use ~21 cores (64 cores / 3 models), allowing all three to run simultaneously.

**Expected time with 64 cores**: 30-45 minutes (vs 90-135 minutes sequential)

## Option 2: Use GNU Parallel

```bash
parallel -j 3 Rscript analysis/fit_single_model.R ::: pvl_delta vse orl
```

## Option 3: Submit as Separate Batch Jobs (HPC)

```bash
# Submit three jobs to SLURM/PBS
sbatch --cpus-per-task=21 --wrap="Rscript analysis/fit_single_model.R pvl_delta"
sbatch --cpus-per-task=21 --wrap="Rscript analysis/fit_single_model.R vse"
sbatch --cpus-per-task=21 --wrap="Rscript analysis/fit_single_model.R orl"
```

## Core Allocation

The pipeline automatically uses all available cores efficiently:
- 64 cores, 4 chains → 16 cores per chain
- 32 cores, 4 chains → 8 cores per chain
- 16 cores, 4 chains → 4 cores per chain
- 8 cores, 4 chains → 2 cores per chain

## After Parallel Fitting

Once all models finish, gather results:

```r
source("analysis/utils/diagnostics.R")
run_full_diagnostics("analysis/outputs")
```
