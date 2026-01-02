# Analysis Pipeline Documentation

This directory contains the complete fitting and validation pipeline for IGT models.

## Core Scripts

**fit_models.R**
Main pipeline that fits all three models sequentially. Includes data caching and automatic resume capability. Models are fit with 4 chains, 2000 adaptation iterations, 2000 burn-in, and 5000 sampling iterations.

**fit_single_model.R**
Fit individual models. Useful for parallel execution across separate processes or debugging specific models.

Usage:
```bash
Rscript analysis/fit_single_model.R pvl_delta
Rscript analysis/fit_single_model.R vse
Rscript analysis/fit_single_model.R orl
```

**quick_test.R**
Runs a quick validation using Ahn2014_HC data (48 subjects). Verifies that models compile and basic convergence is achieved. Takes about 2 minutes.

## JAGS Models

Located in `models/` directory. Each model has theoretical documentation in its header.

**pvl_delta_v2.jags**
- 4 parameters: A (learning rate), alpha (outcome sensitivity), cons (consistency), lambda (loss aversion)
- Uses delta-rule updating: EV[t+1] = EV[t] + A * (u[t] - EV[t])
- Utility function: u(x) = x^alpha for gains, -lambda * |x|^alpha for losses

**vse_v2.jags**
- 8 parameters: PVL-Delta parameters plus epP, epN (perseverance boosts), K (decay), w (weight)
- Combines expected value learning with sequential exploration tendencies
- Choice based on weighted average: w * EV + (1-w) * perseverance

**orl_v2.jags**
- 5 parameters: Arew, Apun (learning rates), K (decay), betaF, betaP (weights)
- Tracks expected value (EV) and expected frequency (EF) separately
- Includes fictive updating for unchosen options

## Utility Functions

**load_data.R**
Harmonizes data from three sources into consistent format. Creates unique subject IDs across studies and validates choice ranges.

**prepare_jags_data.R**
Converts harmonized data into JAGS-compatible format with proper matrix dimensions and variable-length trial handling.

**diagnostics.R**
Computes convergence diagnostics (R-hat, ESS), generates trace plots, and creates diagnostic reports for all fitted models.

**ppc.R**
Posterior predictive checks. Simulates data from fitted models and compares to observed choice proportions.

**parameter_recovery.R**
Tests parameter identifiability by fitting models to simulated data with known parameters.

**model_comparison.R**
Computes model comparison metrics (DIC, WAIC) for selecting best-fitting model.

**visualization.R**
Creates publication-quality figures including trace plots and posterior densities.

**reporting.R**
Generates APA-formatted tables and result summaries for manuscripts.

## Workflow

1. Data loading and validation
2. Model fitting (sequential or parallel)
3. Convergence diagnostics
4. Posterior predictive checks
5. Parameter recovery (optional)
6. Model comparison
7. Publication outputs

## Outputs

Results are saved in `outputs/` directory:
- `cached_data.rds` - Preprocessed data
- `pvl_delta_fit.rds`, `vse_fit.rds`, `orl_fit.rds` - Fitted models
- `*_trace_plots.pdf` - MCMC diagnostics
- `all_diagnostics.rds` - Convergence metrics
- `ppc_results.rds` - Posterior predictive checks

## Configuration

Default MCMC settings in `fit_models.R`:
- 4 chains
- 2000 adaptation iterations
- 2000 burn-in iterations
- 5000 sampling iterations
- Parallel execution enabled

Adjust these based on your computational resources and convergence requirements.

## Troubleshooting

See `TROUBLESHOOTING.md` for common issues and solutions.

For parallel fitting on high-core systems, see `PARALLEL_FITTING.md`.

For preparing publication outputs, see `PAPER_WORKFLOW.md`.
