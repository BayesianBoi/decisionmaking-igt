# Analysis Pipeline Documentation

This directory contains the complete fitting and validation pipeline for IGT models.

## Core Scripts

### Fitting Scripts (in `scripts/`)

**fit_eef_clinical.R**
Fits the EEF (Exploitation-Exploration with Forgetting) model to clinical populations. This is the primary model for testing memory decay in substance users.

**fit_pvl_delta.R**
Fits the PVL-Delta model as a baseline for model comparison.

**fit_vse.R**
Fits the VSE (Value + Sequential Exploration) model for model comparison.

**parameter_recovery_eef.R**
Tests parameter identifiability by simulating data with known parameters and fitting the model.

**compare_groups.R**
Compares forgetting rates (lambda) between clinical groups (HC vs substance users).

**compare_models.R**
Compares convergence diagnostics and parameter estimates across models.

**create_figures.R**
Generates publication-ready figures for the paper.

## JAGS Models

Located in `models/` directory.

**eef_clinical.jags**
- 4 parameters: theta, lambda_forget, phi, cons
- Key innovation: Memory decay applied to both exploitation and exploration
- Uses first-choice priors (w_ini) from data

**pvl_delta.jags**
- 4 parameters: A (learning rate), alpha (outcome sensitivity), cons (consistency), lambda (loss aversion)
- Uses delta-rule updating
- Baseline model without forgetting mechanism

**vse.jags**
- 8 parameters: PVL-Delta parameters plus epP, epN (perseverance boosts), K (decay), w (weight)
- Combines value learning with sequential exploration tendencies

**orl.jags**
- 5 parameters: Arew, Apun (learning rates), K (decay), betaF, betaP (weights)
- Tracks expected value and expected frequency separately

## Utility Functions

**load_data.R**
Loads and harmonizes IGT data from Ahn 2014, Fridberg 2010, and Steingroever 2014 studies.

**prepare_jags_data.R**
Converts harmonized data into JAGS-compatible format for PVL-Delta, VSE, and ORL models.

**prepare_eef_data.R**
Prepares data specifically for EEF model with first-choice priors calculation and group structure.

**diagnostics.R**
Computes convergence diagnostics (R-hat, ESS), generates trace plots.

**ppc.R**
Posterior predictive checks.

**parameter_recovery.R**
Parameter recovery utility functions.

**model_comparison.R**
Model comparison metrics (DIC, WAIC).

**visualization.R**
Creates publication-quality figures.

**reporting.R**
Generates APA-formatted tables.

## Workflow

1. **Parameter Recovery** (optional but recommended)
   ```bash
   Rscript analysis/scripts/parameter_recovery_eef.R
   ```

2. **Fit Models**
   ```bash
   Rscript analysis/scripts/fit_eef_clinical.R
   Rscript analysis/scripts/fit_pvl_delta.R
   Rscript analysis/scripts/fit_vse.R
   ```

3. **Group Comparisons**
   ```bash
   Rscript analysis/scripts/compare_groups.R
   ```

4. **Model Comparison**
   ```bash
   Rscript analysis/scripts/compare_models.R
   ```

5. **Generate Figures**
   ```bash
   Rscript analysis/scripts/create_figures.R
   ```

## Outputs

Results are saved in `results/` directory:
- `results/eef_clinical/` - EEF model outputs
- `results/pvl_delta/` - PVL-Delta model outputs
- `results/vse/` - VSE model outputs
- `results/group_comparison/` - Group comparison results
- `results/model_comparison/` - Model comparison results
- `results/parameter_recovery/` - Parameter recovery results
- `results/figures/` - Publication figures

## MCMC Configuration

Default MCMC settings:
- 4 chains
- 5000 adaptation iterations
- 10000 burn-in iterations
- 20000 sampling iterations
- Thinning = 2

Estimated runtime: 4-8 hours per model.

## Requirements

- R >= 4.0
- JAGS >= 4.3
- R packages: rjags, coda, ggplot2, loo (optional)

## Troubleshooting

See `TROUBLESHOOTING.md` for common issues and solutions.
