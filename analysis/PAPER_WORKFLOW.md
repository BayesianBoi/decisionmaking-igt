# Paper Workflow Guide

## Overview

This guide walks through the complete workflow from model fitting to manuscript submission.

## Step-by-Step Workflow

### 1. Parameter Recovery (Validation)

Before fitting real data, verify the model can recover known parameters:

```bash
Rscript analysis/scripts/parameter_recovery_eef.R
```

**Output:** `results/parameter_recovery/`
- `parameter_recovery.pdf` - Scatter plots of true vs recovered values
- Correlation > 0.8 indicates good recovery

### 2. Fit EEF Model

```bash
Rscript analysis/scripts/fit_eef_clinical.R
```

**Runtime:** 4-8 hours on M1 Mac

**Output:** `results/eef_clinical/`
- `mcmc_samples.rds` - Posterior samples
- `parameter_summary.rds` - Mean, SD, quantiles
- `diagnostics.rds` - R-hat, ESS
- `trace_plots.pdf` - Convergence diagnostics
- `density_plots.pdf` - Posterior distributions

### 3. Fit Comparison Models

```bash
Rscript analysis/scripts/fit_pvl_delta.R
Rscript analysis/scripts/fit_vse.R
```

### 4. Check Convergence

Verify in diagnostic output:
- R-hat < 1.1 for all parameters
- Effective sample size > 1000
- Trace plots show good mixing

### 5. Group Comparisons

```bash
Rscript analysis/scripts/compare_groups.R
```

**Key outputs:**
- Posterior probability that λ_substance > λ_HC
- 95% credible intervals for group differences
- Violin plots by group

### 6. Model Comparison

```bash
Rscript analysis/scripts/compare_models.R
```

Compares convergence and parameter estimates across EEF, PVL-Delta, and VSE.

### 7. Generate Figures

```bash
Rscript analysis/scripts/create_figures.R
```

**Output:** `results/figures/`
- `figure1_model_comparison.pdf`
- `figure2_group_differences.pdf`
- `figure3_parameter_recovery.pdf`
- `figure4_posterior_predictive.pdf`

## Manuscript Integration

### Methods Section

```
We employed the Exploitation-Exploration with Forgetting (EEF) model
(Yang et al., 2025) to analyze Iowa Gambling Task performance. The model
was fit using hierarchical Bayesian estimation in JAGS with 4 chains,
5000 adaptation iterations, 10000 burn-in, and 20000 sampling iterations.

Prior specifications followed Yang et al. (2025) with modifications for
clinical populations (see Supplementary Materials for full prior justification).
```

### Results Section

1. **Convergence:** Report from diagnostics.rds
2. **Group Comparisons:** Report P(λ_substance > λ_HC) and 95% CIs
3. **Model Comparison:** Report which model best fits the data

Example:
```
Substance users exhibited higher forgetting rates than healthy controls
(λ_substance = 0.52, 95% CI [0.45, 0.59] vs λ_HC = 0.38, 95% CI [0.31, 0.45];
posterior probability P(λ_substance > λ_HC) = 0.94).
```

### Figures

- **Figure 1:** Group differences in forgetting rate
- **Figure 2:** Parameter recovery validation
- **Figure 3:** Posterior predictive checks

### Supplementary Materials

Include:
- `docs/PRIOR_JUSTIFICATIONS.md` - Full prior specifications
- Trace plots for all parameters
- Complete parameter estimates table

## Publication Checklist

Before submission:

- [ ] All models converged (R-hat < 1.1)
- [ ] Parameter recovery validated (r > 0.8)
- [ ] Group comparisons computed
- [ ] Figures generated at publication quality
- [ ] Prior justifications documented
- [ ] MCMC settings reported in methods
- [ ] Data and code availability statement
