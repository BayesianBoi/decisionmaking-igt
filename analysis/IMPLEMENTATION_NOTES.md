# Implementation Notes

## Summary

This directory contains a complete R + JAGS analysis pipeline for fitting decision-making models to Iowa Gambling Task (IGT) data. The implementation follows the detailed plan specified in `coding_agent_prompt.md`.

## Implementation Status

All planned deliverables have been completed:

### 1. Data Harmonization ✓
- **load_data.R**: Functions to load and standardize data from 3 studies (8 datasets total)
  - Ahn et al. (2014): 3 groups
  - Fridberg et al. (2010): 2 groups
  - Steingroever et al. (2014): 3 trial lengths
- **prepare_jags_data.R**: Converts harmonized data to JAGS format with validation

### 2. JAGS Models ✓
Three decision-making models implemented:
- **pvl_delta.jags**: PVL-Delta (4 parameters: A, alpha, cons, lambda)
- **vse.jags**: VSE/VPP (8 parameters: adds epP, epN, K, w to PVL-Delta)
- **orl.jags**: ORL (5 parameters: Arew, Apun, K, betaF, betaP)

All models use hierarchical Bayesian structure with group-level and subject-level parameters.

### 3. Fitting Pipeline ✓
- **fit_models.R**: Main script to fit all models
  - Configurable MCMC settings (chains, burn-in, iterations)
  - Can fit to combined dataset or individual studies
  - Saves results as .rds files

### 4. Diagnostics ✓
- **diagnostics.R**: MCMC convergence checking
  - R-hat (Gelman-Rubin diagnostic)
  - Effective sample size (ESS)
  - Trace plots for visual inspection
  - Automated checks for all fitted models

### 5. Posterior Predictive Checks ✓
- **ppc.R**: Model validation
  - Simulates choice data from posterior
  - Compares observed vs predicted choice proportions
  - Credible interval checks
  - Implemented for PVL-Delta and ORL

### 6. Model Comparison ✓
- **model_comparison.R**: Compare model fit
  - Parameter counts
  - Framework for DIC/WAIC/LOO (requires additional computation)
  - Generates markdown report

### 7. Parameter Recovery ✓
- **parameter_recovery.R**: Validate model identifiability
  - Sample parameters from priors
  - Simulate data
  - Re-fit and check recovery correlations
  - Implemented for PVL-Delta (extensible to other models)

## Quick Start

1. **Test data loading**:
   ```bash
   Rscript analysis/test_data_loading.R
   ```

2. **Fit models** (warning: computationally intensive):
   ```bash
   Rscript analysis/fit_models.R
   ```

3. **Run diagnostics**:
   ```R
   source("analysis/utils/diagnostics.R")
   diagnostics <- run_full_diagnostics()
   ```

4. **Run posterior predictive checks**:
   ```R
   source("analysis/utils/ppc.R")
   ppc_results <- run_all_ppc()
   ```

5. **Compare models**:
   ```R
   source("analysis/utils/model_comparison.R")
   comparison <- compare_models()
   ```

6. **Test parameter recovery**:
   ```R
   source("analysis/utils/parameter_recovery.R")
   recovery <- run_parameter_recovery("pvl_delta", n_subj = 10)
   ```

## Requirements

Required R packages:
- `rjags` (JAGS interface)
- `coda` (MCMC diagnostics)
- `dplyr` (data manipulation)

JAGS must be installed separately: https://mcmc-jags.sourceforge.io/

## Known Issues and Limitations

1. **Data Loading Warnings**: The Steingroever dataset has some trial numbering inconsistencies across subjects (different trial counts). This is handled correctly by the code but generates warnings.

2. **Model Comparison Metrics**: Full implementation of WAIC/LOO requires computing pointwise log-likelihoods for each observation. The framework is in place but requires additional computation during model fitting.

3. **VSE Model Simulation**: Posterior predictive checks are not yet implemented for the VSE model (more complex due to perseverance dynamics).

4. **Computational Cost**: Fitting all models to the full dataset (173 subjects, 81831 trials) is computationally intensive. For testing, use `fit_all_studies = FALSE` in the config.

5. **JAGS Model Syntax**: The models use `ifelse()` statements which work in JAGS but may be less efficient than explicit loops. This tradeoff favors code clarity.

## Model Comparison Notes

When comparing the three models:

- **PVL-Delta** is the simplest (4 parameters)
- **VSE/VPP** adds sequential exploration/perseverance (8 parameters)
- **ORL** uses different learning mechanisms (5 parameters)

The VSE model has the most parameters and may overfit. Use cross-validation or information criteria to assess out-of-sample prediction.

## Future Enhancements

Potential improvements not in the initial scope:

1. Implement WAIC/LOO-CV using the `loo` package
2. Add parallel chain execution for faster fitting
3. Implement simulation-based calibration (SBC) for model validation
4. Add visualization functions for parameter estimates
5. Implement cross-validation for model selection
6. Add model recovery analysis (fit wrong model to simulated data)

## File Structure Summary

```
analysis/
├── README.md                       # Main documentation
├── IMPLEMENTATION_NOTES.md         # This file
├── fit_models.R                    # Main fitting script
├── test_data_loading.R             # Data loading test
├── .gitignore                      # Ignore outputs
├── models/                         # JAGS model definitions
│   ├── pvl_delta.jags
│   ├── vse.jags
│   └── orl.jags
├── utils/                          # Helper functions
│   ├── load_data.R
│   ├── prepare_jags_data.R
│   ├── diagnostics.R
│   ├── ppc.R
│   ├── model_comparison.R
│   └── parameter_recovery.R
└── outputs/                        # Results (gitignored)
    └── .gitkeep
```

## Testing Status

- ✓ Data loading tested across all 8 datasets
- ✓ JAGS data preparation validated
- ✓ Model syntax verified (JAGS parseable)
- ⚠ Full model fitting not tested (requires ~hours of computation)
- ⚠ Diagnostics/PPC require fitted models

## References

Model implementations based on:
- **PVL-Delta**: Ahn et al. (2014). Computational Psychiatry
- **VSE/VPP**: Haines et al. (2018). Psychonomic Bulletin & Review
- **ORL**: Haines et al. (2018). Psychonomic Bulletin & Review

Data sources:
- Ahn, W. Y., et al. (2014). Frontiers in Psychology
- Fridberg, D. J., et al. (2010). Behavioral Neuroscience
- Steingroever, H., et al. (2014). Journal of Problem Gambling

Generated: 2026-01-02
