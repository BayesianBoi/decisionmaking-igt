# Analysis Pipeline

## Scripts

| Script | Purpose |
|--------|---------|
| `fit_eef_clinical.R` | Fit EEF model to clinical populations |
| `fit_pvl_delta.R` | Fit PVL-Delta baseline model |
| `fit_vse.R` | Fit VSE model with perseverance |
| `parameter_recovery_eef.R` | Validate model identifiability |
| `compare_groups.R` | Compare forgetting rates between groups |
| `compare_models.R` | Compare model fit across EEF, PVL-Delta, VSE |
| `create_figures.R` | Generate publication figures |

## Models

| File | Parameters | Description |
|------|-----------|-------------|
| `eef_clinical.jags` | theta, lambda_forget, phi, cons | Exploitation-Exploration with Forgetting |
| `pvl_delta.jags` | A, alpha, cons, lambda | Prospect Valence Learning with delta rule |
| `vse.jags` | A, alpha, cons, lambda, epP, epN, K, w | Value plus Sequential Exploration |

## Utilities

| File | Purpose |
|------|---------|
| `load_data.R` | Load and harmonize IGT datasets |
| `prepare_eef_data.R` | Prepare data for EEF model (includes first-choice priors) |
| `prepare_jags_data.R` | Prepare data for PVL-Delta and VSE models |
| `diagnostics.R` | Compute convergence diagnostics |

## MCMC Configuration

All models use identical settings for comparability:

- Adaptation: 5,000 iterations
- Burn-in: 10,000 iterations
- Sampling: 20,000 iterations per chain
- Chains: 4
- Thinning: 2

Total effective samples: ~40,000

## Output

Results are saved to `results/` subdirectories:

- `results/eef_clinical/`
- `results/pvl_delta/`
- `results/vse/`
- `results/group_comparison/`
- `results/model_comparison/`
- `results/figures/`
