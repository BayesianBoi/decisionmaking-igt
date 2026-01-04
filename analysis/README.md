# Analysis Pipeline

## Scripts

| Script | Purpose |
|--------|---------|
| `fit_eef.R` | Fit EEF model (forgetting mechanism) |
| `fit_pvl_delta.R` | Fit PVL-Delta baseline model |
| `fit_orl.R` | Fit ORL model (asymmetric learning) |
| `parameter_recovery.R` | Validate model identifiability (Simulate & Recover) |
| `compare_groups.R` | Compare forgetting rates between groups |
| `compare_models.R` | Compare model fit across EEF, PVL-Delta, ORL |
| `create_figures.R` | Generate publication figures |

## Models

| File | Parameters | Description |
|------|-----------|-------------|
| `eef.jags` | θ, λ, φ, cons | Exploitation-Exploration with Forgetting |
| `pvl_delta.jags` | A, α, cons, λ | Prospect Valence Learning with delta rule |
| `orl.jags` | A_rew, A_pun, K, βF, βP | Outcome Representation Learning |

## Utilities

| File | Purpose |
|------|---------|
| `load_data.R` | Load and harmonise IGT datasets |
| `prepare_eef_data.R` | Prepare data for EEF model (includes first-choice priors) |
| `prepare_jags_data.R` | Prepare data for PVL-Delta and ORL models |

## MCMC Configuration

All models use identical settings for comparability:

- Adaptation: 5,000 iterations
- Burn-in: 10,000 iterations
- Sampling: 20,000 iterations per chain
- Chains: 4 (parallel)
- Thinning: 2

## Output

Results are saved to `results/` subdirectories:

- `results/eef/`
- `results/pvl_delta/`
- `results/orl/`
- `results/group_comparison/`
- `results/model_comparison/`
- `results/figures/`
