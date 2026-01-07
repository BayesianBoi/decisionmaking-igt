# Memory Decay in Reinforcement Learning

A Hierarchical Bayesian analysis investigating whether substance use disorders are associated with elevated memory decay rates during reinforcement learning in the Iowa Gambling Task (IGT).

## Research Questions

**Primary:** Do individuals with substance use disorders exhibit elevated memory decay rates in reinforcement learning compared to healthy controls?

**Secondary:**
1. Can impaired IGT performance in substance use disorders be explained by memory retention deficits (high λ) rather than learning acquisition deficits (low A)?
2. Does the EEF model (incorporating forgetting) provide a superior fit compared to models assuming perfect memory (PVL-Delta) or asymmetric reward/punishment processing (ORL)?
3. Do different substance classes (stimulants vs. opioids) demonstrate distinct forgetting profiles?

## Computational Models

Three reinforcement learning models were evaluated:

| Model | Parameters | Key Mechanism | Reference |
|-------|-----------|-------------|-----------|
| **EEF** | 4 | Exploration and Exploitation with Forgetting (λ) | Yang et al. (2025) |
| **PVL-Delta** | 4 | Prospect Theory + Delta Learning (Perfect Memory) | Ahn et al. (2008) |
| **ORL** | 5 | Outcome-Representation Learning (Asymmetric Weights) | Haines et al. (2018) |

## Data

This project utilizes the dataset from Ahn et al. (2014), comprising decision-making data from three distinct groups:

| Group | N | Description |
|-------|---|-------------|
| **Healthy Controls** | 48 | No history of substance dependence |
| **Amphetamine** | 38 | History of amphetamine dependence (abstinent > 1 month) |
| **Heroin** | 43 | History of heroin dependence (abstinent > 1 month) |
| **Total** | **129** | |

## Requirements

- R >= 4.0
- JAGS >= 4.3
- R packages: `R2jags`, `rjags`, `coda`, `ggplot2`, `gridExtra`, `pacman`

## Usage

### 1. Model Fitting
Models are fit using hierarchical Bayesian estimation via JAGS. Scripts are located in `scripts/fitting/`.

```bash
# Example: Fit ORL model for Healthy Controls
nohup Rscript scripts/fitting/fit_orl.R HC > logs/orl_HC.log 2>&1 &

# Fit other groups/models
nohup Rscript scripts/fitting/fit_eef.R Amph > logs/eef_Amph.log 2>&1 &
```

### 2. Posterior Predictive Checks (PPC)
After fitting, verify model performance using the MPD (Maximum Posterior Density) approach through the scripts in `scripts/ppc/`.

```bash
Rscript scripts/ppc/run_ppc_mpd.R orl HC
Rscript scripts/ppc/plot_ppc_combined.R  # Generate summary plots
```

### 3. Analysis and Visualization
Generate parameter recovery plots and group comparisons:

```bash
Rscript scripts/recovery/recovery_orl.R
Rscript scripts/plotting/plot_group_posteriors.R
```

## Repository Structure

```
decision_making/
├── models/                  # JAGS model definitions (.txt)
├── scripts/
│   ├── fitting/             # Hierarchical model fitting scripts
│   ├── recovery/            # Parameter recovery simulations
│   ├── ppc/                 # Posterior predictive checks
│   ├── group_comparison/    # Group difference analysis
│   └── plotting/            # Visualization utilities
├── utils/                   # Shared R functions (data loading, etc.)
├── outputs/                 # Model fits and results (not tracked)
├── plots/                   # Generated figures
└── data/                    # Raw IGT data
```

## References

- Ahn, W.Y., et al. (2008). Comparison of decision learning models. *Cognitive Science*.
- Ahn, W.Y., et al. (2014). Decision-making in stimulant and opiate addicts. *Neuropsychopharmacology*.
- Haines, N., Vassileva, J., & Ahn, W.Y. (2018). The Outcome-Representation Learning model. *Cognitive Science*.
- Yang, X., et al. (2025). Exploitation and Exploration with Forgetting. *Frontiers in Psychology*.
