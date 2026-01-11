# Memory Decay in Reinforcement Learning

Hierarchical Bayesian analysis of whether substance users show elevated memory decay during reinforcement learning on the Iowa Gambling Task.

## Models

| Model | Params | Key Mechanism | Reference |
|-------|--------|---------------|-----------|
| **EEF** | 4 | Forgetting rate (λ) + exploration bonus | Yang et al. (2025) |
| **PVL-Delta** | 4 | Prospect theory + delta learning | Steingroever et al. (2013) |
| **ORL** | 5 | Separate reward/punishment learning | Haines et al. (2018) |

## Data

IGT data from Ahn et al. (2014): 48 healthy controls, 38 amphetamine users, and 43 heroin users (all abstinent >1 month).

## Results

EEF fits clinical groups best, while ORL is slightly better for healthy controls. No group differences in forgetting rate.

![Model Comparison](figures/paper/dic_model_comparison.png)

## Quick Start

```bash
# fit a model
Rscript scripts/fitting/fit_eef.R HC

# posterior predictive check
Rscript scripts/ppc/run_ppc.R eef HC

# parameter recovery
Rscript scripts/recovery/recovery_eef.R
```

## Structure

```
├── models/           # JAGS model files
├── scripts/
│   ├── fitting/      # model estimation
│   ├── recovery/     # parameter recovery
│   ├── ppc/          # predictive checks
│   └── plotting/     # figures
├── data/
│   ├── raw/          # Ahn 2014 data
│   └── processed/    # fits, ppc, recovery
└── figures/paper/    # publication plots
```

## Requirements

R ≥ 4.0, JAGS ≥ 4.3, and packages: `R2jags`, `coda`, `ggplot2`, `dplyr`, `HDInterval`

## References

- Ahn et al. (2014). Decision-making in stimulant and opiate addicts. *Neuropsychopharmacology*
- Haines et al. (2018). The ORL model. *Cognitive Science*
- Steingroever et al. (2013). Comparing RL models. *J Math Psychol*
- Yang et al. (2025). Exploitation and Exploration with Forgetting. *Front Psychol*
