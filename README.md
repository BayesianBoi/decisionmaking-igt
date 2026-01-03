# Memory Decay in Reinforcement Learning: A Hierarchical Bayesian Analysis of the Iowa Gambling Task

This repository provides the analysis pipeline for investigating whether substance users exhibit elevated memory decay rates compared to healthy controls in the Iowa Gambling Task.

## Research Questions

**Primary RQ:** Do individuals with substance use disorders exhibit elevated memory decay rates in reinforcement learning compared to healthy controls?

**Secondary RQs:**
1. Can poor IGT performance be explained by memory deficits (high λ) rather than learning deficits (low A)?
2. Does EEF (with forgetting) outperform PVL-Delta (perfect memory) and VSE (strategic perseverance)?
3. Do different substance classes (stimulants, opioids, cannabis) show differential forgetting profiles?

## Models

The EEF (Exploitation-Exploration with Forgetting) model is the primary model, based on Yang et al. (2025). We also fit PVL-Delta and VSE for model comparison.

**EEF** (4 parameters)
- theta: Outcome sensitivity
- lambda_forget: Memory decay rate (key parameter of interest)
- phi: Exploration incentive
- cons: Choice consistency

**PVL-Delta** (4 parameters)
Baseline model assuming perfect memory retention.

**VSE** (8 parameters)
Alternative model with perseverance mechanisms.

## Data

The pipeline analyzes 146 subjects from clinical populations:
- Ahn et al. (2014): Healthy controls (n=48), amphetamine users (n=38), heroin users (n=43)
- Fridberg et al. (2010): Cannabis users (n=17)

Total: 14,600 trials across 4 groups.

## Quick Start

**1. Parameter Recovery (verify model identifiability)**
```bash
Rscript analysis/scripts/parameter_recovery_eef.R
```

**2. Fit EEF Model (4-8 hours)**
```bash
Rscript analysis/scripts/fit_eef_clinical.R
```

**3. Group Comparisons**
```bash
Rscript analysis/scripts/compare_groups.R
```

**4. Generate Figures**
```bash
Rscript analysis/scripts/create_figures.R
```

## Requirements

- R >= 4.0
- JAGS >= 4.3
- R packages: rjags, coda, ggplot2

Installation:
```bash
# macOS
brew install jags
R -e "install.packages(c('rjags', 'coda', 'ggplot2'))"

# Ubuntu
sudo apt-get install jags
R -e "install.packages(c('rjags', 'coda', 'ggplot2'))"
```

## Repository Structure

```
├── analysis/
│   ├── scripts/                 # Main analysis scripts
│   │   ├── fit_eef_clinical.R   # Fit EEF model
│   │   ├── fit_pvl_delta.R      # Fit PVL-Delta baseline
│   │   ├── fit_vse.R            # Fit VSE model
│   │   ├── compare_groups.R     # Group comparisons
│   │   ├── compare_models.R     # Model comparisons
│   │   ├── create_figures.R     # Publication figures
│   │   └── parameter_recovery_eef.R
│   ├── models/                  # JAGS model definitions
│   │   ├── eef_clinical.jags    # Primary model
│   │   ├── pvl_delta.jags
│   │   ├── vse.jags
│   │   └── orl.jags
│   └── utils/                   # Helper functions
│       ├── load_data.R
│       ├── prepare_eef_data.R
│       └── prepare_jags_data.R
├── data/raw/                    # Source data files
├── docs/                        # Reference materials
│   └── PRIOR_JUSTIFICATIONS.md  # Prior specification documentation
└── results/                     # Output directory
```

## Prior Specifications

All prior choices are documented with justifications in `docs/PRIOR_JUSTIFICATIONS.md`. Key priors:

| Parameter | Prior | Justification |
|-----------|-------|---------------|
| lambda_forget | Beta(2, 3) | Centered on Yang's empirical range (0.4-0.5) |
| theta | Beta(1.5, 3) | Allows risk-seeking (typical in IGT) |
| phi | Normal(0, 2) | Allows negative values (exploration aversion) |
| cons | Normal(1, 2)T(0, 5) | Standard IGT range |

## MCMC Settings

- 4 chains
- 5000 adaptation iterations
- 10000 burn-in iterations
- 20000 sampling iterations
- Thinning = 2

Total effective samples: ~40,000

## Citations

**Model:**
- Yang, X., et al. (2025). Exploitation and Exploration with Forgetting.

**Data:**
- Ahn, W. Y., et al. (2014). Frontiers in Psychology
- Fridberg, D. J., et al. (2010). Behavioral Neuroscience

## License

This project is for academic research purposes.
