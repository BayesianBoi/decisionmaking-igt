# Memory Decay in Reinforcement Learning

Hierarchical Bayesian analysis of Iowa Gambling Task data investigating whether substance use disorders are associated with elevated memory decay rates during reinforcement learning.

## Research Questions

**Primary:** Do individuals with substance use disorders exhibit elevated memory decay rates in reinforcement learning compared to healthy controls?

**Secondary:**
1. Can poor IGT performance in substance users be explained by memory retention deficits (high λ) rather than learning acquisition deficits (low A)?
2. Does the EEF model (with forgetting) provide better fit than the PVL-Delta model (assumes perfect memory)?
3. Do different substance classes (stimulants, opioids, cannabis) show differential forgetting profiles?

## Models

Three computational models of reinforcement learning are compared:

| Model | Parameters | Key Feature | Reference |
|-------|-----------|-------------|-----------|
| **EEF** | 4 | Memory decay (λ) | Yang et al. (2025) |
| **PVL-Delta** | 4 | Perfect memory | Ahn et al. (2008) |
| **VSE** | 8 | Perseverance | Worthy et al. (2013) |

## Data

Clinical populations from published studies:

| Study | Group | N |
|-------|-------|---|
| Ahn et al. (2014) | Healthy Controls | 48 |
| Ahn et al. (2014) | Amphetamine | 38 |
| Ahn et al. (2014) | Heroin | 43 |
| Fridberg et al. (2010) | Healthy Controls | 15 |
| Fridberg et al. (2010) | Cannabis | 17 |
| **Total** | | **161** |

## Requirements

- R >= 4.0
- JAGS >= 4.3
- R packages: rjags, coda, ggplot2

## Usage

```bash
# Parameter recovery (validate model identifiability)
Rscript analysis/scripts/parameter_recovery_eef.R

# Fit models
Rscript analysis/scripts/fit_eef_clinical.R
Rscript analysis/scripts/fit_pvl_delta.R
Rscript analysis/scripts/fit_vse.R

# Analysis
Rscript analysis/scripts/compare_groups.R
Rscript analysis/scripts/compare_models.R
Rscript analysis/scripts/create_figures.R
```

## Repository Structure

```
├── analysis/
│   ├── scripts/           # Analysis scripts
│   ├── models/            # JAGS model definitions
│   │   ├── eef_clinical.jags
│   │   ├── pvl_delta.jags
│   │   └── vse.jags
│   └── utils/             # Helper functions
├── data/raw/              # Source data
└── results/               # Output (not tracked)
```

## References

- Ahn, W.Y., Busemeyer, J.R., & Wagenmakers, E.J. (2008). Comparison of decision learning models. *Cognitive Science*.
- Ahn, W.Y., et al. (2014). Decision-making in stimulant and opiate addicts. *Neuropsychopharmacology*.
- Fridberg, D.J., et al. (2010). Cognitive mechanisms underlying risky decision-making in chronic cannabis users. *Journal of Mathematical Psychology*.
- Worthy, D.A., et al. (2013). Decomposing the roles of perseveration and expected value. *Frontiers in Psychology*.
- Yang, X., et al. (2025). Exploitation and Exploration with Forgetting. *Frontiers in Psychology*.
