# Hierarchical Bayesian Models for Iowa Gambling Task

This repository provides a complete pipeline for fitting reinforcement learning models to Iowa Gambling Task data using hierarchical Bayesian estimation in JAGS.

## Models

Three models are implemented, each capturing different aspects of decision-making:

**PVL-Delta** (4 parameters)
Prospect-valence learning with delta-rule updating. Participants learn expected values through prediction errors, with utility shaped by outcome sensitivity and loss aversion.

**VSE** (8 parameters)
Extends PVL-Delta by adding perseverance mechanisms. Choices depend on both learned values and recent action tendencies, weighted by parameter w.

**ORL** (5 parameters)
Tracks both outcome magnitude (EV) and frequency (EF) with separate learning rates for wins and losses. Includes fictive learning for unchosen options.

## Data

The pipeline analyzes 173 subjects (81,831 trials) from three published studies:
- Ahn et al. (2014): Healthy controls, amphetamine users, heroin users
- Fridberg et al. (2010): Healthy controls, cannabis users
- Steingroever et al. (2014): Multiple trial lengths (95, 100, 150)

## Quick Start

Test the installation:
```bash
Rscript analysis/quick_test.R
```

Fit all models (3-5 hours on 8 cores):
```bash
Rscript analysis/fit_models.R
```

For systems with 32+ cores, you can fit models in parallel:
```bash
# In separate terminals
Rscript analysis/fit_single_model.R pvl_delta
Rscript analysis/fit_single_model.R vse
Rscript analysis/fit_single_model.R orl
```

See `analysis/PARALLEL_FITTING.md` for details on high-performance computing setups.

## Requirements

- R >= 4.0
- JAGS >= 4.3
- R packages: `rjags`, `coda`, `dplyr`, `parallel`

Installation on Ubuntu:
```bash
sudo apt-get install r-base jags
R -e "install.packages(c('rjags', 'coda', 'dplyr'))"
```

## Repository Structure

```
├── analysis/
│   ├── fit_models.R           # Main fitting script
│   ├── fit_single_model.R     # Fit individual models
│   ├── quick_test.R           # Validation test
│   ├── models/                # JAGS model definitions
│   │   ├── pvl_delta_v2.jags
│   │   ├── vse_v2.jags
│   │   └── orl_v2.jags
│   └── utils/                 # Analysis functions
│       ├── load_data.R
│       ├── diagnostics.R
│       ├── ppc.R
│       └── parameter_recovery.R
├── data/raw/                  # Source data files
└── docs/                      # Reference papers
```

## Validation

The pipeline includes comprehensive validation tools:

**Convergence diagnostics**: R-hat, effective sample size, trace plots
```r
source("analysis/utils/diagnostics.R")
run_full_diagnostics()
```

**Parameter recovery**: Verify identifiability with simulated data
```r
source("analysis/utils/parameter_recovery.R")
run_parameter_recovery("pvl_delta")
```

**Posterior predictive checks**: Test model fit to observed data
```r
source("analysis/utils/ppc.R")
run_all_ppc()
```

## Model Specifications

All models use weakly informative priors based on typical IGT parameter ranges. The v2 versions include numerical stability improvements:
- Outcome scaling (÷100) to prevent overflow
- Bounded parameters with appropriate truncation
- Beta distributions for rates, truncated normals for unbounded parameters

See individual model files for detailed prior specifications and theoretical documentation.

## Citations

Model implementations based on:
- Ahn, W. Y., et al. (2008). Journal of Neuroscience, Psychology, and Economics
- Worthy, D. A., et al. (2013). Psychonomic Bulletin & Review
- Haines, N., et al. (2018). Journal of Experimental Psychology: General

Data from:
- Ahn, W. Y., et al. (2014). Frontiers in Psychology
- Fridberg, D. J., et al. (2010). Behavioral Neuroscience
- Steingroever, H., et al. (2014). Journal of Problem Gambling

## Further Documentation

- `analysis/README.md` - Detailed pipeline documentation
- `analysis/PAPER_WORKFLOW.md` - Publication preparation guide
- `analysis/TROUBLESHOOTING.md` - Common issues and solutions
- `analysis/PARALLEL_FITTING.md` - High-performance computing guide
