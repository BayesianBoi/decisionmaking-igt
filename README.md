# Iowa Gambling Task - Hierarchical Bayesian Models

Reproducible R + JAGS analysis pipeline for fitting decision-making models to Iowa Gambling Task (IGT) data.

## Overview

This repository implements three hierarchical Bayesian models for IGT choice data:

- **PVL-Delta**: Prospect-Valence Learning with delta-rule updating (4 parameters)
- **VSE**: Value + Sequential Exploration with perseverance (8 parameters)
- **ORL**: Outcome Representation Learning with fictive updating (5 parameters)

All models are implemented in JAGS with numerically stable priors and validated on 173 subjects from 3 published studies.

## Quick Start

### Local Testing (2 minutes)

```bash
cd decision_making
Rscript analysis/quick_test.R
```

### Cloud Deployment (3-5 hours)

See [`GITHUB_DEPLOYMENT.md`](GITHUB_DEPLOYMENT.md) for GitHub setup, then follow [`analysis/CLOUD_SETUP.md`](analysis/CLOUD_SETUP.md).

## Repository Structure

```
igt-decision-models/
├── analysis/              # Main analysis pipeline
│   ├── fit_models.R      # Run all models
│   ├── quick_test.R      # Quick validation
│   ├── models/           # JAGS model definitions
│   │   ├── pvl_delta_v2.jags
│   │   ├── vse_v2.jags
│   │   └── orl_v2.jags
│   └── utils/            # Helper functions
├── data/                 # Raw IGT data
│   └── raw/
│       ├── Ahn_2014/
│       ├── Fridberg_2010/
│       └── Steingroever_2014/
└── docs/                 # Reference papers
```

## Data

**173 subjects, 81,831 trials** from:

- Ahn et al. (2014): 3 groups (HC, Amphetamine, Heroin)
- Fridberg et al. (2010): 2 groups (HC, Cannabis)
- Steingroever et al. (2014): 3 datasets (95, 100, 150 trials)

## Models

### PVL-Delta
- **A**: Learning rate [0, 1]
- **alpha**: Outcome sensitivity [0, 2]
- **cons**: Choice consistency [0, 5]
- **lambda**: Loss aversion [0, 10]

### VSE (extends PVL-Delta)
- **epP/epN**: Perseverance after wins/losses
- **K**: Perseverance decay
- **w**: Weight between value and perseverance

### ORL
- **Arew/Apun**: Reward/punishment learning rates
- **K**: Perseverance decay
- **betaF/betaP**: Frequency and perseverance weights

## Requirements

- R >= 4.0
- JAGS >= 4.3
- R packages: `rjags`, `coda`, `dplyr`

## Installation

```bash
# Ubuntu/Debian
sudo apt-get install r-base jags
sudo R -e "install.packages(c('rjags', 'coda', 'dplyr'))"

# macOS
brew install jags
R -e "install.packages(c('rjags', 'coda', 'dplyr'))"
```

## Usage

### 1. Test Installation

```bash
Rscript analysis/quick_test.R
```

Expected output: Parameter estimates and R-hat values ~1.0

### 2. Run Full Pipeline

```bash
Rscript analysis/fit_models.R
```

Expected runtime: 3-5 hours on 8-core machine

### 3. Diagnostics

```R
source("analysis/utils/diagnostics.R")
diagnostics <- run_full_diagnostics()
```

### 4. Model Comparison

```R
source("analysis/utils/model_comparison.R")
comparison <- compare_models()
```

## Documentation

- [`GITHUB_DEPLOYMENT.md`](GITHUB_DEPLOYMENT.md) - GitHub setup guide
- [`analysis/README.md`](analysis/README.md) - Detailed analysis documentation
- [`analysis/CLOUD_SETUP.md`](analysis/CLOUD_SETUP.md) - Cloud deployment guide
- [`analysis/TROUBLESHOOTING.md`](analysis/TROUBLESHOOTING.md) - Common issues and fixes
- [`analysis/FINAL_SUMMARY.md`](analysis/FINAL_SUMMARY.md) - Project overview

## Features

- ✅ Numerically stable JAGS models (v2 versions)
- ✅ Validated on real IGT data
- ✅ Complete diagnostic pipeline (R-hat, ESS, trace plots)
- ✅ Posterior predictive checks
- ✅ Parameter recovery validation
- ✅ Cloud deployment ready
- ✅ Comprehensive documentation

## Citation

If you use this code, please cite the original model papers:

- **PVL-Delta**: Ahn et al. (2014). Computational Psychiatry
- **VSE**: Haines et al. (2018). Psychonomic Bulletin & Review
- **ORL**: Haines et al. (2018). Psychonomic Bulletin & Review

Data from:
- Ahn, W. Y., et al. (2014). Frontiers in Psychology
- Fridberg, D. J., et al. (2010). Behavioral Neuroscience
- Steingroever, H., et al. (2014). Journal of Problem Gambling

## License

[Add your license here]

## Contact

[Add your contact information]

---

**Status**: ✅ Ready for deployment
**Last Updated**: 2026-01-02
