# Iowa Gambling Task Analysis Pipeline

This directory contains a reproducible R + JAGS analysis pipeline for fitting decision-making models to Iowa Gambling Task (IGT) data.

## Quick Start

**Test the pipeline (2 minutes):**
```bash
Rscript analysis/quick_test.R
```

**Run full pipeline (3-5 hours on cloud):**
```bash
Rscript analysis/fit_models.R
```

See `CLOUD_SETUP.md` for running on cloud compute.

## Overview

This pipeline fits three reinforcement learning models (PVL-Delta, VSE, and ORL) to IGT data from three published studies. The models are implemented in JAGS with hierarchical Bayesian estimation.

## Data Sources

### 1. Ahn et al. (2014)
- **Location**: `data/raw/Ahn_2014/`
- **Files**:
  - `IGTdata_HC.txt` - Healthy controls
  - `IGTdata_Amph.txt` - Amphetamine group
  - `IGTdata_Hero.txt` - Heroin group
- **Format**: Tab-delimited with columns: `trial`, `deck`, `gain`, `loss`, `subjID`
- **Structure**: Long format with one row per trial, deck choices 1-4

### 2. Fridberg et al. (2010)
- **Location**: `data/raw/Fridberg_2010/`
- **Files**:
  - `IGTdata_HC.txt` - Healthy controls
  - `IGTdata_Cbis.txt` - Cannabis users
- **Format**: Tab-delimited with columns: `trial`, `deck`, `deckCopy`, `gain`, `loss`, `subjID`
- **Structure**: Long format with one row per trial, deck choices 1-4

### 3. Steingroever et al. (2014)
- **Location**: `data/raw/Steingroever_2014/`
- **Files**: Separate files for choice, wins (wi), losses (lo), and subject indices
  - `choice_95.txt/csv`, `choice_100.txt/csv`, `choice_150.txt/csv`
  - `wi_95.txt/csv`, `wi_100.txt/csv`, `wi_150.txt/csv`
  - `lo_95.txt/csv`, `lo_100.txt/csv`, `lo_150.txt/csv`
  - `index_95.txt/csv`, `index_100.txt/csv`, `index_150.txt/csv`
- **Format**: Wide format with subjects as rows and trials as columns
- **Structure**: Three datasets with 95, 100, and 150 trials per subject

## Models

### PVL-Delta Model
- **Parameters**: A (learning rate), alpha (outcome sensitivity), cons (choice consistency), lambda (loss aversion)
- **Description**: Prospect-valence learning with delta-rule updating and power utility function

### VSE Model (Value + Sequential Exploration)
- **Parameters**: A, alpha, cons, lambda, epP (positive perseverance), epN (negative perseverance), K (perseverance decay), w (weight)
- **Description**: PVL-Delta extended with sequential exploration (perseverance) tendency

### ORL Model (Outcome Representation Learning)
- **Parameters**: Arew (reward learning rate), Apun (punishment learning rate), K (perseverance decay), betaF (frequency weight), betaP (perseverance weight)
- **Description**: Separate learning for outcome valence and frequency, with fictive updating

## Directory Structure

```
analysis/
├── README.md                  # This file
├── fit_models.R              # Main fitting pipeline
├── models/                   # JAGS model files
│   ├── pvl_delta.jags
│   ├── vse.jags
│   └── orl.jags
├── utils/                    # Helper functions
│   ├── load_data.R          # Data loading and harmonization
│   ├── prepare_jags_data.R  # JAGS data preparation
│   ├── diagnostics.R        # R-hat, ESS, trace plots
│   ├── ppc.R                # Posterior predictive checks
│   └── model_comparison.R   # DIC/WAIC/LOO comparison
└── outputs/                 # Results (not tracked in git)
    ├── pvl_delta_fit.rds
    ├── vse_fit.rds
    ├── orl_fit.rds
    ├── diagnostics.rds
    ├── ppc_results.rds
    └── model_comparison.md
```

## Requirements

- R >= 4.0
- rjags package
- dplyr, tidyr for data manipulation
- coda for MCMC diagnostics

## Notes

- All models use hierarchical Bayesian estimation with group-level and subject-level parameters
- Default MCMC settings: 4 chains, 2000 iterations burn-in, 5000 sampling iterations
- Data preprocessing maintains consistency with original publications
