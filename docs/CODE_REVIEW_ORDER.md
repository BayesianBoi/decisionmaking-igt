# Human Code Review Order

This document provides the recommended order for manually reviewing the codebase. Start with the simplest/most foundational scripts and work up to the more complex ones.

---

## Tier 1: Data Foundation (Start Here)

These scripts define how data enters the pipeline. Review these first.

| Order | Script | Lines | Purpose |
|-------|--------|-------|---------|
| 1 | `analysis/utils/load_data.R` | 146 | Loads raw IGT data from Ahn (2014) and Steingroever (2014). Creates standardized format. |
| 2 | `analysis/utils/prepare_jags_data.R` | 148 | Converts harmonized data into matrices for JAGS. Handles subject indexing and scaling. |
| 3 | `analysis/0_preprocess/0_format_raw_data.R` | 36 | Entry point script that calls the above utilities. |

**Key things to verify:**
- [ ] Choice values are correctly coded as 1-4
- [ ] Outcomes are scaled by 100 for numerical stability
- [ ] Subject indices are contiguous

---

## Tier 2: JAGS Models (Core Math)

These are the probabilistic models. Most critical to understand correctly.

| Order | Script | Lines | Purpose |
|-------|--------|-------|---------|
| 4 | `analysis/models/pvl_delta.jags` | 113 | Prospect Valence Learning model. Baseline comparison. |
| 5 | `analysis/models/orl.jags` | 152 | Outcome-Representation Learning model. Tracks EV, EF, Pers. |
| 6 | `analysis/models/eef.jags` | 149 | Exploitation-Exploration with Forgetting. **Primary model of interest.** |

**Key things to verify:**
- [ ] Prior distributions match `docs/PRIOR_JUSTIFICATIONS.md`
- [ ] Update equations match referenced papers
- [ ] log_lik is computed correctly for model comparison

---

## Tier 3: Model Fitting

These scripts run MCMC inference. All use the same parallel pattern.

| Order | Script | Lines | Purpose |
|-------|--------|-------|---------|
| 7 | `analysis/1_analysis/0_Fitting/fit_pvl_delta.R` | 193 | Fits PVL-Delta to clinical data. |
| 8 | `analysis/1_analysis/0_Fitting/fit_orl.R` | 191 | Fits ORL to clinical data. |
| 9 | `analysis/1_analysis/0_Fitting/fit_eef.R` | 214 | Fits EEF to clinical data. Uses first-choice priors. |

**Key things to verify:**
- [ ] JAGS RNG seeds are explicitly set per chain (lines ~89-92)
- [ ] Burn-in and thinning are appropriate (10k burn, thin=5)
- [ ] log_lik is included for WAIC/LOO

---

## Tier 4: Utilities (Simulation & Diagnostics)

These support recovery, PPC, and model comparison.

| Order | Script | Lines | Purpose |
|-------|--------|-------|---------|
| 10 | `analysis/utils/plotting_utils.R` | 160 | Helper functions for plots (scale_rel, HDI, themes). |
| 11 | `analysis/utils/ppc.R` | 444 | Simulation functions for each model (simulate_pvl_delta, etc.). |
| 12 | `analysis/utils/prepare_eef_data.R` | ~300 | EEF-specific data prep including first-choice analysis. |

**Key things to verify:**
- [ ] Simulation functions match JAGS model equations
- [ ] HDI computation is standard (shortest interval)

---

## Tier 5: Plotting Scripts

Final outputs. Review after understanding the models.

| Order | Script | Lines | Purpose |
|-------|--------|-------|---------|
| 13 | `analysis/2_plotting/plot_figures.R` | 282 | Parameter recovery, learning curves, densities. |
| 14 | `analysis/2_plotting/plot_ppc.R` | 236 | Per-subject prediction accuracy plots. |
| 15 | `analysis/2_plotting/plot_group_compare.R` | 222 | Group comparisons (HC vs Amph vs Hero). |

**Key things to verify:**
- [ ] R² is computed correctly for recovery
- [ ] Group labels match metadata

---

## Tier 6: Analysis Scripts

Higher-level analyses that combine multiple components.

| Order | Script | Lines | Purpose |
|-------|--------|-------|---------|
| 16 | `analysis/1_analysis/2_Recovery/parameter_recovery.R` | ~300 | Full parameter recovery simulation. |
| 17 | `analysis/1_analysis/1_Simulation/0_simulate_all_models.R` | ~150 | Simulation for mechanism exploration. |

---

## Quick Reference: What Each Model Computes

| Model | Key Equation | Clinical Relevance |
|-------|--------------|-------------------|
| **PVL-Delta** | `u(x) = x^α` (gains), `-λ|x|^α` (losses) | Loss aversion, outcome sensitivity |
| **ORL** | EV + βF·EF + βP·Pers | Asymmetric learning, frequency sensitivity |
| **EEF** | exploit + explore, both decay by λ | **Memory decay hypothesis** |

---

## Files You Can Skip

These are secondary or not critical for understanding:
- `analysis/utils/diagnostics.R` - Standard convergence checks
- `analysis/utils/model_comparison.R` - WAIC/LOO wrappers
- `analysis/utils/reporting.R` - Table formatting
- `analysis/utils/visualization.R` - Legacy plotting utilities
