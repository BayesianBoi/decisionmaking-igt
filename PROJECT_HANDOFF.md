# Project Handoff: IGT Behavioral Modeling & Recovery

## 1. Project Overview
This project involves hierarchical Bayesian modeling of Iowa Gambling Task (IGT) data using three models:
- **PVL-Delta**: Prospect Valence Learning with Delta rule.
- **ORL**: Outcome-Representation Learning.
- **EEF**: Explore-Exploit-Forget (custom model).

The goal is to fit these models to clinical data (Ahn et al., 2014) and assess parameter recovery and group differences.

## 2. Current Status (As of Jan 6, 2026)
- **Parameter Recovery**: Scripts (`recovery_*.R`) are optimized (N=24, MCMC=100) and ready for/running on cloud.
- **Model Fitting**: Comparison models (`*_compare.txt`) are aligned with base models.
- **Verification**: The entire pipeline (Diagnostics -> PPC -> Plotting) has been verified.
- **Visualization**: Plotting scripts produce publication-quality figures.

## 3. Recent Achievements (Last Session)

### A. Model Consistency
- Aligned `*_compare.txt` models with base `*.txt` models.
- **PVL-Delta**: Fixed `a[s]` truncation (added `T(0,1)`) and removed bounds on `w` and `theta` to match base.
- **ORL**: Removed bounds on subject-level parameters to match base.

### B. Pipeline Verification
- **Diagnostics**: Validated `analysis/utils/diagnostics.R` (R-hat, ESS, trace plots).
- **PPC**: Validated `analysis/utils/ppc.R` simulation functions.
  - *Fix Applied:* Added missing closing brace to `simulate_pvl_delta`.
- **Plotting**: Created pseudo-data generators to verify plotting scripts without waiting for long runs:
  - `analysis/utils/generate_pseudo_ppc.R`
  - `analysis/utils/generate_pseudo_param_estimation.R` (inline in test script)

### C. Aesthetic Refinements
Refined all posterior density and parameter plots to meet publication standards (Example2 style):
- **Visuals**: Black line densities (no fill), minimal theme.
- **Annotations**: `N = ... Bandwidth = ...` added below plots.
- **Labels**: Greek letters used for parameters (e.g., $\mu_w$).
- **Verified Files**:
  - `analysis/2_plotting/test_ppc_plots.R`
  - `analysis/2_plotting/test_param_estimation_plots.R`

## 4. Key Files
- **Recovery Scripts**: `analysis/1_analysis/2_Recovery/recovery_*.R`
- **Plotting Utilities**: `analysis/utils/plotting_publication.R`
- **Verification Scripts**: 
  - `analysis/2_plotting/test_ppc_plots.R`
  - `analysis/2_plotting/test_param_estimation_plots.R`

## 5. Next Steps
1. **Run Full Analysis**: Execute the recovery and fitting scripts on the full dataset/cloud infrastructure.
2. **Generate Final Figures**: Once `results/` are populated with real data, run the plotting scripts to generate the final PDFs.
3. **Analyze Group Differences**: Use the validated `plot_group_compare.R` logic to assess clinical group differences.
