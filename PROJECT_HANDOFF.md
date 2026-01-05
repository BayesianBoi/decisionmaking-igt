# Project Status & Handoff (Decision Making Pipeline)

**Date:** 2026-01-05
**Current Phase:** Phase 1 (Parameter Recovery Validation) - **RUNNING ON CLOUD**

## 1. Project Goal
Develop and validate a robust behavioral modeling pipeline for the Iowa Gambling Task (IGT), focusing on three models:
1.  **PVL-Delta** (Prospect Valence Learning with Delta rule)
2.  **ORL** (Outcome-Representation Learning)
3.  **EEF** (Explore-Exploit-Forget) - *Custom/Derived implementation*

The ultimate goal is to fit these models to empirical data (Ahn et al., 2014) and compare their performance.

## 2. Methodology & Key Decisions
We have strictly aligned the pipeline with **Ahn et al. (2014)** practices:

### A. Payoff Scheme (`analysis/utils/payoff_scheme.R`)
-   **Decision:** Implemented the **EXACT** hardcoded payoff schedule from Ahn 2014 (Modified IGT).
-   **Details:** Fixed loss magnitudes (e.g., Deck B has -2500 loss at trial 10 within block) rather than probabilistic. Match is precise.
-   **Indexing:** **Deck-based**, meaning outcomes depend on how many times a *specific deck* was chosen, not the global trial count.

### B. Model Specifications (`analysis/models/`)
-   **ORL (`orl.txt`):** Perseverance (`PS`) initialized to **0**. (Matches inspiration code).
-   **PVL-Delta (`pvl_delta.txt`):** Hierarchical implementation.
-   **EEF (`eef.txt`):** Custom model. Note that strict "Exploit" probability logic is used.

### C. Parameter Recovery
-   **Sample Size:** N=48 (Simulated). **Reason:** Matches the "Healthy Control" group size in Ahn 2014.
-   **Parallelization:** Expert implementation using `mclapply` (outer loop) and sequential JAGS chains. Optimizes cloud core usage.
-   **Metrics:** We check Pearson correlation (`r`) between True and Recovered (mean posterior) parameters.

## 3. Repository Structure
The code was refactored for clarity. **Do not use the old `analysis/scripts` folder.**

| Folder | Contents |
| :--- | :--- |
| `analysis/1_analysis/0_Fitting/` | Scripts to fit models to **Empirical Data** (Next Phase). |
| `analysis/1_analysis/2_Recovery/` | **(Current Focus)** Scripts to simulate & recover parameters (`recovery_*.R`). |
| `analysis/models/` | JAGS model definitions (`.txt`). |
| `analysis/utils/` | Helpers: `payoff_scheme.R`, `plotting_utils.R`. |
| `data/raw/Ahn_2014/` | The empirical dataset. |

**Inspiration Code:** Located in `R/` (Legacy/Reference files from other researchers).

## 4. Current Status (What is running?)
We launched the recovery jobs on the cloud:
-   `recovery_pvl_delta.R`
-   `recovery_orl.R`
-   `recovery_eef.R`

These are running in the background (PIDs were ~8777).

## 5. Next Steps (For the Next Agent)
1.  **Check Results:**
    -   Wait for cloud jobs to finish (~6-12 hours).
    -   Download plots from `analysis/plots/recovery/`.
    -   **Criteria:** Look for correlations $r > 0.7$ (ideally >0.8) and points along the $y=x$ line.
2.  **Phase 2: Empirical Fitting:**
    -   Once recovery is valid, run the scripts in `analysis/1_analysis/0_Fitting/`.
    -   Ensure they load the **Ahn 2014** data correctly using the `payoff_scheme` logic.
3.  **Model Comparison:**
    -   Compare LOOIC / WAIC scores between models.
    -   Generate Posterior Predictive Checks (PPC).

## 6. Known "Gotchas"
-   **JAGS Installation:** On cloud, `apt-get install jags` often fails. **Solution:** Use Conda (`conda install -c conda-forge jags`). (See `CLOUD_DEPLOYMENT_GUIDE.md`).
-   **Git:** Recovery scripts are in `analysis/1_analysis/2_Recovery`. If missing, pull `main`.

---
*Created by Antigravity Agent, Jan 2026.*
