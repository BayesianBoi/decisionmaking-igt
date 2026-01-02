# Coding Agent Prompt: Decision-Making Models in R + JAGS

## Role and Goal
You are a coding agent working in a manuscript-first repo. Your goal is to implement a clean, reproducible R/JAGS analysis pipeline for Iowa Gambling Task data that fits three models: **PVL-Delta**, **VSE**, and **ORL**. The code should be modular, well-documented, and easy to run without a build system.

## Constraints and Project Rules
- Treat `data/` and `docs/` as read-only. Write any derived files to a new folder (e.g., `analysis/` or `outputs/`).
- Use clear scientific code style and consistent terminology (IGT, PVL-Delta, VSE, ORL).
- If you introduce scripts, document the exact run command at the top of each script (e.g., `Rscript analysis/fit_models.R`).
- Keep outputs out of version control unless explicitly requested.
- We have the code from the papers that we have taken the data from (Found in `R/`, and subfolders in here). These are written in STAN (And we are writing in JAGS), but you can use those as inspiration

## Deliverables
1. An `analysis/` directory containing:
   - `analysis/README.md` with run instructions.
   - `analysis/fit_models.R` that loads data, prepares inputs, and runs JAGS for all models.
   - `analysis/models/` containing JAGS model files:
     - `pvl_delta.jags`
     - `vse.jags`
     - `orl.jags`
   - `analysis/utils/` with helper functions for data prep, parameter recovery, PPC, and diagnostics.
2. Minimal example output (e.g., saved summaries as `.rds`), unless told to skip outputs.

## Detailed Plan (Follow This Exactly)
1. **Repo Audit and Data Understanding**
   - List files in `data/raw` and identify IGT datasets (Ahn_2014, Fridberg_2010, Steingroever_2014).
   - Inspect data formats and required variables (subject id, trial, deck, outcome, gain/loss, etc.).
   - Write a short summary of each dataset’s structure in `analysis/README.md`.

2. **Project Skeleton**
   - Create `analysis/`, `analysis/models/`, `analysis/utils/`.
   - Add a top-level `analysis/README.md` with the run command and a one-paragraph overview.

3. **Data Harmonization**
   - Implement `analysis/utils/load_data.R` to load each dataset and standardize columns:
     - `subj`, `trial`, `choice`, `reward`, `loss`, `net` (or equivalent), and `deck` if needed.
   - Implement `analysis/utils/prepare_jags_data.R` to create JAGS-ready arrays:
     - Subject-level trial counts, choice indices (1–4), outcomes per trial.
   - Include a check function that validates ranges and missing values.

4. **Model Definitions (JAGS)**
   - Write three model files:
     - **PVL-Delta**: utility function, delta-rule updates, softmax choice.
     - **VSE**: value learning plus sequential exploration tendency.
     - **ORL**: outcome-representation learning as per Haines et al. (2018).
   - Use parameter naming consistent with literature.
   - Provide reasonable priors and hierarchical structure (group-level + subject-level).

5. **Fitting Pipeline**
   - In `analysis/fit_models.R`, build a unified pipeline:
     - Load data, prepare JAGS inputs.
     - Fit each model with the same settings (chains, burn-in, thinning).
     - Save posterior samples as `.rds` to `analysis/outputs/`.

6. **Diagnostics and Model Checking**
   - Implement `analysis/utils/diagnostics.R` to compute R-hat, ESS, and trace plots.
   - Implement `analysis/utils/ppc.R` to run posterior predictive checks.
   - Save diagnostic outputs as `.rds` and optional plots to `analysis/outputs/`.

7. **Model Comparison**
   - Implement `analysis/utils/model_comparison.R` for DIC/WAIC or PSIS-LOO.
   - Write a short results summary to `analysis/outputs/model_comparison.md`.

8. **Parameter Recovery (Optional but Preferred)**
   - Implement a simulation routine:
     - Sample parameters from priors.
     - Simulate choices and outcomes.
     - Re-fit models and check recovery correlations.
   - Save results to `analysis/outputs/parameter_recovery.rds`.

9. **Final QA**
   - Verify scripts run end-to-end.
   - Add a short note in `analysis/README.md` describing any unresolved issues.

## Coding Quality Checklist
- Modular functions, clear naming, minimal global state.
- Deterministic behavior via set seeds.
- Comments only where needed (non-obvious math or design).
- Avoid hard-coded file paths; use relative paths.

## What to Ask the User If Anything Is Ambiguous
- Confirm data column names if not documented.
- Confirm whether to include parameter recovery and PPC plots.
- Confirm preferred model comparison metric (DIC/WAIC/LOO).

## Execution Start
Begin by listing files in `data/raw` and drafting `analysis/README.md` based on your findings.
