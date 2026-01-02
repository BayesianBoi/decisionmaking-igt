# Paper/Manuscript Workflow Guide

Complete workflow from data to publication-ready outputs.

## Overview

This pipeline generates all figures, tables, and formatted results needed for manuscript preparation:

- ✅ **Trace plots** for convergence assessment
- ✅ **Posterior density plots** for parameter distributions
- ✅ **Forest plots** for model comparison
- ✅ **CSV tables** ready for manuscript
- ✅ **APA-formatted results** for copy-paste
- ✅ **LaTeX-ready tables** (optional)

## Workflow Steps

### 1. Fit Models (Cloud)

```bash
# On cloud machine
Rscript analysis/fit_models.R > fitting.log 2>&1
```

**Output:** 3 model fit files (~3-5 hours)
- `pvl_delta_fit.rds`
- `vse_fit.rds`
- `orl_fit.rds`

### 2. Download Results (Local)

```bash
# On cloud
tar -czf model_fits.tar.gz analysis/outputs/*_fit.rds

# Locally
scp user@cloud-ip:~/igt-decision-models/analysis/outputs/model_fits.tar.gz .
tar -xzf model_fits.tar.gz -C analysis/outputs/
```

### 3. Generate Publication Outputs (Local)

```bash
# This generates ALL figures and tables
Rscript analysis/generate_paper_outputs.R
```

**Outputs:**
```
analysis/outputs/publication/
├── figures/
│   ├── pvl_delta_traces.pdf       # Convergence checking
│   ├── pvl_delta_posteriors.pdf   # Parameter distributions
│   ├── vse_traces.pdf
│   ├── vse_posteriors.pdf
│   ├── orl_traces.pdf
│   ├── orl_posteriors.pdf
│   └── model_comparison_forest.pdf # Compare all models
├── tables/
│   ├── table1_parameter_estimates.csv
│   ├── table2_model_comparison.csv
│   └── table3_convergence.csv
├── RESULTS_SUMMARY.md            # Quick overview
├── APA_FORMATTED_RESULTS.txt     # Copy-paste to manuscript
└── MODEL_COMPARISON.md           # Comparison notes
```

### 4. Review Convergence

Check `tables/table3_convergence.csv`:

| Model      | Max_Rhat | Mean_Rhat | Min_ESS | Mean_ESS | Converged |
|------------|----------|-----------|---------|----------|-----------|
| pvl_delta  | 1.05     | 1.02      | 450     | 2500     | TRUE      |
| vse        | 1.08     | 1.03      | 380     | 2200     | TRUE      |
| orl        | 1.06     | 1.02      | 420     | 2400     | TRUE      |

**Criteria:**
- ✅ R-hat < 1.1 for all parameters
- ✅ ESS > 400 (ideally > 1000)
- ✅ Converged = TRUE

**If not converged:**
- Re-run with more iterations
- Check trace plots for issues
- See TROUBLESHOOTING.md

### 5. Review Figures

**Trace Plots** (`*_traces.pdf`):
- Should show "fuzzy caterpillar" pattern
- All chains should overlap
- No trends or patterns

**Posterior Plots** (`*_posteriors.pdf`):
- Smooth distributions (not multimodal)
- Reasonable ranges
- Consistent with literature

**Forest Plot** (`model_comparison_forest.pdf`):
- Compare parameter estimates across models
- Check credible interval overlap
- Identify key differences

### 6. Prepare Manuscript

#### Method A: Copy APA-Formatted Results

Open `APA_FORMATTED_RESULTS.txt`:

```
PVL_DELTA Model
--------------------------------------------------

mu_A: M = 0.046, SD = 0.023, 95% CI [0.008, 0.096]
mu_alpha: M = 0.231, SD = 0.146, 95% CI [0.012, 0.533]
mu_cons: M = 2.422, SD = 0.524, 95% CI [1.396, 3.431]
mu_lambda: M = 0.328, SD = 0.230, 95% CI [0.011, 0.828]
```

Copy directly into manuscript Results section.

#### Method B: Import CSV Tables

Import `tables/table1_parameter_estimates.csv` into Word/LaTeX:

**Word:**
1. Insert → Table → Insert Table from File
2. Select CSV file
3. Format as needed

**LaTeX:**
```latex
\begin{table}
\input{tables/table1_parameter_estimates.tex}
\caption{Parameter estimates with 95\% credible intervals}
\end{table}
```

#### Method C: Insert Figures

```latex
\begin{figure}
\includegraphics[width=\textwidth]{figures/model_comparison_forest.pdf}
\caption{Parameter estimates across models}
\label{fig:forest}
\end{figure}
```

## Advanced: Custom Visualizations

### Generate Specific Figure

```R
source("analysis/utils/visualization.R")

# Load fitted model
fit <- readRDS("analysis/outputs/pvl_delta_fit.rds")

# Custom trace plot for specific parameters
create_trace_plots(
  fit$samples,
  params = c("mu_A", "mu_alpha"),
  output_file = "my_custom_traces.pdf"
)

# Custom density plot
create_density_plots(
  fit$samples,
  params = c("mu_A", "mu_cons"),
  output_file = "my_custom_densities.pdf"
)
```

### Compare Specific Parameters

```R
source("analysis/utils/visualization.R")

# Load all models
pvl <- readRDS("analysis/outputs/pvl_delta_fit.rds")
vse <- readRDS("analysis/outputs/vse_fit.rds")
orl <- readRDS("analysis/outputs/orl_fit.rds")

fits <- list(pvl_delta = pvl, vse = vse, orl = orl)

# Compare learning rate across models
compare_parameter_across_models(
  fits,
  param = "mu_A",
  output_file = "learning_rate_comparison.pdf"
)
```

### Create Custom Tables

```R
source("analysis/utils/reporting.R")

# Generate specific table
fit <- readRDS("analysis/outputs/pvl_delta_fit.rds")
param_table <- generate_parameter_table(fit, "PVL-Delta")

# Export for manuscript
write.csv(param_table, "my_param_table.csv", row.names = FALSE)

# Generate LaTeX version
latex_code <- generate_latex_table(
  param_table,
  caption = "PVL-Delta parameter estimates",
  label = "tab:pvl_params"
)
writeLines(latex_code, "my_table.tex")
```

## Publication Checklist

Before submitting:

- [ ] All models converged (R-hat < 1.1)
- [ ] Trace plots reviewed and acceptable
- [ ] Posterior distributions reasonable
- [ ] Parameter estimates match expectations
- [ ] Figures generated and formatted
- [ ] Tables exported and verified
- [ ] Results formatted for manuscript
- [ ] Model comparison conducted
- [ ] Convergence diagnostics reported
- [ ] Supplementary materials prepared

## Manuscript Sections

### Methods

Report from `table2_model_comparison.csv`:

```
"We fit three hierarchical Bayesian models to the IGT data:
PVL-Delta (4 parameters), VSE (8 parameters), and ORL (5 parameters).
All models were implemented in JAGS with 4 chains, 2000 burn-in iterations,
and 5000 sampling iterations."
```

### Results - Convergence

Report from `table3_convergence.csv`:

```
"All models showed excellent convergence (R-hat < 1.1 for all parameters).
Effective sample sizes ranged from XXX to XXX (M = XXX), indicating
adequate sampling from the posterior distributions."
```

### Results - Parameter Estimates

Use `APA_FORMATTED_RESULTS.txt`:

```
"For the PVL-Delta model, the group-level learning rate was
M = 0.046 (SD = 0.023, 95% CI [0.008, 0.096]), indicating..."
```

### Figures

- **Figure 1**: Posterior distributions (`pvl_delta_posteriors.pdf`)
- **Figure 2**: Model comparison (`model_comparison_forest.pdf`)
- **Figure 3**: Convergence diagnostics (`pvl_delta_traces.pdf` - supplementary)

### Tables

- **Table 1**: Parameter estimates (`table1_parameter_estimates.csv`)
- **Table 2**: Model comparison (`table2_model_comparison.csv`)
- **Table S1**: Convergence diagnostics (`table3_convergence.csv` - supplementary)

## Tips for Publication

1. **Report Exact Values**: Use values from CSV tables, not figures
2. **Show Uncertainty**: Always report 95% CIs
3. **Document Priors**: Note that models use weakly informative truncated priors
4. **Convergence**: Include trace plots in supplementary materials
5. **Reproducibility**: Mention code/data availability

## Example Results Section

```
Results

Model Convergence

All three models converged successfully (R̂ < 1.1 for all parameters;
see Table S1). Effective sample sizes exceeded 400 for all parameters,
with a mean ESS of 2300 across all models.

Parameter Estimates

PVL-Delta Model. The group-level learning rate showed moderate updating
(μ_A = 0.046, 95% CI [0.008, 0.096]), with outcome sensitivity
μ_α = 0.231 [0.012, 0.533], choice consistency μ_c = 2.422 [1.396, 3.431],
and loss aversion μ_λ = 0.328 [0.011, 0.828]. See Figure 1 for full
posterior distributions.

[Continue for other models...]

Model Comparison

Forest plot comparison (Figure 2) revealed substantial differences in
learning rates across models [describe key findings].
```

## Need More?

- **Custom analysis**: Modify `visualization.R` or `reporting.R`
- **Additional metrics**: Add to `model_comparison.R`
- **Different priors**: Edit model files and re-run
- **Parameter recovery**: Run `parameter_recovery.R` and visualize

## Troubleshooting

**Issue**: Figures look wrong

- Check that fitted models loaded correctly
- Verify parameter names match expectations
- Try regenerating with `generate_paper_outputs.R`

**Issue**: Tables have missing values

- Check convergence (some parameters may not have converged)
- Verify all models completed successfully
- Check for NA values in fitted models

**Issue**: Want different figure format

- Modify output_file extension (.png, .pdf, .svg)
- Adjust width/height in visualization functions
- Export from R manually with custom settings

---

**Ready for your paper?** Run `generate_paper_outputs.R` after fitting completes!
