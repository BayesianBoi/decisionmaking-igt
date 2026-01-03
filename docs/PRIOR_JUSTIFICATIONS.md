# Prior Specifications for EEF Clinical Model

**Document Purpose:** Justify all hyperprior choices for hierarchical Bayesian modeling of Iowa Gambling Task data from substance use populations.

**Model:** Exploitation-Exploration with Forgetting (EEF), adapted from Yang et al. (2025)
**Application:** Clinical populations (Ahn 2014, Fridberg 2010)
**Sample:** N=146 (48 HC, 38 Amphetamine, 43 Heroin, 17 Cannabis)

---

## Executive Summary

**CRITICAL ISSUE:** Yang et al. (2025) used Variational Bayes Approximation (VBA toolbox), NOT MCMC. Their paper does NOT explicitly state priors.

**Our approach:** Design weakly informative priors based on:
1. Yang's empirical posterior means (reverse-engineering reasonable priors)
2. IGT computational modeling literature (Ahn et al., 2008, 2014)
3. Clinical neuroscience evidence (memory deficits in addiction)
4. Parameter identifiability constraints

---

## 1. θ (Outcome Sensitivity) - **NEW IN YANG'S MODEL**

### Parameter Definition
- **Range:** [0, 1]
- **Role:** Exponent in utility function: `V(t) = Gain^θ - Loss^θ`
- **Interpretation:**
  - θ = 0.5: Risk-seeking (concave utility)
  - θ = 1.0: Risk-neutral (linear utility)
  - θ > 1.0: Not used in Yang's model (bounded at 1)

### Difference from PVL-Delta
PVL-Delta uses:
- `alpha` (outcome sensitivity) ∈ [0, 2]
- `lambda` (loss aversion) ∈ [0, 10]
- Utility: `u(x) = x^alpha` if x ≥ 0, else `-lambda * |x|^alpha`

Yang's EEF uses:
- **Symmetric** utility (no separate loss aversion)
- θ ∈ [0, 1] only

**Rationale:** Simplifies model (4 params instead of 5), focuses on forgetting.

### Prior Specification

```jags
mu_theta ~ dbeta(1.5, 3)
```

**Justification:**

1. **Empirical evidence (Yang Fig. 9, page 13):**
   - Young adults: θ ≈ 0.35
   - Older adults: θ ≈ 0.30
   - Both groups show θ < 0.5 (risk-seeking in IGT context)

2. **Beta(1.5, 3) properties:**
   - Mean = 1.5/(1.5+3) = 0.33
   - Mode = (1.5-1)/(1.5+3-2) = 0.2
   - SD = 0.19
   - 95% credible interval: [0.07, 0.70]

3. **Why NOT uniform?**
   - Uniform Beta(1,1) is too diffuse
   - IGT literature shows θ rarely > 0.7
   - Weakly informative prior centers on empirically plausible range

4. **Clinical prediction:**
   - Substance users may show **lower θ** (blunted reward sensitivity)
   - Eppinger et al. (2012): Reduced striatal reward processing in addiction

### Group-Level Variability

```jags
sigma_theta ~ dunif(0, 0.3)
```

**Justification:**
- Constrains individual differences to prevent over-dispersion
- Uniform(0, 0.3) allows SD up to 30% of parameter range
- Small clinical samples (N=38-48 per group) need regularization
- If sigma_theta → 0.3, individual θ values span [0.03, 0.63] (95% CI)

---

## 2. λ (Forgetting Rate) - **KEY PARAMETER FOR RESEARCH QUESTIONS**

### Parameter Definition
- **Range:** [0, 1]
- **Role:** Memory decay per trial in both exploitation and exploration modules
- **Interpretation:**
  - λ = 0: Perfect memory (EEF reduces to PVL-Delta-like)
  - λ = 0.2: Lose 20% of learned value per trial (slow forgetting)
  - λ = 0.5: Lose 50% per trial (fast forgetting)
  - λ = 1.0: Complete amnesia (only current trial matters)

### Theoretical Background

**Memory decay in decision-making:**
- Ebbinghaus (1885): Exponential forgetting curves
- ACT-R (Anderson & Schooler, 1991): Power law of forgetting
- IGT context: 100 trials over ~10 minutes, so λ should be modest (0.1-0.5)

**Clinical prediction (PRIMARY RESEARCH HYPOTHESIS):**
- **Substance users should show higher λ** due to:
  - Working memory deficits (Ersche et al., 2006)
  - Impaired memory consolidation (Verdejo-García et al., 2008)
  - Hippocampal damage in opioid users (Tramullas et al., 2008)

### Prior Specification

```jags
mu_lambda_forget ~ dbeta(2, 3)
```

**Justification:**

1. **Empirical evidence (Yang Fig. 9, page 13):**
   - Young adults: λ ≈ 0.46
   - Older adults: λ ≈ 0.54
   - **Higher in older adults** (consistent with age-related memory decline)

2. **Beta(2, 3) properties:**
   - Mean = 2/(2+3) = 0.4
   - Mode = (2-1)/(2+3-2) = 0.33
   - SD = 0.19
   - 95% credible interval: [0.09, 0.75]

3. **Why these hyperparameters?**
   - Centers prior on Yang's empirical range (0.4-0.5)
   - Mildly skeptical of extreme forgetting (λ > 0.7)
   - Still allows data to dominate if λ is very low/high
   - Beta(2,3) is "weakly informative" - regularizes but doesn't constrain strongly

4. **Alternative priors considered:**
   - Beta(1, 1): Uniform, too diffuse, no information from Yang
   - Beta(5, 5): Too concentrated around 0.5, over-regularizes
   - Beta(2, 3): Balances prior information with flexibility

### Group-Level Variability

```jags
sigma_lambda_forget ~ dunif(0, 0.3)
```

**Justification:**
- Same rationale as sigma_theta
- Clinical samples need moderate regularization
- Allows substantial individual differences (SD up to 30% of range)

### Expected Effect Size

Based on Yang's age comparison (0.46 vs 0.54):
- Effect size d = (0.54 - 0.46) / 0.19 ≈ 0.42 (small-to-medium)

For substance users, we expect **medium-to-large effects** (d = 0.5-0.8):
- λ_HC ≈ 0.40
- λ_substance ≈ 0.55-0.60 (predicted)

---

## 3. φ (Exploration Incentive)

### Parameter Definition
- **Range:** [-5, 5] (can be NEGATIVE!)
- **Role:** Attraction/repulsion to unchosen decks in exploration module
- **Interpretation:**
  - φ > 0: Positive exploration incentive (prefer unexplored options)
  - φ = 0: Neutral (no exploration bias)
  - φ < 0: Exploration aversion (avoid unexplored options, stick to known)

### Update Rule (Yang Eq. 5)

For UNCHOSEN decks:
```
Exploration_d(t+1) = λ × Exploration_d(t) + (1-λ) × φ
```

- If φ > 0 and deck unchosen for many trials → exploration weight increases
- If φ < 0 → exploration weight becomes negative (repels choice)

### Prior Specification

```jags
mu_phi ~ dnorm(0, 2)T(-5, 5)
```

**Justification:**

1. **Empirical evidence (Yang Fig. 9):**
   - Young adults: φ ≈ 0.69
   - Older adults: φ ≈ 0.50
   - Both positive, but moderate (not extreme)

2. **Normal(0, 2) properties (before truncation):**
   - Mean = 0
   - Precision = 2 → SD = 1/sqrt(2) = 0.71
   - 95% interval: [-1.4, 1.4]
   - Truncation at [-5, 5] has minimal effect (tails are tiny)

3. **Why centered at 0?**
   - **Theoretically neutral:** No a priori assumption about exploration direction
   - **Clinical prediction:** Substance users may show **negative φ**:
     - Habitual behavior in addiction (Everitt & Robbins, 2016)
     - Reduced novelty-seeking after chronic use (Belin et al., 2008)
     - Impulsivity → stick with known options, avoid uncertain exploration

4. **Why allow negative values?**
   - PVL-Delta and VSE models don't have this parameter
   - Exploration aversion (φ < 0) is psychologically meaningful
   - Substance users might AVOID exploration due to cognitive load

### Group-Level Variability

```jags
sigma_phi ~ dunif(0, 2)
```

**Justification:**
- Wider variability than other parameters (SD up to 2)
- φ has large range [-5, 5], so SD=2 is ~20% of range
- Exploration strategies show high individual differences

---

## 4. Consistency (Inverse Temperature)

### Parameter Definition
- **Range:** [0, 5]
- **Role:** Softmax determinism parameter (β in most models)
- **Interpretation:**
  - cons = 0: Completely random choices
  - cons = 1: Moderate sensitivity to value differences
  - cons = 5: Very deterministic (always choose highest-value deck)

### Softmax Choice Rule

```
P(choice = d) = exp(cons × weight_d) / Σ exp(cons × weight_i)
```

### Prior Specification

```jags
mu_cons ~ dnorm(1, 2)T(0, 5)
```

**Justification:**

1. **Standard in IGT literature:**
   - Ahn et al. (2008): cons typically 0.5-2.0
   - Ahn et al. (2014): Used Normal(1, 2)T(0,5) as hyperprior
   - **We use EXACTLY the same prior as Ahn (2014)**

2. **Normal(1, 2) properties:**
   - Mean = 1
   - Precision = 2 → SD = 0.71
   - 95% interval: [0, 2.4] after truncation at [0,5]

3. **Clinical prediction:**
   - **No group difference expected**
   - Consistency reflects choice noise, not learning/memory
   - Substance users may have slightly lower cons (more noisy), but effect is typically small

4. **Why NOT Yang's transformation?**
   - Yang uses: C = 3β - 1, where β ∈ [0, 5]
   - This gives C ∈ [-1, 14] (unusual range)
   - **Reason:** VBA toolbox compatibility, NOT cognitive
   - **We use standard β directly** for interpretability

### Group-Level Variability

```jags
sigma_cons ~ dunif(0, 1.5)
```

**Justification:**
- Same as Ahn et al. (2014)
- Allows moderate individual differences
- SD=1.5 means 95% of individuals span [0, 3.9] if mu_cons=1

---

## 5. First-Choice Priors (w_ini) - **DATA-DRIVEN**

### Background

Yang et al. (2025, Section 2.4, page 5) initialize deck weights using **empirical first-choice frequencies** from Steingroever et al. (2015) data (N=504):

- Deck A: 38.29% (most popular) → w_ini[A] = 0.0184
- Deck B: 25.00% → w_ini[B] = 0.0020
- Deck C: 20.83% → w_ini[C] = -0.0050
- Deck D: 15.87% (least popular) → w_ini[D] = -0.0155

**Calculation:** Reverse softmax transformation
```r
freq <- c(0.3829, 0.25, 0.2083, 0.1587)
w_ini <- log(freq) / 3  # beta = 3 (Yang's choice)
```

### Our Approach: Calculate from YOUR Data

**CRITICAL DECISION:** Do NOT use Yang's w_ini values. Calculate from YOUR 173 subjects.

**Rationale:**
2. **Data-driven approach:** w_ini should reflect YOUR sample, not Steingroever's healthy sample
3. **Validation:** Check if groups show different first-choice patterns (chi-square test)

### Implementation

```r
# In prepare_eef_data.R:
w_ini <- calculate_first_choice_priors(dat, beta = 3)
```

**If groups differ significantly** in first-choice patterns:
```r
w_ini_by_group <- calculate_first_choice_priors_by_group(dat, group_var = "group")
# Then use group-specific w_ini in JAGS model
```

### Why β = 3?

Yang doesn't justify this choice. Likely options:
- β = 3 is arbitrary but reasonable
- Smaller β (e.g., 1) would amplify deck differences
- Larger β (e.g., 5) would compress differences

**We keep β = 3** to match Yang, but this is a modeling choice, not a psychological parameter.

---

## 6. Subject-Level Parameter Distributions

### Beta Distribution Parameterization

For parameters bounded [0, 1] (θ, λ):

```jags
theta[s] ~ dbeta(
  mu_theta * (1/sigma_theta^2 - 1),
  (1 - mu_theta) * (1/sigma_theta^2 - 1)
)T(0.01, 0.99)
```

**Method of moments:** Convert (mu, sigma) → (alpha, beta) for Beta distribution.

**Truncation T(0.01, 0.99):**
- Prevents extreme values (0 or 1) that cause numerical issues
- Still allows >99% of probability mass

### Normal Distribution Parameterization

For unbounded parameters (φ, cons):

```jags
phi[s] ~ dnorm(mu_phi, 1/sigma_phi^2)T(-5, 5)
cons[s] ~ dnorm(mu_cons, 1/sigma_cons^2)T(0.01, 5)
```

**JAGS precision parameterization:** `precision = 1/variance = 1/sigma^2`

**Truncation:**
- Enforces parameter bounds
- Prevents MCMC from sampling impossible values

---

## 7. MCMC Settings - **CRITICAL FOR HIERARCHICAL MODELS**

### Configuration

```r
n_adapt = 5000      # Adaptation phase
n_burnin = 10000    # Burn-in (discarded)
n_iter = 20000      # Sampling iterations per chain
n_chains = 4        # Number of chains
thin = 2            # Store every 2nd sample
```

### Justifications

**1. Why n_adapt = 5000 (not default 1000)?**
- Hierarchical models need more adaptation
- 4 group-level means + 4 group-level SDs + 173×4 individual parameters = 700+ parameters
- JAGS needs to tune proposal distributions for efficient sampling

**2. Why n_burnin = 10000 (not 2000)?**
- **Small N problem:** N=38-48 per group → less data → harder inference
- Need longer burn-in to ensure convergence from arbitrary starting points
- Hierarchical structure: Individual parameters depend on group parameters → slower mixing

**3. Why n_iter = 20000 (not 5000)?**
- **Effective sample size:** Want N_eff > 1000 per parameter for stable inference
- Autocorrelation in MCMC → collect more samples
- 4 chains × 20K iter / 2 thin = 40K samples → likely N_eff = 5-15K

**4. Why n_chains = 4 (not 1)?**
- Gelman-Rubin R-hat diagnostic requires multiple chains
- R-hat < 1.1 = good convergence
- 4 chains is standard in Bayesian literature (Gelman et al., 2013)

**5. Why thin = 2 (not 1)?**
- Reduces autocorrelation in stored samples
- Saves memory (40K samples instead of 80K)
- Minimal information loss (samples are highly correlated anyway)

### Expected Runtime

**Estimate:** 4-8 hours on M1 Mac

**Calculation:**
- 173 subjects × 100 trials = 17,300 choice observations
- 700+ parameters
- (5K adapt + 10K burnin + 20K sample) × 4 chains = 140K total iterations
- Each iteration evaluates likelihood for all 17,300 observations
- Total: ~2.4 billion likelihood evaluations

**Variance:**
- Faster if JAGS parallelizes chains (likely)
- Slower if sigma parameters explore extreme values (triggers resampling)

---

## 8. Summary Table of All Priors

| Parameter | Range | Distribution | Hyperparameters | Mean | SD | Justification |
|-----------|-------|--------------|-----------------|------|----|--------------|
| **mu_theta** | [0,1] | Beta | (1.5, 3) | 0.33 | 0.19 | Yang's empirical: ~0.3-0.35 |
| **mu_lambda_forget** | [0,1] | Beta | (2, 3) | 0.40 | 0.19 | Yang's empirical: 0.46-0.54 |
| **mu_phi** | [-5,5] | Normal (trunc) | (0, 2) | 0 | 0.71 | Neutral, allows negative |
| **mu_cons** | [0,5] | Normal (trunc) | (1, 2) | 1 | 0.71 | Ahn et al. (2014) standard |
| **sigma_theta** | [0,∞) | Uniform | (0, 0.3) | 0.15 | 0.087 | Regularization for small N |
| **sigma_lambda_forget** | [0,∞) | Uniform | (0, 0.3) | 0.15 | 0.087 | Regularization for small N |
| **sigma_phi** | [0,∞) | Uniform | (0, 2) | 1.0 | 0.58 | Allow wide individual diffs |
| **sigma_cons** | [0,∞) | Uniform | (0, 1.5) | 0.75 | 0.43 | Ahn et al. (2014) standard |

---

## 9. Sensitivity Analysis (Future Work)

**Question:** How sensitive are results to prior choices?

**Plan:**
1. Refit model with alternative priors:
   - Beta(1,1) for theta, lambda (uniform)
   - Beta(5,5) for theta, lambda (concentrated at 0.5)
   - Normal(1,1) for phi (wider)

2. Compare posterior estimates:
   - If posteriors change < 10%, priors are "weakly informative" ✓
   - If posteriors change > 20%, priors are "too strong" ✗

3. Report in supplementary materials:
   - "Results were robust to alternative prior specifications"
   - OR "Prior choice affected X parameter, future work needs more data"

---

## 10. Comparison to Yang et al. (2025)

| Aspect | Yang et al. (2025) | Our Implementation |
|--------|--------------------|--------------------|
| **Inference method** | Variational Bayes (VBA toolbox) | MCMC (JAGS) |
| **Priors stated?** | No (not in paper) | Yes (this document) |
| **Sample size** | N=504 healthy | N=146 (48 HC, 98 substance users) |
| **w_ini** | Calculated from their N=504 | Calculated from sample data |
| **Consistency param** | C = 3β - 1 (VBA-specific) | β directly (standard) |
| **Research question** | Does forgetting matter? (methods) | Memory deficits in addiction (clinical) |

This implementation applies Yang's validated model to investigate memory deficits in addiction populations.

---

## References

- Ahn et al. (2008). Comparison of decision learning models using the generalization criterion method. *Cognitive Science, 32*(8), 1376-1402.
- Ahn et al. (2014). Decision-making in stimulant and opiate addicts in protracted abstinence. *Neuropsychopharmacology, 39*(5), 1244-1253.
- Ersche et al. (2006). Profile of executive and memory function associated with amphetamine and opiate dependence. *Neuropsychopharmacology, 31*, 1036-1047.
- Gelman et al. (2013). *Bayesian Data Analysis* (3rd ed.). CRC Press.
- Verdejo-García et al. (2008). Neurocognitive processes in addiction. In *Addiction Neuroethics* (pp. 75-100).
- Yang et al. (2025). Forgetting phenomena in the Iowa Gambling Task. *Frontiers in Psychology, 16*, 1510151.

---

**Document Version:** 1.0  
**Date:** January 2026  
**Status:** Ready for supplementary materials

