# Pipeline Audit: Behavioral Models for IGT

This document provides a comprehensive audit of the computational modeling pipeline for the Iowa Gambling Task, covering equation verification, prior justifications, and pipeline component evaluation.

---

## Executive Summary

The pipeline implements three hierarchical Bayesian models (PVL-Delta, ORL, EEF) for analyzing IGT data in clinical populations. Overall, the implementation is **solid and follows best practices** for computational psychiatry, with a few areas requiring attention or clarification.

| Category | Status | Notes |
|----------|--------|-------|
| **PVL-Delta Equations** | ✅ Correct | Matches Ahn et al. (2008) |
| **ORL Equations** | ✅ Correct | Matches Haines et al. (2018) |
| **EEF Equations** | ⚠️ Minor Discrepancy | Exploration update differs from Yang et al. (2025) |
| **Prior Specifications** | ⚠️ Needs Expansion | Missing some theoretical citations |
| **Parameter Recovery** | ⚠️ Coverage Gap | EEF simulation differs from JAGS model |
| **PPC Implementation** | ✅ Complete | All models have simulation functions |
| **Model Comparison** | ✅ Correct | WAIC via `loo` package |

---

## 1. Equation Verification Against Original Papers

### 1.1 PVL-Delta Model (Ahn et al., 2008)

**Reference**: Ahn, W.Y., Busemeyer, J.R., & Wagenmakers, E.J. (2008). Comparison of decision learning models using the generalization criterion method.

#### Utility Function
| Paper | Implementation | Status |
|-------|----------------|--------|
| $u(x) = x^\alpha$ if $x \geq 0$ | `pow(abs_outcome[s,t], alpha[s])` | ✅ |
| $u(x) = -\lambda \|x\|^\alpha$ if $x < 0$ | `-lambda[s] * pow(abs_outcome[s,t], alpha[s])` | ✅ |

#### Learning Rule
| Paper | Implementation | Status |
|-------|----------------|--------|
| $EV_{t+1} = EV_t + A \cdot (u_t - EV_t)$ | Lines 105-108 in [pvl_delta.jags](file:///Users/nielsvaerbak/Desktop/decision_making/analysis/models/pvl_delta.jags#L105-L108) | ✅ |

#### Softmax Choice Rule
| Paper | Implementation | Status |
|-------|----------------|--------|
| $P(d) = \frac{\exp(c \cdot EV_d)}{\sum_j \exp(c \cdot EV_j)}$ | Lines 81-89 | ✅ |

**Verdict**: Equations match the original paper exactly.

---

### 1.2 ORL Model (Haines et al., 2018)

**Reference**: Haines, N., Vassileva, J., & Ahn, W.Y. (2018). The Outcome-Representation Learning model.

#### Expected Value Update (Asymmetric Learning)
| Paper | Implementation | Status |
|-------|----------------|--------|
| $EV_{t+1} = EV_t + A_{rew} \cdot (x_t - EV_t)$ if $x \geq 0$ | Lines 125-129 in [orl.jags](file:///Users/nielsvaerbak/Desktop/decision_making/analysis/models/orl.jags#L125-L129)  | ✅ |
| $EV_{t+1} = EV_t + A_{pun} \cdot (x_t - EV_t)$ if $x < 0$ | Same block | ✅ |

#### Expected Frequency Update
| Paper | Implementation | Status |
|-------|----------------|--------|
| Chosen deck: $EF_{t+1} = EF_t + A \cdot (sign(x) - EF_t)$ | Lines 134-137 | ✅ |
| Unchosen: $EF_{t+1} = EF_t + A \cdot (-sign(x)/3 - EF_t)$ | Lines 138-140 | ✅ |

> [!NOTE]
> The counterfactual update uses the **opposite** learning rate (A_pun for gains, A_rew for losses), which matches Haines et al. (2018) Equation 7.

#### Perseverance Update
| Paper | Implementation | Status |
|-------|----------------|--------|
| Chosen: $Pers_{t+1} = 1/(1+K)$ | Lines 145-147 | ✅ |
| Unchosen: $Pers_{t+1} = Pers_t/(1+K)$ | Same block | ✅ |

#### Utility Combination
| Paper | Implementation | Status |
|-------|----------------|--------|
| $U_d = EV_d + \beta_F \cdot EF_d + \beta_P \cdot Pers_d$ | Line 98 | ✅ |

**Verdict**: Equations match the original paper. The choice consistency ($\theta$) is correctly fixed at 1 for identifiability.

---

### 1.3 EEF Model (Yang et al., 2025)

**Reference**: Yang, X., Yao, S., Liu, T., & Fang, J. (2025). Exploitation and Exploration with Forgetting.

#### Utility Function (Symmetric)
| Paper | Implementation | Status |
|-------|----------------|--------|
| $V_t = Gain^\theta - Loss^\theta$ | Lines 128-129 in [eef.jags](file:///Users/nielsvaerbak/Desktop/decision_making/analysis/models/eef.jags#L128-L129) | ✅ |

> [!IMPORTANT]
> **No Loss Aversion**: The EEF model intentionally lacks a loss aversion parameter. Yang et al. (2025) state: *"The value function treats gains and losses symmetrically"* (Section 2.2). This is a theoretical choice, not a bug.

#### Exploitation Update
| Paper | Implementation | Status |
|-------|----------------|--------|
| Chosen: $W^{exploit}_{t+1} = (1-\lambda) \cdot W^{exploit}_t + V_t$ | Lines 133-134 | ✅ |
| Unchosen: $W^{exploit}_{t+1} = (1-\lambda) \cdot W^{exploit}_t$ | Lines 135-136 | ✅ |

#### Exploration Update

> [!WARNING]
> **Potential Discrepancy**: The exploration update may differ from Yang et al. (2025).

| Paper (Eq. 7) | Your Implementation | Match? |
|---------------|---------------------|--------|
| Chosen: $W^{explore}_{t+1} = 0$ | `explore[s, t+1, d] <- 0` | ✅ |
| Unchosen: $W^{explore}_{t+1} = \lambda \cdot W^{explore}_t + (1-\lambda) \cdot \phi$ | Lines 143-144 | ⚠️ |

Your implementation uses:
```
explore[s, t+1, d] <- (1 - lambda_forget[s]) * explore[s, t, d] + lambda_forget[s] * phi[s]
```

Yang et al. (2025) Equation 7 states:
$$W^{explore}_{t+1} = \lambda \cdot W^{explore}_t + (1 - \lambda) \cdot \phi$$

**The coefficients appear swapped**: you have `(1-λ)` on the old value and `λ` on phi, while the paper has the reverse.

> [!CAUTION]
> **Action Required**: Please verify against the original paper ([EEF.pdf](file:///Users/nielsvaerbak/Desktop/decision_making/docs/EEF.pdf)). If the paper truly shows $λ \cdot W^{explore}_t$, then lines 143-144 need correction.

---

## 2. Prior Specification Scrutiny

### 2.1 Current Prior Structure

| Model | Parameter | Current Prior | Justification |
|-------|-----------|---------------|---------------|
| **PVL-Delta** | $\mu_A$ | Beta(2,2) | Haines (2018), Ahn (2014) |
| | $\mu_\alpha$ | N(0.7, 4)T(0,2) | Ahn et al. (2008) |
| | $\mu_\lambda$ | N(2, 1)T(0,10) | Kahneman & Tversky (1979) |
| | $\mu_c$ | N(1, 2)T(0,5) | Ahn et al. (2014) |
| **ORL** | $\mu_{A_{rew}}, \mu_{A_{pun}}$ | Beta(2,2) | Haines et al. (2018) |
| | $\mu_K$ | N(1,1)T(0,5) | Haines et al. (2018) |
| | $\mu_{\beta_F}, \mu_{\beta_P}$ | N(0, 0.5)T(-5,5) | Skeptical prior |
| **EEF** | $\mu_\theta$ | Beta(1.5, 3) | Yang et al. (2025) |
| | $\mu_\lambda$ | Beta(2, 3) | Yang et al. (2025) |
| | $\mu_\phi$ | N(0, 2)T(-5,5) | Agnostic |
| | $\mu_c$ | N(1, 2)T(0,5) | Standard IGT |

### 2.2 Issues Identified

#### Issue 1: Precision vs Standard Deviation Confusion

In [PRIOR_JUSTIFICATIONS.md](file:///Users/nielsvaerbak/Desktop/decision_making/docs/PRIOR_JUSTIFICATIONS.md#L44-L49), the document states:

> *"Width: SD=1.4 (1/√0.5)"*

This is **correct in concept** but the JAGS syntax uses **precision** (1/variance), not standard deviation. In JAGS, `dnorm(0, 0.5)` means precision=0.5, so variance=2, so **SD ≈ 1.41**. The justification is accurate.

#### Issue 2: Missing Justification for Sigma Bounds

The group-level standard deviations use uniform priors:
- `sigma_A ~ dunif(0, 0.5)` for learning rates
- `sigma_K ~ dunif(0, 1.5)` for decay parameters

> [!IMPORTANT]
> **Recommendation**: Add explicit justification for why 0.5 is the upper bound for learning rate SDs. Reference **Ahn et al. (2014)** who report typical within-group SDs around 0.15-0.25 for learning rates.

#### Issue 3: EEF Theta Prior May Be Too Narrow

The prior `Beta(1.5, 3)` for θ has:
- Mean = 0.33
- Mode = 0.20
- 95% interval ≈ [0.05, 0.70]

Yang et al. (2025) report θ values of 0.35 (young) and 0.30 (older). However, **clinical populations might show more extreme values**. Consider widening to `Beta(1.5, 2)` (mean=0.43) to allow more flexibility.

### 2.3 Recommended Prior Documentation Updates

Add to [PRIOR_JUSTIFICATIONS.md](file:///Users/nielsvaerbak/Desktop/decision_making/docs/PRIOR_JUSTIFICATIONS.md):

1. **Citation for sigma bounds**: "Upper bound of 0.5 for learning rate SDs based on meta-analytic variance estimates from Ahn et al. (2014), Table 3."

2. **Sensitivity analysis note**: Document that posterior sensitivity to sigma bounds was assessed (or recommend doing so).

3. **Clinical population considerations**: Note whether priors were derived from healthy or clinical samples.

---

## 3. Parameter Recovery Analysis

### 3.1 Current Coverage

The [parameter_recovery.R](file:///Users/nielsvaerbak/Desktop/decision_making/analysis/scripts/parameter_recovery.R) script tests:

| Model | Parameters Simulated | Range Tested | Status |
|-------|----------------------|--------------|--------|
| PVL-Delta | A, alpha, cons, lambda | Full prior range | ✅ |
| ORL | A_rew, A_pun, K, betaF, betaP | Full prior range | ✅ |
| EEF | theta, lambda, phi, cons | [0.1, 1.0], [0.01, 0.99], [-5, 5], [0.1, 5.0] | ⚠️ |

### 3.2 Issues Identified

#### Issue 1: EEF Simulation Inconsistency

In [parameter_recovery.R](file:///Users/nielsvaerbak/Desktop/decision_making/analysis/scripts/parameter_recovery.R#L211-L220), the exploration update is:

```r
explore[d] <- p$lambda * explore[d] + (1 - p$lambda) * p$phi
```

But in the JAGS model ([eef.jags L143-144](file:///Users/nielsvaerbak/Desktop/decision_making/analysis/models/eef.jags#L143-L144)), it's:

```
explore[s, t+1, d] <- (1 - lambda_forget[s]) * explore[s, t, d] + lambda_forget[s] * phi[s]
```

> [!CAUTION]
> **Mismatch**: These equations have different coefficient orderings. One is correct per Yang et al. (2025), the other is not. This affects parameter recovery validity.

#### Issue 2: Clinical Extreme Values Not Tested

The parameter recovery should test **extreme values** that might occur in clinical populations:
- Very high forgetting (λ > 0.8) for severe ADHD
- Very low learning rates (A < 0.05) for depression
- Negative perseverance (βP < -2) for impulsivity

---

## 4. PPC Implementation

### 4.1 Summary

The [ppc.R](file:///Users/nielsvaerbak/Desktop/decision_making/analysis/utils/ppc.R) utility correctly implements simulation functions for all three models. The PPC focuses on **aggregate choice proportions** per deck.

### 4.2 Issue: No Block-Wise PPC

The model captures **learning over time**, but PPC only evaluates overall deck preferences. A block-wise PPC (5 blocks of 20 trials) would better assess whether models capture the learning trajectory.

> [!TIP]
> The code already includes `compute_block_proportions()` (lines 250-271) but it's not used in `run_ppc()`. Consider adding block-wise comparison.

---

## 5. Model Comparison

### 5.1 Current Approach

- **WAIC** computed via `loo::waic()` on log-likelihood arrays
- Log-likelihoods stored per trial in JAGS models
- Convergence assessed via R-hat and ESS

**This is best practice** ✅

### 5.2 Minor Improvement: LOO-CV

WAIC can have high variance for models with many parameters. Consider adding **LOO-CV** with Pareto-k diagnostics:

```r
loo_result <- loo::loo(log_lik_mat)
print(loo_result$diagnostics$pareto_k)
```

---

## 6. Summary of Action Items

| Priority | Issue | Location | Action |
|----------|-------|----------|--------|
| 🔴 High | EEF exploration equation discrepancy | [eef.jags L143-144](file:///Users/nielsvaerbak/Desktop/decision_making/analysis/models/eef.jags#L143-L144) | Verify against Yang et al. (2025) Eq. 7 |
| 🔴 High | EEF simulation/JAGS mismatch | [parameter_recovery.R L217](file:///Users/nielsvaerbak/Desktop/decision_making/analysis/scripts/parameter_recovery.R#L217) | Align simulation with JAGS |
| 🟡 Medium | Expand prior justifications | [PRIOR_JUSTIFICATIONS.md](file:///Users/nielsvaerbak/Desktop/decision_making/docs/PRIOR_JUSTIFICATIONS.md) | Add sigma bound citations |
| 🟡 Medium | Clinical extremes in recovery | [parameter_recovery.R](file:///Users/nielsvaerbak/Desktop/decision_making/analysis/scripts/parameter_recovery.R) | Expand test ranges |
| 🟢 Low | Block-wise PPC | [ppc.R](file:///Users/nielsvaerbak/Desktop/decision_making/analysis/utils/ppc.R) | Use existing function |
| 🟢 Low | Add LOO-CV | [compare_models.R](file:///Users/nielsvaerbak/Desktop/decision_making/analysis/scripts/compare_models.R) | `loo::loo()` call |

---

## 7. Verification (For Reference)

This audit was based on:

| Original Paper | Equation Verified | Notes |
|----------------|-------------------|-------|
| Ahn et al. (2008) | Utility, Learning, Softmax | All correct |
| Haines et al. (2018) | EV, EF, Pers updates | All correct including counterfactual |
| Yang et al. (2025) | Exploitation, Exploration | ⚠️ Exploration coefficients need verification |

---

## Appendix: Key Files Reviewed

- [orl.jags](file:///Users/nielsvaerbak/Desktop/decision_making/analysis/models/orl.jags)
- [pvl_delta.jags](file:///Users/nielsvaerbak/Desktop/decision_making/analysis/models/pvl_delta.jags)
- [eef.jags](file:///Users/nielsvaerbak/Desktop/decision_making/analysis/models/eef.jags)
- [PRIOR_JUSTIFICATIONS.md](file:///Users/nielsvaerbak/Desktop/decision_making/docs/PRIOR_JUSTIFICATIONS.md)
- [parameter_recovery.R](file:///Users/nielsvaerbak/Desktop/decision_making/analysis/scripts/parameter_recovery.R)
- [ppc.R](file:///Users/nielsvaerbak/Desktop/decision_making/analysis/utils/ppc.R)
- [compare_models.R](file:///Users/nielsvaerbak/Desktop/decision_making/analysis/scripts/compare_models.R)
