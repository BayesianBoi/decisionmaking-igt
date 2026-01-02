# Troubleshooting Guide

## Issue: JAGS Model Initialization Errors

### Problem
The original JAGS models (`pvl_delta.jags`, `vse.jags`, `orl.jags`) encountered initialization errors:
1. "Attempt to redefine node" - Variable `sens[s,t]` was being assigned inside a deck loop
2. "Invalid parent values" - Extreme parameter values during initialization caused numerical issues

### Root Causes

1. **Redefinition Error**: The sensitivity parameter was computed inside the deck loop, causing it to be redefined 4 times per trial.

2. **Numerical Instability**:
   - Unbounded priors (e.g., `alpha ~ exp(alpha_raw)` where `alpha_raw ~ dnorm(...)`) can produce extreme values
   - `pow(negative_number, fractional_alpha)` is undefined in real numbers
   - Softmax with extreme values causes overflow/underflow
   - Initial EV values with high sensitivity can create invalid probabilities

### Solution

Created `pvl_delta_v2.jags` with:

1. **Truncated Priors**: Direct truncated distributions instead of transformations
   ```jags
   A[s] ~ dbeta(...)T(0.01, 0.99)      # bounded [0.01, 0.99]
   alpha[s] ~ dnorm(...)T(0.01, 2)      # bounded [0.01, 2]
   cons[s] ~ dnorm(...)T(0.01, 5)       # bounded [0.01, 5]
   lambda[s] ~ dnorm(...)T(0.01, 10)    # bounded [0.01, 10]
   ```

2. **Numerical Safety**:
   ```jags
   abs_outcome[s, t] <- abs(outcome[s, t]) + 0.001  # avoid pow(0, alpha)
   ```

3. **Cleaner Structure**: Used `equals()` function instead of nested `ifelse()` for updates

### Test Results

The fixed model (`pvl_delta_v2.jags`) successfully fits data:
- 48 subjects, 100 trials each
- Converges in ~2 minutes (500 adapt, 500 burn-in, 1000 samples)
- R-hat values ~1.0-1.05 (good convergence)
- Reasonable parameter estimates:
  - mu_A: 0.046 (learning rate)
  - mu_alpha: 0.231 (outcome sensitivity)
  - mu_cons: 2.42 (choice consistency)
  - mu_lambda: 0.328 (loss aversion)

## How to Proceed

### Option 1: Use the Working v2 Model (Recommended)

The `pvl_delta_v2.jags` model works and can be used immediately:

```bash
# Test it
Rscript analysis/quick_test.R

# Use it in the main pipeline (modify fit_models.R to use pvl_delta_v2.jags)
```

### Option 2: Fix Original Models

Apply the same principles to fix `pvl_delta.jags`, `vse.jags`, and `orl.jags`:

1. Replace probit/log transformations with direct truncated priors
2. Add small constants to prevent `pow(0, alpha)`
3. Consider clamping extreme values in softmax

### Option 3: Simplify Models Further

For even better stability:
- Use simpler priors (beta/gamma instead of normal+transform)
- Pre-center and scale data
- Use fixed reasonable bounds on all parameters

## Known Limitations

1. **Unused Variables Warning**: The warnings about unused variables (`T`, `reward`, `loss`) are harmless - the data preparation creates these but the models only use `outcome` (which combines reward and loss).

2. **Computational Cost**: With all 173 subjects:
   - Expect 30-120 minutes per model
   - Consider using `fit_all_studies = FALSE` for testing

3. **VSE and ORL Models**: These still need the same fixes applied. Use the v2 approach as a template.

## Quick Fix Checklist

To fix any model encountering similar issues:

- [ ] Move variable assignments outside loops where they're defined
- [ ] Add truncation bounds to all priors: `T(lower, upper)`
- [ ] Add small constants before `pow()` operations
- [ ] Test with small dataset first (`quick_test.R`)
- [ ] Check R-hat < 1.1 for all parameters
- [ ] Verify parameter estimates are reasonable

## Next Steps

1. ✅ PVL-Delta v2 model works
2. ⏳ Apply same fixes to VSE model
3. ⏳ Apply same fixes to ORL model
4. ⏳ Update `fit_models.R` to use v2 models
5. ⏳ Run full pipeline on complete dataset

## References

- JAGS manual: https://mcmc-jags.sourceforge.io/
- Truncated distributions in JAGS: Use `T(lower, upper)` or `T(lower,)` or `T(,upper)`
- Common JAGS errors: https://sourceforge.net/p/mcmc-jags/discussion/610037/
