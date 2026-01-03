# Troubleshooting Guide

## Common Issues

### JAGS Model Initialization Errors

**Problem:** "Attempt to redefine node" or "Invalid parent values"

**Root Causes:**
1. Variable assignments inside loops being redefined
2. Numerical instability from extreme parameter values
3. `pow(0, alpha)` undefined operations

**Solution:**
The JAGS models use truncated priors and numerical safeguards:
```jags
# Truncated priors prevent extreme values
A[s] ~ dbeta(...)T(0.01, 0.99)
theta[s] ~ dbeta(...)T(0.01, 0.99)

# Small constants avoid pow(0, alpha) issues
util[s, t] <- pow(gain[s, t] + 0.001, theta[s])
```

### Model Not Converging

**Symptoms:** R-hat > 1.1, trace plots showing trends

**Solutions:**
1. Increase adaptation: `n_adapt = 10000`
2. Increase burn-in: `n_burnin = 20000`
3. Increase samples: `n_iter = 40000`
4. Check for data issues (extreme outcomes, missing values)

### Memory Issues

**Symptoms:** R crashes during fitting

**Solutions:**
1. Reduce number of chains: `n_chains = 2`
2. Increase thinning: `thin = 5`
3. Fit on a cloud server with more RAM
4. Split subjects into batches

### Slow Performance

**Expected runtimes:**
- Parameter recovery (50 subjects): 2-4 hours
- Full EEF model (146 subjects): 4-8 hours
- All three models: 12-24 hours

**Speed optimizations:**
1. Use cloud computing for full dataset
2. Run models in parallel on separate machines
3. Start with shorter runs for debugging

## Quick Diagnostic Checklist

After fitting, verify:

- [ ] R-hat < 1.1 for all parameters
- [ ] Effective sample size > 1000
- [ ] Trace plots show mixing (no trends)
- [ ] Posterior distributions are unimodal
- [ ] Parameter estimates are in expected ranges

## Parameter Ranges

Expected values for EEF model:

| Parameter | Range | Typical Values |
|-----------|-------|----------------|
| theta | [0, 1] | 0.2-0.5 |
| lambda_forget | [0, 1] | 0.3-0.6 |
| phi | [-5, 5] | 0-2 |
| cons | [0, 5] | 0.5-3 |

## Getting Help

1. Check JAGS error messages carefully
2. Run `analysis/quick_test.R` to isolate issues
3. Enable verbose output: `quiet = FALSE` in jags.model()

## References

- JAGS manual: https://mcmc-jags.sourceforge.io/
- Truncated distributions: Use `T(lower, upper)` syntax
