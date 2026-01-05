# Prior Specification and Justification

This document provides transparent justification for all prior distributions used in the hierarchical Bayesian models. We distinguish between priors derived directly from published values and priors that represent our own reasonable choices informed by the literature.

## General Approach

We used informative priors where published empirical estimates exist, and weakly informative priors elsewhere. All individual-level parameters were drawn from group-level distributions estimated simultaneously with the data.

For bounded parameters (learning rates, sensitivity parameters), we used Beta distributions that place mass away from boundary values. For unbounded parameters (frequency weights, exploration incentives), we used Normal distributions centred at zero to remain agnostic about effect direction.

---

## PVL-Delta Model

### Learning Rate (A)
- **Prior**: Beta(2, 2), bounds [0, 1]
- **Justification**: Weakly informative prior symmetric around 0.5. This is a standard choice for learning rates in reinforcement learning models, not derived from any specific paper. The beta shape places mass away from extreme values (0 or 1) which would indicate no learning or complete replacement of expectations.

### Outcome Sensitivity (α)
- **Prior**: Normal(0.7, precision=4), truncated [0, 2]
- **Justification**: Centred on empirical estimates from Ahn et al. (2008), who reported mean α values in the range 0.6-0.8 across multiple studies. The specific prior distribution is our choice; Ahn et al. (2008) did not specify this prior.

### Loss Aversion (λ)
- **Prior**: Normal(2, precision=1), truncated [0, 10]
- **Justification**: Centred on the well-established finding from Kahneman and Tversky (1979) that losses are weighted approximately 2-2.5 times more heavily than equivalent gains. This is one of the most robust findings in behavioural economics.

### Consistency (c)
- **Prior**: Normal(1, precision=2), truncated [0, 5]
- **Justification**: Weakly informative prior. Values around 1-2 produce reasonable choice stochasticity. This is our choice, not derived from a specific source.

---

## ORL Model

### Learning Rates (Arew, Apun)
- **Prior**: Beta(2, 2), bounds [0, 1]
- **Justification**: Same rationale as PVL-Delta learning rate. Haines et al. (2018) used probit-transformed priors with normal(0,1) group means, which produces approximately uniform coverage over (0,1). Our Beta(2,2) prior is similar in spirit but not identical to their implementation.

### Decay (K)
- **Prior**: Normal(1, precision=1), truncated [0, 5]
- **Justification**: Centred on typical values observed in Haines et al. (2018). The specific prior distribution is our choice.

### Frequency Weight (βF) and Perseverance Weight (βP)
- **Prior**: Normal(0, precision=2), truncated [-5, 5]
- **Justification**: Sceptical priors centred at zero, allowing the data to determine whether these effects are positive or negative. Haines et al. (2018) used similar sceptical priors with half-Cauchy(0,1) on the group-level standard deviation.

---

## EEF Model

### Outcome Sensitivity (θ)
- **Prior**: Beta(1.5, 3), bounds [0, 1]
- **Justification**: Places most prior mass below 0.5, consistent with the expectation that outcome sensitivity is typically moderate. We followed the model structure from Yang et al. (2025), though the specific prior values were chosen by us to match their reported empirical estimates.

### Forgetting Rate (λ)
- **Prior**: Beta(2, 3), bounds [0, 1]
- **Justification**: Centred near 0.4, reflecting values reported in Yang et al. (2025) for healthy adults. The prior allows for both low forgetting (memory retention) and high forgetting (myopic responding) depending on the data.

### Exploration Incentive (φ)
- **Prior**: Normal(0, precision=2), truncated [-5, 5]
- **Justification**: Sceptical prior centred at zero. Positive values indicate attraction to unchosen decks; negative values indicate avoidance of uncertainty.

### Consistency (c)
- **Prior**: Normal(1, precision=2), truncated [0, 5]
- **Justification**: Same as PVL-Delta.

---

## Group-Level Standard Deviations

- **Prior**: Uniform(0, upper)
- **Upper bounds**: 0.5 for Beta-distributed parameters; 1.5-2.0 for Normal-distributed parameters
- **Justification**: Weakly informative priors that prevent implausibly large individual variation. Upper bounds were chosen based on the sample size and expected heterogeneity in clinical populations.

---

## References

- Ahn, W. Y., Busemeyer, J. R., & Wagenmakers, E. J. (2008). Comparison of decision learning models using the generalization criterion method. Cognitive Science, 32(8), 1376-1402.
- Haines, N., Vassileva, J., & Ahn, W. Y. (2018). The Outcome-Representation Learning model: A novel reinforcement learning model of the Iowa Gambling Task. Cognitive Science, 42(8), 2534-2561.
- Kahneman, D., & Tversky, A. (1979). Prospect theory: An analysis of decision under risk. Econometrica, 47(2), 263-291.
- Yang, T., Xie, C., & Wang, X. (2025). Forgetting phenomena in the Iowa Gambling Task: A new computational model among diverse participants. Frontiers in Psychology.
