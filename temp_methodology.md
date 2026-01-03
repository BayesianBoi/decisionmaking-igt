# 2. Methodology

## 2.1 Data Sources

We drew trial-level data from two published studies that examined IGT performance in substance-using populations and healthy controls.

Ahn et al. (2014) contributed 129 participants across three groups. Forty-eight healthy controls, 38 individuals with a history of amphetamine dependence, and 43 with a history of heroin dependence. Substance-using participants had maintained abstinence for at least one month following acute withdrawal. Each participant completed 100 IGT trials.

Fridberg et al. (2010) contributed 32 participants. Fifteen healthy controls and 17 chronic cannabis users who had consumed regularly for at least two years. Abstinence was verified by urinalysis. Each completed 100 trials.

We pooled healthy controls from both studies into a single reference group of 63 participants. The final sample comprised 161 individuals across four groups for comparison.

## 2.2 Computational Models

The introduction described three models in detail. Here we summarise their implementation for hierarchical Bayesian estimation.

### 2.2.1 Model Parameterisation

**PVL-Delta** was implemented with four parameters following Ahn et al. (2008). Learning rate A governs the speed of expected value updating. Outcome sensitivity α shapes the curvature of the prospect-theoretic utility function. Loss aversion λ weights losses against equivalent gains. Choice consistency c governs how deterministically participants select higher-valued decks via the softmax rule. Expected values for each deck initialise at zero.

**VSE** extends PVL-Delta with a perseverance mechanism following Worthy et al. (2013). Beyond the four PVL-Delta parameters, it adds perseverance boosts after positive outcomes (ε_P) and negative outcomes (ε_N), a decay rate K governing how quickly perseverance fades, and a weight w that balances value-based learning against perseverance. Choice probabilities depend on a weighted combination of expected values and perseverance traces.

**EEF** differs from the other two models in structure. Following Yang et al. (2025), it separates choice into exploitation and exploration components, both subject to memory decay governed by a single forgetting parameter λ. The model uses four parameters. Outcome sensitivity θ shapes subjective utility through a symmetric power function applied equally to gains and losses without separate loss aversion. Forgetting rate λ controls decay of both exploitation weights and exploration tendencies. Exploration incentive φ determines the attractiveness of unchosen decks. Choice consistency c governs the softmax temperature. Exploitation weights initialise from empirically-derived first-choice frequencies rather than zero.

### 2.2.2 Hierarchical Prior Specification

All models were estimated hierarchically, with group-level distributions and individual-level parameters estimated simultaneously. This approach borrows strength across participants while preserving individual differences (Ahn et al., 2011).

For PVL-Delta and VSE, learning rate means received Beta(2, 2) priors. Outcome sensitivity means received Normal(0.7, 0.5) priors truncated at zero and two. Loss aversion means received Normal(2, 1) priors truncated at zero and ten. Consistency means received Normal(1, 0.5) priors truncated at zero and five. VSE-specific perseverance parameters received Normal(0.5, 0.5) priors truncated at zero, with decay and weight parameters receiving Beta(2, 2) priors.

For EEF, outcome sensitivity means received Beta(1.5, 3) priors centred near 0.33 to reflect values reported in Yang et al. (2025). Forgetting rate means received Beta(2, 3) priors centred at 0.40, allowing the posterior to identify elevated decay in substance-using groups if present. Exploration incentive means received Normal(0, 0.5) priors truncated between -5 and 5, remaining neutral about whether participants would show novelty-seeking or exploration aversion. Consistency priors matched those in PVL-Delta.

Group-level standard deviations came from uniform priors with upper bounds chosen to constrain individual variation within plausible ranges given the sample sizes.

## 2.3 Estimation

We implemented all models in JAGS 4.3.2 using Markov Chain Monte Carlo sampling. Four independent chains ran in parallel, each beginning with 5,000 adaptation iterations for sampler tuning. We discarded 10,000 burn-in iterations before collecting 20,000 posterior samples per chain. A thinning interval of 2 reduced autocorrelation, yielding 40,000 retained samples.

Convergence was assessed via the Gelman-Rubin statistic (R-hat below 1.1 for all parameters) and effective sample size (above 1,000 for all parameters). We inspected trace plots visually to confirm adequate mixing.

## 2.4 Model Validation

**Parameter recovery.** Before fitting empirical data, we simulated synthetic datasets from known EEF parameter values drawn from plausible ranges. Re-estimating parameters from these synthetic data allowed us to verify that the model can reliably recover its parameters given the structure and sample size of the observed data.

**Convergence diagnostics.** All parameters were required to meet the R-hat and effective sample size thresholds before we accepted posterior estimates.

**Posterior predictive checks.** We generated synthetic choice sequences by sampling from posterior predictive distributions, then compared summary statistics of simulated and observed data to assess model fit.

## 2.5 Statistical Analysis

For group comparisons, we extracted subject-level posterior distributions of the forgetting parameter λ from the EEF model. We computed the posterior probability that substance-using groups exhibit higher forgetting rates than healthy controls by calculating the proportion of posterior samples where the group difference exceeded zero. This Bayesian approach avoids the dichotomous logic of null hypothesis testing while providing direct probability statements about the magnitude and direction of group differences.

Model comparison relied on convergence quality and parameter recovery performance across the three models.
