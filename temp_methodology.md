# 2. Methodology

## 2.1 Data Sources

We drew trial-level data from two published studies examining IGT performance in substance-using populations and healthy controls.

Ahn et al. (2014) contributed 129 participants across three groups. Forty-eight healthy controls, 38 individuals with a history of amphetamine dependence, and 43 with a history of heroin dependence. Substance-using participants had maintained abstinence for at least one month following acute withdrawal. Each participant completed 100 IGT trials.

Fridberg et al. (2010) contributed 32 participants. Fifteen healthy controls and 17 chronic cannabis users who had consumed regularly for at least two years. Abstinence was verified by urinalysis. Each completed 100 trials.

We pooled healthy controls from both studies into a single reference group of 63 participants. The final sample comprised 161 individuals across four groups for comparison.

## 2.2 Computational Models

The introduction described three models in detail. Here we summarise their implementation for hierarchical Bayesian estimation.

### 2.2.1 Model Parameterisation

**PVL-Delta** was implemented with four parameters following Ahn et al. (2008). Learning rate A governs the speed of expected value updating. Outcome sensitivity α shapes the curvature of the prospect-theoretic utility function. Loss aversion λ weights losses against equivalent gains. Choice consistency c governs how deterministically participants select higher-valued decks via the softmax rule. Expected values for each deck initialise at zero.

**ORL** extends the basic reinforcement learning framework by allowing asymmetric learning from rewards versus punishments, following Haines et al. (2018). Reward learning rate A_rew controls updating after gains. Punishment learning rate A_pun controls updating after losses. This asymmetry can capture optimistic or pessimistic learning biases observed in clinical populations. Win frequency sensitivity β_F determines how strongly the rate of winning (as opposed to magnitude) influences choice. Perseverance tendency β_P captures the inclination to repeat recent choices. Decay parameter K governs how quickly perseverance fades. Choice consistency was fixed at 1 for identifiability, following Haines et al. (2018).

**EEF** differs from the other two models in structure. Following Yang et al. (2025), it separates choice into exploitation and exploration components, both subject to memory decay governed by a single forgetting parameter λ. Outcome sensitivity θ shapes subjective utility through a symmetric power function applied equally to gains and losses without separate loss aversion. Forgetting rate λ controls decay of both exploitation weights and exploration tendencies. Exploration incentive φ determines the attractiveness of unchosen decks. Choice consistency c governs the softmax temperature. Exploitation weights initialise from empirically-derived first-choice frequencies rather than zero.

### 2.2.2 Hierarchical Prior Specification

All models were estimated hierarchically, with group-level distributions and individual-level parameters estimated simultaneously. This approach borrows strength across participants while preserving individual differences (Ahn et al., 2011).

For PVL-Delta, learning rate means received Beta(2, 2) priors. Outcome sensitivity means received Normal(0.7, 0.5) priors truncated at zero and two. Loss aversion means received Normal(2, 1) priors truncated at zero and ten. Consistency means received Normal(1, 0.5) priors truncated at zero and five.

For ORL, both learning rate means received Beta(2, 2) priors centred at 0.5, allowing the data to identify asymmetries between reward and punishment learning. Decay means received Normal(1, 1) priors truncated at zero and five, reflecting typical values observed in Haines et al. (2018). Frequency weight and perseverance weight means each received Normal(0, 0.5) priors truncated between -5 and 5, remaining neutral about the direction of these influences.

For EEF, outcome sensitivity means received Beta(1.5, 3) priors centred near 0.33 to reflect values reported in Yang et al. (2025). Forgetting rate means received Beta(2, 3) priors centred at 0.40, allowing the posterior to identify elevated decay in substance-using groups if present. Exploration incentive means received Normal(0, 0.5) priors truncated between -5 and 5. Consistency priors matched those in PVL-Delta.

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

Model comparison focused on convergence quality, parameter recovery, and theoretical interpretability. We compared the three models to determine whether forgetting (EEF), asymmetric learning (ORL), or baseline reinforcement learning (PVL-Delta) best accounts for substance-related IGT deficits.
