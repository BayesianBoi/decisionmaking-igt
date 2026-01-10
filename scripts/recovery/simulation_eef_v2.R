# =============================================================================
# Simulation function for EEF v2 (Yang et al. 2025)
# =============================================================================
# Generates choice data from known parameter values
# Matches the exact specification from the original paper
# =============================================================================

simulation_eef_v2 <- function(payoff_struct, nsubs, ntrials,
                              mu_theta, mu_lambda, mu_phi, mu_cons,
                              sigma_theta, sigma_lambda, sigma_phi, sigma_cons) {
    # payoff_struct should have $gain and $loss matrices (ntrials_max x 4)
    # ntrials is a vector of length nsubs

    ntrials_max <- max(ntrials)

    # Output arrays
    x <- matrix(NA, nrow = nsubs, ncol = ntrials_max)
    X <- matrix(NA, nrow = nsubs, ncol = ntrials_max) # net outcome for reference

    # Store individual parameters for recovery check
    true_theta <- numeric(nsubs)
    true_lambda <- numeric(nsubs)
    true_phi <- numeric(nsubs)
    true_cons <- numeric(nsubs)

    for (s in 1:nsubs) {
        # Sample subject-level parameters from truncated normal
        theta_s <- rtruncnorm(1, a = 0, b = 1, mean = mu_theta, sd = sigma_theta)
        lambda_s <- rtruncnorm(1, a = 0, b = 1, mean = mu_lambda, sd = sigma_lambda)
        phi_s <- rtruncnorm(1, a = -5, b = 5, mean = mu_phi, sd = sigma_phi)
        cons_s <- rtruncnorm(1, a = 0, b = 5, mean = mu_cons, sd = sigma_cons)

        true_theta[s] <- theta_s
        true_lambda[s] <- lambda_s
        true_phi[s] <- phi_s
        true_cons[s] <- cons_s

        # Consistency transformation
        C_s <- 3 * cons_s - 1

        # Initialize exploitation and exploration weights
        Exploit <- rep(0, 4)
        Explore <- rep(0, 4)

        # Track deck-specific positions for payoff lookup
        deck_pos <- rep(1, 4)

        # Trial 1: random choice (no prior info)
        x[s, 1] <- sample(1:4, 1)
        deck_chosen <- x[s, 1]
        gain_t1 <- payoff_struct$gain[deck_pos[deck_chosen], deck_chosen]
        loss_t1 <- payoff_struct$loss[deck_pos[deck_chosen], deck_chosen]
        X[s, 1] <- gain_t1 + loss_t1
        deck_pos[deck_chosen] <- deck_pos[deck_chosen] + 1

        for (t in 2:ntrials[s]) {
            # Get outcome from previous trial
            prev_choice <- x[s, t - 1]
            gain_prev <- payoff_struct$gain[deck_pos[prev_choice] - 1, prev_choice]
            loss_prev <- abs(payoff_struct$loss[deck_pos[prev_choice] - 1, prev_choice])

            # Value function: V(t) = Gain^theta - Loss^theta
            V_t <- gain_prev^theta_s - loss_prev^theta_s

            # Update exploitation weights
            for (d in 1:4) {
                if (d == prev_choice) {
                    # Chosen deck: (1-lambda)*Exploit + V(t)
                    Exploit[d] <- (1 - lambda_s) * Exploit[d] + V_t
                } else {
                    # Unchosen deck: (1-lambda)*Exploit
                    Exploit[d] <- (1 - lambda_s) * Exploit[d]
                }
            }

            # Update exploration weights
            for (d in 1:4) {
                if (d == prev_choice) {
                    # Chosen deck: reset to 0
                    Explore[d] <- 0
                } else {
                    # Unchosen deck: lambda*Explore + (1-lambda)*phi
                    Explore[d] <- lambda_s * Explore[d] + (1 - lambda_s) * phi_s
                }
            }

            # Compute choice probabilities via softmax
            total_val <- Exploit + Explore
            a <- C_s * total_val
            a <- a - max(a) # stable softmax shift
            exp_vals <- exp(a)
            probs <- exp_vals / sum(exp_vals)

            # Sample choice
            x[s, t] <- sample(1:4, 1, prob = probs)

            # Get outcome for current trial
            deck_chosen <- x[s, t]
            gain_t <- payoff_struct$gain[deck_pos[deck_chosen], deck_chosen]
            loss_t <- payoff_struct$loss[deck_pos[deck_chosen], deck_chosen]
            X[s, t] <- gain_t + loss_t
            deck_pos[deck_chosen] <- deck_pos[deck_chosen] + 1
        }
    }

    # Build Gain and Loss matrices for JAGS (trial t gets outcome from that trial)
    Gain_mat <- matrix(0, nrow = nsubs, ncol = ntrials_max)
    Loss_mat <- matrix(0, nrow = nsubs, ncol = ntrials_max)

    for (s in 1:nsubs) {
        deck_pos <- rep(1, 4)
        for (t in 1:ntrials[s]) {
            deck_chosen <- x[s, t]
            Gain_mat[s, t] <- payoff_struct$gain[deck_pos[deck_chosen], deck_chosen]
            Loss_mat[s, t] <- payoff_struct$loss[deck_pos[deck_chosen], deck_chosen]
            deck_pos[deck_chosen] <- deck_pos[deck_chosen] + 1
        }
    }

    return(list(
        x = x,
        X = X,
        Gain = Gain_mat,
        Loss = Loss_mat,
        true_theta = true_theta,
        true_lambda = true_lambda,
        true_phi = true_phi,
        true_cons = true_cons
    ))
}
