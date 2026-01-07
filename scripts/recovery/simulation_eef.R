# ==============================================================================
# Dependencies
# Dependencies
required_packages <- c("extraDistr", "truncnorm")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")
lapply(required_packages, library, character.only = TRUE)
# ==============================================================================
# EEF (Explore-Exploit with Forgetting) Model Simulation
# ==============================================================================
#
# Reference: Yang et al. (2025)
#
# The EEF model conceptualizes decision-making as a balance between
# exploitation (using known values) and exploration (seeking information).
#
# IMPLEMENTATION NOTES (matching other simulations):
# 1. Uses DECK-BASED indexing: payoff position depends on how many times
#    that specific deck has been chosen, not the trial number
# 2. Payoff structure: list with $gain and $loss matrices
# ==============================================================================

simulation_eef <- function(payoff_struct, nsubs, ntrials,
                           mu_theta, mu_lambda, mu_phi, mu_cons,
                           sigma_theta, sigma_lambda, sigma_phi, sigma_cons) {
    # Output arrays
    x <- array(NA, c(nsubs, max(ntrials))) # Choices (1-4)
    X <- array(NA, c(nsubs, max(ntrials))) # Net outcomes

    for (s in 1:nsubs) {
        # -------------------------------------------------------------------------
        # Sample subject-level parameters from group distributions
        # -------------------------------------------------------------------------
        theta <- rtruncnorm(1, a = 0, b = 1, mean = mu_theta, sd = sigma_theta)
        lambda <- rtruncnorm(1, a = 0, b = 1, mean = mu_lambda, sd = sigma_lambda)
        phi <- rtruncnorm(1, a = -5, b = 5, mean = mu_phi, sd = sigma_phi)
        cons <- rtruncnorm(1, a = 0, b = 5, mean = mu_cons, sd = sigma_cons)

        # Consistency transformation: C = 3^cons - 1
        C <- 3^cons - 1

        # -------------------------------------------------------------------------
        # Initialize state variables
        # -------------------------------------------------------------------------
        Exploit <- c(0, 0, 0, 0)
        Explore <- c(0, 0, 0, 0)
        deckCount <- c(0, 0, 0, 0)

        # Get number of trials for this subject
        n_t <- ifelse(length(ntrials) == 1, ntrials, ntrials[s])

        # -------------------------------------------------------------------------
        # Trial loop
        # -------------------------------------------------------------------------
        for (t in 1:n_t) {
            # --- Compute choice probabilities ---
            V <- Exploit + Explore
            exp_p <- exp(C * V)

            # Handle numerical overflow
            if (any(is.infinite(exp_p))) {
                exp_p <- exp(C * V - max(C * V))
            }

            pChoose <- exp_p / sum(exp_p)

            if (any(is.na(pChoose))) {
                pChoose <- c(0.25, 0.25, 0.25, 0.25)
            }

            # --- Make choice ---
            choice <- rcat(1, pChoose)
            x[s, t] <- choice

            # --- Get outcome using DECK-BASED indexing ---
            deckCount[choice] <- deckCount[choice] + 1
            deck_position <- deckCount[choice]

            if (deck_position > nrow(payoff_struct$gain)) {
                deck_position <- nrow(payoff_struct$gain)
            }

            gain <- payoff_struct$gain[deck_position, choice]
            loss <- payoff_struct$loss[deck_position, choice]
            outcome <- gain + loss
            X[s, t] <- outcome

            # --- Calculate Utility ---
            if (outcome >= 0) {
                u_t <- abs(outcome)^theta
            } else {
                u_t <- -1 * abs(outcome)^theta
            }

            # --- Update Exploit and Explore ---
            for (d in 1:4) {
                if (d == choice) {
                    # Chosen deck: update exploit, reset explore
                    Exploit[d] <- (1 - lambda) * Exploit[d] + lambda * u_t
                    Explore[d] <- 0
                } else {
                    # Unchosen decks: decay exploit, accumulate explore
                    Exploit[d] <- (1 - lambda) * Exploit[d]
                    Explore[d] <- (1 - lambda) * Explore[d] + lambda * phi
                }
            }
        }
    }

    return(list(x = x, X = X))
}
