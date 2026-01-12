# ORL simulation - generate fake data with known params
# deck-based indexing, PS starts at 1, omega weights unbounded

required_packages <- c("extraDistr", "truncnorm")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")
invisible(lapply(required_packages, library, character.only = TRUE))

hier_ORL_sim <- function(payoff_struct, nsubs, ntrials,
                         mu_a_rew, mu_a_pun, mu_K, mu_omega_f, mu_omega_p,
                         sigma_a_rew, sigma_a_pun, sigma_K,
                         sigma_omega_f, sigma_omega_p) {
    # payoff_struct should be a list with $gain and $loss matrices
    # Each matrix is (ntrials x 4 decks)

    # Output arrays
    x <- array(NA, c(nsubs, max(ntrials))) # Choices (1-4)
    X <- array(NA, c(nsubs, max(ntrials))) # Net outcomes

    for (s in 1:nsubs) {
        # Sample subject-level parameters from group distributions
        a_rew <- rtruncnorm(1, a = 0, b = 1, mean = mu_a_rew, sd = sigma_a_rew)
        a_pun <- rtruncnorm(1, a = 0, b = 1, mean = mu_a_pun, sd = sigma_a_pun)
        K <- rtruncnorm(1, a = 0, b = Inf, mean = mu_K, sd = sigma_K)
        theta <- 1 # Fixed to 1

        # Omega weights are UNBOUNDED (can be negative)
        omega_f <- rnorm(1, mean = mu_omega_f, sd = sigma_omega_f)
        omega_p <- rnorm(1, mean = mu_omega_p, sd = sigma_omega_p)

        # Initialize state variables (matching inspiration code)
        ev <- c(0, 0, 0, 0) # Expected values per deck
        ef <- c(0, 0, 0, 0) # Expected frequencies per deck
        pers <- c(1, 1, 1, 1) # Perseverance initialized to 1 (matching JAGS model)
        pChoose <- c(0.25, 0.25, 0.25, 0.25) # Choice probabilities
        deckCount <- c(0, 0, 0, 0) # How many times each deck has been chosen

        # Trial loop
        for (t in 1:ntrials[s]) {
            # --- Compute choice probabilities ---
            V <- ev + ef * omega_f + pers * omega_p
            exp_p <- exp(theta * V)
            pChoose <- exp_p / sum(exp_p)

            # Handle numerical issues
            if (any(is.na(pChoose)) || any(is.infinite(pChoose))) {
                pChoose <- c(0.25, 0.25, 0.25, 0.25)
            }

            # --- Make choice ---
            choice <- rcat(1, pChoose)
            x[s, t] <- choice

            # --- Get outcome using DECK-BASED indexing ---
            deckCount[choice] <- deckCount[choice] + 1
            deck_position <- deckCount[choice]

            # Ensure we don't exceed payoff matrix dimensions
            if (deck_position > nrow(payoff_struct$gain)) {
                deck_position <- nrow(payoff_struct$gain) # Cap at max
            }

            gain <- payoff_struct$gain[deck_position, choice]
            loss <- payoff_struct$loss[deck_position, choice]
            outcome <- gain + loss # loss is already negative
            X[s, t] <- outcome

            # --- Compute binary outcome sign ---
            if (outcome > 0) {
                binOutcome <- 1
            } else if (outcome < 0) {
                binOutcome <- -1
            } else {
                binOutcome <- 0
            }

            # --- Update Expected Value (Ev) ---
            PEval <- outcome - ev[choice]
            if (outcome >= 0) {
                ev[choice] <- ev[choice] + a_rew * PEval
            } else {
                ev[choice] <- ev[choice] + a_pun * PEval
            }

            # --- Update Expected Frequency (Ef) ---
            PEfreq <- binOutcome - ef[choice]
            PEfreq_fic <- -binOutcome / 3 - ef[-choice] # Fictive update for unchosen

            if (outcome >= 0) {
                ef[choice] <- ef[choice] + a_rew * PEfreq
                ef[-choice] <- ef[-choice] + a_pun * PEfreq_fic
            } else {
                ef[choice] <- ef[choice] + a_pun * PEfreq
                ef[-choice] <- ef[-choice] + a_rew * PEfreq_fic
            }

            # --- Update Perseverance ---
            # Chosen deck gets reset to 1, then all decay
            pers[choice] <- 1
            pers <- pers / (1 + K)
        }
    }

    return(list(x = x, X = X))
}
