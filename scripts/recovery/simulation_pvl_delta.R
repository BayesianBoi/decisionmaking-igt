# PVL-Delta simulation - generate fake data with known params
# uses deck-based indexing (position = how many times that deck was picked)

required_packages <- c("extraDistr", "truncnorm")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")
invisible(lapply(required_packages, library, character.only = TRUE))

hier_PVL_sim <- function(payoff_struct, nsubs, ntrials,
                         mu_w, mu_A, mu_a, mu_theta,
                         sigma_w, sigma_A, sigma_a, sigma_theta) {
    # Output arrays
    x <- array(NA, c(nsubs, max(ntrials))) # Choices (1-4)
    X <- array(NA, c(nsubs, max(ntrials))) # Net outcomes

    for (s in 1:nsubs) {
        # Sample subject-level parameters from group distributions
        # Loss aversion w >= 0: typically 1.5-2.5 in healthy adults
        w <- rtruncnorm(1, a = 0, b = Inf, mean = mu_w, sd = sigma_w)

        # Outcome sensitivity A: controls curvature of utility function
        A <- rtruncnorm(1, a = 0, b = Inf, mean = mu_A, sd = sigma_A)

        # Inverse temperature theta >= 0
        theta <- rtruncnorm(1, a = 0, b = Inf, mean = mu_theta, sd = sigma_theta)

        # Learning rate a (0-1)
        a <- rtruncnorm(1, a = 0, b = 1, mean = mu_a, sd = sigma_a)

        # Initialize state variables
        ev <- c(0, 0, 0, 0) # Expected values per deck
        pChoose <- c(0.25, 0.25, 0.25, 0.25) # Choice probabilities
        deckCount <- c(0, 0, 0, 0) # How many times each deck has been chosen

        # Trial loop
        for (t in 1:ntrials[s]) {
            # --- Compute choice probabilities with softmax ---
            exp_p <- exp(theta * ev)
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
                deck_position <- nrow(payoff_struct$gain)
            }

            gain <- payoff_struct$gain[deck_position, choice]
            loss <- payoff_struct$loss[deck_position, choice]
            outcome <- gain + loss # loss is already negative
            X[s, t] <- outcome

            # --- Calculate Utility using Prospect Theory ---
            # For gains: u = X^A
            # For losses: u = -w * |X|^A
            if (outcome >= 0) {
                curUtil <- abs(outcome)^A
            } else {
                curUtil <- -w * abs(outcome)^A
            }

            # --- Delta Rule Learning ---
            # Ev(t+1) = Ev(t) + a * [u - Ev(t)]
            PE <- curUtil - ev[choice]
            ev[choice] <- ev[choice] + a * PE
        }
    }

    return(list(x = x, X = X))
}
