# ==============================================================================
# Modified IGT Payoff Schedule - Ahn et al. (2014) Exact Schedule
# ==============================================================================
#
# This file provides the exact payoff schedule used in the Ahn et al. (2014)
# datasets (HC, Amphetamine, Heroin groups). The schedule is taken directly
# from the inspiration code (R/utils/set_version_igt.R).
#
# Key Features of Modified IGT (vs Original):
# - Variable gains (not constant $100/$50)
# - Increasing loss magnitudes for Decks B and D across positions
# - 60 trials per deck (but we extend to 100 by cycling)
#
# Reference: Code adapted from inspiration implementation for methodological
# consistency with published recovery studies.
# ==============================================================================

#' Generate Modified IGT Payoff - Exact Ahn 2014 Schedule
#'
#' Returns the exact payoff structure used in Ahn et al. (2014) datasets.
#' This uses hardcoded schedules matching the original experimental task.
#'
#' @param ntrials Number of trials (default 100)
#' @param scale Logical, if TRUE divide all values by 100
#' @return List with $gain and $loss matrices (ntrials x 4)

generate_modified_igt_payoff <- function(ntrials = 100, scale = FALSE) {
    # =========================================================================
    # GAINS - Exact schedule from Ahn 2014
    # =========================================================================
    # Note: Original schedule has 60 values, we repeat for 100 trials

    gain_A_60 <- c(
        100, 120, 80, 90, 110, 100, 80, 120, 110, 90,
        110, 130, 90, 100, 120, 110, 90, 130, 120, 100,
        120, 140, 110, 110, 100, 120, 130, 110, 140, 120,
        130, 120, 140, 130, 110, 150, 140, 120, 150, 110,
        100, 120, 80, 90, 110, 100, 80, 120, 110, 90,
        110, 130, 90, 100, 120, 110, 90, 130, 120, 100
    )

    gain_B_60 <- c(
        100, 80, 110, 120, 90, 100, 90, 120, 110, 80,
        110, 100, 90, 130, 120, 130, 110, 90, 100, 120,
        120, 110, 140, 130, 100, 110, 120, 120, 140, 110,
        130, 140, 120, 110, 130, 150, 110, 150, 120, 140,
        140, 150, 130, 120, 140, 160, 120, 160, 130, 150,
        150, 160, 140, 130, 150, 170, 130, 170, 140, 160
    )

    gain_C_60 <- c(
        50, 60, 40, 55, 55, 45, 50, 45, 60, 40,
        55, 55, 65, 45, 70, 40, 50, 60, 70, 40,
        60, 65, 55, 80, 40, 60, 55, 65, 40, 80,
        65, 75, 55, 60, 70, 65, 55, 75, 45, 85,
        70, 80, 60, 65, 75, 70, 60, 80, 50, 90,
        75, 85, 65, 70, 80, 75, 65, 85, 55, 95
    )

    gain_D_60 <- c(
        50, 40, 45, 45, 55, 60, 40, 55, 50, 60,
        55, 40, 60, 40, 45, 55, 65, 70, 50, 70,
        60, 55, 65, 80, 40, 80, 40, 65, 55, 60,
        65, 75, 60, 65, 75, 85, 45, 55, 70, 55,
        70, 80, 65, 70, 80, 90, 50, 60, 75, 60,
        75, 85, 70, 75, 85, 95, 55, 65, 80, 65
    )

    # =========================================================================
    # LOSSES - Exact schedule from Ahn 2014
    # =========================================================================
    # Deck A: Frequent small losses (50% of trials)
    loss_A_60 <- c(
        0, 0, -150, 0, -300, 0, -200, 0, -250, -350,
        0, -350, 0, -250, -200, 0, -300, -150, -250, 0,
        -250, -300, 0, -350, 0, -200, -250, -150, -250, 0,
        -350, -200, -250, -250, -150, 0, -150, -300, -350, 0,
        0, 0, -150, 0, -300, 0, -200, 0, -250, -350,
        0, -350, 0, -250, -200, 0, -300, -150, -250, 0
    )

    # Deck B: Infrequent large losses (increasing magnitude)
    # Position 9: -1250, Position 14: -1500, Position 21: -1750,
    # Position 32: -2000, Position 46: -2250, Position 58: -2500
    loss_B_60 <- c(
        0, 0, 0, 0, 0, 0, 0, 0, -1250, 0,
        0, 0, 0, -1500, 0, 0, 0, 0, 0, 0,
        -1750, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, -2000, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, -2250, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, -2500, 0, 0
    )

    # Deck C: Frequent small losses
    loss_C_60 <- c(
        0, 0, -50, 0, -50, 0, -50, 0, -50, -50,
        0, -25, -75, 0, -25, 0, -25, -75, 0, -50,
        0, -25, 0, -50, -25, -50, 0, -25, -75, -50,
        -25, 0, -25, -25, -25, 0, -75, -25, -50, -75,
        -25, 0, -25, -25, -25, -25, -75, -25, -50, -75,
        -25, -25, -25, -25, -25, -25, -75, -25, -50, -75
    )

    # Deck D: Infrequent small losses (increasing magnitude)
    loss_D_60 <- c(
        0, 0, 0, 0, 0, 0, 0, 0, 0, -250,
        0, 0, 0, 0, 0, 0, 0, 0, 0, -275,
        0, 0, 0, 0, 0, 0, 0, 0, -300, 0,
        0, 0, 0, 0, -325, 0, 0, 0, 0, 0,
        0, 0, 0, 0, -350, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, -375, 0, 0
    )

    # =========================================================================
    # Extend to ntrials by cycling (if needed)
    # =========================================================================
    extend_schedule <- function(vec, n) {
        if (n <= length(vec)) {
            return(vec[1:n])
        } else {
            # Cycle the schedule
            return(rep(vec, ceiling(n / length(vec)))[1:n])
        }
    }

    gain <- matrix(0, nrow = ntrials, ncol = 4)
    loss <- matrix(0, nrow = ntrials, ncol = 4)
    colnames(gain) <- colnames(loss) <- c("A", "B", "C", "D")

    gain[, 1] <- extend_schedule(gain_A_60, ntrials)
    gain[, 2] <- extend_schedule(gain_B_60, ntrials)
    gain[, 3] <- extend_schedule(gain_C_60, ntrials)
    gain[, 4] <- extend_schedule(gain_D_60, ntrials)

    loss[, 1] <- extend_schedule(loss_A_60, ntrials)
    loss[, 2] <- extend_schedule(loss_B_60, ntrials)
    loss[, 3] <- extend_schedule(loss_C_60, ntrials)
    loss[, 4] <- extend_schedule(loss_D_60, ntrials)

    if (scale) {
        gain <- gain / 100
        loss <- loss / 100
    }

    return(list(gain = gain, loss = loss))
}


#' Get net payoff matrix (gain + loss)
#'
#' Convenience function for models that use net outcomes
get_net_payoff <- function(ntrials = 100, scale = FALSE) {
    payoff <- generate_modified_igt_payoff(ntrials, scale = scale)
    return(payoff$gain + payoff$loss)
}


#' Verify payoff structure
verify_payoff <- function(ntrials = 100) {
    payoff <- generate_modified_igt_payoff(ntrials)
    net <- payoff$gain + payoff$loss

    cat("=== Ahn 2014 Exact Payoff Schedule ===\n\n")

    cat("Gains per deck (first 10 positions):\n")
    cat("  Deck A:", payoff$gain[1:10, 1], "\n")
    cat("  Deck B:", payoff$gain[1:10, 2], "\n")
    cat("  Deck C:", payoff$gain[1:10, 3], "\n")
    cat("  Deck D:", payoff$gain[1:10, 4], "\n\n")

    cat("Deck B losses (non-zero):\n")
    loss_B_pos <- which(payoff$loss[, 2] != 0)
    cat("  Positions:", loss_B_pos, "\n")
    cat("  Values:", payoff$loss[loss_B_pos, 2], "\n\n")

    cat("Net outcomes per block (10 trials):\n")
    cat("Block |    A    |    B    |    C    |    D    |\n")
    cat("------|---------|---------|---------|---------|")
    for (b in 1:min(6, ceiling(ntrials / 10))) {
        start <- (b - 1) * 10 + 1
        end <- min(b * 10, ntrials)
        cat(sprintf(
            "\n  %d   | %6.0f  | %6.0f  | %6.0f  | %6.0f  |",
            b, sum(net[start:end, 1]), sum(net[start:end, 2]),
            sum(net[start:end, 3]), sum(net[start:end, 4])
        ))
    }
    cat("\n")
}
