# Simulate EEF Mechanism (Generative)
# Generates Figure 3: Probability of Choosing Disadvantageous Deck B
# Based on 1000 simulations of the EEF model.

library(ggplot2)
library(dplyr)
library(tidyr)

# Set seed
set.seed(123)

# ==============================================================================
# 1. Define IGT Payoff Structure (Standard Payoff Scheme)
# ==============================================================================
# Deck A: -250/10 trials (Freq -250 loss every 2)
# Deck B: -250/10 trials (Freq -1250 loss every 10)
# Deck C: +250/10 trials (Freq -50 loss every 2)
# Deck D: +250/10 trials (Freq -250 loss every 10)

cards <- list(
    A = list(win = 100, loss = -250, freq = 0.5),
    B = list(win = 100, loss = -1250, freq = 0.1),
    C = list(win = 50, loss = -50, freq = 0.5),
    D = list(win = 50, loss = -250, freq = 0.1)
)

get_feedback <- function(deck_idx) {
    # 1=A, 2=B, 3=C, 4=D
    deck_char <- LETTERS[deck_idx]
    pars <- cards[[deck_char]]

    # Determine loss
    is_loss <- runif(1) < pars$freq
    loss_val <- ifelse(is_loss, pars$loss, 0)

    net <- pars$win + loss_val
    return(net)
}

# ==============================================================================
# 2. Generative Agent (EEF)
# ==============================================================================
# Using Equation from eef.jags:
# Q(t+1) = (1 - lambda) * Q(t) + Utility

simulate_agent <- function(n_trials, lambda, theta = 0.5, cons = 1.5) {
    q_values <- numeric(4) # A, B, C, D (start at 0)
    p_history_B <- numeric(n_trials) # Track Prob(Deck B)

    for (t in 1:n_trials) {
        # 1. Calculate Probabilities (Softmax)
        exp_v <- exp(cons * q_values)
        probs <- exp_v / sum(exp_v)
        p_history_B[t] <- probs[2] # Store Prob(B)

        # 2. Make Choice
        choice <- sample(1:4, 1, prob = probs)

        # 3. Get Outcome
        net_outcome <- get_feedback(choice)

        # 4. Utility
        # util = sign(x) * |x|^theta
        util <- sign(net_outcome) * (abs(net_outcome)^theta)

        # 5. Update Q-Values (EEF Rule)
        # Chosen: Q_new = (1-lambda)*Q + Utility
        # Unchosen: Q_new = (1-lambda)*Q

        # Apply decay to ALL
        q_values <- (1 - lambda) * q_values

        # Add utility to CHOSEN only
        q_values[choice] <- q_values[choice] + util
    }
    return(p_history_B)
}

# ==============================================================================
# 3. Run Monte Carlo Simulation (Average over agents)
# ==============================================================================
n_sims <- 1000
n_trials <- 100

# Condition 1: Low Forgetting (Control-like)
# lambda = 0.1 (Retains 90% per trial)
p_matrix_low <- replicate(n_sims, simulate_agent(n_trials, lambda = 0.1))
mean_p_low <- rowMeans(p_matrix_low)

# Condition 2: High Forgetting (Substance-like)
# lambda = 0.5 (Retains 50% per trial)
p_matrix_high <- replicate(n_sims, simulate_agent(n_trials, lambda = 0.6))
mean_p_high <- rowMeans(p_matrix_high)

# ==============================================================================
# 4. Prepare Plot Data
# ==============================================================================
df_plot <- data.frame(
    Trial = 1:n_trials,
    Control = mean_p_low,
    Substance = mean_p_high
) %>%
    pivot_longer(
        cols = c("Control", "Substance"),
        names_to = "Group",
        values_to = "Prob_B"
    ) %>%
    mutate(Group = factor(Group,
        levels = c("Control", "Substance"),
        labels = c("Retention Agent (lambda=0.1)", "Decay Agent (lambda=0.6)")
    ))

# ==============================================================================
# 5. Plot
# ==============================================================================
p <- ggplot(df_plot, aes(x = Trial, y = Prob_B, color = Group)) +
    # Chance Line
    geom_hline(yintercept = 0.25, linetype = "dashed", color = "grey50") +

    # Smooth Traces
    geom_line(linewidth = 1.2) +

    # Colors
    scale_color_manual(values = c("Retention Agent (lambda=0.1)" = "#0072B2", "Decay Agent (lambda=0.6)" = "#D55E00")) +

    # Theme
    theme_classic(base_size = 14) +
    scale_y_continuous(limits = c(0, 0.6), expand = c(0, 0)) +
    labs(
        x = "Trial Number",
        y = "Probability of Choosing Deck B",
        color = NULL
    ) +
    theme(
        legend.position = c(0.70, 0.85),
        legend.background = element_rect(fill = "white", color = NA),
        axis.line = element_line(linewidth = 0.5),
        panel.grid.major.y = element_line(color = "grey95")
    )

ggsave("results/figures/mechanism_simulation.png", p, width = 7, height = 4.5, dpi = 300)
print("Saved generative simulation plot (Prob Deck B).")
