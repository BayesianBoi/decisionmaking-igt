# Install required packages for the IGT Analysis Pipeline

# Package Management
if (!require("pacman")) install.packages("pacman")

# Core Bayesian Modeling
pacman::p_load(rjags, R2jags, coda)

# Parallel Computing
pacman::p_load(parallel)

# Data Manipulation & Visualization
pacman::p_load(ggplot2, gridExtra, ggpubr)

# Model Comparison & Simulation
pacman::p_load(loo, extraDistr, truncnorm)

message("All requirements checked/installed.")
