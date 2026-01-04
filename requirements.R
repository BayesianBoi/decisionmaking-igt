# Install required packages for the IGT Analysis Pipeline

# Core Bayesian Modeling
if (!require("rjags")) install.packages("rjags")
if (!require("coda")) install.packages("coda")

# Parallel Computing
if (!require("parallel")) install.packages("parallel")

# Data Manipulation & Visualization
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("gridExtra")) install.packages("gridExtra")

# Model Comparison
if (!require("loo")) install.packages("loo")

message("All requirements checked/installed.")
