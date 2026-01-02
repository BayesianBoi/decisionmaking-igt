rm(list=ls())

library(rstan)

dataList <- readRDS("Data/Simulated/orl_originalIGT_N48_stan_ready.rds")

##### Fit stan models #####
m <- stan_model("Stan/igt_orl.stan")

set.seed(43211)
fit <- sampling(m, data = dataList, iter = 4000, warmup = 1500, chains = 4, cores = 4)
saveRDS(fit, file = "Data/Recovered/orl_originalIGT_N48_recovered.rds")
