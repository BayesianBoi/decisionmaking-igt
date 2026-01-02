rm(list=ls())

library(rstan)
library(foreach)
library(doSNOW)

dataList <- readRDS("Data/Preprocessed/dataList_12datasets_stan_ready.rds")

##### Fit stan models #####
m <- stan_model("Stan/igt_vpp.stan")

seeds <- c(rep(43210, 4), 43218, rep(43210, 6), 43211)
deltas <- c(rep(.8, 4), .9, rep(.8, 6), .9)

cl <- makeCluster(12)
registerDoSNOW(cl)
fits <- foreach(i=1:12, .export = c("m", "dataList")) %dopar% {
  library(rstan)
  set.seed(seeds[i])
  sampling(m, data = dataList[[names(dataList)[i]]], iter = 4000, warmup = 1500, chains = 4, cores = 4, control = list(adapt_delta=deltas[i]))
}
stopCluster(cl)

saveRDS(fits, file = "Data/Fitted/fits_all_vpp.rds")
