rm(list=ls())

library(rstan)
library(foreach)
library(doSNOW)

dataList <- readRDS("Data/Preprocessed/dataList_12datasets_stan_ready.rds")

##### Fit stan models #####
m <- stan_model("Stan/igt_pvl_delta.stan")

seeds <- c(rep(43210, 6), 43212, rep(43210, 5))

cl <- makeCluster(12)
registerDoSNOW(cl)
fits <- foreach(i=1:12, .export = c("m", "dataList")) %dopar% {
  library(rstan)
  set.seed(seeds[i])
  sampling(m, data = dataList[[names(dataList)[i]]], iter = 4000, warmup = 1500, chains = 4, cores = 4)
}
stopCluster(cl)

saveRDS(fits, file = "Data/Fitted/fits_all_pvl_delta.rds")
