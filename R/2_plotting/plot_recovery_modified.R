rm(list = ls())

# Libraries
library(ggplot2)
library(hBayesDM)
library(tidyr)


# Function for relative scaling
scale_rel <- function(x,by) {
  (x - mean(by)) / sd(by)
}
point_size <- 2

# ORL
orl_AhnReal <- read.table("Data/Simulated/orl_simulated_indPars.txt", header = T)
fit_orl <- readRDS("Data/Recovered/orl_modifiedIGT_N48_recovered.rds") # 
# PVL Delta
pvl_AhnReal <- read.table("Data/Simulated/pvl_delta_simulated_indPars.txt", header = T)
fit_pvl <- readRDS("Data/Recovered/pvl_delta_modifiedIGT_N48_recovered.rds") # 
# VPP
vpp_AhnReal <- read.table("Data/Simulated/vpp_simulated_indPars.txt", header = T)
fit_vpp <- readRDS("Data/Recovered/vpp_modifiedIGT_N48_recovered.rds") # 

## Extract parameters
pars_orl <- rstan::extract(fit_orl, permuted=T)
pars_pvl <- rstan::extract(fit_pvl, permuted=T)
pars_vpp <- rstan::extract(fit_vpp, permuted=T)

# Individual parameters (e.g., individual posterior means)
indPars_orl <- array(NA, c(48, 5))
indPars_orl <- as.data.frame(indPars_orl)
indPars_pvl <- array(NA, c(48, 4))
indPars_pvl <- as.data.frame(indPars_pvl)
indPars_vpp <- array(NA, c(48, 8))
indPars_vpp <- as.data.frame(indPars_vpp)

for (i in 1:48) {
    indPars_orl[i, ] <- c( mean(pars_orl$Arew[, i]), 
                          mean(pars_orl$Apun[, i]),
                          mean(pars_orl$K[, i]),
                          mean(pars_orl$betaF[, i]),
                          mean(pars_orl$betaP[, i]))
    indPars_pvl[i, ] <- c( mean(pars_pvl$A[, i]), 
                           mean(pars_pvl$alpha[, i]),
                           mean(pars_pvl$cons[, i]),
                           mean(pars_pvl$lambda[, i]))
    indPars_vpp[i, ] <- c( mean(pars_vpp$A[, i]), 
                           mean(pars_vpp$alpha[, i]),
                           mean(pars_vpp$cons[, i]),
                           mean(pars_vpp$lambda[, i]),
                           mean(pars_vpp$epP[, i]), 
                           mean(pars_vpp$epN[, i]),
                           mean(pars_vpp$K[, i]),
                           mean(pars_vpp$w[, i]))
}
indPars_orl           <- cbind(indPars_orl, 1:48)
indPars_pvl           <- cbind(indPars_pvl, 1:48)
indPars_vpp           <- cbind(indPars_vpp, 1:48)
colnames(indPars_orl) <- c("Arew", "Apun", "K", "betaF", "betaP", "subjID")
colnames(indPars_pvl) <- c("A", "alpha", "cons", "lambda", "subjID")
colnames(indPars_vpp) <- c("A", "alpha", "cons", "lambda", "epP", "epN", "K", "w", "subjID")

# OVL
indPars_orl$Arew_real  <- orl_AhnReal$Arew
indPars_orl$Apun_real  <- orl_AhnReal$Apun
indPars_orl$K_real     <- orl_AhnReal$K
indPars_orl$betaF_real <- orl_AhnReal$betaF
indPars_orl$betaP_real <- orl_AhnReal$betaP
# PVL Delta
indPars_pvl$A_real      <- pvl_AhnReal$A
indPars_pvl$alpha_real  <- pvl_AhnReal$alpha
indPars_pvl$cons_real   <- pvl_AhnReal$cons
indPars_pvl$lambda_real <- pvl_AhnReal$lambda
# VPP
indPars_vpp$A_real      <- vpp_AhnReal$A
indPars_vpp$alpha_real  <- vpp_AhnReal$alpha
indPars_vpp$cons_real   <- vpp_AhnReal$cons
indPars_vpp$lambda_real <- vpp_AhnReal$lambda
indPars_vpp$epP_real    <- vpp_AhnReal$epP
indPars_vpp$epN_real    <- vpp_AhnReal$epN
indPars_vpp$K_real      <- vpp_AhnReal$K
indPars_vpp$w_real      <- vpp_AhnReal$w


# Plotting 
# PVL Delta
h1 <- ggplot() +
  geom_point(data = indPars_pvl, 
             aes(x = scale(A_real), y = scale_rel(A,A_real)), 
             size = point_size, alpha = 0.4, color = I("#e41a1c")) + 
  geom_point(data = indPars_pvl, 
             aes(x = scale(alpha_real), y = scale_rel(alpha,alpha_real)), 
             size = point_size, alpha = 0.4, color = I("#377eb8")) + 
  geom_point(data = indPars_pvl, 
             aes(x = scale(cons_real), y = scale_rel(cons,cons_real)), 
             size = point_size, alpha = 0.4, color = I("#4daf4a")) + 
  geom_point(data = indPars_pvl, 
             aes(x = scale(lambda_real), y = scale_rel(lambda,lambda_real)), 
             size = point_size, alpha = 0.4, color = I("#984ea3")) +  
  coord_cartesian(xlim = c(-4,4), ylim = c(-4,4)) +
  geom_abline(intercept = 0, slope = 1) + 
  geom_abline(intercept = 1, slope = 1, lty = 2, color = "black") + 
  geom_abline(intercept = 2, slope = 1, lty = 3, color = "black") + 
  geom_abline(intercept = -1, slope = 1, lty = 2, color = "black") +
  geom_abline(intercept = -2, slope = 1, lty = 3, color = "black") +
  theme(legend.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  lazerhawk::theme_trueMinimal()

h2 <- ggplot() +
  geom_point(data = indPars_vpp, 
             aes(x = scale(A_real), y = scale_rel(A,A_real)), 
             size = point_size, alpha = 0.4, color = I("#e41a1c")) + 
  geom_point(data = indPars_vpp, 
             aes(x = scale(alpha_real), y = scale_rel(alpha,alpha_real)), 
             size = point_size, alpha = 0.4, color = I("#377eb8")) + 
  geom_point(data = indPars_vpp, 
             aes(x = scale(cons_real), y = scale_rel(cons,cons_real)), 
             size = point_size, alpha = 0.4, color = I("#4daf4a")) + 
  geom_point(data = indPars_vpp, 
             aes(x = scale(lambda_real), y = scale_rel(lambda,lambda_real)), 
             size = point_size, alpha = 0.4, color = I("#984ea3")) +  
  geom_point(data = indPars_vpp, 
             aes(x = scale(epP_real), y = scale_rel(epP,epP_real)), 
             size = point_size, alpha = 0.4, color = I("#ff7f00")) + 
  geom_point(data = indPars_vpp, 
             aes(x = scale(epN_real), y = scale_rel(epN,epN_real)), 
             size = point_size, alpha = 0.4, color = I("#ffff33")) + 
  geom_point(data = indPars_vpp, 
             aes(x = scale(K_real), y = scale_rel(K,K_real)), 
             size = point_size, alpha = 0.4, color = I("#a65628")) + 
  geom_point(data = indPars_vpp, 
             aes(x = scale(w_real), y = scale_rel(w,w_real)), 
             size = point_size, alpha = 0.4, color = I("#f781bf")) +
  coord_cartesian(xlim = c(-4,4), ylim = c(-4,4)) +
  geom_abline(intercept = 0, slope = 1) + 
  geom_abline(intercept = 1, slope = 1, lty = 2, color = "black") + 
  geom_abline(intercept = 2, slope = 1, lty = 3, color = "black") + 
  geom_abline(intercept = -1, slope = 1, lty = 2, color = "black") +
  geom_abline(intercept = -2, slope = 1, lty = 3, color = "black") +
  theme(legend.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  lazerhawk::theme_trueMinimal()

h3 <- ggplot() +
  geom_point(data = indPars_orl, 
             aes(x = scale(Arew_real), y = scale_rel(Arew,Arew_real)), 
             size = point_size, alpha = 0.4, color = I("#e41a1c")) + 
  geom_point(data = indPars_orl, 
             aes(x = scale(Apun_real), y = scale_rel(Apun,Apun_real)), 
             size = point_size, alpha = 0.4, color = I("#377eb8")) + 
  geom_point(data = indPars_orl, 
             aes(x = scale(K_real), y = scale_rel(K,K_real)), 
             size = point_size, alpha = 0.4, color = I("#4daf4a")) + 
  geom_point(data = indPars_orl, 
             aes(x = scale(betaF_real), y = scale_rel(betaF,betaF_real)), 
             size = point_size, alpha = 0.4, color = I("#984ea3")) + 
  geom_point(data = indPars_orl, 
             aes(x = scale(betaP_real), y = scale_rel(betaP,betaP_real)), 
             size = point_size, alpha = 0.4, color = I("#ff7f00")) + 
  coord_cartesian(xlim = c(-4,4), ylim = c(-4,4)) +
  geom_abline(intercept = 0, slope = 1) + 
  geom_abline(intercept = 1, slope = 1, lty = 2, color = "black") + 
  geom_abline(intercept = 2, slope = 1, lty = 3, color = "black") + 
  geom_abline(intercept = -1, slope = 1, lty = 2, color = "black") +
  geom_abline(intercept = -2, slope = 1, lty = 3, color = "black") +
  theme(legend.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  lazerhawk::theme_trueMinimal()

# ORL
cor_method <- "spearman"
orl_cor <- data.frame(Arew = apply(pars_orl$Arew, 1, 
                                   function(x) {cor(x,orl_AhnReal$Arew, method = cor_method)}),
                      Apun = apply(pars_orl$Apun, 1, 
                                   function(x) {cor(x,orl_AhnReal$Apun, method = cor_method)}),
                      K = apply(pars_orl$K, 1, 
                                function(x) {cor(x,orl_AhnReal$K, method = cor_method)}),
                      betaF = apply(pars_orl$betaF, 1, 
                                    function(x) {cor(x,orl_AhnReal$betaF, method = cor_method)}),
                      betaP = apply(pars_orl$betaP, 1, 
                                    function(x) {cor(x,orl_AhnReal$betaP, method = cor_method)}))
orl_cor <- gather(orl_cor)
orl_cor$key <- factor(orl_cor$key, 
                      levels = c("Arew", "Apun", "K", "betaF", "betaP"), 
                      labels = c("Arew", "Apun", "K", "betaF", "betaP"))

# PVL
pvl_cor <- data.frame(A = apply(pars_pvl$A, 1, 
                                function(x) {cor(x,pvl_AhnReal$A, method = cor_method)}),
                      alpha = apply(pars_pvl$alpha, 1, 
                                    function(x) {cor(x,pvl_AhnReal$alpha, method = cor_method)}),
                      c = apply(pars_pvl$cons, 1, 
                                function(x) {cor(x,pvl_AhnReal$cons, method = cor_method)}),
                      lambda = apply(pars_pvl$lambda, 1, 
                                     function(x) {cor(x,pvl_AhnReal$lambda, method = cor_method)}))
pvl_cor <- gather(pvl_cor)
pvl_cor$key <- factor(pvl_cor$key, 
                      levels = c("A", "alpha", "c", "lambda"), 
                      labels = c("A", "alpha", "c", "lambda"))
# Correlation for vpp
vpp_cor <- data.frame(A = apply(pars_vpp$A, 1, 
                                function(x) {cor(x,vpp_AhnReal$A)}),
                      alpha = apply(pars_vpp$alpha, 1, 
                                    function(x) {cor(x,vpp_AhnReal$alpha, method = cor_method)}),
                      c = apply(pars_vpp$cons, 1, 
                                function(x) {cor(x,vpp_AhnReal$cons, method = cor_method)}),
                      lambda = apply(pars_vpp$lambda, 1, 
                                     function(x) {cor(x,vpp_AhnReal$lambda, method = cor_method)}),
                      epP = apply(pars_vpp$epP, 1, 
                                  function(x) {cor(x,vpp_AhnReal$epP, method = cor_method)}),
                      epN = apply(pars_vpp$epN, 1, 
                                  function(x) {cor(x,vpp_AhnReal$epN, method = cor_method)}),
                      K = apply(pars_vpp$K, 1, 
                                function(x) {cor(x,vpp_AhnReal$K, method = cor_method)}),
                      w = apply(pars_vpp$w, 1, 
                                function(x) {cor(x,vpp_AhnReal$w, method = cor_method)}))
vpp_cor <- gather(vpp_cor)
vpp_cor$key <- factor(vpp_cor$key, 
                      levels = c("A", "alpha", "c", "lambda", "epP", "epN", "K", "w"), 
                      labels = c("A", "alpha", "c", "lambda", "epP", "epN", "K", "w"))

# Plots
p1 <- ggplot(pvl_cor, aes(x = value, group = key, fill = key)) + 
  geom_density(alpha = 0.4) +
  scale_fill_brewer(type = "qual", palette = 6, 
                    labels = c(expression(paste(A)), expression(paste(phantom(.),alpha)), 
                               expression(paste(c)), expression(paste(lambda)))) + 
  geom_vline(color = I("black"), linetype = 3, xintercept = quantile(pvl_cor$value, probs = 0.05)) +
  geom_vline(color = I("black"), linetype = 2, xintercept = quantile(pvl_cor$value, probs = 0.25)) +
  geom_vline(color = I("black"), linetype = 1, xintercept = quantile(pvl_cor$value, probs = 0.5)) +
  geom_vline(color = I("black"), linetype = 2, xintercept = quantile(pvl_cor$value, probs = 0.75)) +
  geom_vline(color = I("black"), linetype = 3, xintercept = quantile(pvl_cor$value, probs = 0.95)) +
  coord_cartesian(xlim = c(-0.5,1.1), ylim = c(0,22)) + 
  theme(legend.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  lazerhawk::theme_trueMinimal()
p2 <- ggplot(vpp_cor, aes(x = value, group = key, fill = key)) + 
  geom_density(alpha = 0.4) +
  scale_fill_brewer(type = "qual", palette = 6, 
                    labels = c(expression(A), expression(paste(alpha)), 
                               expression(paste(c)), expression(paste(lambda)),
                               expression(paste(epsilon["P"])), expression(paste(epsilon["N"])), 
                               expression(paste(K,phantom(.))), expression(paste(omega)))) + 
  geom_vline(color = I("black"), linetype = 3, xintercept = quantile(vpp_cor$value, probs = 0.05)) +
  geom_vline(color = I("black"), linetype = 2, xintercept = quantile(vpp_cor$value, probs = 0.25)) +
  geom_vline(color = I("black"), linetype = 1, xintercept = quantile(vpp_cor$value, probs = 0.5)) +
  geom_vline(color = I("black"), linetype = 2, xintercept = quantile(vpp_cor$value, probs = 0.75)) +
  geom_vline(color = I("black"), linetype = 3, xintercept = quantile(vpp_cor$value, probs = 0.95)) +
  coord_cartesian(xlim = c(-0.5,1.1), ylim = c(0,22)) + 
  theme(legend.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  lazerhawk::theme_trueMinimal()
p3 <- ggplot(orl_cor, aes(x = value, group = key, fill = key)) + 
  geom_density(alpha = 0.4) +
  scale_fill_brewer(type = "qual", palette = 6,
                    labels = c(expression(paste(A["+"])), expression(paste(A["-"])), 
                               expression(paste(K,phantom(.))),
                               expression(paste(beta["F"])), expression(paste(beta["P"])))) + 
  geom_vline(color = I("black"), linetype = 3, xintercept = quantile(orl_cor$value, probs = 0.05)) +
  geom_vline(color = I("black"), linetype = 2, xintercept = quantile(orl_cor$value, probs = 0.25)) +
  geom_vline(color = I("black"), linetype = 1, xintercept = quantile(orl_cor$value, probs = 0.5)) +
  geom_vline(color = I("black"), linetype = 2, xintercept = quantile(orl_cor$value, probs = 0.75)) +
  geom_vline(color = I("black"), linetype = 3, xintercept = quantile(orl_cor$value, probs = 0.95)) +
  coord_cartesian(xlim = c(-0.5,1.1), ylim = c(0,22)) + 
  theme(legend.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  lazerhawk::theme_trueMinimal()

# Plot everything
all_plots <- cowplot::plot_grid(h3,p3,h1,p1,h2,p2, ncol = 2, rel_widths = c(1,2))
ggsave(filename="~/Desktop/recover_mod.pdf", plot = all_plots, units = "in",
       height = 7.97, width = 8.92)
