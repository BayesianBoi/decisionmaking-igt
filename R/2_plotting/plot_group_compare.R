rm(list=ls())

library(rstan)
library(reshape2)
library(hBayesDM)
library(cowplot)

source("R/utils/stan_rhat.R")
source("R/utils/plot_diff_HDI_info.R")
source("R/utils/plot_diff.R")

# Load in all model fits
orl_fits <- readRDS("Data/Fitted/fits_all_orl.rds")

# Parameter values
orl_Ahn_HC_parVals    <- rstan::extract(orl_fits[[8]])[13:17]
orl_Ahn_Amph_parVals  <- rstan::extract(orl_fits[[9]])[13:17]
orl_Ahn_Hero_parVals  <- rstan::extract(orl_fits[[10]])[13:17]
orl_Frid_HC_parVals   <- rstan::extract(orl_fits[[11]])[13:17]
orl_Frid_Cbis_parVals <- rstan::extract(orl_fits[[12]])[13:17]


title <- c(expression(A["+"]), expression(A["-"]), expression(K), 
           expression(beta[F]), expression(beta[P]))
# HC - Amph
p1 <- plot_diff_info(orl_Ahn_HC_parVals, orl_Ahn_Amph_parVals, title = title, cols = 1)
# HC - Hero
p2 <- plot_diff_info(orl_Ahn_HC_parVals, orl_Ahn_Hero_parVals, title = title, cols = 1)
# Amph - Hero
p3 <- plot_diff_info(orl_Ahn_Amph_parVals, orl_Ahn_Hero_parVals, title = title, cols = 1)
# Frid HC - Cbis
p4 <- plot_diff_info(orl_Frid_HC_parVals, orl_Frid_Cbis_parVals, title = title, cols = 1)

p_all <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 4)
ggsave(filename="~/Desktop/group_compare.pdf", plot = p_all)

mu_pars <- c("mu_Arew", "mu_Apun", "mu_K", "mu_betaF", "mu_betaP")
groups  <- c("Control", "Amphetamine", "Heroin")
mu_Arew <- mu_Apun <- mu_Adiff <- mu_K <- mu_betaF <- mu_betaP <- data.frame(Arew = rep(NA, 10000), 
                                                                             Apun = rep(NA, 10000),
                                                                             Adiff = rep(NA, 10000),
                                                                             K = rep(NA, 10000),
                                                                             betaF = rep(NA, 10000),
                                                                             betaP = rep(NA, 10000))
all_dat <- list()
all_fits <- list(Control = orl_Ahn_HC_parVals, 
                 Amphetamine = orl_Ahn_Amph_parVals, 
                 Heroin = orl_Ahn_Hero_parVals)
for (i in mu_pars) {
  for (j in groups) {
    all_dat[[i]][[j]] <- all_fits[[j]][[i]]
  }
}

Arew_m <- melt(all_dat$mu_Arew); names(Arew_m) <- c("iter", "mu_Arew", "Group")
Apun_m <- melt(all_dat$mu_Apun); names(Apun_m) <- c("iter", "mu_Apun", "Group")
Adiff_m <- melt(list(Control = all_dat$mu_Arew$Control - all_dat$mu_Apun$Control,
                     Amphetamine = all_dat$mu_Arew$Amphetamine - all_dat$mu_Apun$Amphetamine,
                     Heroin = all_dat$mu_Arew$Heroin - all_dat$mu_Apun$Heroin)); names(Adiff_m) <- c("iter", "mu_Adiff", "Group")
K_m <- melt(all_dat$mu_K); names(K_m) <- c("iter", "mu_K", "Group")
betaF_m <- melt(all_dat$mu_betaF); names(betaF_m) <- c("iter", "mu_betaF", "Group")
betaP_m <- melt(all_dat$mu_betaP); names(betaP_m) <- c("iter", "mu_betaP", "Group")
Arew_m$Group <- Apun_m$Group <- K_m$Group <- betaF_m$Group <- betaP_m$Group <-  Adiff_m$Group <-
  factor(Arew_m$Group, 
         levels = c("Control", "Amphetamine", "Heroin"), 
         labels = c("Control", "Amphetamine", "Heroin"))

h1 <- ggplot(Arew_m, aes(x = mu_Arew, group = Group, color = Group, fill = Group)) +
  scale_fill_manual(values = c("Control" = "#fb9a99", "Amphetamine" = "#1f78b4", "Heroin" = "#542788")) + 
  scale_color_manual(values = c("Control" = "#fb9a99", "Amphetamine" = "#1f78b4", "Heroin" = "#542788")) + 
  geom_density(alpha = 0.5) + 
  ggtitle(expression(A["+"])) + 
  xlab("") +
  ylab("") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 20))
h2 <- ggplot(Apun_m, aes(x = mu_Apun, group = Group, color = Group, fill = Group)) +
  scale_fill_manual(values = c("Control" = "#fb9a99", "Amphetamine" = "#1f78b4", "Heroin" = "#542788")) + 
  scale_color_manual(values = c("Control" = "#fb9a99", "Amphetamine" = "#1f78b4", "Heroin" = "#542788")) + 
  geom_density(alpha = 0.5) + 
  ggtitle(expression(A["-"])) + 
  xlab("") +
  ylab("") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 20))
h3 <- ggplot(K_m, aes(x = mu_K, group = Group, color = Group, fill = Group)) +
  scale_fill_manual(values = c("Control" = "#fb9a99", "Amphetamine" = "#1f78b4", "Heroin" = "#542788")) + 
  scale_color_manual(values = c("Control" = "#fb9a99", "Amphetamine" = "#1f78b4", "Heroin" = "#542788")) + 
  geom_density(alpha = 0.5) + 
  ggtitle(expression(K)) + 
  xlab("") +
  ylab("") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 20))
h4 <- ggplot(betaF_m, aes(x = mu_betaF, group = Group, color = Group, fill = Group)) +
  scale_fill_manual(values = c("Control" = "#fb9a99", "Amphetamine" = "#1f78b4", "Heroin" = "#542788")) + 
  scale_color_manual(values = c("Control" = "#fb9a99", "Amphetamine" = "#1f78b4", "Heroin" = "#542788")) + 
  geom_density(alpha = 0.5) + 
  ggtitle(expression(beta[F])) + 
  xlab("") +
  ylab("") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 20))
h5 <- ggplot(betaP_m, aes(x = mu_betaP, group = Group, color = Group, fill = Group)) +
  scale_fill_manual(values = c("Control" = "#fb9a99", "Amphetamine" = "#1f78b4", "Heroin" = "#542788")) + 
  scale_color_manual(values = c("Control" = "#fb9a99", "Amphetamine" = "#1f78b4", "Heroin" = "#542788")) + 
  geom_density(alpha = 0.5) + 
  ggtitle(expression(beta[P])) + 
  xlab("") +
  ylab("") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 20))
ggplot(Adiff_m, aes(x = mu_Adiff, group = Group, color = Group, fill = Group)) +
  scale_fill_manual(values = c("Control" = "#fb9a99", "Amphetamine" = "#1f78b4", "Heroin" = "#542788")) + 
  scale_color_manual(values = c("Control" = "#fb9a99", "Amphetamine" = "#1f78b4", "Heroin" = "#542788")) + 
  geom_density(alpha = 0.5) + 
  ggtitle(expression(A[diff])) + 
  xlab("") +
  ylab("") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 20))
h_key <- ggplot(betaP_m, aes(x = mu_betaP, group = Group, color = Group, fill = Group)) +
  scale_fill_manual(values = c("Control" = "#fb9a99", "Amphetamine" = "#1f78b4", "Heroin" = "#542788")) + 
  scale_color_manual(values = c("Control" = "#fb9a99", "Amphetamine" = "#1f78b4", "Heroin" = "#542788")) + 
  geom_density(alpha = 0.5) + 
  ggtitle(expression(beta[P])) + 
  xlab("") +
  ylab("") +
  theme(plot.title = element_text(size = 20))
  

h_all <- plot_grid(h1, h2, h3, h4, h5, ncol = 3)
ggsave(filename="~/Desktop/group_pars_ahn.pdf", plot = h_all)
ggsave(filename="~/Desktop/group_pars_ahn_key.pdf", plot = h_key)

### For original IGT
mu_pars <- c("mu_Arew", "mu_Apun", "mu_K", "mu_betaF", "mu_betaP")
groups  <- c("Control", "Cannabis")
mu_Arew <- mu_Apun <- mu_K <- mu_betaF <- mu_betaP <- data.frame(Arew = rep(NA, 10000), 
                                                                 Apun = rep(NA, 10000),
                                                                 K = rep(NA, 10000),
                                                                 betaF = rep(NA, 10000),
                                                                 betaP = rep(NA, 10000))
all_dat <- list()
all_fits <- list(Control = orl_Frid_HC_parVals, 
                 Cannabis = orl_Frid_Cbis_parVals)
for (i in mu_pars) {
  for (j in groups) {
    all_dat[[i]][[j]] <- all_fits[[j]][[i]]
  }
}

Arew_m <- melt(all_dat$mu_Arew); names(Arew_m) <- c("iter", "mu_Arew", "Group")
Apun_m <- melt(all_dat$mu_Apun); names(Apun_m) <- c("iter", "mu_Apun", "Group")
Adiff_m <- melt(list(Control = all_dat$mu_Arew$Control - all_dat$mu_Apun$Control,
                     Cannabis = all_dat$mu_Arew$Cannabis - all_dat$mu_Apun$Cannabis)); names(Adiff_m) <- c("iter", "mu_Adiff", "Group")
K_m <- melt(all_dat$mu_K); names(K_m) <- c("iter", "mu_K", "Group")
betaF_m <- melt(all_dat$mu_betaF); names(betaF_m) <- c("iter", "mu_betaF", "Group")
betaP_m <- melt(all_dat$mu_betaP); names(betaP_m) <- c("iter", "mu_betaP", "Group")
Arew_m$Group <- Apun_m$Group <- K_m$Group <- betaF_m$Group <- betaP_m$Group <- Adiff_m$Group <-
  factor(Arew_m$Group, 
         levels = c("Control", "Cannabis"), 
         labels = c("Control", "Cannabis"))


h1 <- ggplot(Arew_m, aes(x = mu_Arew, group = Group, color = Group, fill = Group)) +
  scale_fill_manual(values = c("Control" = "#fb9a99", "Cannabis" = "#33a02c")) + 
  scale_color_manual(values = c("Control" = "#fb9a99", "Cannabis" = "#33a02c")) + 
  geom_density(alpha = 0.5) + 
  ggtitle(expression(A["+"])) + 
  xlab("") +
  ylab("") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 20))
h2 <- ggplot(Apun_m, aes(x = mu_Apun, group = Group, color = Group, fill = Group)) +
  scale_fill_manual(values = c("Control" = "#fb9a99", "Cannabis" = "#33a02c")) + 
  scale_color_manual(values = c("Control" = "#fb9a99", "Cannabis" = "#33a02c")) + 
  geom_density(alpha = 0.5) + 
  ggtitle(expression(A["-"])) + 
  xlab("") +
  ylab("") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 20))
h3 <- ggplot(K_m, aes(x = mu_K, group = Group, color = Group, fill = Group)) +
  scale_fill_manual(values = c("Control" = "#fb9a99", "Cannabis" = "#33a02c")) + 
  scale_color_manual(values = c("Control" = "#fb9a99", "Cannabis" = "#33a02c")) + 
  geom_density(alpha = 0.5) + 
  ggtitle(expression(K)) + 
  xlab("") +
  ylab("") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 20))
h4 <- ggplot(betaF_m, aes(x = mu_betaF, group = Group, color = Group, fill = Group)) +
  scale_fill_manual(values = c("Control" = "#fb9a99", "Cannabis" = "#33a02c")) + 
  scale_color_manual(values = c("Control" = "#fb9a99", "Cannabis" = "#33a02c")) + 
  geom_density(alpha = 0.5) + 
  ggtitle(expression(beta[F])) + 
  xlab("") +
  ylab("") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 20))
h5 <- ggplot(betaP_m, aes(x = mu_betaP, group = Group, color = Group, fill = Group)) +
  scale_fill_manual(values = c("Control" = "#fb9a99", "Cannabis" = "#33a02c")) + 
  scale_color_manual(values = c("Control" = "#fb9a99", "Cannabis" = "#33a02c")) + 
  geom_density(alpha = 0.5) + 
  ggtitle(expression(beta[P])) + 
  xlab("") +
  ylab("") +
  theme(legend.position = "none", 
        plot.title = element_text(size = 20))
h_key <- ggplot(betaP_m, aes(x = mu_betaP, group = Group, color = Group, fill = Group)) +
  scale_fill_manual(values = c("Control" = "#fb9a99", "Cannabis" = "#33a02c")) + 
  scale_color_manual(values = c("Control" = "#fb9a99", "Cannabis" = "#33a02c")) + 
  geom_density(alpha = 0.5) + 
  ggtitle(expression(beta[P])) + 
  xlab("") +
  ylab("") +
  theme(plot.title = element_text(size = 20))


h_all <- plot_grid(h1, h2, h3, h4, h5, ncol = 3)
ggsave(filename="~/Desktop/group_pars_frid.pdf", plot = h_all)
ggsave(filename="~/Desktop/group_pars_frid_key.pdf", plot = h_key)


