rm(list=ls())

library(rethinking)
library(rstan)
library(foreach)
library(loo)
library(cowplot)
source("R/utils/compare.R")


# Healthy controls
orl_fits <- readRDS("Data/Fitted/fits_all_orl.rds")
pvl_fits <- readRDS("Data/Fitted/fits_all_pvl_delta.rds")
vpp_fits <- readRDS("Data/Fitted/fits_all_vpp.rds")

# Remove excluded datasets
orl_fits[[1]] <- pvl_fits[[1]] <- vpp_fits[[1]] <- NULL
orl_fits[[3]] <- pvl_fits[[3]] <- vpp_fits[[3]] <- NULL
orl_fits[[5]] <- pvl_fits[[5]] <- vpp_fits[[5]] <- NULL

plots <- foreach(i=1:9) %do% {
  x <- compare(vpp_fits[[i]], pvl_fits[[i]], orl_fits[[i]])
  y <- data.frame(Model = factor(x      = c("VPP", "PVL-Delta", "ORL"), 
                                 levels = c("VPP", "PVL-Delta", "ORL")),
                  LOOIC = x@output$LOO - min(x@output$LOO),
                  sem = x@output$dSE)
  if (i %in% c(1, 4, 7)) {
    axis.y = NULL
  } else {
    axis.y = element_blank()
  }
  z <- ggplot(data = y, aes(x = Model, y = LOOIC, fill = Model)) +
    geom_bar(stat="identity") + 
    geom_errorbar(aes(ymin = LOOIC - sem*2, ymax = LOOIC + sem*2), width = 0.2) + 
    geom_hline(yintercept = 0, lty = 2) + 
    scale_fill_manual(values = c("blue", "red", "green")) + 
    ylab("") +
    xlab("") + 
    theme_minimal() + 
    theme(text = element_text(size=20),
          axis.text.y = axis.y,
          legend.position = "none") + 
    coord_flip()
}

# All plots
all_plots <- plot_grid(plotlist = plots, ncol = 3, 
                       rel_widths = c(rep(c(1.2,1,1), 3)))

ggsave("~/Desktop/loo.pdf", plot = all_plots, units = "in", 
       height = 9.79 * .75, width = 16.61 * .75)


# Make tables
tab_mod <- data.frame(Model = c("VPP", "ORL", "PVL-Delta"),
                      LOOIC = round(comp_dat_mod@output$LOO, 2),
                      dLOOIC = round(comp_dat_mod@output$dLOO, 2),
                      dSE = round(comp_dat_mod@output$dSE))
tab_org <- data.frame(Model = c("VPP", "ORL", "PVL-Delta"),
                      LOOIC = round(comp_dat_org@output$LOO, 2),
                      dLOOIC = round(comp_dat_org@output$dLOO, 2),
                      dSE = round(comp_dat_org@output$dSE))
xtable::xtable(tab_mod)
xtable::xtable(tab_org)
