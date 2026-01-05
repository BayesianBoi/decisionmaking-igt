rm(list=ls())

library(ggplot2)
library(cowplot)
library(dplyr)
library(doParallel)
library(foreach)
library(abind)
library(zoo)

source("R/utils/simulate_igt_master.R")
source("R/utils/set_version_igt.R")
source("R/utils/plot_IGT_dat.R")

# Run parameters
set.seed(43210)
models <- c("orl", "pvl_delta", "vpp")
iters  <- 1000
cores  <- 40

# Behavioral data
behav_dat <- readRDS("Data/Preprocessed/all_dat_12_studies_IGT_Data_stan_ready.rds")

# Max trials for each dataset
all_max_trials <- c(60, 60, 60, 100, 60, 60, 60, 100, 100)

results <- foreach(m=models) %do% {
  if (m=="orl") {
    fits <- readRDS("Data/Fitted/fits_all_orl.rds")
  } else if (m=="pvl_delta") {
    fits <- readRDS("Data/Fitted/fits_all_pvl_delta.rds")
  } else if (m=="vpp") {
    fits <- readRDS("Data/Fitted/fits_all_vpp.rds")
  }
  # Rename for later (indexing was wrong before)
  names(fits) <- unique(behav_dat$Study)
  names(fits)[4] <- "Premkumar"
  names(fits)[5] <- "SteingroverInPrep"
    
  # Remove studies with randomized card order
  fits$Horstmann <- NULL
  fits$SteingroverInPrep <- NULL
    
  # Remove study where subjects asked to introspect
  fits$Maia <- NULL
  
  # For storage
  plot_behav <- list()
  plot_sim   <- list()
  MSD <- vector(length=length(fits))
  
  for (i in seq_along(names(fits))) {
    plot_behav[[i]] <- plot_IGT_dat(data=behav_dat %>% 
                                      filter(Study==names(fits)[i]), 
                                    ylim = c(0,0.5))
    plot_sim[[i]] <- simulate_igt(dat = fits[[names(fits)[i]]], 
                                  model = paste0(m, "_maxDeck"),
                                  title = m,
                                  max_trials = all_max_trials[i],
                                  dataset = names(fits)[i],
                                  iters = iters,
                                  n_cores = cores,
                                  pred_den = NULL)
    avg <- behav_dat %>% 
      filter(Study==names(fits)[i]) %>%
      group_by(trial) %>% 
      summarize(A = mean(deck==1),
                B = mean(deck==2),
                C = mean(deck==3),
                D = mean(deck==4)) %>% 
      mutate(A = rollapply(A, FUN = mean, width = 7, partial = T, align = "center"),
             B = rollapply(B, FUN = mean, width = 7, partial = T, align = "center"),
             C = rollapply(C, FUN = mean, width = 7, partial = T, align = "center"),
             D = rollapply(D, FUN = mean, width = 7, partial = T, align = "center"))
    sim_A <- (plot_sim[[i]]$data$mean[plot_sim[[i]]$data$Deck==1][1:length(avg$A)]*100 - avg$A*100)^2
    sim_B <- (plot_sim[[i]]$data$mean[plot_sim[[i]]$data$Deck==2][1:length(avg$A)]*100 - avg$B*100)^2
    sim_C <- (plot_sim[[i]]$data$mean[plot_sim[[i]]$data$Deck==3][1:length(avg$A)]*100 - avg$C*100)^2
    sim_D <- (plot_sim[[i]]$data$mean[plot_sim[[i]]$data$Deck==4][1:length(avg$A)]*100 - avg$D*100)^2
    MSD[i] <- sum(sim_A, sim_B, sim_C, sim_D) / (4 * length(avg$A))
  }
  all_data <- list(behavior   = plot_behav,
                   simulation = plot_sim,
                   MSD        = MSD)
  return(all_data)
}

plot_MSD <- data.frame(Model = c(rep(models[1], 9),
                                 rep(models[2], 9),
                                 rep(models[3], 9)),
                       Study = rep(names(fits), 3),
                       MSD   = c(results[[1]]$MSD,
                                 results[[2]]$MSD,
                                 results[[3]]$MSD))

# MSD performance
ggplot(plot_MSD, aes(x = as.factor(Study), y = MSD, fill = Model)) +
  geom_bar(position = position_dodge(.5), stat = "identity", width = .5) +
  scale_fill_manual(values = c("green", "red", "blue")) + 
  theme(axis.text.x = element_text(angle = 25, vjust = .7))

for (i in 1:9) {
  results[[1]]$behavior[[i]]   <- results[[1]]$behavior[[i]] + xlab("") + ggtitle("")
  results[[1]]$simulation[[i]] <- results[[1]]$simulation[[i]] + xlab("") + ggtitle("")
  results[[2]]$simulation[[i]] <- results[[2]]$simulation[[i]] + xlab("") + ggtitle("")
  results[[3]]$simulation[[i]] <- results[[3]]$simulation[[i]] + xlab("") + ggtitle("")
}

results[[1]]$behavior[[8]]   <- results[[1]]$behavior[[8]] + coord_cartesian(ylim = c(0,.6))
results[[1]]$simulation[[8]] <- results[[1]]$simulation[[8]] + coord_cartesian(ylim = c(0,.6))
results[[2]]$simulation[[8]] <- results[[2]]$simulation[[8]] + coord_cartesian(ylim = c(0,.6))
results[[3]]$simulation[[8]] <- results[[3]]$simulation[[8]] + coord_cartesian(ylim = c(0,.6))

# Behavioral versus simulated plots
p1 <- plot_grid(plotlist = results[[1]]$behavior, ncol = 1)
p2 <- plot_grid(plotlist = results[[1]]$simulation, ncol = 1)
p3 <- plot_grid(plotlist = results[[2]]$simulation, ncol = 1)
p4 <- plot_grid(plotlist = results[[3]]$simulation, ncol = 1)

all_plots <- plot_grid(p1, p2, p3, p4, ncol = 4)
ggsave(all_plots, filename = "~/Desktop/all_sims_maxDeck.pdf", units = "in", 
       height = 30, width = 15)
