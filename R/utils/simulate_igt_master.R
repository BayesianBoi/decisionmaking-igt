simulate_igt <- function(dat        = NULL, 
                         file       = NULL, 
                         model      = "orl",
                         dataset    = NULL,
                         title      = NULL,
                         n_trials   = 100,
                         max_trials = 100,
                         n_cores    = 4,
                         pred_den   = "50", 
                         iters      = 100,
                         y_lim      = c(0, 0.5)) {
  
  # Find class of object
  if (class(dat)[1] == "stanfit") {
    simul_pars <- rstan::extract(dat)
    num_subjs  <- dim(simul_pars$log_lik)[2]
    samp_vec   <- 1:dim(simul_pars$log_lik)[1]
  } else if (class(dat)[1] == "hBayesDM") {
    simul_pars <- dat$parVals
    num_subjs  <- dim(simul_pars$log_lik)[2]
    samp_vec   <- 1:dim(simul_pars$log_lik)[1]
  } 
  
  # For using multiple cores and combining into array
  acomb <- function(...) abind(..., along=4)
  registerDoParallel(cores=n_cores)
  
  pDeck <- foreach(subj=1:num_subjs, .combine = "acomb", .multicombine = T) %dopar% {
    # Probability of choosing deck array
    pDeck_subj <- array(.25, c(n_trials, 4, iters))
    
    # Loop through all iterations
    for (iter in 1:iters) {
      all_decks <- set_version_igt(max_trials = max_trials, 
                                   dataset = dataset,
                                   payscale = 100)
      # Run one iteration of model
      source(paste0("R/1_analysis/Simulation/Models/", model, ".R"), local = T)
    }
    return(pDeck_subj)
  }
  stopImplicitCluster()
  
  # ggplot2 version
  meanIter = apply(pDeck, c(1,2,4), mean)
  meanSubj = apply(meanIter, c(1,2), mean)
  sdSubj = apply(meanIter, c(1,2), sd) 
  q_05  = apply(meanIter, c(1,2), quantile, probs = 0.05) 
  q_1   = apply(meanIter, c(1,2), quantile, probs = 0.1) 
  q_25  = apply(meanIter, c(1,2), quantile, probs = 0.25) 
  q_5   = apply(meanIter, c(1,2), quantile, probs = 0.5) 
  q_75  = apply(meanIter, c(1,2), quantile, probs = 0.75) 
  q_9   = apply(meanIter, c(1,2), quantile, probs = 0.9) 
  q_95 = apply(meanIter, c(1,2), quantile, probs = 0.95) 
  semSubj = apply(meanIter, c(1,2), sd) / sqrt( num_subjs)
  pCorrPlot = data.frame(
    Trial=c(1:n_trials, 1:n_trials, 1:n_trials, 1:n_trials),
    mean = c( meanSubj[, 1], meanSubj[, 2], meanSubj[, 3], meanSubj[, 4] ),
    sd = c( sdSubj[, 1], sdSubj[, 2], sdSubj[, 3], sdSubj[, 4]),
    sem = c( semSubj[,1], semSubj[,2], semSubj[,3], semSubj[,4] ),
    q_05 = c( q_05[,1], q_05[,2], q_05[,3], q_05[,4]),
    q_1 = c( q_1[,1], q_1[,2], q_1[,3], q_1[,4]),
    q_25 = c( q_25[,1], q_25[,2], q_25[,3], q_25[,4]),
    q_5 = c( q_5[,1], q_5[,2], q_5[,3], q_5[,4]),
    q_75 = c( q_75[,1], q_75[,2], q_75[,3], q_75[,4]),
    q_9 = c( q_9[,1], q_9[,2], q_9[,3], q_9[,4]),
    q_95 = c( q_95[,1], q_95[,2], q_95[,3], q_95[,4]),
    Deck = as.factor( c(rep(1, n_trials), rep(2, n_trials), rep(3, n_trials), rep(4, n_trials)) )
  )
  
  if (!is.null(pred_den)) {
    eb <- switch(pred_den,
                 "50" = aes(x = Trial, y = mean, group = Deck, colour = Deck, 
                            ymax = q_75, ymin = q_25),
                 "80" = aes(x = Trial, y = mean, group = Deck, fill = Deck, 
                            ymax = q_9, ymin = q_1),
                 "90" = aes(x = Trial, y = mean, group = Deck, fill = Deck, 
                            ymax = q_95, ymin = q_05))
  } else {
    eb <- aes(x = Trial, y = mean, group = Deck, colour = Deck, 
              ymax = mean + 2*sem, ymin = mean - 2*sem)
  }
  
  f1 <- ggplot(data = pCorrPlot, aes(x = Trial, y = mean, group=Deck, colour=Deck) ) + 
    scale_colour_manual(values = c("red", "orange", "blue", "green" ), 
                        labels=c("Deck A", "Deck B", "Deck C", "Deck D") ) +
    scale_fill_manual(values = c("red", "orange", "blue", "green" ), 
                      labels=c("Deck A", "Deck B", "Deck C", "Deck D") ) +
    geom_line(size = 2) +
    geom_ribbon(eb, alpha = 0.25, colour = NA) +
    coord_cartesian(ylim = y_lim) + 
    scale_y_continuous(breaks = seq(0, 0.6, by=0.1)) +
    scale_x_discrete(breaks = seq(0, n_trials, by=20)) +
    theme_minimal() + 
    ggtitle(title) + 
    xlab("Trial\n") + 
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 15),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 15),
          legend.position = "none")
  return(f1)
}