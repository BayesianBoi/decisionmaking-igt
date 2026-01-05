# Function that plots IGT performance across subjects

plot_IGT_dat <- function(file      = NULL,
                         data      = NULL,
                         nBlocks   = 5, 
                         pointSize = 2,
                         lineSize  = 1,
                         title     = "IGT Choices",
                         ylim      = c(0,0.5),
                         xlab      = "Block",
                         ylab      = "Proportion of Choices from Each Deck",
                         xyAxis    = c(1,1),
                         dodge     = 1, 
                         simDat    = NULL) {
  
  source("~/Dropbox/Box/CCSL/Projects/IGT/manuscript/Code/R/utils/IGTbehav_v2.R")
  if (is.null(file)) {
    dat <- data
  } else {
    dat <- read.table(file = file, header = T)
  }
  numSubjs <- length(unique(dat$subjID))
  
  dat <- IGTbehav_v2(dat)
  
  # ggplot2 version
  meanSubj = apply(dat, c(1,2), mean)
  sdSubj = apply(dat, c(1,2), sd) 
  semSubj = apply(dat, c(1,2), sd) / sqrt(numSubjs)
  pCorrPlot = data.frame(
    Block=c(seq(10,90,20), seq(10,90,20), seq(10,90,20), seq(10,90,20)),
    mean = c( meanSubj[, 1], meanSubj[, 2], meanSubj[, 3], meanSubj[, 4] ),
    sd = c( sdSubj[, 1], sdSubj[, 2], sdSubj[, 3], sdSubj[, 4]),
    sem = c( semSubj[,1], semSubj[,2], semSubj[,3], semSubj[,4] ),
    Deck = as.factor( c(rep(1, nBlocks), rep(2, nBlocks), rep(3, nBlocks), rep(4, nBlocks)) )
  )
  
  if (xyAxis[1]==0) {
    xAx <- element_blank()
  } else {
    xAx <- element_text(size = 20)
  }
  if (xyAxis[2]==0) {
    yAx <- element_blank()
  } else {
    yAx <- element_text(size = 20)
  }
  
  # Plotting parameters
  eb <- aes(ymin = mean- 1*sem, ymax = mean + 1*sem)
  
  f1 = ggplot(data = pCorrPlot, aes(x = Block, y = mean, group=Deck, colour=Deck) ) + 
    theme_bw() +
    scale_colour_manual(values = c("red", "orange", "blue", "green" ), 
                        labels=c("Deck A", "Deck B", "Deck C", "Deck D") ) +
    geom_line(size=lineSize) +
    geom_point(size=pointSize) +
    geom_ribbon(eb, alpha = 0.25, colour=NA ) +
    ggtitle("Choice Behavior") +
    xlab("Block\n") + 
    coord_cartesian(xlim = c(0,100), ylim = ylim) + 
    scale_y_continuous(breaks = seq(0, 0.6, by=0.1)) +
    scale_x_discrete(breaks = seq(0, 100, by=20)) +
    theme_minimal() + 
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          axis.title.x = element_text(size = 12),
          axis.ticks.x = element_line(inherit.blank = F),
          axis.text.x = element_text(size = 15),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 15),
          legend.position = "none")
  return(f1)
}